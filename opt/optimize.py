"""
Optimize the LoopGap antenna object for the following:
Cost equation:
1. Maximum B1 per sqrt(W) at the diamond surface
Hard thresholds:
2. S11 <= -10dB across +- 0.1 GHz centered on 2.87 Ghz
"""

import os, sys, ctypes, json, hashlib, argparse
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import cma
import lib.simulate as simulate, lib.postprocess as postprocess
from lib.knobs import design, KNOBS
from designs import rev2 as base
from settings import Runs, Data, sim_time

F0, HALF_BAND, S11_MAX = 2.87e9, 0.1e9, -10.0
PLATE_T = 0.5   # mm: B1 is read at the plate's top face, the weakest point

# search box in KNOBS order, mm; D2's floor is optical access, g's the fab's minimum spacing
LO = np.array([ 6.0, 0.5, 0.0, 0.1, 0.3])
HI = np.array([20.0, 2.0, 6.0, 1.0, 4.0])
THREADS = 0   # openEMS threads per solve (0 = all cores); set_threads sets it in each worker


def solve(geo) -> str:
    """Solve into a folder named by a hash of the design's params, so parallel workers never
    collide and a repeat is read back; openEMS's output goes to openems.log there.  Code isn't
    in the hash: clear runs/opt after changing the mesh or the dumps."""
    p = json.dumps(simulate.params(geo, sim_time), sort_keys=True)
    path = os.path.join(Runs, 'opt', hashlib.sha1(p.encode()).hexdigest()[:12])
    if not simulate.done(geo, path, sim_time):
        os.makedirs(path, exist_ok=True)
        # openEMS prints from C++, beyond redirect_stdout: swap the file descriptors themselves
        with open(os.path.join(path, 'openems.log'), 'w') as log:
            saved = os.dup(1), os.dup(2)
            os.dup2(log.fileno(), 1); os.dup2(log.fileno(), 2)
            try:
                simulate.run(geo, path, sim_time, THREADS)
            finally:   # flush C's buffers too, or openEMS's last lines land in the job log
                sys.stdout.flush(); sys.stderr.flush(); ctypes.CDLL(None).fflush(None)
                os.dup2(saved[0], 1); os.dup2(saved[1], 2)
                for fd in saved:
                    os.close(fd)
    return path


def cost_function(D1, D2, yc, g, tap_depth) -> float:
    """Feasible (spec 2 met): -B1.  Infeasible: how far off, positive, so every feasible design
    outranks every infeasible one with no penalty weight; fine for a rank-based optimiser."""
    try:
        geo = design(D1, D2, yc, g, tap_depth)
    except AssertionError:
        return np.inf
    path = solve(geo)
    band = np.linspace(F0 - HALF_BAND, F0 + HALF_BAND, 41)
    worst = postprocess.s11_dB(geo, path, band).max()
    if worst <= S11_MAX:
        return -postprocess.b1_at(geo, path, geo.subs_h + geo.mask_t + PLATE_T)
    # off resonance, worst sits flat near -0.5 dB: the f_res term gives the optimiser a slope
    f = np.linspace(geo.f0 - geo.fc, geo.f0 + geo.fc, 401)
    f_res = f[np.argmin(postprocess.s11_dB(geo, path, f))]
    # each term is 1 for a full unit off: total reflection, or f_res a band edge away
    return (1 - worst / S11_MAX) + abs(f_res - F0) / HALF_BAND


def knobs(u) -> np.ndarray:
    """Unit cube to knob values in mm: CMA-ES searches the cube, so one step size suits all knobs."""
    return LO + np.asarray(u) * (HI - LO)


def legal(u) -> bool:
    try:
        design(*knobs(u))
        return True
    except AssertionError:
        return False


def cost_u(u) -> float:
    return cost_function(*knobs(u))


def set_threads(n) -> None:
    global THREADS
    THREADS = n


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--popsize', type=int, default=8,
                    help='designs per generation, all solved at once (8 = the CMA-ES default for 5 knobs)')
    ap.add_argument('--threads', type=int, default=4, help='openEMS threads per solve; cores used = popsize x threads')
    ap.add_argument('--seed', type=int, default=1,
                    help='same seed and settings = same run, so rerunning a killed job replays it from the solve cache')
    ap.add_argument('--maxfevals', type=int, default=800, help='solve budget')
    a = ap.parse_args()

    u0 = (np.array([getattr(base, k) for k in KNOBS]) - LO) / (HI - LO)
    es = cma.CMAEvolutionStrategy(u0, 0.25, dict(   # a quarter of the box: the usual bounded start
        bounds=[0, 1], popsize=a.popsize, seed=a.seed, maxfevals=a.maxfevals,
        tolx=0.01,      # 1 % of each range: below that the mesh, not the design, moves the cost
        verb_log=0))    # cmaes.csv logs every design instead
    log = os.path.join(Data, 'cmaes.csv')   # rewritten: a replayed run re-logs every generation
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, 'w') as f, ProcessPoolExecutor(a.popsize, initializer=set_threads, initargs=(a.threads,)) as pool:
        f.write('gen,' + ','.join(f'{k}_mm' for k in KNOBS) + ',cost\n')
        while not es.stop():
            U = es.ask()
            for i in range(len(U)):   # resample illegal designs rather than spend a slot on them
                while not legal(U[i]):
                    U[i] = es.ask(1)[0]
            F = list(pool.map(cost_u, U))
            es.tell(U, F)
            for u, c in zip(U, F):
                f.write(f'{es.countiter},' + ','.join(f'{v:.4f}' for v in knobs(u)) + f',{c:.4f}\n')
            f.flush()
            print(f'gen {es.countiter:3d}   best {es.result.fbest:9.3f}   '
                  f'{sum(c < 0 for c in F)}/{len(F)} meet spec 2   sigma {es.sigma:.3f}', flush=True)

    best = es.result.fbest
    print(f'\nstopped: {es.stop()}')
    print('best design: ' + '   '.join(f'{k} {v:.3f}' for k, v in zip(KNOBS, knobs(es.result.xbest))) + ' mm')
    print(f'B1 {-best:.1f} uT/sqrt(W) at the plate face' if best < 0 else
          f'no design met spec 2 (best cost {best:.3f}): the band may be out of reach for this layout')
    print(f'every design: {log}')
