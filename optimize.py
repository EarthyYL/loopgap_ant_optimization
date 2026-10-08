"""
Optimize the LoopGap antenna object for the following:
Cost equation:
1. Maximum B1 per sqrt(W) at the diamond surface
Hard thresholds:
2. S11 <= -10dB across +- 0.1 GHz centered on 2.87 Ghz
"""

import os, json, hashlib, argparse
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import cma
import src.simulate as simulate, src.postprocess as postprocess
from sweep import design, KNOBS
from loopgap import geo as base, Runs, sim_time

F0, HALF_BAND, S11_MAX = 2.87e9, 0.1e9, -10.0
# mm: diamond plate resting on the board, not chosen yet (0.3-0.5 seen).  B1 is read at its far
# face: the weakest point for NVs through the whole plate.  0 reads at the board top instead,
# right for a shallow NV layer face down.
PLATE_T = 0.5

# Search box in KNOBS order, mm: D1, D2, yc, g, tap_depth.  The asserts in LoopGap reject the
# illegal combinations inside it.  D2's floor is optical access / the plate; g's is the fab's
# minimum spacing.
LO = np.array([ 6.0, 0.5, 0.0, 0.1, 0.3])
HI = np.array([20.0, 2.0, 6.0, 1.0, 4.0])
THREADS = 0   # openEMS threads per solve, 0 = all cores; the driver sets it in each worker


def solve(geo) -> str:
    """Run folder named by a hash of everything that decides the output, so parallel workers
    never collide and a design asked for twice is read back, not re-solved.  Code changes
    (mesh, dumps) are not in the hash: clear runs/opt after one of those."""
    p = json.dumps(simulate.params(geo, sim_time), sort_keys=True)
    path = os.path.join(Runs, 'opt', hashlib.sha1(p.encode()).hexdigest()[:12])   # absolute: Run chdirs
    if not simulate.done(geo, path, sim_time):
        simulate.run(geo, path, sim_time, THREADS)
    return path


def cost_function(D1, D2, yc, g, tap_depth) -> float:
    """Feasible (spec 2 met): -B1, negative.  Infeasible: how far off, positive.  The signs
    alone rank every feasible design above every infeasible one, so there is no penalty
    weight to tune -- but the jump at 0 suits rank-based optimisers (Nelder-Mead, CMA-ES,
    DE), not a Gaussian process, which wants B1 and worst-S11 as separate outputs."""
    try:
        geo = design(D1, D2, yc, g, tap_depth)
    except AssertionError:      # illegal geometry: worse than anything, and no solve spent
        return np.inf
    path = solve(geo)
    band = np.linspace(F0 - HALF_BAND, F0 + HALF_BAND, 41)    # 5 MHz bins, both edges included
    worst = postprocess.s11_dB(geo, path, band).max()
    if worst <= S11_MAX:
        return -postprocess.b1_at(geo, path, geo.subs_h + geo.mask_t + PLATE_T)   # board top + plate
    # While the resonance is outside the band, worst sits flat near -0.5 dB whatever the
    # knobs do; the f_res term is what gives the optimiser a slope there.  Ranking only: a
    # low-Q design can meet the band with its S11 minimum outside it, and that one is feasible.
    f = np.linspace(geo.f0 - geo.fc, geo.f0 + geo.fc, 401)
    f_res = f[np.argmin(postprocess.s11_dB(geo, path, f))]
    # Both terms dimensionless, 1 = a full unit of "off": the match term runs 0 at the spec to
    # 1 at total reflection (a passive port can't exceed 0 dB); the f_res term is 1 at the band edge.
    return (1 - worst / S11_MAX) + abs(f_res - F0) / HALF_BAND


def knobs(u) -> np.ndarray:
    """Unit cube -> knob values in mm.  CMA-ES searches the cube, so one step size suits
    knobs that span 0.1 to 20 mm."""
    return LO + np.asarray(u) * (HI - LO)


def legal(u) -> bool:
    try:
        design(*knobs(u))
        return True
    except AssertionError:
        return False


def cost_u(u) -> float:
    return cost_function(*knobs(u))


def set_threads(n) -> None:   # pool initializer: each worker process has its own THREADS
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

    u0 = (np.array([getattr(base, k) for k in KNOBS]) - LO) / (HI - LO)   # start at loopgap.py's design
    es = cma.CMAEvolutionStrategy(u0, 0.25, dict(   # 0.25 of the box: the usual start for a bounded search
        bounds=[0, 1], popsize=a.popsize, seed=a.seed, maxfevals=a.maxfevals,
        tolx=0.01,      # 1% of each range: below that the mesh, not the design, moves the cost
        verb_log=0))    # data/cmaes.csv logs every design instead
    # absolute, and rewritten from the top: a replayed run re-logs every generation anyway
    log = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'cmaes.csv')
    with open(log, 'w') as f, ProcessPoolExecutor(a.popsize, initializer=set_threads, initargs=(a.threads,)) as pool:
        f.write('gen,' + ','.join(f'{k}_mm' for k in KNOBS) + ',cost\n')
        while not es.stop():
            U = es.ask()
            for i in range(len(U)):   # illegal geometry costs no solve: resample it instead of spending a slot
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
