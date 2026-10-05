import os
from dataclasses import asdict
import numpy as np
import matplotlib.pyplot as plt
from src.geometry import LoopGap
from scipy.optimize import minimize
from functools import partial 
import src.simulate as simulate, src.postprocess as postprocess
from loopgap import geo as base, Sim_Path, sim_time

# The design knobs, in x-vector order.  yc (disc centre above the hole) replaces l1/l2:
# l1 + l2 == D1, so they are one knob, not two.  mask_t is a knob too, but 0-or-0.03, so it
# stays in loopgap.py.  Everything not a knob is taken from loopgap.py's geo.
KNOBS = ('D1', 'D2', 'yc', 'g', 'tap_depth')
FIXED = {k: v for k, v in asdict(base).items() if k not in ('D1', 'D2', 'l1', 'l2', 'g', 'tap_depth')}
RESONANCE_TARGET = 2.87e9

def design(D1, D2, yc, g, tap_depth) -> LoopGap:
    """LoopGap from the design knobs alone; call as design(**dict(zip(KNOBS, x)))."""
    return LoopGap(**FIXED, D1=D1, D2=D2, l1=D1 / 2 + yc, l2=D1 / 2 - yc, g=g, tap_depth=tap_depth)


def run_sim(geo: LoopGap, sim_time: float, sim_path: str, res_f: float) -> float:
    """Run one simulation and return S11 at the desired resonance, for optimization."""
    simulate.run(geo, sim_path, sim_time)
    return postprocess.s11_dB(geo, sim_path, res_f)



if __name__ == '__main__':
    # The swept knob, one of KNOBS; the others stay at loopgap.py's current design.
    knob, values = 'g', np.geomspace(0.05, 0.8, 6)
    dsgn_parameterized = partial(design, **{k: getattr(base, k) for k in KNOBS})
    s11 = np.full(len(values), np.nan)   # nan = design rejected, never solved
    # absolute, and fixed before the loop: after the first solve the cwd is that run's folder
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', f'sweep_{knob}.csv')
    # ndmin=2: a one-row file would otherwise load 1-D; nan rows were never solved, so re-run them
    saved_results = np.loadtxt(out, delimiter=',', ndmin=2) if os.path.exists(out) else np.empty((0, 2))
    saved_results = saved_results[~np.isnan(saved_results[:, 1])]
    for i, v in enumerate(values):
        try:
            geo = dsgn_parameterized(**{knob: v})
        except AssertionError as e:   # illegal geometry: skip it instead of dying mid-sweep
            print(f'{knob} {v:.3f} mm   skipped: {e}', flush=True)
            continue
        hit = np.isclose(saved_results[:, 0], v, rtol=0, atol=1e-4)   # file keeps 4 decimals, so no ==
        if np.any(hit):
            s11[i] = saved_results[hit, 1][0]
            print(f'{knob} {v:.3f} mm   S11 {s11[i]:.2f} dB   (saved in {out})', flush=True)
            continue
        path = os.path.join(Sim_Path, f'{knob}_{v:.3f}')   # absolute: FDTD.Run chdirs into it
        os.makedirs(path, exist_ok=True)                   # Run only mkdirs the last level
        s11[i] = run_sim(geo=geo, sim_time=sim_time, sim_path=path, res_f=RESONANCE_TARGET)
        print(f'{knob} {v:.3f} mm   S11 at {RESONANCE_TARGET/1e9:.2f} GHz {s11[i]:.2f} dB   ({path})', flush=True)
        # rewritten after every solve, so a crash mid-sweep keeps what finished
        np.savetxt(out, np.column_stack([values, s11]), delimiter=',', fmt='%.4f',
                   header=f'S11 at {RESONANCE_TARGET/1e9:.2f} GHz vs {knob}; nan = rejected or not run yet\n'
                          f'{knob}_mm,s11_dB')
    print(f'\n{knob} (mm)   S11 (dB)')
    for v, s in zip(values, s11):
        print(f'{v:7.3f}   {s:8.2f}')
    print(f'saved {out}')

    fig, ax = plt.subplots(num=f'S11 vs {knob}', tight_layout=True)
    ax.plot(values, s11, 'ko-', lw=2)
    ax.set_xlabel(f'{knob} (mm)')
    ax.set_ylabel(f'S11 at {RESONANCE_TARGET/1e9:.2f} GHz (dB)')
    ax.grid()
    plt.show()

