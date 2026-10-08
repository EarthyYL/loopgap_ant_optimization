import os
from dataclasses import asdict
import numpy as np
import matplotlib.pyplot as plt
from src.geometry import LoopGap
from functools import partial 
import src.simulate as simulate, src.postprocess as postprocess
from loopgap import geo as base, Runs, sim_time

# The design knobs, in x-vector order.  yc (disc centre above the hole) replaces l1/l2:
# l1 + l2 == D1, so they are one knob, not two.  mask_t is a knob too, but 0-or-0.03, so it
# stays in loopgap.py.  Everything not a knob is taken from loopgap.py's geo.
KNOBS = ('D1', 'D2', 'yc', 'g', 'tap_depth')
FIXED = {k: v for k, v in asdict(base).items() if k not in ('D1', 'D2', 'l1', 'l2', 'g', 'tap_depth')}
RESONANCE_TARGET = 2.87e9

def design(D1, D2, yc, g, tap_depth) -> LoopGap:
    """LoopGap from the design knobs alone; call as design(**dict(zip(KNOBS, x)))."""
    return LoopGap(**FIXED, D1=D1, D2=D2, l1=D1 / 2 + yc, l2=D1 / 2 - yc, g=g, tap_depth=tap_depth)


def run_sim(geo: LoopGap, sim_time: float, sim_path: str, res_f: float) -> tuple[float, float]:
    """Solve, unless sim_path already holds a finished run of this exact design, and return
    (S11 in dB at res_f, |E| at the hole centre in V/m per sqrt(W) at geo.f0), for optimization."""
    if not simulate.done(geo, sim_path, sim_time):
        simulate.run(geo, sim_path, sim_time)
    return postprocess.s11_dB(geo, sim_path, res_f)[0], postprocess.e_center(geo, sim_path)



if __name__ == '__main__':
    # The swept knob, one of KNOBS; the others stay at loopgap.py's current design.
    knob, values = 'g', np.geomspace(0.05, 0.8, 6)
    dsgn_parameterized = partial(design, **{k: getattr(base, k) for k in KNOBS})
    res = np.full((len(values), 2), np.nan)   # S11 dB, |E|; nan = design rejected, never solved
    # absolute, and fixed before the loop: after the first solve the cwd is that run's folder
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', f'sweep_{knob}.csv')
    for i, v in enumerate(values):
        try:
            geo = dsgn_parameterized(**{knob: v})
        except AssertionError as e:   # illegal geometry: skip it instead of dying mid-sweep
            print(f'{knob} {v:.3f} mm   skipped: {e}', flush=True)
            continue
        # ponytail: named by knob value, so a sweep after changing loopgap.py's design
        # re-solves into (overwrites) these; name by a hash of params.json to keep every run
        path = os.path.join(Runs, f'{knob}_{v:.3f}')   # absolute: FDTD.Run chdirs into it
        res[i] = run_sim(geo=geo, sim_time=sim_time, sim_path=path, res_f=RESONANCE_TARGET)
        print(f'{knob} {v:.3f} mm   S11 at {RESONANCE_TARGET/1e9:.2f} GHz {res[i, 0]:.2f} dB   '
              f'|E| {res[i, 1]:.0f} V/m/sqrt(W)   ({path})', flush=True)
        # rewritten after every point, so a crash mid-sweep keeps what finished
        np.savetxt(out, np.column_stack([values, res]), delimiter=',', fmt='%.4f',
                   header=f'S11 and |E| at the hole centre, at {RESONANCE_TARGET/1e9:.2f} GHz, vs {knob}; '
                          f'nan = rejected or not run yet\n{knob}_mm,s11_dB,E_V_per_m_sqrtW')
    print(f'\n{knob} (mm)   S11 (dB)   |E| (V/m/sqrt(W))')
    for v, (s, e) in zip(values, res):
        print(f'{v:7.3f}   {s:8.2f}   {e:10.0f}')
    print(f'saved {out}')

    fig, (a1, a2) = plt.subplots(2, 1, sharex=True, num=f'S11 and E vs {knob}', tight_layout=True)
    a1.plot(values, res[:, 0], 'ko-', lw=2)
    a1.set_ylabel(f'S11 at {RESONANCE_TARGET/1e9:.2f} GHz (dB)')
    a2.plot(values, res[:, 1], 'ko-', lw=2)
    a2.set_ylabel('|E| at hole centre (V/m/√W)')
    a2.set_xlabel(f'{knob} (mm)')
    a1.grid(); a2.grid()
    plt.show()

