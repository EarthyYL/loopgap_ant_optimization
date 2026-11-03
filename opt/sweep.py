import os
import numpy as np
import matplotlib.pyplot as plt
from lib.geometry import LoopGap
from functools import partial
import lib.simulate as simulate, lib.postprocess as postprocess
from designs import rev2 as base
from settings import Runs, Data, sim_time
from lib.knobs import KNOBS, design

RESONANCE_TARGET = 2.87e9


def run_sim(geo: LoopGap, sim_time: float, sim_path: str, res_f: float) -> tuple[float, float]:
    """Solve unless already done; return (S11 in dB at res_f, |E| in V/m/sqrt(W) at the hole centre)."""
    if not simulate.done(geo, sim_path, sim_time):
        simulate.run(geo, sim_path, sim_time)
    return postprocess.s11_dB(geo, sim_path, res_f)[0], postprocess.e_center(geo, sim_path)



if __name__ == '__main__':
    knob, values = 'g', np.geomspace(0.05, 0.8, 6)
    dsgn_parameterized = partial(design, **{k: getattr(base, k) for k in KNOBS})
    res = np.full((len(values), 2), np.nan)   # S11 dB, |E|; nan = rejected, never solved
    out = os.path.join(Data, f'sweep_{knob}.csv')
    for i, v in enumerate(values):
        try:
            geo = dsgn_parameterized(**{knob: v})
        except AssertionError as e:
            print(f'{knob} {v:.3f} mm   skipped: {e}', flush=True)
            continue
        # ponytail: named by knob value, so a sweep after changing designs.py's rev2
        # re-solves into (overwrites) these; name by a hash of params.json to keep every run
        path = os.path.join(Runs, f'{knob}_{v:.3f}')
        res[i] = run_sim(geo=geo, sim_time=sim_time, sim_path=path, res_f=RESONANCE_TARGET)
        print(f'{knob} {v:.3f} mm   S11 at {RESONANCE_TARGET/1e9:.2f} GHz {res[i, 0]:.2f} dB   '
              f'|E| {res[i, 1]:.0f} V/m/sqrt(W)   ({path})', flush=True)
        # written after every point, so a crash keeps what finished
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

