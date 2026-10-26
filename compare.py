"""
rev0 vs rev1 vs rev2: S11 and B1 per sqrt(W) against frequency, one figure.
    rev0  Sewani et al.'s board: 0.254 mm FR4, straight feed
    rev1  run 1's CMA-ES best on that board
    rev2  Sasaki et al.'s loop on 1.6 mm FR4 with a tapered feed (loopgap.py's geo)
Solves any revision not yet in runs/<rev> (5-10 min each), then plots into data/compare/.
Run with /home/yluo/opt/openEMS/venv/bin/python.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import src.simulate as simulate, src.postprocess as postprocess
from loopgap import rev0, rev1, rev2, Runs, sim_time

PLATE_T = 0.5   # mm: B1 at the top face of a diamond plate resting on the board, where the NVs are
REVS = [('rev0', rev0), ('rev1', rev1), ('rev2', rev2)]
COLORS = ['#2a78d6', '#eb6834', '#1baf7a']   # categorical slots 1-3, in fixed order
INK, INK2, MUTED, GRID, AXIS, SURFACE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
# Rabi frequency per uT of |B1|: gamma/sqrt(2) * sin(54.7 deg), (100) plate with B1 normal to it
RABI_PER_UT = 28.025e-3 / np.sqrt(2) * np.sqrt(2 / 3)

plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
                     'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2, 'text.color': INK,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'font.size': 10})

if __name__ == '__main__':
    fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(7.5, 7.5), tight_layout=True)
    print('rev    f_res (GHz)  S11 at f_res  B1 at 2.87 GHz  B1 min over 2.6-3.15 GHz  (uT/sqrt(W), plate top)')
    for (name, geo), c in zip(REVS, COLORS):
        path = os.path.join(Runs, name)   # absolute: simulate.run chdirs into it
        if not simulate.done(geo, path, sim_time):
            print(f'solving {name} ...', flush=True)
            simulate.run(geo, path, sim_time)
        f, B = postprocess.b1_spectrum(geo, path, geo.subs_h + geo.mask_t + PLATE_T)
        fs = np.linspace(f[0], f[-1], 501)
        s = postprocess.s11_dB(geo, path, fs)
        a1.plot(fs / 1e9, s, color=c, lw=2, label=name)
        a2.plot(f / 1e9, B, color=c, lw=2, marker='o', ms=5, mec=SURFACE, mew=1, label=name)
        i, j = np.argmin(s), np.argmax(B)
        # direct labels in ink, beside each curve's dip / peak; the legend carries the colour key
        ha, dx = {'rev1': ('left', 6), 'rev2': ('right', -6)}.get(name, ('center', 0))   # rev1/rev2 dips sit close
        a1.annotate(f'{name}  {fs[i]/1e9:.2f} GHz', (fs[i] / 1e9, s[i]), xytext=(dx, -14), textcoords='offset points',
                    ha=ha, color=INK2, fontsize=9)
        a2.annotate(name, (f[j] / 1e9, B[j]), xytext=(dx, 8), textcoords='offset points', ha=ha, color=INK2, fontsize=9)
        band = (f >= 2.6e9 - 1) & (f <= 3.15e9 + 1)
        print(f'{name}   {fs[i]/1e9:.3f}        {s[i]:6.1f} dB     {np.interp(2.87e9, f, B):6.0f}          {B[band].min():6.0f}')

    for ax in (a1, a2):
        ax.axvline(2.87, color=MUTED, ls=':', lw=1)
    a1.text(2.875, 0, 'D = 2.87 GHz', color=MUTED, fontsize=8, va='top')
    a1.set_title('S11 at the SMA', loc='left', color=INK, fontsize=11)
    a1.set_ylabel('S11 (dB)')
    a1.legend(frameon=False, loc='lower right', labelcolor=INK2)
    a2.set_title(f'|B1| at the NVs: hole axis, top of a {PLATE_T} mm plate on the board', loc='left', color=INK, fontsize=11)
    a2.set_ylabel('|B1| (µT per √W incident)')
    a1.set_ylim(a1.get_ylim()[0] - 0.8, 0.3)   # room for the dip labels
    a2.set_ylim(0, a2.get_ylim()[1] * 1.08)    # and the peak labels
    a2.set_xlabel('Frequency (GHz)')
    a2.text(0.0, -0.2, f'Rabi frequency (MHz) ≈ {RABI_PER_UT:.4f} × |B1| × √(power in W), (100) plate.  '
            'Model: openEMS, cu_res 0.25 mm, no diamond.', transform=a2.transAxes, color=MUTED, fontsize=8)
    a2.set_xmargin(0)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'compare')
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, 'rev0_rev1_rev2.png'), dpi=150)
    print('saved', os.path.join(out, 'rev0_rev1_rev2.png'))
    plt.show()
