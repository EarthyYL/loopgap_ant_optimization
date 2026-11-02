"""
Rev 0 vs Rev 1 vs Rev 2: S11, B1 per sqrt(W) and Zin against frequency, one figure each.
    rev0  Sewani et al.'s board: 0.254 mm FR4, straight feed
    rev1  run 1's CMA-ES best on that board
    rev2  Sasaki et al.'s loop on 1.6 mm FR4 with a tapered feed, g 0.25 (loopgap.py's geo)
Solves any revision not yet in runs/<rev> (5-15 min each), then writes data/compare/s11,
data/compare/b1 and data/compare/zin as .pdf and .png; --pgf adds a .pgf of each for \\input in LaTeX (needs pdflatex).
B1 is on the hole axis at the top face of a 0.5 mm diamond plate resting on the board.  For a (100)
plate the Rabi frequency is ~0.0162 MHz per uT of B1, times sqrt(power in W).
Run with /home/yluo/opt/openEMS/venv/bin/python.
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import src.simulate as simulate, src.postprocess as postprocess
from loopgap import rev0, rev1, rev2, Runs, sim_time

PLATE_T = 0.5     # mm above the board: where the NVs are
WIDTH = 3.4       # inches: one column; ~7 for full width
HEIGHT = 0.7 * WIDTH
# colour plus line style, so the curves survive greyscale printing
STYLE = {'rev0': dict(color='#2a78d6', ls='--'), 'rev1': dict(color='#eb6834', ls='-.'),
         'rev2': dict(color='#1baf7a', ls='-')}
# where each curve's name goes, in data coordinates (x GHz, y, alignment): S11, B1, then Re Zin.
# Placed by hand beside each dip and peak -- move them if the curves change.
LABEL = {'rev0': ((2.48, -9.6, 'left'),  (2.40, 255, 'left'),  (2.42, 150, 'left')),
         'rev1': ((2.80, -6.9, 'right'), (3.02, 165, 'right'), (2.64, 290, 'right')),
         'rev2': ((3.05, -3.0, 'left'),  (2.97, 268, 'left'),  (2.95, 140, 'left'))}

# LaTeX look without LaTeX: matplotlib's bundled Computer Modern.  With --pgf the document's own
# fonts take over (pgf.rcfonts False), so every label is valid LaTeX.
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
                     'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False, 'font.size': 9,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'lines.linewidth': 1.2,
                     'pgf.texsystem': 'pdflatex', 'pgf.rcfonts': False})

if __name__ == '__main__':
    pgf = '--pgf' in sys.argv
    (f1, a1), (f2, a2) = (plt.subplots(figsize=(WIDTH, HEIGHT), layout='constrained') for _ in range(2))
    # Zin: Re over Im, as tall as the S11 and B1 figures stacked, so the three tile a 2x2 grid
    f3, (a3, a4) = plt.subplots(2, 1, sharex=True, figsize=(WIDTH, 2 * HEIGHT), layout='constrained')
    print('rev    f_res (GHz)  S11 at f_res  B1 at 2.87 GHz  B1 min over 2.6-3.15 GHz  (uT/sqrt(W), plate top)')
    for name, geo in [('rev0', rev0), ('rev1', rev1), ('rev2', rev2)]:
        path = os.path.join(Runs, name)   # absolute: simulate.run chdirs into it
        if not simulate.done(geo, path, sim_time):
            print(f'solving {name} ...', flush=True)
            simulate.run(geo, path, sim_time)
        f, B = postprocess.b1_spectrum(geo, path, geo.subs_h + geo.mask_t + PLATE_T)
        fs = np.linspace(f[0], f[-1], 501)
        s = postprocess.s11_dB(geo, path, fs)
        Z = postprocess.zin(geo, path, fs)
        a1.plot(fs / 1e9, s, **STYLE[name])
        a2.plot(f / 1e9, B, **STYLE[name])
        a3.plot(fs / 1e9, Z.real, **STYLE[name])
        a4.plot(fs / 1e9, Z.imag, **STYLE[name])
        for ax, (x, y, ha) in zip((a1, a2, a3), LABEL[name]):
            ax.text(x, y, name.replace('rev', 'Rev '), color=STYLE[name]['color'], ha=ha, va='center', fontsize=8)
        i, band = np.argmin(s), (f >= 2.6e9 - 1) & (f <= 3.15e9 + 1)
        print(f'{name}   {fs[i]/1e9:.3f}        {s[i]:6.1f} dB     {np.interp(2.87e9, f, B):6.0f}          {B[band].min():6.0f}')

    for ax in (a1, a2, a3, a4):
        ax.axvline(2.87, color='0.6', ls=':', lw=0.8, zorder=0)   # zero-field splitting D
        ax.set_xmargin(0)
    for ax in (a1, a2, a4):
        ax.set_xlabel(r'$f$ (GHz)')
    a3.axhline(50, color='0.6', ls=':', lw=0.8, zorder=0)   # the 50 ohm source
    a4.axhline(0, color='0.6', ls=':', lw=0.8, zorder=0)
    a3.set_ylabel(r'$\mathrm{Re}\,Z_\mathrm{in}$ ($\Omega$)')
    a4.set_ylabel(r'$\mathrm{Im}\,Z_\mathrm{in}$ ($\Omega$)')
    a3.set_ylim(bottom=0)
    a1.set_ylabel(r'$S_{11}$ (dB)')
    a2.set_ylabel(r'$B_1$ ($\mu$T/$\sqrt{\mathrm{W}}$)')
    a2.set_ylim(bottom=0)

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'compare')
    os.makedirs(outdir, exist_ok=True)
    for fig, stem in [(f1, 's11'), (f2, 'b1'), (f3, 'zin')]:
        out = os.path.join(outdir, stem)
        for ext in ['pdf', 'png'] + (['pgf'] if pgf else []):
            fig.savefig(f'{out}.{ext}', dpi=300)
            print('saved', f'{out}.{ext}')
        if pgf:   # matplotlib's .pgf expects \mathdefault from the preamble; define it in the file
            body = open(f'{out}.pgf').read()
            open(f'{out}.pgf', 'w').write('\\providecommand\\mathdefault[1]{#1}%\n' + body)
    plt.show()
