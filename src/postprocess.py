# -*- coding: utf-8 -*-
"""
Read a finished run from sim_path: S11, Zin, B1 and E, printed and plotted.  Can work
on an old run too (python loopgap.py --post).
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from openEMS.physical_constants import MUE0, C0


def field(sim_path, name):
    """(x, y, z in mm, |F| as (Nx, Ny, Nz)) at f0 from the frequency-domain dump `name`.
    File layout: FieldData/FD/f0 is (3, Nx, Ny, Nz) complex, Mesh/x,y,z in METRES at cell
    centres (dump_mode=2).  The dump holds f0 only (see geometry.py), not the resonance."""
    with h5py.File(os.path.join(sim_path, f'{name}.h5'), 'r') as h:
        x, y, z = (h['Mesh'][c][:] * 1e3 for c in 'xyz')
        F = h['FieldData/FD/f0'][:]
    return x, y, z, np.sqrt((np.abs(F)**2).sum(0))


def analyse(geo, sim_path) -> None:
    f0, fc, subs_h = geo.f0, geo.fc, geo.subs_h
    # The port only carries its probe file names and 50 ohm reference, so a rebuilt one
    # reads exactly the files the run wrote.  Rebuilding costs milliseconds, no solve.
    _, port = geo.build()

    f = np.linspace(f0 - fc, f0 + fc, 401)
    port.CalcPort(sim_path, f)
    s11_dB = 20 * np.log10(np.abs(port.uf_ref / port.uf_inc))
    Zin = port.uf_tot / port.if_tot

    i = np.argmin(s11_dB)
    print(f"resonance {f[i]/1e9:.4f} GHz   S11 {s11_dB[i]:.2f} dB   "
          f"Zin {Zin[i].real:.1f}{Zin[i].imag:+.1f}j ohm")

    fig, ax = plt.subplots(num="S11", tight_layout=True)
    ax.plot(f / 1e9, s11_dB, 'k-', lw=2)
    ax.axvline(2.87, color='r', ls=':', label='NV 2.87 GHz')
    ax.set_xlabel('Frequency (GHz)'); ax.set_ylabel('S11 (dB)')
    ax.grid(); ax.legend(); ax.set_xmargin(0)

    fig, ax = plt.subplots(num="Zin", tight_layout=True)
    ax.plot(f / 1e9, Zin.real, 'k-',  lw=2, label=r'$\Re\{Z_{in}\}$')
    ax.plot(f / 1e9, Zin.imag, 'r--', lw=2, label=r'$\Im\{Z_{in}\}$')
    ax.set_xlabel('Frequency (GHz)'); ax.set_ylabel('Zin (Ohm)'); ax.set_title('Input Impedance')
    ax.grid(); ax.legend(); ax.set_xmargin(0)

    # |B1| and |E| on the hole axis at f0, per sqrt(W) incident from a 50 ohm source.
    # Absolute, not relative: the FD dump and CalcPort share the same DFT scaling (x2, x dt)
    # -- checked against Ampere's law around a lumped port, 0.4% off.
    bx, by, bz, H = field(sim_path, 'B1')
    k0 = np.argmin(np.abs(f - f0))
    B = MUE0 * H / np.sqrt(port.P_inc[k0]) * 1e6                  # uT/sqrt(W), peak
    E = field(sim_path, 'E')[3] / np.sqrt(port.P_inc[k0])         # V/m/sqrt(W), peak; same box, same cells
    ix, iy, iz = np.argmin(np.abs(bx)), np.argmin(np.abs(by)), np.argmin(np.abs(bz - subs_h))
    print(f"B1 at top of hole (z={bz[iz]:.2f} mm): {B[ix, iy, iz]:.1f} uT/sqrt(W) incident, "
          f"f0 {f0/1e9:.2f} GHz, |S11| there {s11_dB[k0]:.1f} dB")
    # E/cB is 1 in a plane wave; the loop-gap keeps E in the slot, so here it should be << 1.
    # Unlike |E| per incident watt, it does not drop just because the match got worse.
    print(f"|E| there: {E[ix, iy, iz]:.0f} V/m/sqrt(W) incident, "
          f"E/cB1 {E[ix, iy, iz] / (C0 * B[ix, iy, iz] * 1e-6):.3f}")

    for name, F, unit in [('B1', B, 'µT'), ('E', E, 'V/m')]:
        fig, ax = plt.subplots(num=name, tight_layout=True)
        ax.plot(bz, F[ix, iy], 'k.-', lw=2)
        ax.axvspan(0, subs_h, color='0.9', label='board thickness (hole is air)')
        ax.set_xlabel('z (mm)'); ax.set_ylabel(f'|{name}| ({unit}/√W, peak)')
        ax.set_title(f'On the hole axis (x={bx[ix]:.3f}, y={by[iy]:.3f} mm) at {f0/1e9:.2f} GHz')
        ax.grid(); ax.legend(); ax.set_xmargin(0)

def on_axis(geo, sim_path, name, z) -> float:
    """|F| of the dump `name` at geo.f0 on the hole axis at height z (mm), per sqrt(W)
    incident, from a finished run: no plots, for an optimiser.  Linear in z between cell
    centres: above the board they are 0.1-0.2 mm apart, so the nearest one can be far off.
    The axis is air all the way up (the hole is a through-hole), so no interface to straddle."""
    _, port = geo.build()
    port.CalcPort(sim_path, geo.f0)   # the dumps hold f0 only, so normalise there
    x, y, zz, F = field(sim_path, name)
    ix, iy = np.argmin(np.abs(x)), np.argmin(np.abs(y))
    return float(np.interp(z, zz, F[ix, iy]) / np.sqrt(port.P_inc[0]))

def e_center(geo, sim_path) -> float:
    """|E| in V/m per sqrt(W) incident at the top of the hole -- the B1 readout point."""
    return on_axis(geo, sim_path, 'E', geo.subs_h)

def b1_at(geo, sim_path, z) -> float:
    """|B1| in uT per sqrt(W) incident, peak, on the hole axis at height z (mm)."""
    return MUE0 * on_axis(geo, sim_path, 'B1', z) * 1e6

def s11_dB(geo, sim_path, f) -> np.ndarray:
    """S11 in dB at exactly the frequencies f (Hz, array) from a finished run: no plots, for
    an optimiser.  A scalar f still gives a length-1 array."""
    _, port = geo.build()
    port.CalcPort(sim_path, f)   # transform evaluated at f itself, no frequency grid
    return 20 * np.log10(np.abs(port.uf_ref / port.uf_inc))