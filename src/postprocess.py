# -*- coding: utf-8 -*-
"""
Read a finished run from sim_path: S11, Zin and B1, printed and plotted.  Needs no
solver, so it works on an old run too (python loopgap.py --post).
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from openEMS.physical_constants import MUE0


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

    # |B1| on the hole axis at f0, per sqrt(W) incident from a 50 ohm source. Absolute, not
    # relative: the FD dump and CalcPort share the same DFT scaling (x2, x dt) -- checked
    # against Ampere's law around a lumped port, 0.4% off. File layout: FieldData/FD/f0 is
    # (3, Nx, Ny, Nz) complex, Mesh/x,y,z in METRES at cell centres (dump_mode=2).
    # The dump holds f0 only (see geometry.py), not the resonance.
    with h5py.File(os.path.join(sim_path, 'B1.h5'), 'r') as h:
        bx, by, bz = (h['Mesh'][c][:] * 1e3 for c in 'xyz')
        H = h['FieldData/FD/f0'][:]
    k0 = np.argmin(np.abs(f - f0))
    B = MUE0 * np.sqrt((np.abs(H)**2).sum(0)) / np.sqrt(port.P_inc[k0]) * 1e6   # uT/sqrt(W), peak
    ix, iy, iz = np.argmin(np.abs(bx)), np.argmin(np.abs(by)), np.argmin(np.abs(bz - subs_h))
    print(f"B1 at top of hole (z={bz[iz]:.2f} mm): {B[ix, iy, iz]:.1f} uT/sqrt(W) incident, "
          f"f0 {f0/1e9:.2f} GHz, |S11| there {s11_dB[k0]:.1f} dB")

    fig, ax = plt.subplots(num="B1", tight_layout=True)
    ax.plot(bz, B[ix, iy], 'k.-', lw=2)
    ax.axvspan(0, subs_h, color='0.9', label='board thickness (hole is air)')
    ax.set_xlabel('z (mm)'); ax.set_ylabel('|B1| (µT/√W, peak)')
    ax.set_title(f'On the hole axis (x={bx[ix]:.3f}, y={by[iy]:.3f} mm) at {f0/1e9:.2f} GHz')
    ax.grid(); ax.legend(); ax.set_xmargin(0)

def s11_dB(geo, sim_path, f) -> float:
    """S11 in dB at exactly f (Hz) from a finished run: no plots, for an optimiser."""
    _, port = geo.build()
    port.CalcPort(sim_path, f)   # transform evaluated at f itself, no frequency grid
    return float(20 * np.log10(abs(port.uf_ref[0] / port.uf_inc[0])))
