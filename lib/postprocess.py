# -*- coding: utf-8 -*-
"""
Read a finished run from sim_path: S11, Zin, B1 and E.  Each reader rebuilds geo's port, which
costs no solve: the port only names the probe files the run wrote.
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from openEMS.physical_constants import MUE0, C0


def field(sim_path, name, f):
    """(x, y, z in mm, |F| as (Nx, Ny, Nz)) from the frequency-domain dump `name`, at its
    dumped frequency nearest f.  In the file, FieldData/FD/f<k> is (3, Nx, Ny, Nz) complex with
    a 'frequency' attribute, and Mesh/x,y,z is in metres."""
    with h5py.File(os.path.join(sim_path, f'{name}.h5'), 'r') as h:
        x, y, z = (h['Mesh'][c][:] * 1e3 for c in 'xyz')
        fd = h['FieldData/FD']
        F = fd[min(fd, key=lambda k: abs(fd[k].attrs['frequency'] - f))][:]
    return x, y, z, np.sqrt((np.abs(F)**2).sum(0))


def analyse(geo, sim_path) -> None:
    f0, fc, subs_h = geo.f0, geo.fc, geo.subs_h
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

    # absolute: the field dumps and CalcPort share one DFT scaling (checked by Ampere's law, 0.4 %)
    bx, by, bz, H = field(sim_path, 'B1', f0)
    k0 = np.argmin(np.abs(f - f0))
    B = MUE0 * H / np.sqrt(port.P_inc[k0]) * 1e6                  # uT/sqrt(W), peak
    E = field(sim_path, 'E', f0)[3] / np.sqrt(port.P_inc[k0])     # V/m/sqrt(W), peak
    ix, iy, iz = np.argmin(np.abs(bx)), np.argmin(np.abs(by)), np.argmin(np.abs(bz - subs_h))
    print(f"B1 at top of hole (z={bz[iz]:.2f} mm): {B[ix, iy, iz]:.1f} uT/sqrt(W) incident, "
          f"f0 {f0/1e9:.2f} GHz, |S11| there {s11_dB[k0]:.1f} dB")
    # E/cB1 is 1 in a plane wave; << 1 here if the slot holds the E field
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
    incident, interpolated in z between cell centres."""
    _, port = geo.build()
    port.CalcPort(sim_path, geo.f0)
    x, y, zz, F = field(sim_path, name, geo.f0)
    ix, iy = np.argmin(np.abs(x)), np.argmin(np.abs(y))
    return float(np.interp(z, zz, F[ix, iy]) / np.sqrt(port.P_inc[0]))

def b1_spectrum(geo, sim_path, z) -> tuple[np.ndarray, np.ndarray]:
    """(f in Hz, |B1| in uT/sqrt(W) incident, peak) on the hole axis at height z (mm), at every
    frequency the B1 dump holds.  Each frequency is normalised by its own incident power."""
    with h5py.File(os.path.join(sim_path, 'B1.h5'), 'r') as h:
        x, y, zz = (h['Mesh'][c][:] * 1e3 for c in 'xyz')
        fd = h['FieldData/FD']
        f = np.array([fd[k].attrs['frequency'] for k in fd])
        H = [np.sqrt((np.abs(fd[k][:])**2).sum(0)) for k in fd]
    ix, iy = np.argmin(np.abs(x)), np.argmin(np.abs(y))
    _, port = geo.build()
    port.CalcPort(sim_path, f)
    B = MUE0 * np.array([np.interp(z, zz, Hk[ix, iy]) for Hk in H]) / np.sqrt(port.P_inc) * 1e6
    o = np.argsort(f)
    return f[o], B[o]

def e_center(geo, sim_path) -> float:
    """|E| in V/m per sqrt(W) incident at the top of the hole -- the B1 readout point."""
    return on_axis(geo, sim_path, 'E', geo.subs_h)

def b1_at(geo, sim_path, z) -> float:
    """|B1| in uT per sqrt(W) incident, peak, on the hole axis at height z (mm)."""
    return MUE0 * on_axis(geo, sim_path, 'B1', z) * 1e6

def s11_dB(geo, sim_path, f) -> np.ndarray:
    """S11 in dB at exactly the frequencies f (Hz); a scalar f gives a length-1 array."""
    _, port = geo.build()
    port.CalcPort(sim_path, f)
    return 20 * np.log10(np.abs(port.uf_ref / port.uf_inc))

def zin(geo, sim_path, f) -> np.ndarray:
    """Complex input impedance (ohm) at the port at the frequencies f (Hz, array)."""
    _, port = geo.build()
    port.CalcPort(sim_path, f)
    return port.uf_tot / port.if_tot