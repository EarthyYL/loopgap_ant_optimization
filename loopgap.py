# -*- coding: utf-8 -*-
"""
Planar ring (loop-gap) antenna for NV- ODMR: Sasaki et al.'s design on standard 1.6 mm FR4
(Rev. Sci. Instrum. 87, 053904 (2016); arXiv:1605.04627), measured at 2.79 GHz with ~400 MHz
bandwidth.  Same loop as the Sewani et al. board (Am. J. Phys. 88, 1156 (2020);
Sewani_Files/NV_PCB_Antenna.brd), which moved it to 0.254 mm FR4.

This file names the values and dispatches.  The work is in:
    geometry.py      LoopGap: parameters -> CSX + port
    check.py         verify the built outline and mesh
    meshview.py      plot the mesh and the copper the solver sees
    simulate.py      run openEMS, results to Sim_Path
    postprocess.py   S11, Zin, B1, E from the files in Sim_Path

Run with /home/yluo/opt/openEMS/venv/bin/python:
    python loopgap.py --check    # verify geometry + mesh, no solve
    python loopgap.py --mesh     # xy mesh at the top copper, no solve
    python loopgap.py --view     # 3D geometry only, no solve
    python loopgap.py            # solve + S11 + B1 readout
    python loopgap.py --post     # S11 + B1 readout of the last solve, no solve
"""

import os, sys
from dataclasses import replace
import matplotlib.pyplot as plt
from src.geometry import LoopGap
import src.check as check, src.meshview as meshview, src.simulate as simulate, src.postprocess as postprocess

# Kept, not temp: one folder per run, each with the params.json that made it.
# optimize.py puts its sweep points beside this design's own folder.
Runs = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs')
Sim_Path = os.path.join(Runs, 'rev2')   # runs/ is shared with main, whose base design is runs/loopgap

# UNITS: every length below is mm -- the grid unit is 1 mm.
# The one exception is cu_t, which openEMS wants in METRES: the conducting-sheet
# operator computes sigma*t as a sheet conductance in S/square, so it is an SI
# quantity independent of the grid unit.  cu_t must never appear in a coordinate.
# Everything else non-length is SI: S/m, Hz, dimensionless.
geo = LoopGap(
    D1=14.0, D2=1.0,
    l1=10.9, l2=3.1,
    g=0.1,
    subs_eps=4.3, subs_tand=0.03,      # FR4
    subs_h=1.6,                        # standard FR4 (Sasaki et al.)
    subs_cells=4,
    cu_res=0.25,                       # 0.125 for final numbers: f0 +0.6%, 2.3x the runtime
    blk_w=26.0, blk_l=26.0,
    neck_w=20.0, neck_l=22.7,          # fits Sasaki's 23 mm feed: 13 + 22.7 - 1.8 (sma_inset) - 10.9 (disc top) = 23
    mnt_d=4.0, mnt_inset=1.0,          # plated 3 mm vias; 4 mm with their annular ring
    mnt_metal=False,
    trace_w=0.54, trace_w0=3.3,        # Sasaki: tapered 3.3 mm at the SMA -> 0.54 mm at the ring, for 50 ohm
    tap_depth=1.0,                     # .brd
    sma_inset=1.8,                     # .brd
    cu_sigma=58e6, cu_t=18e-6,         # 18 um per Sasaki (JLCPCB's thinnest is 35 um; loss only) -- METRES
    # Solder mask over the whole top, as on the board (its top mask opens only the SMA
    # pad).  0 = keepout over the antenna: mask inside the 0.1 mm slot is an uncontrolled
    # fab variable sitting on exactly what sets f0.  Nonzero costs a mask_t-sized z-cell,
    # so the timestep drops with it.
    mask_t=0.03, mask_eps=4.3, mask_tand=0.025,
    air_xy=60.0, air_z=80.0,
    f0=2.87e9, fc=1.5e9,               # NV zero-field splitting is 2.87 GHz
)
# The revisions, for compare.py.  geo above is rev2.
# rev0: Sewani et al.'s board as built -- 0.254 mm FR4, straight 0.476 mm (50 ohm) feed, 36 um Cu.
rev0 = replace(geo, subs_h=0.254, neck_l=5.5, trace_w=0.476, trace_w0=0.476, cu_t=36e-6)
# rev1: run 1's CMA-ES best on rev0's board (data/loopgap_ant_optimization, runs/opt/e57b2d7dd263)
rev1 = replace(rev0, D1=11.959, D2=0.831, l1=11.959 / 2 + 3.390, l2=11.959 / 2 - 3.390, g=0.105, tap_depth=0.342)
rev2 = geo

sim_time = 30e-9   # ~ Q of 60 at 2.87 GHz; raise if EndCriteria is not reached

if __name__ == "__main__":
    if '--check' in sys.argv:
        check.check(geo)
    elif '--mesh' in sys.argv:
        meshview.plot(geo)
        plt.show()
    elif '--view' in sys.argv:
        os.makedirs(Sim_Path, exist_ok=True)
        xml = os.path.join(Sim_Path, 'loopgap.xml')
        csx, _ = geo.build()
        csx.Write2XML(xml)
        nx, ny, nz = (csx.GetGrid().GetQtyLines(d) for d in range(3))
        print(f"grid {nx} x {ny} x {nz} = {nx*ny*nz:,} cells")
        os.environ['PATH'] += os.pathsep + '/home/yluo/opt/openEMS/bin'
        from CSXCAD import AppCSXCAD_BIN
        raise SystemExit(os.system(AppCSXCAD_BIN + ' "{}"'.format(xml)))
    else:
        if '--post' not in sys.argv:
            simulate.run(geo, Sim_Path, sim_time)
        postprocess.analyse(geo, Sim_Path)
        plt.show()
