# -*- coding: utf-8 -*-
"""
Planar loop-gap resonator for NV- spin manipulation: the Sewani et al. board
(Am. J. Phys. 88, 1156 (2020); Sewani_Files/NV_PCB_Antenna.brd).

This file names the values and dispatches.  The work is in:
    geometry.py      LoopGap: parameters -> CSX + port
    check.py         verify the built outline and mesh
    meshview.py      plot the mesh and the copper the solver sees
    simulate.py      run openEMS, results to Sim_Path
    postprocess.py   S11, Zin, B1 from the files in Sim_Path

Run with /home/yluo/opt/openEMS/venv/bin/python:
    python loopgap.py --check    # verify geometry + mesh, no solve
    python loopgap.py --mesh     # xy mesh at the top copper, no solve
    python loopgap.py --view     # 3D geometry only, no solve
    python loopgap.py            # solve + S11 + B1 readout
    python loopgap.py --post     # S11 + B1 readout of the last solve, no solve
"""

import os, sys, tempfile
import matplotlib.pyplot as plt
from src.geometry import LoopGap
import src.check as check, src.meshview as meshview, src.simulate as simulate, src.postprocess as postprocess

Sim_Path = os.path.join(tempfile.gettempdir(), 'LoopGap')

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
    subs_h=0.254,                      # 10 mil FR4 (Sewani et al. Fig. 2, .brd)
    subs_cells=4,
    cu_res=0.25,                       # 0.125 for final numbers: f0 +0.6%, 2.3x the runtime
    blk_w=26.0, blk_l=26.0,
    neck_w=20.0, neck_l=5.5,
    mnt_d=4.0, mnt_inset=1.0,          # plated 3 mm vias; 4 mm with their annular ring
    mnt_metal=False,
    trace_w=0.476,                     # from the .brd: 50 ohm on the 10 mil core
    tap_depth=1.0,                     # .brd
    sma_inset=1.8,                     # .brd
    cu_sigma=58e6, cu_t=36e-6,         # 36 um per the paper -- METRES
    # Solder mask over the whole top, as on the board (its top mask opens only the SMA
    # pad).  0 = keepout over the antenna: mask inside the 0.1 mm slot is an uncontrolled
    # fab variable sitting on exactly what sets f0.  Nonzero costs a mask_t-sized z-cell,
    # so the timestep drops with it.
    mask_t=0.03, mask_eps=4.3, mask_tand=0.025,
    air_xy=60.0, air_z=80.0,
    f0=2.87e9, fc=1.5e9,               # NV zero-field splitting is 2.87 GHz
)
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
