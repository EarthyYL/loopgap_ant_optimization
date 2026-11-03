# -*- coding: utf-8 -*-
"""
Check, view or solve Rev 2.  Run from the repo root with /home/yluo/opt/openEMS/venv/bin/python:
    python -m sim.loopgap --check    # verify geometry + mesh, no solve
    python -m sim.loopgap --mesh     # xy mesh at the top copper, no solve
    python -m sim.loopgap --view     # 3D geometry only, no solve
    python -m sim.loopgap            # solve + S11 + B1 readout
    python -m sim.loopgap --post     # S11 + B1 readout of the last solve, no solve
"""

import os, sys
import matplotlib.pyplot as plt
from designs import rev2 as geo
from settings import Runs, sim_time
import lib.check as check, lib.meshview as meshview, lib.simulate as simulate, lib.postprocess as postprocess

Sim_Path = os.path.join(Runs, 'rev2')

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
        os.environ['PATH'] += os.pathsep + '/home/yluo/opt/openEMS/bin'   # AppCSXCAD on this machine; elsewhere, put it on PATH
        from CSXCAD import AppCSXCAD_BIN
        raise SystemExit(os.system(AppCSXCAD_BIN + ' "{}"'.format(xml)))
    else:
        if '--post' not in sys.argv:
            simulate.run(geo, Sim_Path, sim_time)
        postprocess.analyse(geo, Sim_Path)
        plt.show()
