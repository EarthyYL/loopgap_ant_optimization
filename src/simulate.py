# -*- coding: utf-8 -*-
"""
Run openEMS on a LoopGap.  Solver settings only: the geometry comes from geo.build(),
and the results are left as files in sim_path for postprocess.py to read.
"""

from openEMS import openEMS
from src.geometry import LoopGap


def run(geo: LoopGap, sim_path, sim_time) -> None:
    # NrTS is a STEP COUNT, so refining the mesh shrinks dt and silently shortens the
    # simulated time.  Express the run length in seconds instead: openEMS takes
    # min(NrTS, MaxTime/dt), so NrTS stays a safety ceiling and MaxTime is the real wall.
    # EndCriteria (1e-4 = -40 dB energy) should be what actually stops the run.
    # Energy decays as exp(-w*t/Q), so reaching -40 dB needs t ~ 1.47*Q/f0.
    FDTD = openEMS(NrTS=int(1e7), EndCriteria=1e-4, MaxTime=sim_time)
    FDTD.SetGaussExcite(geo.f0, geo.fc)   # from geo, so excitation and mesh share one band
    FDTD.SetBoundaryCond(['MUR'] * 6)
    csx, _ = geo.build()
    FDTD.SetCSX(csx)
    FDTD.Run(sim_path, cleanup=True)      # note: chdirs into sim_path
