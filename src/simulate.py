# -*- coding: utf-8 -*-
"""
Run openEMS on a LoopGap.  Solver settings only: the geometry comes from geo.build(),
and the results are left as files in sim_path for postprocess.py to read, with
params.json saying what made them.
"""

import os, json
from dataclasses import asdict
from openEMS import openEMS
from src.geometry import LoopGap


def params(geo: LoopGap, sim_time) -> dict:
    """Everything that decides a run's output; JSON-safe (floats, ints, bools)."""
    return {**asdict(geo), 'sim_time': sim_time}


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
    # absolute: Run chdirs into sim_path.  Its cleanup only deletes openEMS's own files, so
    # drop the old params.json too, or a crashed re-run would look finished to done().
    record = os.path.join(os.path.abspath(sim_path), 'params.json')
    os.makedirs(sim_path, exist_ok=True)   # Run only mkdirs the last level
    if os.path.exists(record):
        os.remove(record)
    FDTD.Run(sim_path, cleanup=True)      # note: chdirs into sim_path
    with open(record, 'w') as f:          # last, so it only exists for a finished run
        json.dump(params(geo, sim_time), f, indent=1)


def done(geo: LoopGap, sim_path, sim_time) -> bool:
    """True if sim_path already holds a finished run of exactly this geo and sim_time.
    Code changes (mesh, dumps) are not in params.json: clear runs/ after one of those."""
    record = os.path.join(sim_path, 'params.json')
    if not os.path.exists(record):
        return False
    with open(record) as f:
        return json.load(f) == params(geo, sim_time)
