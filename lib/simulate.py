# -*- coding: utf-8 -*-
"""Run openEMS on a LoopGap, and tell whether a run is already done."""

import os, json
from dataclasses import asdict
from openEMS import openEMS
from lib.geometry import LoopGap


def params(geo: LoopGap, sim_time) -> dict:
    """Everything that decides a run's output, JSON-safe."""
    return {**asdict(geo), 'sim_time': sim_time}


def run(geo: LoopGap, sim_path, sim_time, threads=0) -> None:
    """Solve geo into sim_path on `threads` engine threads (0 = all cores)."""
    # NrTS counts steps, so it is only a ceiling: MaxTime (seconds) and EndCriteria (-40 dB) end the run
    FDTD = openEMS(NrTS=int(1e7), EndCriteria=1e-4, MaxTime=sim_time)
    FDTD.SetGaussExcite(geo.f0, geo.fc)
    FDTD.SetBoundaryCond(['MUR'] * 6)
    csx, _ = geo.build()
    FDTD.SetCSX(csx)
    record = os.path.join(os.path.abspath(sim_path), 'params.json')   # absolute: Run chdirs
    os.makedirs(sim_path, exist_ok=True)
    if os.path.exists(record):
        os.remove(record)   # or a crashed re-run would look done()
    FDTD.Run(sim_path, cleanup=True, numThreads=threads)
    with open(record, 'w') as f:   # last, so it exists only for a finished run
        json.dump(params(geo, sim_time), f, indent=1)


def done(geo: LoopGap, sim_path, sim_time) -> bool:
    """True if sim_path holds a finished run of exactly this geo and sim_time.  Code isn't in
    params.json: clear runs/ after changing the mesh or the dumps."""
    record = os.path.join(sim_path, 'params.json')
    if not os.path.exists(record):
        return False
    with open(record) as f:
        return json.load(f) == params(geo, sim_time)
