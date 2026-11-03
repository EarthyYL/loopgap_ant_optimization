# -*- coding: utf-8 -*-
"""Run 1's search knobs, and design(): a LoopGap from the knobs, everything else from rev2."""

from dataclasses import asdict
from lib.geometry import LoopGap
from designs import rev2 as base

KNOBS = ('D1', 'D2', 'yc', 'g', 'tap_depth')   # yc stands for l1 and l2, since l1 + l2 == D1
FIXED = {k: v for k, v in asdict(base).items() if k not in ('D1', 'D2', 'l1', 'l2', 'g', 'tap_depth')}

def design(D1, D2, yc, g, tap_depth) -> LoopGap:
    """Call as design(**dict(zip(KNOBS, x)))."""
    return LoopGap(**FIXED, D1=D1, D2=D2, l1=D1 / 2 + yc, l2=D1 / 2 - yc, g=g, tap_depth=tap_depth)
