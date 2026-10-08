# -*- coding: utf-8 -*-
"""
Plot the xy mesh at the top copper.  No solve.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

SOLVER_MASK = 0x002 | 0x004   # MATERIAL | METAL: the property types openEMS queries per edge


def plot(geo) -> None:
    # geo is a LoopGap object
    R1, R2, yc, z = geo.R1, geo.R2, geo.yc, geo.subs_h
    csx, _ = geo.build()
    csx.Update()   # openEMS does this at setup; without it polygons and cylinders are never "inside"
    X, Y = (np.asarray(csx.GetGrid().GetLines(d)) for d in range(2))

    def copper(x, y):
        p = csx.GetPropertyByCoordPriority([x, y, z], SOLVER_MASK)
        return p is not None and p.GetName() == 'loop'

    segs = [[(a, y), (b, y)] for a, b in zip(X[:-1], X[1:]) for y in Y if copper((a + b) / 2, y)] + \
           [[(x, a), (x, b)] for a, b in zip(Y[:-1], Y[1:]) for x in X if copper(x, (a + b) / 2)]
    xs, ys = geo.outline()

    fig, A = plt.subplots(2, 2, figsize=(12, 12), num='Mesh', tight_layout=True)
    fig.suptitle(f"grid {len(X)} x {len(Y)} x {csx.GetGrid().GetQtyLines(2)}, "
                 f"cu_res {geo.cu_res} mm over the C, {R2 / 5:g} mm over the hole")
    for row, (lx, ly) in enumerate([((-R1 - 1, R1 + 1), (yc - R1 - 1, yc + R1 + 1)),   # whole C
                                    ((-2 * R2, 2 * R2), (-2 * R2, 2.4 * R2))]):         # hole + slot root
        a = A[row, 0]
        a.fill(xs, ys, color='#e8b070')
        a.vlines(X, *ly, lw=0.4); a.hlines(Y, *lx, lw=0.4)
        a.set_title('outline + grid lines')
        a = A[row, 1]
        a.plot(xs, ys, 'k', lw=0.6)
        a.add_collection(LineCollection(segs, color='#c05000', lw=0.8 + row))
        a.set_title('copper edges the solver uses')
        for a in A[row]:
            a.set_xlim(lx); a.set_ylim(ly); a.set_aspect('equal')
            a.set_xlabel('x (mm)'); a.set_ylabel('y (mm)')
