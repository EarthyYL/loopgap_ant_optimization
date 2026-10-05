# -*- coding: utf-8 -*-
"""
Verify that a LoopGap's geometry came OUT right -- the LoopGap asserts only check that
the inputs were legal.  A wrong arc direction or a self-intersecting outline still
builds, still runs, and gives a plausible-looking wrong answer.  No solve.
"""

import numpy as np
from matplotlib.path import Path
from openEMS.physical_constants import C0, MUE0


def check(geo) -> None:
    R1, R2, yc, yi, g, f0 = geo.R1, geo.R2, geo.yc, geo.yi, geo.g, geo.f0
    csx, _ = geo.build()
    xs, ys = (np.asarray(a) for a in geo.outline())
    poly = Path(np.column_stack([xs, ys]))

    # outline area vs analytic catches a flipped arc or a self-intersecting outline
    area = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(np.roll(xs, -1), ys))
    want = np.pi * R1**2 - np.pi * R2**2 - g * (yc + R1 - yi)
    assert abs(area - want) / want < 0.01, f"outline area {area:.3f} vs {want:.3f} mm^2"

    for px, py, inside, what in [
            (0.0, (yi + yc + R1) / 2, False, 'slot interior'),
            (0.0, 0.0,                False, 'hole centre'),
            (R1 * 0.8, yc,            True,  'right arm'),
            (-R1 * 0.8, yc,           True,  'left arm'),
            (0.0, yc - R1 * 0.8,      True,  'copper below the hole')]:
        assert poly.contains_point((px, py)) == inside, f"{what} should be inside={inside}"

    assert min(geo.trace_x1, g / 2) - max(geo.trace_x0, -g / 2) <= 0, "trace bridges the slot"

    L = [np.asarray(csx.GetGrid().GetLines(d)) for d in range(3)]
    cells_g = int(((L[0] >= -g / 2) & (L[0] <= g / 2)).sum()) - 1
    lam4 = C0 / f0 / 4 * 1e3
    print(f"outline    area {area:.3f} mm^2, {abs(area-want)/want*100:.2f}% off analytic (faceting)")
    print(f"grid       {L[0].size} x {L[1].size} x {L[2].size} = "
          f"{L[0].size*L[1].size*L[2].size:,} cells")
    print(f"min cell   {min(np.diff(l).min() for l in L)*1000:.1f} um   "
          f"(floor {geo.mesh_floor*1000:.1f} um)")
    print(f"slot       g = {g} mm across {cells_g} cells" + ("  <-- want >= 3" if cells_g < 3 else ""))
    print(f"substrate  {geo.subs_cells} cells through {geo.subs_h} mm, "
          f"eps_r {geo.subs_eps} tand {geo.subs_tand}")
    sd = np.sqrt(2 / (2 * np.pi * f0 * MUE0 * geo.cu_sigma))
    print(f"copper     sheets at z=0 and z={geo.subs_h}, {geo.cu_t*1e6:.0f} um = "
          f"{geo.cu_t/sd:.0f} skin depths, Rs {1e3/(geo.cu_sigma*sd):.1f} mOhm/sq")
    print(f"mask       " + (f"{geo.mask_t} mm, eps_r {geo.mask_eps} tand {geo.mask_tand}"
                            if geo.mask_t > 0 else "none (keepout over the antenna)"))
    for d, ax in enumerate('xyz'):
        lo = dict(x=-geo.blk_w/2, y=-geo.blk_l/2, z=0.0)[ax]
        hi = dict(x=geo.blk_w/2,  y=geo.y_sma,    z=geo.subs_h)[ax]
        print(f"air gap {ax}  {lo-L[d][0]:5.1f} / {L[d][-1]-hi:5.1f} mm to boundary"
              + ("   <-- under lambda/4" if min(lo-L[d][0], L[d][-1]-hi) < lam4 else ""))
    print(f"           lambda/4 at {f0/1e9:.2f} GHz = {lam4:.1f} mm")
