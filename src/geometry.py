# -*- coding: utf-8 -*-
"""
Planar loop-gap resonator geometry: parameters in, CSXCAD structure out.

Geometry (top copper), all dims mm, hole D2 centred on the origin:
  - outer copper disc  diameter D1, centre (0, (l1-l2)/2)   [non-concentric with D2]
  - central hole       diameter D2, centre (0, 0)
  - radial slot        width g, from top of D2 out to top of D1  -> copper is a "C"
  - feed trace         microstrip from the SMA pin down onto the left arm of the loop
  - bottom layer       solid ground

No solver and no values here: the values are named in loopgap.py, and simulate.py
hands the CSX from build() to openEMS.
"""

from dataclasses import dataclass
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS.ports import LumpedPort
from openEMS.automesh import mesh_hint_from_primitive
from openEMS.physical_constants import C0, EPS0


def graded(L, max_res, ratio, regions=()):
    """Mesh lines through every line in L, cells up to ~max_res, each cell within ~ratio
    of its neighbours (up to ~1.7x where pinned lines force a whole number of cells).
    regions: (lo, hi, h) caps the cells at h inside [lo, hi], tapering out at the same ratio.
    Replaces CSXCAD's SmoothMeshLines, which only splits gaps WIDER than max_res: it left
    a 25 um slot cell against a 450 um one, and 35x at the trace edge."""
    L = np.unique(L)
    d = np.diff(L)
    # wanted cell size at x: each fixed cell's own size, growing ln(ratio) per unit distance
    # away from it (the continuous form of a geometric taper), capped at max_res
    x = np.unique(np.r_[L, np.linspace(L[0], L[-1], 20001)])
    away = np.maximum(L[:-1, None] - x, x - L[1:, None]).clip(0)
    h = np.minimum(max_res, (d[:, None] + np.log(ratio) * away).min(0))
    for lo, hi, hr in regions:
        h = np.minimum(h, hr + np.log(ratio) * np.maximum(lo - x, x - hi).clip(0))
    n = np.r_[0, np.cumsum(np.diff(x) * (1 / h[1:] + 1 / h[:-1]) / 2)]   # cells needed so far
    out = [L[0]]
    for a, b in zip(L[:-1], L[1:]):
        na, nb = np.interp([a, b], x, n)
        N = max(1, int(np.ceil(nb - na - 0.25)))   # -0.25: don't split a cell for a 1.1x excess
        out += list(np.interp(na + (nb - na) * np.arange(1, N) / N, n, x)) + [b]
    return np.array(out)


@dataclass(frozen=True)  # frozen: the ass below only run on construction, so
class LoopGap:            # changes must go through dataclasses.replace(), which reruns them
    # UNITS: every length is mm (build() sets the grid unit to 1 mm), except cu_t,
    # which is in METRES.  Everything else is SI.  See loopgap.py.
    D1: float           # outer copper diameter
    D2: float           # central hole diameter (sample / optical access)
    l1: float           # hole centre -> top of D1
    l2: float           # hole centre -> bottom of D1
    g: float            # slot width  <-- dominates the resonant frequency

    subs_eps: float     # substrate relative permittivity
    subs_tand: float    # substrate loss tangent: this sets the achievable Q
    subs_h: float       # substrate thickness: between top copper and ground
    subs_cells: int     # mesh cells through that thickness
    cu_res: float       # largest mesh cell over the copper C, in x and y (the hole gets D2/10)

    # Board outline is two rectangles: a lower block the hole D2 is centred on, plus a
    # neck at the top carrying the SMA.  (Named blk/neck, not R2 -- R2 is the hole radius.)
    blk_w: float
    blk_l: float
    neck_w: float
    neck_l: float
    mnt_d: float        # corner through-holes in the block; 0 disables
    mnt_inset: float    # corner hole edge to block edge
    mnt_metal: bool     # True = PEC posts (screws/standoffs), worst case

    trace_w: float      # feed trace width
    tap_depth: float    # how far the trace reaches down onto the left arm
    sma_inset: float    # SMA centre pin to board edge: trace ends and port sits here
    cu_sigma: float     # copper conductivity (S/m): finite, not PEC
    cu_t: float         # copper thickness in METRES -- must never appear in a coordinate
    mask_t: float       # solder mask over the whole top; 0 = keepout
    mask_eps: float
    mask_tand: float

    air_xy: float       # air added to the board's x and y extent (total, both sides)
    air_z: float        # air box height, 1/3 below the board and 2/3 above
    f0: float           # centre frequency: FR4 loss, mesh resolution and the B1/E dumps use it
    fc: float           # excitation half-bandwidth: f0 +- fc

    # ---------------- derived ----------------
    R1 = property(lambda s: s.D1 / 2)
    R2 = property(lambda s: s.D2 / 2)
    yc = property(lambda s: (s.l1 - s.l2) / 2)      # centre of the outer disc, from l1 + l2 == D1
    yi = property(lambda s: np.sqrt(s.R2**2 - (s.g / 2)**2))          # slot edge on the hole
    yo = property(lambda s: s.yc + np.sqrt(s.R1**2 - (s.g / 2)**2))  # slot edge on the outer disc
    y_blk = property(lambda s: s.blk_l / 2)
    y_sma = property(lambda s: s.blk_l / 2 + s.neck_l)      # the board edge
    y_feed = property(lambda s: s.y_sma - s.sma_inset)      # SMA centre pin
    trace_x1 = property(lambda s: -s.g / 2)                 # right edge, flush with the slot's left edge
    trace_x0 = property(lambda s: -s.g / 2 - s.trace_w)     # left edge
    trace_y0 = property(lambda s: s.yc + s.R1 - s.tap_depth)  # bottom end, down on the arm
    # floor for dropping near-duplicate mesh lines; see build()
    mesh_floor = property(lambda s: min(0.025, s.g / 8, s.mask_t / 2 if s.mask_t > 0 else np.inf))

    def __post_init__(self):
        # constraints: violate one and the model is silently wrong
        R1, R2, yc, g = self.R1, self.R2, self.yc, self.g
        assert abs((self.l1 + self.l2) - self.D1) < 1e-9, "l1 + l2 must equal D1"
        assert g < self.D2, "slot must be narrower than the hole, else its edge on D2 is undefined"
        assert abs(yc) + R2 < R1, "hole breaks out through the disc: l1/l2 too asymmetric for D2"
        assert R1 <= self.blk_w / 2 and yc + R1 <= self.blk_l / 2 and yc - R1 >= -self.blk_l / 2, \
            "disc does not fit inside the block"
        assert self.trace_w + g / 2 <= self.neck_w / 2 + self.blk_w / 2, "trace wider than the board allows"
        if self.mnt_d > 0:   # corner holes must clear the disc
            mx = self.blk_w / 2 - self.mnt_inset - self.mnt_d / 2
            my = self.blk_l / 2 - self.mnt_inset - self.mnt_d / 2
            assert np.hypot(mx, my - yc) > R1 + self.mnt_d / 2, "corner hole overlaps the disc"

        # the arm narrows towards the top of the disc; the tap has to land fully on copper
        arm_x = -np.sqrt(R1**2 - (self.trace_y0 - yc)**2)
        assert arm_x <= self.trace_x0, (
            f"trace overhangs the disc: arm reaches x={arm_x:.2f} at the tap but trace "
            f"needs {self.trace_x0:.2f} -- raise tap_depth or cut trace_w")
        assert self.trace_x0 >= -self.neck_w / 2, "trace runs off the side of the neck -- widen neck_w"
        assert self.trace_y0 > self.yi, "tap reaches past the slot down to the hole"
        assert self.trace_y0 < self.y_feed <= self.y_sma, "SMA pin must sit on the board, beyond the tap"

    def outline(self):
        """The "C" as [xs, ys]: around the hole skipping the slot, out along the slot
        edge, then the outer disc back around."""
        R1, R2, yc, g = self.R1, self.R2, self.yc, self.g
        th_i = np.arctan2(self.yi, g / 2)                       # ~90 deg
        th_o = np.arctan2(self.yo - yc, g / 2)
        a_in  = np.linspace(th_i, -np.pi - th_i, 181)           # hole, clockwise, skipping slot
        a_out = np.linspace(np.pi - th_o, 2 * np.pi + th_o, 181)  # outer disc, back around
        # AddPolygon wants [xs, ys], not an (N,2) array -- its docstring says otherwise
        return [np.r_[R2 * np.cos(a_in),  R1 * np.cos(a_out)],
                np.r_[R2 * np.sin(a_in),  yc + R1 * np.sin(a_out)]]

    def build(self) -> tuple[ContinuousStructure, LumpedPort]:
        """Return (csx, port), new on every call: FDTD.SetCSX takes ownership and the
        C++ structure is destroyed with that FDTD, so one CSX can never serve two runs."""
        R1, R2, yc, g = self.R1, self.R2, self.yc, self.g
        subs_h, f0 = self.subs_h, self.f0
        y_blk, y_sma, y_feed = self.y_blk, self.y_sma, self.y_feed
        trace_x0, trace_x1, trace_y0 = self.trace_x0, self.trace_x1, self.trace_y0
        bw, nw = self.blk_w / 2, self.neck_w / 2

        csx = ContinuousStructure()
        mesh = csx.GetGrid()
        mesh.SetDeltaUnit(1e-3)

        ### ---------------- the "C": one polygon,erts no boolean ops needed ----------------
        # CSXCAD has no CSG; either overlap primitives by priority or, as here, trace the
        # outline directly -- exact, one primitive, and no priority bookkeeping.
        loop = csx.AddConductingSheet('loop', conductivity=self.cu_sigma, thickness=self.cu_t)
        loop.AddPolygon(self.outline(), norm_dir='z', elevation=subs_h, priority=10)

        # Feed the LEFT arm only: the trace's right edge sits exactly on the left edge of the
        # slot, so it contacts one side of g and never reaches the right arm.  A trace centred
        # on x=0 bridges the gap over its bottom few tenths of a mm, shorting the capacitor
        # that sets the resonance.
        loop.AddBox(priority=10, start=[trace_x0, trace_y0, subs_h],
                                 stop =[trace_x1, y_feed,   subs_h])

        ### ---------------- substrate: block + neck ----------------
        kappa = 2 * np.pi * f0 * EPS0 * self.subs_eps * self.subs_tand   # only valid near f0
        subs = csx.AddMaterial('FR4', epsilon=self.subs_eps, kappa=kappa)
        subs.AddBox(priority=0, start=[-bw, -y_blk, 0], stop=[bw, y_blk, subs_h])
        subs.AddBox(priority=0, start=[-nw,  y_blk, 0], stop=[nw, y_sma, subs_h])

        ### ---------------- ground: solid, same outline as the substrate ----------------
        # Openings come only from the through-holes below (D2 and the four corners), which cut
        # copper and FR4 together at higher priority.
        gnd = csx.AddConductingSheet('gnd', conductivity=self.cu_sigma, thickness=self.cu_t)
        gnd.AddBox(priority=10, start=[-bw, -y_blk, 0], stop=[bw, y_blk, 0])
        gnd.AddBox(priority=10, start=[-nw,  y_blk, 0], stop=[nw, y_sma, 0])

        ### ---------------- solder mask (optional) ----------------
        if self.mask_t > 0:
            top = subs_h + self.mask_t
            mask = csx.AddMaterial('mask', epsilon=self.mask_eps,
                                   kappa=2 * np.pi * f0 * EPS0 * self.mask_eps * self.mask_tand)
            # priority 1: above the FR4, below the copper sheets, and it fills the slot -- which
            # is where it actually changes the answer, not over the open copper.
            mask.AddBox(priority=1, start=[-bw, -y_blk, subs_h], stop=[bw, y_blk, top])
            mask.AddBox(priority=1, start=[-nw,  y_blk, subs_h], stop=[nw, y_sma, top])
            mesh.AddLine('z', [top])

        ### ---------------- through-holes: cut FR4, not just copper ----------------
        # D2 is a real through-hole, so the plug of lossy FR4 has to go -- it sits exactly
        # where B1 peaks.  CSXCAD has no CSG, so overlap at higher priority instead.
        air = csx.AddMaterial('air', epsilon=1.0)
        air.AddCylinder([0, 0, -1], [0, 0, subs_h + 1], R2, priority=20)

        if self.mnt_d > 0:   # four corner holes in the block
            m = self.mnt_d / 2
            posts = csx.AddMetal('posts') if self.mnt_metal else air
            for sx in (-1, 1):
                for sy in (-1, 1):
                    cx, cy = sx * (bw - self.mnt_inset - m), sy * (y_blk - self.mnt_inset - m)
                    posts.AddCylinder([cx, cy, -1], [cx, cy, subs_h + 1], m, priority=20)
                    mesh.AddLine('x', [cx - m, cx + m]); mesh.AddLine('y', [cy - m, cy + m])

        ### ---------------- mesh: air box, then pin every critical edge ----------------
        box_x = self.blk_w + self.air_xy
        box_y = self.blk_l + self.neck_l + self.air_xy
        mesh.AddLine('x', [-box_x / 2, box_x / 2])
        # board is not y-symmetric (block below, neck above), so centre the air box
        # on the board, not on the origin, or the +y boundary ends up inside lambda/4
        y_mid = (y_sma - y_blk) / 2
        mesh.AddLine('y', [y_mid - box_y / 2, y_mid + box_y / 2])
        mesh.AddLine('z', [-self.air_z / 3, self.air_z * 2 / 3])

        # what FDTD.AddEdges2Grid(dirs='xy', ...) does, without needing the FDTD object
        for prim in subs.GetAllPrimitives() + gnd.GetAllPrimitives():
            hint = mesh_hint_from_primitive(prim, 'xy')
            for n in range(3):
                if hint is not None and hint[n] is not None:
                    mesh.AddLine(n, hint[n])

        # the slot is the capacitor -- resolve it or the resonance is simply wrong
        # 4 cells: once graded, 2 and 8 gave the same f0 (+-0.3%) and Zin (+-3 ohm). Over 7
        # cells and they drop under mesh_floor (g/8), which silently deletes them.
        mesh.AddLine('x', np.linspace(-g / 2, g / 2, 5))   # odd count -> x=0 lands on a line
        mesh.AddLine('x', [-R1, -R2, 0, R2, R1, trace_x0, trace_x1])
        mesh.AddLine('y', [yc - R1, -R2, 0, R2, self.yi, self.yo, yc + R1, trace_y0])
        # NOTE: do not pin lines along the slot's LENGTH. The slot is uniform that way, so the
        # lines buy nothing, and they collide with other pinned features -- one landed 19.4 um
        # from a corner-hole edge, which forced a 106x grading jump and cost ~1.6x on dt.
        # Only the slot's WIDTH needs resolving, which the x-linspace above does.
        mesh.AddLine('z', np.linspace(0, subs_h, self.subs_cells + 1))

        ### ---------------- feed port (through-hole SMA centre pin modelled as a lumped port) -----
        port = LumpedPort(csx, 1, 50, [trace_x0, y_feed, 0], [trace_x1, y_feed, subs_h],
                          'z', 1.0, priority=5)
        mesh.AddLine('x', [trace_x0, trace_x1]); mesh.AddLine('y', y_feed)   # = edges2grid='xy'

        ### ---------------- B1 in the hole: the number that actually matters ----------------
        b1 = csx.AddDump('B1', dump_type=11, dump_mode=2, file_type=1, frequency=[f0])  # H, freq domain, HDF5
        b1.AddBox(start=[-R2, -R2, 0], stop=[R2, R2, subs_h + 2])
        # ...and E over the same box (so the same cells), which has to stay small there
        e = csx.AddDump('E', dump_type=10, dump_mode=2, file_type=1, frequency=[f0])    # E, freq domain, HDF5
        e.AddBox(start=[-R2, -R2, 0], stop=[R2, R2, subs_h + 2])

        # LAST: everything that adds mesh lines has to come before this.
        # Drop near-duplicate mesh lines. Pinning both the slot edge on the outer circle
        # and the top of that circle leaves them ~2 um apart, which cuts the Courant
        # timestep ~20x for zero accuracy gain -- and it silently comes back every time
        # D1/D2/g change.
        # The floor must scale with g or a narrow slot gets its OWN lines deleted, which
        # silently destroys the capacitor that sets the resonance.
        # Then grade, not SmoothMeshLines: its 35x jumps put f0 7.6% low, Zin off by ~15 ohm,
        # and cost 1.7x on dt (2026-10-03 runs).
        # Cap the cells over the copper, not only at its edges. Graded out to lambda/20 alone,
        # the C was 1-2 mm cells and the 1 mm hole 2 cells tall -- a square, to the solver --
        # which put f0 7% low: 2.285 GHz, vs 2.435 at cu_res 0.25 and 2.450 at 0.125.  dt was
        # 60.5 fs in all of them; the cost is cells only (2026-10-04 runs).
        max_res = C0 / (f0 + self.fc) / 1e-3 / 20   # lambda/20 at the top of the band
        fine = dict(x=[(-R2, R2, R2 / 5), (-R1, R1, self.cu_res)],
                    y=[(-R2, R2, R2 / 5), (yc - R1, yc + R1, self.cu_res)], z=[])
        for d, ax in enumerate('xyz'):
            L = np.unique(mesh.GetLines(d))
            mesh.SetLines(ax, graded(L[np.r_[True, np.diff(L) > self.mesh_floor]], max_res, 1.4, fine[ax]))
        return csx, port
