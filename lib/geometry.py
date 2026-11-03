# -*- coding: utf-8 -*-
"""
LoopGap: the antenna's parameters and their constraints, its copper shapes, and build(), the
openEMS model.  Lengths in mm, the hole D2 at the origin, +y toward the SMA.
"""

from dataclasses import dataclass
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS.ports import LumpedPort
from openEMS.automesh import mesh_hint_from_primitive
from openEMS.physical_constants import C0, EPS0


def graded(L, max_res, ratio, regions=()):
    """Mesh lines through every line in L, cells up to ~max_res, each within ~ratio of its
    neighbours (up to ~1.7x where pinned lines force a whole number of cells).  regions:
    (lo, hi, cap) limits the cells to cap inside [lo, hi].  CSXCAD's SmoothMeshLines doesn't
    grade: it only splits gaps wider than max_res."""
    L = np.unique(L)
    x = np.unique(np.r_[L, np.linspace(L[0], L[-1], 20001)])
    growth = np.log(ratio)   # per unit distance: the continuous form of a geometric taper
    dist = np.maximum(L[:-1, None] - x, x - L[1:, None]).clip(0)
    size = np.minimum(max_res, (np.diff(L)[:, None] + growth * dist).min(0))
    for lo, hi, cap in regions:
        size = np.minimum(size, cap + growth * np.maximum(lo - x, x - hi).clip(0))
    cells = np.r_[0, np.cumsum(np.diff(x) * (1 / size[1:] + 1 / size[:-1]) / 2)]
    out = [L[0]]
    for a, b in zip(L[:-1], L[1:]):
        ca, cb = np.interp([a, b], x, cells)
        n = max(1, int(np.ceil(cb - ca - 0.25)))   # -0.25: don't split a cell for a 1.1x excess
        out += list(np.interp(ca + (cb - ca) * np.arange(1, n) / n, cells, x)) + [b]
    return np.array(out)


@dataclass(frozen=True)   # the asserts run only on construction: change it via dataclasses.replace()
class LoopGap:
    # Lengths in mm, except cu_t in metres: openEMS takes sigma * t as an SI sheet conductance.
    D1: float           # outer copper diameter
    D2: float           # central hole diameter
    l1: float           # hole centre to the top of D1
    l2: float           # hole centre to the bottom of D1
    g: float            # slot width

    subs_eps: float
    subs_tand: float
    subs_h: float
    subs_cells: int     # mesh cells through subs_h
    cu_res: float       # largest mesh cell over the C, in x and y

    # board outline: a block centred on the hole, and a neck from it to the SMA edge
    blk_w: float
    blk_l: float
    neck_w: float
    neck_l: float
    mnt_d: float        # corner hole diameter; 0 = none
    mnt_inset: float    # corner hole edge to block edge
    mnt_metal: bool     # metal posts in the corner holes instead of air

    trace_w: float      # feed width at the loop
    trace_w0: float     # feed width at the SMA
    tap_depth: float    # feed overlap down onto the left arm
    sma_inset: float    # board edge to the SMA pin, where the port sits
    cu_sigma: float
    cu_t: float         # metres
    mask_t: float       # solder mask over the whole top; 0 = none
    mask_eps: float
    mask_tand: float

    air_xy: float       # air beyond the board in x and y, both sides together
    air_z: float        # air box height, 1/3 below the board and 2/3 above
    f0: float           # centre frequency: sets the FR4 loss, the mesh and the dumps
    fc: float           # excitation half-bandwidth

    R1 = property(lambda s: s.D1 / 2)
    R2 = property(lambda s: s.D2 / 2)
    yc = property(lambda s: (s.l1 - s.l2) / 2)                         # outer disc centre
    yi = property(lambda s: np.sqrt(s.R2**2 - (s.g / 2)**2))           # slot edge on the hole
    yo = property(lambda s: s.yc + np.sqrt(s.R1**2 - (s.g / 2)**2))   # slot edge on the disc
    y_blk = property(lambda s: s.blk_l / 2)
    y_sma = property(lambda s: s.blk_l / 2 + s.neck_l)                 # the board edge
    y_feed = property(lambda s: s.y_sma - s.sma_inset)
    trace_x1 = property(lambda s: -s.g / 2)   # flush with the slot's left edge: feeds one arm only
    trace_x0 = property(lambda s: -s.g / 2 - s.trace_w)
    # the taper widens equally about the tap's centreline
    feed_x0 = property(lambda s: s.trace_x1 - (s.trace_w + s.trace_w0) / 2)
    feed_x1 = property(lambda s: s.trace_x1 - (s.trace_w - s.trace_w0) / 2)
    trace_y0 = property(lambda s: s.yc + s.R1 - s.tap_depth)
    # mesh lines closer than this merge; it scales with g and mask_t, or it deletes their own lines
    mesh_floor = property(lambda s: min(0.025, s.g / 8, s.mask_t / 2 if s.mask_t > 0 else np.inf))

    def __post_init__(self):
        R1, R2, yc, g = self.R1, self.R2, self.yc, self.g
        assert abs((self.l1 + self.l2) - self.D1) < 1e-9, "l1 + l2 must equal D1"
        assert g < self.D2, "slot must be narrower than the hole, else its edge on D2 is undefined"
        assert abs(yc) + R2 < R1, "hole breaks out through the disc: l1/l2 too asymmetric for D2"
        assert R1 <= self.blk_w / 2 and yc + R1 <= self.blk_l / 2 and yc - R1 >= -self.blk_l / 2, \
            "disc does not fit inside the block"
        assert self.trace_w + g / 2 <= self.neck_w / 2 + self.blk_w / 2, "trace wider than the board allows"
        if self.mnt_d > 0:
            mx = self.blk_w / 2 - self.mnt_inset - self.mnt_d / 2
            my = self.blk_l / 2 - self.mnt_inset - self.mnt_d / 2
            assert np.hypot(mx, my - yc) > R1 + self.mnt_d / 2, "corner hole overlaps the disc"
        arm_x = -np.sqrt(R1**2 - (self.trace_y0 - yc)**2)
        assert arm_x <= self.trace_x0, (
            f"trace overhangs the disc: arm reaches x={arm_x:.2f} at the tap but trace "
            f"needs {self.trace_x0:.2f} -- raise tap_depth or cut trace_w")
        assert min(self.trace_x0, self.feed_x0) >= -self.neck_w / 2 and self.feed_x1 <= self.neck_w / 2, \
            "feed runs off the side of the neck -- widen neck_w"
        assert self.trace_y0 > self.yi, "tap reaches past the slot down to the hole"
        assert self.trace_y0 < self.y_feed <= self.y_sma, "SMA pin must sit on the board, beyond the tap"

    def outline(self):
        """The C as [xs, ys] (what AddPolygon takes; its docstring says (N, 2)): around the hole
        skipping the slot, out along the slot edge, then back around the outer disc."""
        R1, R2, yc, g = self.R1, self.R2, self.yc, self.g
        t_hole = np.arctan2(self.yi, g / 2)
        t_disc = np.arctan2(self.yo - yc, g / 2)
        hole = np.linspace(t_hole, -np.pi - t_hole, 181)
        disc = np.linspace(np.pi - t_disc, 2 * np.pi + t_disc, 181)
        return [np.r_[R2 * np.cos(hole), R1 * np.cos(disc)],
                np.r_[R2 * np.sin(hole), yc + R1 * np.sin(disc)]]

    def feed(self):
        """The feed as [xs, ys]: a linear taper (Sasaki et al. give no shape) from trace_w0 at
        the SMA pin to trace_w at the top of the disc, then straight down onto the left arm."""
        y_top = self.yc + self.R1
        return [[self.feed_x0, self.feed_x1, self.trace_x1, self.trace_x1, self.trace_x0, self.trace_x0],
                [self.y_feed,  self.y_feed,  y_top,         self.trace_y0, self.trace_y0, y_top]]

    def build(self) -> tuple[ContinuousStructure, LumpedPort]:
        """(csx, port), built fresh each call: FDTD.SetCSX takes ownership, so a CSX can't serve
        two runs."""
        csx = ContinuousStructure()
        csx.GetGrid().SetDeltaUnit(1e-3)
        self._add_loop(csx)
        board = self._add_board(csx)
        self._add_holes(csx)
        self._pin_mesh_lines(csx, board)
        port = self._add_port(csx)
        self._add_dumps(csx)
        self._grade_mesh(csx)
        return csx, port

    def _kappa(self, eps, tand):
        """The conductivity that gives loss tangent tand, exact at f0 only."""
        return 2 * np.pi * self.f0 * EPS0 * eps * tand

    def _board_boxes(self, z0, z1):
        """The board outline, block then neck, as (start, stop) boxes from z0 to z1."""
        bw, nw = self.blk_w / 2, self.neck_w / 2
        return [([-bw, -self.y_blk, z0], [bw, self.y_blk, z1]),
                ([-nw,  self.y_blk, z0], [nw, self.y_sma, z1])]

    def _add_loop(self, csx):
        loop = csx.AddConductingSheet('loop', conductivity=self.cu_sigma, thickness=self.cu_t)
        for shape in (self.outline(), self.feed()):
            loop.AddPolygon(shape, norm_dir='z', elevation=self.subs_h, priority=10)

    def _add_board(self, csx):
        """FR4, the ground under it and the mask over it.  Returns the FR4 and ground
        primitives: the mesh needs their edges."""
        h = self.subs_h
        subs = csx.AddMaterial('FR4', epsilon=self.subs_eps, kappa=self._kappa(self.subs_eps, self.subs_tand))
        for start, stop in self._board_boxes(0, h):
            subs.AddBox(start, stop, priority=0)
        gnd = csx.AddConductingSheet('gnd', conductivity=self.cu_sigma, thickness=self.cu_t)
        for start, stop in self._board_boxes(0, 0):
            gnd.AddBox(start, stop, priority=10)
        if self.mask_t > 0:
            top = h + self.mask_t
            mask = csx.AddMaterial('mask', epsilon=self.mask_eps, kappa=self._kappa(self.mask_eps, self.mask_tand))
            for start, stop in self._board_boxes(h, top):
                mask.AddBox(start, stop, priority=1)
            csx.GetGrid().AddLine('z', [top])
        return subs.GetAllPrimitives() + gnd.GetAllPrimitives()

    def _add_holes(self, csx):
        """Through-holes as air (or metal) cylinders at a higher priority, cutting the FR4 and
        the ground in one go: CSXCAD has no CSG."""
        air = csx.AddMaterial('air', epsilon=1.0)
        air.AddCylinder([0, 0, -1], [0, 0, self.subs_h + 1], self.R2, priority=20)
        if self.mnt_d > 0:
            m = self.mnt_d / 2
            posts = csx.AddMetal('posts') if self.mnt_metal else air
            mesh = csx.GetGrid()
            for sx in (-1, 1):
                for sy in (-1, 1):
                    cx = sx * (self.blk_w / 2 - self.mnt_inset - m)
                    cy = sy * (self.y_blk - self.mnt_inset - m)
                    posts.AddCylinder([cx, cy, -1], [cx, cy, self.subs_h + 1], m, priority=20)
                    mesh.AddLine('x', [cx - m, cx + m]); mesh.AddLine('y', [cy - m, cy + m])

    def _pin_mesh_lines(self, csx, board):
        """The air box, and a mesh line on every edge that matters."""
        mesh = csx.GetGrid()
        R1, R2, yc, g = self.R1, self.R2, self.yc, self.g
        box_x = self.blk_w + self.air_xy
        box_y = self.blk_l + self.neck_l + self.air_xy
        board_mid_y = (self.y_sma - self.y_blk) / 2
        mesh.AddLine('x', [-box_x / 2, box_x / 2])
        mesh.AddLine('y', [board_mid_y - box_y / 2, board_mid_y + box_y / 2])
        mesh.AddLine('z', [-self.air_z / 3, self.air_z * 2 / 3])
        for prim in board:   # what FDTD.AddEdges2Grid(dirs='xy') does, without an FDTD object
            hint = mesh_hint_from_primitive(prim, 'xy')
            for n in range(3):
                if hint is not None and hint[n] is not None:
                    mesh.AddLine(n, hint[n])
        mesh.AddLine('x', np.linspace(-g / 2, g / 2, 5))   # 4 cells across the slot, a line on x = 0
        mesh.AddLine('x', [-R1, -R2, 0, R2, R1, self.trace_x0, self.trace_x1])
        # none along the slot's length: they collide with other edges and wreck the grading
        mesh.AddLine('y', [yc - R1, -R2, 0, R2, self.yi, self.yo, yc + R1, self.trace_y0])
        mesh.AddLine('z', np.linspace(0, self.subs_h, self.subs_cells + 1))

    def _add_port(self, csx):
        """The SMA centre pin, as a 50 ohm lumped port through the board."""
        port = LumpedPort(csx, 1, 50, [self.feed_x0, self.y_feed, 0], [self.feed_x1, self.y_feed, self.subs_h],
                          'z', 1.0, priority=5)
        mesh = csx.GetGrid()
        mesh.AddLine('x', [self.feed_x0, self.feed_x1]); mesh.AddLine('y', self.y_feed)
        return port

    def _add_dumps(self, csx):
        """B1 over f0 +- 0.5 GHz and E at f0, over the hole."""
        box = dict(start=[-self.R2, -self.R2, 0], stop=[self.R2, self.R2, self.subs_h + 2])
        # dump_type 11 = H, 10 = E, frequency domain; dump_mode 2 = at cell centres; file_type 1 = HDF5
        b1 = csx.AddDump('B1', dump_type=11, dump_mode=2, file_type=1,
                         frequency=list(np.linspace(self.f0 - 0.5e9, self.f0 + 0.5e9, 21)))
        b1.AddBox(**box)
        e = csx.AddDump('E', dump_type=10, dump_mode=2, file_type=1, frequency=[self.f0])
        e.AddBox(**box)

    def _grade_mesh(self, csx):
        """Merge near-duplicate lines, then grade each axis.  Runs last: it replaces every line."""
        mesh = csx.GetGrid()
        R1, R2, yc = self.R1, self.R2, self.yc
        max_res = C0 / (self.f0 + self.fc) / 1e-3 / 20   # lambda/20 at the top of the band
        # a taper's slanted edges need the fine cells all the way up to the SMA
        fine_y_top = self.y_feed if self.trace_w0 != self.trace_w else yc + R1
        fine = dict(x=[(-R2, R2, R2 / 5), (-R1, R1, self.cu_res)],
                    y=[(-R2, R2, R2 / 5), (yc - R1, fine_y_top, self.cu_res)], z=[])
        for d, ax in enumerate('xyz'):
            L = np.unique(mesh.GetLines(d))
            kept = L[np.r_[True, np.diff(L) > self.mesh_floor]]
            mesh.SetLines(ax, graded(kept, max_res, 1.4, fine[ax]))
