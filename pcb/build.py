# -*- coding: utf-8 -*-
"""
Write Rev 2 as a KiCad board, every length from designs.rev2, so the board is the simulated
geometry.  Run from the repo root with /home/yluo/opt/openEMS/venv/bin/python:
    python -m pcb.build    # boards/rev2/rev2.kicad_pcb + rev2.kicad_pro
Then open it in KiCad and plot the Gerbers and drill files per JLCPCB's KiCad guide.  Rebuilding
overwrites both files, so change the design in designs.py, not in KiCad.

KiCad (X, Y) = (X0 - y, Y0 - x) from the model's (x, y): turned 90 degrees, not mirrored, SMA on
the left.  That way the SMA footprint goes in at 0 degrees, unedited (at any other angle KiCad's
board file needs every pad's angle rewritten).
"""

import json, os, re
import numpy as np
from matplotlib.path import Path
from designs import rev2 as geo
from settings import Boards

X0, Y0 = 150.0, 100.0   # where the hole sits on the KiCad page
HOLE_DRILL = 0.6        # unplated; the copper opening stays the model's D2, 0.2 mm out (JLCPCB clears 0.2)
MNT_DRILL = 3.2         # corner holes: M3 clearance (the user, 2026-11-02); their pad is the model's mnt_d
# KiCad's DRC rules (Board Setup > Constraints): JLCPCB's standard 2-layer minimums, from
# jlcpcb.com/capabilities/pcb-capabilities (2026-11-02).  Its same-net spacing (0.25) and smallest
# unplated hole (0.5) aren't KiCad constraints: check() has them.
JLCPCB = dict(min_clearance=0.10, min_track_width=0.10,
              min_through_hole_diameter=0.15, min_via_diameter=0.25,
              min_via_annular_width=0.18,        # their PTH ring minimum (0.25 recommended)
              min_hole_to_hole=0.45,             # pad holes (0.2 between via holes)
              min_hole_clearance=0.2,            # hole to copper; the centre hole sits exactly on it
              min_text_height=1.0, min_text_thickness=0.15,
              min_copper_edge_clearance=0.5)     # JLCPCB allows 0.2; KiCad's 0.5 kept, a margin for routing
SMA_FP = ('/mnt/c/Program Files/KiCad/9.0/share/kicad/footprints/Connector_Coaxial.pretty/'
          'SMA_Molex_73251-1153_EdgeMount_Horizontal.kicad_mod')
SMA_EDGE = -4.76        # footprint x of the connector's face = the board edge (from its F.Fab drawing)
SMA_PAD_END = -4.26     # footprint x of its centre pad's outer end

w = geo.blk_w / 2
board = [[-w, w, w, -w], [-geo.y_blk, -geo.y_blk, geo.y_sma, geo.y_sma]]
sma_x = (geo.feed_x0 + geo.feed_x1) / 2   # the feed's centreline at the SMA: the centre pin goes here
sma_y = geo.y_sma + SMA_EDGE              # footprint origin; the footprint's -x points off the board
y_pad = sma_y - SMA_PAD_END
# the model's port sits sma_inset in from the edge: run the wide end on to the SMA's pad
runon = [[geo.feed_x0, geo.feed_x1, geo.feed_x1, geo.feed_x0], [geo.y_feed, geo.y_feed, y_pad, y_pad]]
copper = [geo.outline(), geo.feed(), runon]   # top copper, all one net
mx = w - geo.mnt_inset - geo.mnt_d / 2        # corner holes, where the model has them
my = geo.y_blk - geo.mnt_inset - geo.mnt_d / 2


def check():
    """KiCad's DRC can't see a gap inside one net, so check the slot here: open and g wide along its
    whole length, copper on both sides (the feed touches one arm only).  And the taper's ends."""
    cu = [Path(np.column_stack(s)) for s in copper]
    on = lambda x, y: any(p.contains_point((x, y)) for p in cu)
    g, e = geo.g, 0.005
    for y in np.linspace(geo.yi + e, geo.yo - e, 100):
        assert not any(on(x, y) for x in (-g / 2 + e, 0, g / 2 - e)), f"copper in the slot at y = {y:.2f}"
        assert on(-g / 2 - e, y) and on(g / 2 + e, y), f"slot wider than g at y = {y:.2f}"
    xs, ys = geo.feed()
    assert np.isclose(xs[1] - xs[0], geo.trace_w0) and np.isclose(xs[2] - xs[4], geo.trace_w), "taper ends"
    assert geo.neck_w == geo.blk_w, "the board is one rectangle: neck_w must equal blk_w"
    assert g >= 0.25 and HOLE_DRILL >= 0.5, "JLCPCB: same-net gaps >= 0.25 mm, unplated holes >= 0.5 mm"
    print(f"slot       open, {g} mm wide, y {geo.yi:.3f} to {geo.yo:.3f}")
    print(f"taper      {xs[1] - xs[0]:.2f} -> {xs[2] - xs[4]:.2f} mm over {ys[0] - ys[2]:.2f} mm")


def xy(shape):
    """A model-frame [xs, ys] as KiCad's (xy X Y) list."""
    return ' '.join(f'(xy {X0 - y:.6f} {Y0 - x:.6f})' for x, y in zip(*shape))


def poly(shape, layer, extra):
    return f'\t(gr_poly (pts {xy(shape)}) (stroke (width {0.05 if layer == "Edge.Cuts" else 0}) ' \
           f'(type default)) {extra} (layer "{layer}"))\n'


def hole(ref, x, y, drill, size, gnd):
    """A one-pad footprint: plated to ground, or unplated."""
    pad = (f'(pad "1" thru_hole circle (at 0 0) (size {size} {size}) (drill {drill}) (layers "*.Cu" "*.Mask") '
           f'(net 1 "GND"))' if gnd else
           f'(pad "" np_thru_hole circle (at 0 0) (size {size} {size}) (drill {drill}) (layers "*.Cu" "*.Mask"))')
    return (f'\t(footprint "MountingHole" (layer "F.Cu") (at {X0 - y:.6f} {Y0 - x:.6f})\n'
            f'\t\t(property "Reference" "{ref}" (at 0 0 0) (layer "F.Fab") (hide yes) '
            f'(effects (font (size 1 1) (thickness 0.15))))\n'
            f'\t\t(attr exclude_from_pos_files exclude_from_bom)\n\t\t{pad})\n')


def sma():
    """KiCad's own footprint for the Molex 73251-1153, placed, with nets on its pads."""
    fp = open(SMA_FP).read()
    fp = re.sub(r'\n\t\((version|generator|generator_version) [^\n]*', '', fp)   # footprint-file header only
    fp = fp.replace('(footprint "', '(footprint "Connector_Coaxial:', 1).replace('"REF**"', '"J1"')
    fp = fp.replace('\n', f'\n\t(at {X0 - sma_y:.6f} {Y0 - sma_x:.6f})\n', 1)
    fp = re.sub(r'\(pad "1" \w+ \w+', r'\g<0> (net 2 "SIG")', fp)
    return re.sub(r'\(pad "2" \w+ \w+', r'\g<0> (net 1 "GND")', fp)


def text(s, x, y, layer='F.SilkS', size=1.0):
    """Silkscreen text, size mm tall (JLCPCB: >= 1.0, stroke >= 0.15).  Bottom text is mirrored so
    it reads right from below."""
    mirror = ' (justify mirror)' if layer.startswith('B.') else ''
    return (f'\t(gr_text "{s}" (at {X0 - y:.6f} {Y0 - x:.6f} 0) (layer "{layer}") '
            f'(effects (font (size {size} {size}) (thickness {0.15 * size:g})){mirror}))\n')


HEADER = '''(kicad_pcb
\t(version 20241229)
\t(generator "pcb.build")
\t(generator_version "9.0")
\t(general (thickness 1.6) (legacy_teardrops no))
\t(paper "A4")
\t(layers
\t\t(0 "F.Cu" signal) (2 "B.Cu" signal)
\t\t(9 "F.Adhes" user "F.Adhesive") (11 "B.Adhes" user "B.Adhesive")
\t\t(13 "F.Paste" user) (15 "B.Paste" user)
\t\t(5 "F.SilkS" user "F.Silkscreen") (7 "B.SilkS" user "B.Silkscreen")
\t\t(1 "F.Mask" user) (3 "B.Mask" user)
\t\t(17 "Dwgs.User" user "User.Drawings") (19 "Cmts.User" user "User.Comments")
\t\t(21 "Eco1.User" user "User.Eco1") (23 "Eco2.User" user "User.Eco2")
\t\t(25 "Edge.Cuts" user) (27 "Margin" user)
\t\t(31 "F.CrtYd" user "F.Courtyard") (29 "B.CrtYd" user "B.Courtyard")
\t\t(35 "F.Fab" user) (33 "B.Fab" user)
\t)
\t(setup (pad_to_mask_clearance 0) (allow_soldermask_bridges_in_footprints no) (tenting front back))
\t(net 0 "")
\t(net 1 "GND")
\t(net 2 "SIG")
'''

if __name__ == "__main__":
    check()
    corners = [(sx * mx, sy * my) for sx in (-1, 1) for sy in (-1, 1)]
    pcb = (HEADER
           + poly(board, 'Edge.Cuts', '(fill no)')
           + ''.join(poly(s, 'F.Cu', '(fill yes) (net 2)') for s in copper)
           # bottom: one ground zone; clearance 0.2 gives the model's 1.0 mm opening at the centre hole
           + f'\t(zone (net 1) (net_name "GND") (layer "B.Cu") (hatch edge 0.5) (connect_pads yes (clearance 0.2))\n'
             f'\t\t(min_thickness 0.25) (filled_areas_thickness no) (fill (thermal_gap 0.5) (thermal_bridge_width 0.5))\n'
             f'\t\t(polygon (pts {xy(board)})))\n'
           + hole('H1', 0, 0, HOLE_DRILL, HOLE_DRILL, False)
           + ''.join(hole(f'H{i + 2}', x, y, MNT_DRILL, geo.mnt_d, True) for i, (x, y) in enumerate(corners))
           + sma()
           + text('Rev 2', -w / 2, 20)   # feed end, clear of the loop
           # bottom ink does nothing electrically (the field is between loop and ground); clear of the holes
           + text('Carolina Chip Fab', -5, (geo.y_sma - geo.y_blk) / 2, 'B.SilkS', 2.0)
           + ')\n')
    out = os.path.join(Boards, 'rev2')
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'rev2.kicad_pcb'), 'w') as f:
        f.write(pcb)
    with open(os.path.join(out, 'rev2.kicad_pro'), 'w') as f:
        json.dump({'board': {'design_settings': {'rules': JLCPCB}},
                   'meta': {'filename': 'rev2.kicad_pro', 'version': 1}}, f, indent=2)
    print(f"SMA pin    model ({sma_x:.3f}, {geo.y_sma}), the board edge")
    print(f"board      {geo.blk_w} x {geo.y_sma + geo.y_blk} mm  ->  {out}")
