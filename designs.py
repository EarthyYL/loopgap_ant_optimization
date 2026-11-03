# -*- coding: utf-8 -*-
"""
The designs, rev by rev: values only.  Plots say "Rev 0"; code and runs/ say rev0.
    rev0  Sewani et al.'s board as built (Am. J. Phys. 88, 1156 (2020); Sewani_Files/NV_PCB_Antenna.brd)
    rev1  run 1's CMA-ES best on rev0's board
    rev2  the board to order: Sasaki et al.'s design (Rev. Sci. Instrum. 87, 053904 (2016)) on 1.6 mm
          FR4, slot widened to 0.25 mm, 1 oz copper, one 26 x 48.7 mm rectangle
"""

from dataclasses import replace
from lib.geometry import LoopGap

# Lengths in mm; cu_t in metres.
rev2 = LoopGap(
    D1=14.0, D2=1.0,
    l1=10.9, l2=3.1,
    g=0.25,                            # Sasaki 0.1; 0.25 moves f_res 2.79 -> 2.86 GHz
    subs_eps=4.3, subs_tand=0.03,      # FR4
    subs_h=1.6,
    subs_cells=4,
    cu_res=0.25,                       # 0.125 for final numbers: f0 +0.6 %, 2.3x the runtime
    blk_w=26.0, blk_l=26.0,
    neck_w=26.0, neck_l=22.7,          # neck_w = blk_w: one rectangle; neck_l fits Sasaki's 23 mm feed
    mnt_d=4.0, mnt_inset=1.0,          # 3 mm plated holes, 4 mm with their annular ring
    mnt_metal=False,
    trace_w=0.54, trace_w0=3.3,        # Sasaki's 50 ohm taper
    tap_depth=1.0, sma_inset=1.8,      # Sewani's .brd
    cu_sigma=58e6, cu_t=35e-6,         # 1 oz
    # no mask: mask_eps=1, mask_tand=0 (mask_t=0 coarsens the z mesh above the board)
    mask_t=0.03, mask_eps=4.3, mask_tand=0.025,
    air_xy=60.0, air_z=80.0,
    f0=2.87e9, fc=1.5e9,
)
rev0 = replace(rev2, subs_h=0.254, neck_w=20.0, neck_l=5.5, g=0.1, trace_w=0.476, trace_w0=0.476, cu_t=36e-6)
rev1 = replace(rev0, D1=11.959, D2=0.831, l1=11.959 / 2 + 3.390, l2=11.959 / 2 - 3.390, g=0.105,
               tap_depth=0.342)   # runs/opt/e57b2d7dd263
