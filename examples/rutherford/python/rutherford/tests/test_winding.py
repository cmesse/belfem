# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

import math
import sys

sys.path.insert(0, ".")

from rutherford.winding import Winding, WindingError


def check(label, ok):
    print("{:s}: {:s}".format("PASS" if ok else "FAIL", label))
    if not ok:
        sys.exit(1)


def make(**kw):
    args = dict(tapeWidth=4, layAngle=30, edgeRadius=0.9, cableWidth=12,
                numTapes=1, numLayers=1, layerGap=200, numPitches=1,
                tapeResolution=0.4, gapResolution=0.8, sense=1)
    args.update(kw)
    return Winding(**args)


# default config of main.py
W = make()
P = 2 * 10.2 + 2 * math.pi * 0.9
check("flat length", abs(W.ws - 10.2) < 1e-14)
check("perimeter", abs(W.layers[0].perimeter - P) < 1e-12)
check("pitch = P / tan 30", abs(W.pitch - P / math.tan(math.radians(30))) < 1e-12)
check("length = one pitch", abs(W.length - W.pitch) < 1e-12)
L0 = W.layers[0]
check("tape span across = wt / cos", abs(L0.tape_span - 4 / math.cos(math.radians(30))) < 1e-12)
check("gap span closes the sector", abs(L0.tape_span + L0.gap_span - P) < 1e-12)
check("tape centered at the top face", abs(L0.tapes[0].center) < 1e-14)
check("tape shift 0 for one pitch", W.tape_shift() == 0)
check("n_s tape = ceil(4.62/0.4) = 12", W.n_s(0, "tape") == 12)
check("n_s gap = ceil(gap/0.8)", W.n_s(0, "gap") == math.ceil(L0.gap_span / 0.8))
dz = 0.4 * math.cos(math.radians(30))
check("rows from the along-tape step", W.z_partition_count() == math.ceil(W.length / dz))
check("bending radius r / sin^2", abs(L0.bending_radius() - 0.9 / 0.25) < 1e-12)
check("coverage", abs(L0.coverage() - L0.tape_span / P) < 1e-14)

# three tapes, a third of a pitch: closes by one tape position
W3 = make(numTapes=3, numPitches=1.0 / 3.0)
check("third pitch with 3 tapes closes", W3.tape_shift() == 1)
check("four thirds too", make(numTapes=3, numPitches=4.0 / 3.0).tape_shift() == 1)
check("sector = P / 3", abs(W3.layers[0].sector - P / 3) < 1e-12)

# fractional pitch without closure: refused
try:
    make(numPitches=0.5)
    ok = False
except WindingError:
    ok = True
check("half pitch refused (no rotation closure)", ok)

# overlap: 6 mm tapes, 5 of them at 30 deg need 34.6 mm > P
try:
    make(tapeWidth=6, numTapes=5)
    ok = False
except WindingError:
    ok = True
check("overlapping tapes refused", ok)

# two layers: common pitch, layer-specific lay angle, counter-wound
W2 = make(numLayers=2, tapeResolution=0.2, gapResolution=0.4)
l0, l1 = W2.layers
check("layer 1 edge radius", abs(l1.radius - 1.1) < 1e-14)
check("layer 1 counter-wound", l1.sigma == -1 and l0.sigma == 1)
check("common pitch: tan theta_1 = P_1 / L",
      abs(l1.tan_theta - l1.perimeter / W2.pitch) < 1e-14)
check("layer 1 steeper", l1.theta > l0.theta)
check("sense flips both", make(sense=-1, numLayers=2, tapeResolution=0.2,
                                gapResolution=0.4).layers[1].sigma == 1)

# thin-slab rule
try:
    make(numLayers=2, tapeResolution=0.5)
    ok = False
except WindingError:
    ok = True
check("coarse tape resolution on a thin slab refused", ok)

# geometry guards
for kw in (dict(edgeRadius=7.0), dict(layAngle=0), dict(layAngle=90),
           dict(numTapes=0), dict(sense=2)):
    try:
        make(**kw)
        ok = False
    except WindingError:
        ok = True
    check("refused: {}".format(kw), ok)

# conductor groups: tape k continues into tape k + q across the planes
check("5 tapes, 1/5 pitch: one conductor",
      make(numTapes=5, numPitches=0.2).conductor_groups() == [[0, 1, 2, 3, 4]])
check("5 tapes, 1 pitch: five conductors",
      make(numTapes=5, numPitches=1).conductor_groups() == [[0], [1], [2], [3], [4]])
check("4 tapes, 1/2 pitch: two conductors",
      make(numTapes=4, numPitches=0.5).conductor_groups() == [[0, 2], [1, 3]])
check("two layers group per layer",
      make(numTapes=3, numPitches=1.0 / 3.0, numLayers=2, tapeResolution=0.2,
           gapResolution=0.4).conductor_groups() == [[0, 1, 2], [3, 4, 5]])

# round core (edge radius = half width) is allowed
R = make(edgeRadius=6.0)
check("round core allowed", abs(R.ws) < 1e-14)
print("WINDING TESTS DONE")
