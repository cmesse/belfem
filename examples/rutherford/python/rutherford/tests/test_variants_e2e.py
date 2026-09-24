# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""End-to-end variants: the 45 degree lay with four tapes, the single
tape of the reference study over a full pitch, and two criss-crossed
layers. Each build runs gmsh; pass a variant name to run one of them."""

import os
import sys
import tempfile

sys.path.insert(0, ".")

from rutherford import Cable


def check(label, ok):
    print("{:s}: {:s}".format("PASS" if ok else "FAIL", label))
    if not ok:
        sys.exit(1)


def run(name):
    C = Cable()
    if name == "lay45":
        C.layAngle = 45
        C.numTapes = 4
        C.numPitches = 1.0 / 4.0
    elif name == "onetape":
        C.numTapes = 1
        C.numPitches = 1
    elif name == "twolayers":
        C.numLayers = 2
        C.tapeResolution = 0.2
        C.gapResolution = 0.4
        C.innerResolution = 0.4
    else:
        raise SystemExit("unknown variant " + name)
    C.build()
    C.print()
    out = tempfile.mkdtemp(prefix="rutherford_{}_".format(name))
    C.save(os.path.join(out, "rutherford.msh"))
    C.write_topology(os.path.join(out, "rutherford_topology.conf"))
    if name == "lay45":
        check("four tape sidesets", C._post.tape_physicals == [3, 4, 5, 6])
    if name == "onetape":
        check("one tape sideset", C._post.tape_physicals == [3])
        check("no tape position advanced", C.winding.tape_shift() == 0)
    if name == "twolayers":
        check("three volumes", C._post.volume_count == 3)
        check("ten tape sidesets", C._post.tape_physicals == list(range(4, 14)))
        check("layer 1 counter-wound", C.winding.layers[1].sigma == -1)
    print("VARIANT {} DONE".format(name))


names = sys.argv[1:] or ["lay45", "onetape", "twolayers"]
for n in names:
    run(n)
