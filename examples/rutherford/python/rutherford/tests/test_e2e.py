# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""End-to-end: the exact config of main.py, then the physical tags and
the periodic block of the written mesh are checked. Calls gmsh."""

import os
import sys
import tempfile

sys.path.insert(0, ".")

from rutherford import Cable


def check(label, ok):
    print("{:s}: {:s}".format("PASS" if ok else "FAIL", label))
    if not ok:
        sys.exit(1)


out = tempfile.mkdtemp(prefix="rutherford_e2e_")
C = Cable()
C.build()
C.print()
msh = os.path.join(out, "rutherford.msh")
conf = os.path.join(out, "rutherford_topology.conf")
C.save(msh)
C.write_topology(conf)

text = open(msh).read()
check("physical names", all(
    n in text for n in ('"Volume_1"', '"Volume_2"', '"Tape_1"', '"Tape_5"',
                        '"DomainBoundary"', '"DomainFront"', '"DomainBack"')))
check("periodic block written", "$Periodic" in text)
check("tape sidesets are 3..7", C._post.tape_physicals == [3, 4, 5, 6, 7])
check("front/back are 9/10",
      C._post.front_physical == 9 and C._post.back_physical == 10)
check("one tape position advanced", C.winding.tape_shift() == 1)
topo = open(conf).read()
check("topology names the tapes", "sidesets : 3:7 ;" in topo)
check("topology periodic triple", "source : 2, 3, 4 ;" in topo)
check("one conductor, one bracket group",
      "input curves  : [1] ;" in topo and "amplitude     : 30 A ;" in topo)
print("E2E DONE ({})".format(out))
