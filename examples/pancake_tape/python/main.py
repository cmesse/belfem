# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

import sys

from pancake import Pancake

C = Pancake()

# in mm. One 4 mm tape wound around the z axis, 200 turns from 7 to 17 mm:
# a turn pitch of 50 um. The solder between the turns is a solid that
# touches itself, and the contact surface is the tape (thin shell).
C.tapeWidth = 4
C.innerRadius = 7
C.outerRadius = 17
C.numTurns = 200

# leads. The inner one leaves the spiral at its inner end on a circle of
# radius leadRadius tangent to the spiral, wraps leadWrap degrees into the
# bore, bends hard-way in the unrolled plane of that small cylinder with
# radius riserRadius (None = 2 x tapeWidth) until it runs along z (+1) or
# -z (-1), and climbs to the cap on a solder bar. The outer one bends
# easy-way with bendRadius and leaves radially on a solder bar. The bars'
# end faces are the terminals.
C.innerLeadDirection = 1
C.leadRadius = 3
C.leadWrap = 90
C.riserRadius = None
C.innerBarThickness = 1.0
C.bendRadius = 5
C.outerBarThickness = None     # None = the wedge gap = tapeResolution

# air cylinder about z
C.domainRadius = 80
C.domainHalfHeight = 50

# resolutions in mm; the pack is one cell of the turn pitch thick, so
# the tets there have an aspect ratio of tapeResolution / pitch
C.tapeResolution = 1.0
C.leadResolution = 1.0
C.domainResolution = 8.0

# A, the current the deck ramps to
C.tapeCurrent = 100

# smoke run: python3 main.py <turns>, keeping the turn pitch
if len(sys.argv) > 1:
    pitch = (C.outerRadius - C.innerRadius) / C.numTurns
    C.numTurns = int(sys.argv[1])
    C.outerRadius = C.innerRadius + C.numTurns * pitch

C.build()
C.print()

# Relative, so both land in the directory the script is run from -- which
# is where ../input.conf looks for them (`mesh { file : pancake.msh ; }`).
C.save("pancake.msh")
C.write_topology("pancake_topology.conf")
