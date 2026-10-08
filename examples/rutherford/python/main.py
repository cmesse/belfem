# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

from rutherford import Cable

C = Cable()

# in mm. The reference study (Romancew 2026) wraps 2, 4 and 6 mm tapes at
# lay angles of 30 and 45 degrees, measured from the cable axis.
C.tapeWidth = 4
C.layAngle = 30

# The core is a flat Rutherford-type copper cable with rounded edges; it
# is not meshed as a conductor here, its section is air. The edge radius
# is half the core thickness. The core width is not documented in the
# reference study: 12 mm is a placeholder.
C.edgeRadius = 0.9
C.cableWidth = 12

# Five tapes on one layer: the most 4 mm tapes that fit on this core at
# 30 degrees (a tape takes tapeWidth / cos theta of perimeter), 89 % of
# the surface covered. The reference study wraps a single tape
# (numTapes = 1). Tapes sit at equal perimeter spacing; odd layers are
# counter-wound (criss-cross) at the same pitch, so their lay angle
# differs slightly from layAngle. numPitches * numTapes must be whole so
# the end planes are translated images of each other: a fifth of a pitch
# closes with every tape continuing into its neighbor.
C.numTapes = 5
C.numLayers = 1
C.layerGap = 200          # um, only used with numLayers > 1
C.numPitches = 1.0 / 5.0

# in mm
C.domainRadius = 50

# resolutions in mm. The tape strip is built from rows at constant z
# spaced so the step along the tape equals tapeResolution; the bare core
# between the tape turns is meshed at gapResolution.
C.tapeResolution = 0.4
C.gapResolution = 0.8
C.innerResolution = 0.5
C.domainResolution = 4

C.build()
C.print()

# Relative, so both land in the directory the script is run from -- which
# is where ../input.conf looks for them (`mesh { file : rutherford.msh ; }`).
# Only the gmsh scratch files stay in Cable.workdir.
C.save("rutherford.msh")
C.write_topology("rutherford_topology.conf")
