"""End-to-end: bent cable along the Lame coil-end curve (single layer,
fractional turns -> rotated cap periodicity)."""

import sys

sys.path.insert(0, ".")

from corc import Cable, Lame

C = Cable()
C.tapeWidth = 4
C.pitch = 10
C.gap = 1
C.tapeThickness = 100
C.solderThickness = 10
C.domainRadius = 20          # kappa_max = 0.0377/mm allows < 26.5 mm
C.numTapesPerLayer = 3
C.numLayers = 1
C.delta = 5
C.centerline = Lame(76.1, 25, 1.55, order=2)
C.workdir = "/tmp"

C.build()
C.print()
C.save("/tmp/corc_lame.msh")
C.write_topology("/tmp/corc_lame_topology.conf")
print("LAME E2E DONE")
