"""End-to-end: straight cable with the exact config of main.py."""

import sys

sys.path.insert(0, ".")

from corc import Cable, StraightLine

C = Cable()
C.tapeWidth = 4
C.pitch = 10
C.numTurns = 1
C.gap = 1
C.tapeThickness = 100
C.solderThickness = 10
C.domainRadius = 50
C.numTapesPerLayer = 3
C.numLayers = 1
C.delta = 5
C.centerline = StraightLine()
C.workdir = "/tmp"

C.build()
C.print()
C.save("/tmp/corc_new.msh")
C.write_topology("/tmp/corc_topology.conf")
print("STRAIGHT E2E DONE")
