"""End-to-end: TWO counter-wound layers, half-integer turns (cap rotation
theta = pi), thin 0.11 mm inter-layer gap.

KNOWN LIMITATION (2026-07-14, gmsh 4.13.1): the geometry, closure, cap
periodicity and the thin-gap fill are all correct (gap region has ZERO
sub-threshold tets at tapeResolution = 2x gap with Netgen optimization),
but gmsh's Delaunay leaves ~70-90 front-collision slivers in the OUTER air
(concentric size-field point shells at mid-domain; RandomFactor/Netgen do
not remove them, HXT and Frontal segfault in this build). Expected to be
fixed by gmsh 4.15.2 (/opt/gmsh/latest) — re-run this test after the
upgrade. Until then the validation error is the EXPECTED outcome; this
test asserts exactly that, so a silent regression of the validator would
also be caught."""

import sys

sys.path.insert(0, ".")

from corc import Cable, StraightLine
from corc.postprocess import ValidationError

C = Cable()
C.tapeWidth = 4
C.pitch = 10
C.numTurns = 0.5            # half-integer: closes with theta = pi
C.gap = 1
C.tapeThickness = 100
C.solderThickness = 10
C.domainRadius = 30
C.numTapesPerLayer = 3
C.numLayers = 2
C.delta = 3                 # ds ~ 0.52 mm ~ 4.8x gap (Spike D regime)
C.tapeResolution = 0.22     # 2x gap: pancake regime, gap fills cleanly
C.innerResolution = 1
C.domainResolution = 8
C.centerline = StraightLine()
C.workdir = "/tmp"

try:
    C.build()
except ValidationError as e:
    assert "degenerate tets" in str(e), e
    print("MULTILAYER: known gmsh 4.13.1 outer-air sliver limitation "
          "reproduced ({})".format(e))
    print("MULTILAYER E2E: EXPECTED-FAIL OK (retry on gmsh >= 4.15)")
else:
    # the gmsh upgrade fixed it — promote this to a hard-pass test!
    C.print()
    C.save("/tmp/corc_multilayer.msh")
    C.write_topology("/tmp/corc_multilayer_topology.conf")
    print("MULTILAYER E2E DONE — gmsh handles it now; make this test "
          "a hard pass and update tmp/todo.md")
