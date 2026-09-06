import math
import sys

sys.path.insert(0, ".")

from corc.winding import Winding, WindingError


def check(label, ok):
    print("{:s}: {:s}".format("PASS" if ok else "FAIL", label))
    if not ok:
        sys.exit(1)


# default config of main.py: 3 tapes x (4 + 1) mm on layer 0
W = Winding(tapeWidth=4, gap=1, tapeThickness=100, solderThickness=10,
            numTapesPerLayer=3, numLayers=1, pitch=10, numTurns=1,
            delta=5, tapeResolution=0.5)
check("layer0 radius = 15/2pi", abs(W.layers[0].radius - 15 / (2 * math.pi)) < 1e-14)
check("length = 2 pi pitch", abs(W.length - 20 * math.pi) < 1e-12)
check("cap rotation 0 for integer turns", abs(W.cap_rotation()) < 1e-12)
check("s partition = 72 rows", W.s_partition_count() == 72)
check("n_phi tape >= 8", W.n_phi(0, "tape") == 8)
check("n_phi gap = 2", W.n_phi(0, "gap") == 2)

# single layer, fractional turns: closure via composite rotation
W2 = Winding(4, 1, 100, 10, 3, 1, 10, 1.25, 5, 0.5)
check("single layer fractional turns allowed",
      abs(W2.cap_rotation() - 0.5 * math.pi) < 1e-12)

# two counter-wound layers, fractional turns: must refuse
try:
    Winding(4, 1, 100, 10, 3, 2, 10, 1.25, 5, 0.1).cap_rotation()
    ok = False
except WindingError:
    ok = True
check("counter-wound fractional turns refused", ok)

# two counter-wound layers, a whole number of tape positions advanced
# (numTurns * numTapesPerLayer whole): closes with zero rotation and the
# tapes permute -- the periodic twist. Argument order:
# (tapeWidth, gap, tapeThickness, solderThickness, numTapesPerLayer,
#  numLayers, pitch, numTurns, delta, tapeResolution)
W4 = Winding(4, 1, 100, 10, 3, 2, 10, 1.0 / 3.0, 3, 0.11)
check("third turn with 3 tapes closes at 0", abs(W4.cap_rotation()) < 1e-12)
check("third turn advances one tape position", W4.tape_shift() == 1)
W5 = Winding(4, 1, 100, 10, 3, 2, 10, 4.0 / 3.0, 3, 0.11)
check("four thirds advance one tape position too",
      W5.tape_shift() == 1 and abs(W5.cap_rotation()) < 1e-12)
W6 = Winding(4, 1, 100, 10, 3, 2, 10, 1.0, 3, 0.11)
check("integer turns advance zero positions", W6.tape_shift() == 0)
check("a sixth of a turn is still refused", Winding(4, 1, 100, 10, 3, 2, 10, 1.0 / 6.0, 3, 0.11).tape_shift() is None)

# two layers, half-integer turns allowed (delta small enough for the
# thin-slab s-spacing law)
W3 = Winding(4, 1, 100, 10, 3, 2, 10, 1.5, 3, 0.11)
check("half-integer closes at pi", abs(W3.cap_rotation() - math.pi) < 1e-9)

# s-spacing law: coarse delta refused for multilayer
try:
    Winding(4, 1, 100, 10, 3, 2, 10, 1.5, 10, 0.11)
    ok = False
except WindingError:
    ok = True
check("thin-slab s-spacing law enforced", ok)

# thin-gap sizing law
try:
    Winding(4, 1, 100, 10, 3, 2, 10, 1, 5, 0.5)
    ok = False
except WindingError:
    ok = True
check("interlayer resolution law enforced", ok)

# zero gap refused (overlap on layer 0 is impossible by construction:
# the packing radius guarantees tape_angle < 2 pi / n whenever gap > 0)
try:
    Winding(4, 0.0, 100, 10, 3, 1, 10, 1, 5, 0.5)
    ok = False
except WindingError:
    ok = True
check("zero gap refused", ok)

print("ALL WINDING TESTS PASSED")
