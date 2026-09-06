"""Unit checks for corc.curve and corc.frame (plain asserts, exit 1 on
failure): jet arithmetic, curvature of known curves, arc-length round
trips, frame orthonormality/continuity, straight-line frame identity."""

import math
import sys

import numpy as np

sys.path.insert(0, ".")

from corc import jet
from corc.jet import Jet
from corc.curve import StraightLine, Lame, Trefoil, Curve
from corc.frame import Frame, FrameError


def check(label, ok):
    print("{:s}: {:s}".format("PASS" if ok else "FAIL", label))
    if not ok:
        sys.exit(1)


# ------------------------------------------------------------------ jets
t = Jet.variable(0.7)
f = jet.sin(3.0 * t) * jet.cos(t) + jet.power(t, 2.5)
# reference values via analytic formulas
x = 0.7
val = math.sin(3 * x) * math.cos(x) + x ** 2.5
d1 = 3 * math.cos(3 * x) * math.cos(x) - math.sin(3 * x) * math.sin(x) + 2.5 * x ** 1.5
d2 = (-9 * math.sin(3 * x) * math.cos(x) - 3 * math.cos(3 * x) * math.sin(x)
      - 3 * math.cos(3 * x) * math.sin(x) - math.sin(3 * x) * math.cos(x)
      + 2.5 * 1.5 * x ** 0.5)
check("jet value", abs(f.value - val) < 1e-14)
check("jet d1", abs(f.d1 - d1) < 1e-13)
check("jet d2", abs(f.d2 - d2) < 1e-12)


# ------------------------------------------------- circle curvature via jets
class _Circle(Curve):
    def __init__(self, R):
        self.R = R
        Curve.__init__(self, 0.0, 2 * math.pi)

    def position(self, tt):
        return [self.R * jet.cos(tt), self.R * jet.sin(tt), Jet(0.0)]


C = _Circle(7.5)
F = Frame(C, nsamples=256)
check("circle kappa = 1/R", abs(F.kappa_max - 1.0 / 7.5) < 1e-10)
check("circle length = 2 pi R", abs(C.length - 2 * math.pi * 7.5) < 1e-9)

# arc length round trip
for s in (0.0, 1.0, 11.7, 2 * math.pi * 7.5):
    check("circle t(s(t)) roundtrip s={:g}".format(s),
          abs(C.s_of_t(C.t_of_s(s)) - s) < 1e-10)

# ------------------------------------------------------------ straight line
S = StraightLine(100.0)
FS = Frame(S, nsamples=64)
R0 = FS.R(37.0)
check("straight frame is identity", np.allclose(R0, np.eye(3), atol=1e-14))
p = FS.transform(12.5, 3.0, -4.0)
check("straight transform = (x, y, s)", np.allclose(p, [3.0, -4.0, 12.5], atol=1e-12))

# ------------------------------------------------------------------- lame
L = Lame(76.1, 25, 1.55, order=2)
check("lame length positive", L.length > 100)
for s in (0.0, 0.3 * L.length, L.length):
    tt = L.t_of_s(s)
    check("lame roundtrip s={:g}".format(s), abs(L.s_of_t(tt) - s) < 1e-9 * L.length)

FL = Frame(L, nsamples=512)
# orthonormality + continuity along the curve
prev_n = None
worst_orth = 0.0
worst_step = 0.0
for s in np.linspace(0, L.length, 200):
    R = FL.R(float(s))
    worst_orth = max(worst_orth, float(np.abs(R.T @ R - np.eye(3)).max()))
    det = np.linalg.det(R)
    assert det > 0.99, "left-handed frame at s={:g}".format(s)
    n = R[:, 0]
    if prev_n is not None:
        worst_step = max(worst_step, float(np.arccos(np.clip(np.dot(prev_n, n), -1, 1))))
    prev_n = n
check("lame frame orthonormal (worst {:.1e})".format(worst_orth), worst_orth < 1e-12)
check("lame frame continuous (worst step {:.3f} rad)".format(worst_step), worst_step < 0.2)
print("lame kappa_max = {:.5f} 1/mm -> max domainRadius = {:.2f} mm".format(
    FL.kappa_max, 1 / FL.kappa_max))

# admissibility: generous radius passes, huge radius fails
try:
    FL.check_domain_admissibility(1.0 / FL.kappa_max * 0.5)
    ok = True
except FrameError:
    ok = False
check("lame admissibility r < 1/kappa", ok)
try:
    FL.check_domain_admissibility(1.0 / FL.kappa_max * 1.5)
    ok = False
except FrameError:
    ok = True
check("lame admissibility rejects folding radius", ok)

# ---------------------------------------------------------------- trefoil
T = Trefoil(8.718490620126541)
check("trefoil closed", T.closed)
p0, pL = T.p(T.tmin), T.p(T.tmax)
check("trefoil endpoints coincide", np.linalg.norm(p0 - pL) < 1e-9)
FT = Frame(T, nsamples=1024)
check("trefoil kappa positive", FT.kappa_max > 0)

# jerk correctness on the trefoil vs the (fixed) analytic formulas
tt = 1.234
_, v, a, j3 = T.derivatives(tt)
r = T.radius
va = np.array([r * (math.cos(tt) + 4 * math.cos(2 * tt)),
               r * (-math.sin(tt) + 4 * math.sin(2 * tt)),
               -3 * r * math.cos(3 * tt)])
ja = np.array([r * (-math.cos(tt) - 16 * math.cos(2 * tt)),
               r * (math.sin(tt) - 16 * math.sin(2 * tt)),
               27 * r * math.cos(3 * tt)])
check("trefoil v matches analytic", np.allclose(v, va, atol=1e-10))
check("trefoil jerk matches analytic", np.allclose(j3, ja, atol=1e-9))

print("ALL CURVE/FRAME TESTS PASSED")
