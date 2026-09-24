# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Stadium map: closure, the four junction points, arc-length parameter-
ization, outward normals and the signed distance."""

import math
import sys

import numpy as np

sys.path.insert(0, ".")

from rutherford.stadium import Stadium


def check(label, ok):
    print("{:s}: {:s}".format("PASS" if ok else "FAIL", label))
    if not ok:
        sys.exit(1)


ws, r = 10.2, 0.9
S = Stadium(ws, r)
P = S.perimeter
check("perimeter", abs(P - (2 * ws + 2 * math.pi * r)) < 1e-14)
check("width", abs(S.width - 12.0) < 1e-14)

# landmarks
for s, exp in [(0.0, (0.0, r)),
               (ws / 2, (ws / 2, r)),
               (ws / 2 + math.pi * r / 2, (ws / 2 + r, 0.0)),
               (ws / 2 + math.pi * r, (ws / 2, -r)),
               (ws / 2 + math.pi * r + ws, (-ws / 2, -r)),
               (1.5 * ws + 1.5 * math.pi * r, (-ws / 2 - r, 0.0)),
               (1.5 * ws + 2 * math.pi * r, (-ws / 2, r)),
               (P, (0.0, r))]:
    x, y = S.point(s)
    check("landmark s={:.4g}".format(s),
          abs(x - exp[0]) < 1e-12 and abs(y - exp[1]) < 1e-12)

# arc-length parameterization: chord length -> ds for fine steps
s = np.linspace(0.0, P, 20001)
x, y = S.points(s)
seg = np.hypot(np.diff(x), np.diff(y))
ds = P / 20000
check("chord/arc within 1e-6", np.abs(seg / ds - 1.0).max() < 1e-6)

# points on the curve, normals outward and unit
d = S.distance(x, y)
check("points on the curve", np.abs(d).max() < 1e-12)
nx, ny = S.normals(s)
check("unit normals", np.abs(np.hypot(nx, ny) - 1.0).max() < 1e-14)
eps = 1e-6
check("normals point outward",
      (S.distance(x + eps * nx, y + eps * ny) > 0.5 * eps).all())
# the walk is clockwise, so the tangent is the normal rotated by -90
# degrees, (ny, -nx); compare chord tangents with midpoint normals
tx, ty = np.diff(x) / ds, np.diff(y) / ds
mx, my = S.normals(0.5 * (s[1:] + s[:-1]))
check("tangent = clockwise rotation of the normal",
      max(np.abs(tx - my).max(), np.abs(ty + mx).max()) < 1e-4)

# negative distance inside, modulo wrap
check("center inside", float(S.distance(0.0, 0.0)) < 0.0)
x1, y1 = S.point(-1.0)
x2, y2 = S.point(P - 1.0)
check("mod P wrap", abs(x1 - x2) < 1e-12 and abs(y1 - y2) < 1e-12)

# degenerate: round core
C = Stadium(0.0, 1.0)
s = np.linspace(0.0, C.perimeter, 9)
x, y = C.points(s)
check("round core is a circle", np.abs(np.hypot(x, y) - 1.0).max() < 1e-12)

try:
    Stadium(-0.1, 1.0)
    ok = False
except ValueError:
    ok = True
check("negative flat length refused", ok)
print("STADIUM TESTS DONE")
