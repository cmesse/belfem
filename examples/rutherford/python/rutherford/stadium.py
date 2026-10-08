# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Cross-section of the flat core: a stadium in the x-y plane, width w
along x, thickness 2 r along y. Two flat faces of length ws = w - 2 r at
y = +-r are joined by semicircles of radius r centered at (+-ws/2, 0).

The arc length s runs clockwise seen from +z, starting at the top-face
center (0, r): top face in +x, right edge down, bottom face in -x, left
edge up, back to the start after one perimeter P = 2 ws + 2 pi r.

The core extruded along z is a developable surface, so a helix of
constant lay angle theta is a straight band in the unrolled (s, z)
plane with ds/dz = tan theta. Every strip node is generated in that
plane and mapped through point(); the map is piecewise analytic and
tangent-continuous, with a curvature jump from 0 to 1/r at the four
flat/round junctions.
"""

import math

import numpy as np


class Stadium:

    def __init__(self, ws: float, r: float):
        if r <= 0.0:
            raise ValueError("edge radius must be positive (got {:g})".format(r))
        if ws < 0.0:
            raise ValueError(
                "flat length ws = w - 2 r is negative ({:g}); the edge "
                "radius exceeds half the cable width".format(ws))
        self.ws = ws
        self.r = r
        self.width = ws + 2.0 * r
        self.perimeter = 2.0 * ws + 2.0 * math.pi * r
        # segment starts: top-right flat, right arc, bottom flat, left arc,
        # top-left flat
        self.s_b = 0.5 * ws
        self.s_c = 0.5 * ws + math.pi * r
        self.s_d = 1.5 * ws + math.pi * r
        self.s_e = 1.5 * ws + 2.0 * math.pi * r

    # --------------------------------------------------------------- map

    def points(self, s):
        """(x, y) of the perimeter points at arc length s (array), taken
        modulo the perimeter."""
        s = np.mod(np.asarray(s, dtype=float), self.perimeter)
        ws, r = self.ws, self.r
        x = np.empty_like(s)
        y = np.empty_like(s)

        m = s < self.s_b
        x[m] = s[m]
        y[m] = r

        m = (s >= self.s_b) & (s < self.s_c)
        phi = (s[m] - self.s_b) / r
        x[m] = 0.5 * ws + r * np.sin(phi)
        y[m] = r * np.cos(phi)

        m = (s >= self.s_c) & (s < self.s_d)
        x[m] = 0.5 * ws - (s[m] - self.s_c)
        y[m] = -r

        m = (s >= self.s_d) & (s < self.s_e)
        phi = (s[m] - self.s_d) / r
        x[m] = -0.5 * ws - r * np.sin(phi)
        y[m] = -r * np.cos(phi)

        m = s >= self.s_e
        x[m] = -0.5 * ws + (s[m] - self.s_e)
        y[m] = r

        return x, y

    def point(self, s: float):
        x, y = self.points(np.array([s]))
        return float(x[0]), float(y[0])

    def normals(self, s):
        """Outward unit normals at arc length s (array)."""
        s = np.mod(np.asarray(s, dtype=float), self.perimeter)
        r = self.r
        nx = np.zeros_like(s)
        ny = np.zeros_like(s)

        m = (s < self.s_b) | (s >= self.s_e)
        ny[m] = 1.0

        m = (s >= self.s_b) & (s < self.s_c)
        phi = (s[m] - self.s_b) / r
        nx[m] = np.sin(phi)
        ny[m] = np.cos(phi)

        m = (s >= self.s_c) & (s < self.s_d)
        ny[m] = -1.0

        m = (s >= self.s_d) & (s < self.s_e)
        phi = (s[m] - self.s_d) / r
        nx[m] = -np.sin(phi)
        ny[m] = -np.cos(phi)

        return nx, ny

    def distance(self, x, y):
        """Signed distance of (x, y) (arrays) to the perimeter: negative
        inside, zero on the curve."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        qx = np.maximum(np.abs(x) - 0.5 * self.ws, 0.0)
        return np.hypot(qx, y) - self.r
