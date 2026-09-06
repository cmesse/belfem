"""Analytic centerline curves for the CORC cable generator.

A Curve provides the position and its first three derivatives at any
parameter t (via Jet arithmetic — see jet.py), an exact-to-quadrature
arc length s(t), and the inverse t(s). All lengths in mm.
"""

import math

import numpy as np

from corc import jet
from corc.jet import Jet


class Curve:

    #: closed curves have p(tmin) == p(tmax) and no end caps
    closed = False

    def __init__(self, tmin: float, tmax: float, nsamples=4096):
        self.tmin = tmin
        self.tmax = tmax
        self._build_arclength_table(nsamples)

    # ------------------------------------------------------------- geometry

    def position(self, t):
        """Return [x, y, z] at parameter t; components are Jets if t is a
        Jet. Subclasses implement this with jet-friendly operations."""
        raise NotImplementedError

    def p(self, t: float) -> np.ndarray:
        return np.array([c.value if isinstance(c, Jet) else float(c)
                         for c in self.position(Jet.variable(t))])

    def derivatives(self, t: float):
        """Return (p, v, a, j) as 3-vectors at parameter t."""
        comps = self.position(Jet.variable(t))
        comps = [c if isinstance(c, Jet) else Jet(c) for c in comps]
        p = np.array([c.value for c in comps])
        v = np.array([c.d1 for c in comps])
        a = np.array([c.d2 for c in comps])
        b = np.array([c.d3 for c in comps])
        return p, v, a, b

    def speed(self, t: float) -> float:
        _, v, _, _ = self.derivatives(t)
        return float(np.linalg.norm(v))

    # ----------------------------------------------------------- arc length

    def _build_arclength_table(self, nsamples: int):
        # composite Gauss-Legendre (5 point) per interval: exact enough that
        # the table itself limits nothing (error ~ h^10)
        xg, wg = np.polynomial.legendre.leggauss(5)
        tgrid = np.linspace(self.tmin, self.tmax, nsamples + 1)
        seg = np.zeros(nsamples)
        for k in range(nsamples):
            t0, t1 = tgrid[k], tgrid[k + 1]
            hm = 0.5 * (t1 - t0)
            mid = 0.5 * (t0 + t1)
            seg[k] = hm * sum(w * self.speed(mid + hm * x) for x, w in zip(xg, wg))
        self._tgrid = tgrid
        self._sgrid = np.concatenate([[0.0], np.cumsum(seg)])
        self.length = float(self._sgrid[-1])

    def s_of_t(self, t: float) -> float:
        """Arc length from tmin to t (table + local Gauss correction)."""
        t = float(t)
        k = int(np.clip(np.searchsorted(self._tgrid, t) - 1, 0, len(self._tgrid) - 2))
        t0 = self._tgrid[k]
        xg, wg = np.polynomial.legendre.leggauss(5)
        hm = 0.5 * (t - t0)
        mid = 0.5 * (t + t0)
        loc = hm * sum(w * self.speed(mid + hm * x) for x, w in zip(xg, wg))
        return float(self._sgrid[k] + loc)

    def t_of_s(self, s: float) -> float:
        """Invert the arc length: Newton on s(t) - s with table start."""
        s = float(np.clip(s, 0.0, self.length))
        k = int(np.clip(np.searchsorted(self._sgrid, s) - 1, 0, len(self._sgrid) - 2))
        # linear initial guess inside the bracketing interval
        s0, s1 = self._sgrid[k], self._sgrid[k + 1]
        t0, t1 = self._tgrid[k], self._tgrid[k + 1]
        t = t0 if s1 == s0 else t0 + (s - s0) / (s1 - s0) * (t1 - t0)
        for _ in range(60):
            f = self.s_of_t(t) - s
            if abs(f) < 1e-13 * max(1.0, self.length):
                break
            sp = self.speed(t)
            step = f / sp
            # keep Newton inside the bracket for robustness
            tn = t - step
            if tn < t0 or tn > t1:
                if f > 0:
                    t1 = t
                else:
                    t0 = t
                tn = 0.5 * (t0 + t1)
            t = tn
        return float(t)


class StraightLine(Curve):
    """Straight centerline along +z starting at the origin. The length is
    set by the cable (pitch * numTurns * 2*pi)."""

    def __init__(self, length: float = 1.0):
        Curve.__init__(self, 0.0, float(length), nsamples=8)

    def set_length(self, length: float):
        Curve.__init__(self, 0.0, float(length), nsamples=8)

    def position(self, t):
        return [Jet(0.0), Jet(0.0), t]


class Lame(Curve):
    """Quarter superellipse (Lame curve) wrapped on a cylinder barrel of
    given radius — the coil-end path from the legacy code.

    x = R cos(alpha), y = R sin(alpha), z = b sin(t)^(2/order)
    alpha = a cos(t)^(2/order) / (0.5 pi R),  t in (0, pi/2)
    a = cos(phi_deg) * 0.5 pi R,  b = a * ratio

    For order > 2 the parametrization is singular AT the endpoints, so the
    domain is clipped by a tiny margin; order = 2 is regular everywhere.
    """

    def __init__(self, radius: float, phi: float, ratio: float, order=2):
        self.radius = float(radius)
        self.a = math.cos(math.radians(phi)) * 0.5 * math.pi * radius
        self.b = self.a * ratio
        self.order = float(order)
        eps = 0.0 if order == 2 else 1e-6
        Curve.__init__(self, eps, 0.5 * math.pi - eps)

    def position(self, t):
        m = 2.0 / self.order
        if self.order == 2:
            u = jet.cos(t)
            w = jet.sin(t)
        else:
            u = jet.power(jet.cos(t), m)
            w = jet.power(jet.sin(t), m)
        alpha = u * (self.a / (0.5 * math.pi * self.radius))
        x = self.radius * jet.cos(alpha)
        y = self.radius * jet.sin(alpha)
        z = self.b * w
        return [x, y, z]


class Trefoil(Curve):
    """Trefoil knot (closed, show-off only)."""

    closed = True

    def __init__(self, radius: float):
        self.radius = float(radius)
        Curve.__init__(self, 0.0, 2.0 * math.pi)

    def position(self, t):
        r = self.radius
        x = r * (jet.sin(t) + 2.0 * jet.sin(2.0 * t))
        y = r * (jet.cos(t) - 2.0 * jet.cos(2.0 * t))
        z = -r * jet.sin(3.0 * t)
        return [x, y, z]
