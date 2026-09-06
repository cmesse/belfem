"""Planar spiral windings built as offset curves of a closed convex base curve.

A tape wound on a mandrel forms turn k+1 by laying it onto turn k, i.e. every
turn is the *parallel (offset) curve* of the previous one at the distance of
one winding pitch p (tape thickness plus insulation).  The continuous version
of this is the spiral

    r(theta) = c(theta) + delta(theta) nu(theta),   delta = p (theta - theta0) / (2 pi)

where c(theta) is the base curve (traversed counter-clockwise for increasing
theta) and nu(theta) its outward unit normal.  Because parallel curves share
their normals, the normal distance between consecutive turns is p up to a
relative error of order (p / (2 pi R))^2.  Outward offsets of convex curves
never develop cusps, so the construction is regular for any number of turns.

All derivatives with respect to theta up to third order are generated
symbolically with sympy; the curve is then reparametrized by arc length s so
that v = T, a = dT/ds and b = d2T/ds2 (unit speed, Russenschuck Sec. 3.1).

The strip frame of the planar spiral is

    T = v,   b = e_z,   n = b x T          (n points inward for CCW winding)

which coincides with the Frenet frame for counter-clockwise winding and gives
kappa_n = kappa (signed), kappa_g = 0, tau = 0.  Unlike the Frenet frame it
stays defined and continuous when the curvature changes sign or vanishes.
"""

import math
import numpy as np
import sympy as sp

from .basecurve import Basecurve, _GL_X, _GL_W

_EZ = np.array([0.0, 0.0, 1.0])


class PlanarSpiral(Basecurve):
    """Offset spiral r = c + delta nu of a symbolic base curve c(theta) = (xexpr, yexpr)."""

    def __init__(self, theta: sp.Symbol, xexpr, yexpr, pitch: float, nturns: float,
                 theta0: float = 0.0, sense: int = +1, nodes_per_turn: int = 64):
        Basecurve.__init__(self)
        if pitch <= 0.0:
            raise ValueError("pitch must be positive")
        if nturns <= 0.0:
            raise ValueError("nturns must be positive")
        if sense not in (+1, -1):
            raise ValueError("sense must be +1 (counter-clockwise) or -1 (clockwise)")

        self.pitch = float(pitch)
        self.nturns = float(nturns)
        self.theta0 = float(theta0)
        self.theta1 = self.theta0 + 2.0 * math.pi * self.nturns
        self.sense = int(sense)
        self._theta_sym = theta

        # --- symbolic definition ------------------------------------------
        dx = sp.diff(xexpr, theta)
        dy = sp.diff(yexpr, theta)
        nrm = sp.sqrt(dx ** 2 + dy ** 2)
        delta = sp.Float(self.pitch) * (theta - sp.Float(self.theta0)) / (2 * sp.pi)
        X = xexpr + delta * dy / nrm      # outward normal nu = ( dy, -dx ) / |c'|
        Y = yexpr - delta * dx / nrm

        exprs = [(X, Y)]
        for _ in range(3):
            Xp, Yp = exprs[-1]
            exprs.append((sp.diff(Xp, theta), sp.diff(Yp, theta)))
        self._fun = [sp.lambdify(theta, [ex, ey], modules="numpy", cse=True) for ex, ey in exprs]
        self._base_fun = sp.lambdify(theta, [xexpr, yexpr], modules="numpy")
        self._normal_fun = sp.lambdify(theta, [dy / nrm, -dx / nrm], modules="numpy")

        # --- arc-length table -------------------------------------------------
        m = max(8, int(math.ceil(self.nturns * nodes_per_turn)))
        self._theta_nodes = np.linspace(self.theta0, self.theta1, m + 1)
        self._s_nodes = np.zeros(m + 1)
        for k in range(m):
            self._s_nodes[k + 1] = self._s_nodes[k] + self._gl_length(self._theta_nodes[k], self._theta_nodes[k + 1])
        self.tmin = 0.0
        self.tmax = float(self._s_nodes[-1])

        self._cache_s = None
        self._cache_val = None

    # ------------------------------------------------------------------
    # theta domain
    # ------------------------------------------------------------------
    def derivs_theta(self, theta: float) -> np.ndarray:
        """Array (4, 3): position and theta-derivatives 1..3 of the spiral at theta."""
        out = np.zeros((4, 3))
        for k in range(4):
            x, y = self._fun[k](theta)
            if not (np.isfinite(x) and np.isfinite(y)):
                # non-even exponents can be singular exactly on the axes; nudge
                eps = 1e-9
                x0, y0 = self._fun[k](theta - eps)
                x1, y1 = self._fun[k](theta + eps)
                x, y = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
            out[k, 0] = x
            out[k, 1] = self.sense * y
        return out

    def base_point(self, theta: float) -> np.ndarray:
        x, y = self._base_fun(theta)
        return np.array([float(x), self.sense * float(y), 0.0])

    def outward_normal(self, theta: float) -> np.ndarray:
        x, y = self._normal_fun(theta)
        return np.array([float(x), self.sense * float(y), 0.0])

    def point_theta(self, theta: float) -> np.ndarray:
        return self.derivs_theta(theta)[0]

    def speed_theta(self, theta: float) -> float:
        d = self.derivs_theta(theta)
        return float(np.linalg.norm(d[1]))

    def _gl_length(self, ta: float, tb: float) -> float:
        hm = 0.5 * (tb - ta)
        tm = 0.5 * (tb + ta)
        acc = 0.0
        for x, w in zip(_GL_X, _GL_W):
            acc += w * self.speed_theta(tm + hm * x)
        return hm * acc

    def s_of_theta(self, theta: float) -> float:
        """Arc length from the start of the spiral to the point at theta."""
        k = int(np.searchsorted(self._theta_nodes, theta, side="right") - 1)
        k = min(max(k, 0), len(self._theta_nodes) - 2)
        return float(self._s_nodes[k] + self._gl_length(self._theta_nodes[k], theta))

    def theta_of_s(self, s: float) -> float:
        """Inverse of s_of_theta by Newton iteration on the arc-length table."""
        if s <= 0.0:
            return self.theta0 + s / self.speed_theta(self.theta0)
        if s >= self.tmax:
            return self.theta1 + (s - self.tmax) / self.speed_theta(self.theta1)
        k = int(np.searchsorted(self._s_nodes, s, side="right") - 1)
        k = min(max(k, 0), len(self._s_nodes) - 2)
        t0, t1 = self._theta_nodes[k], self._theta_nodes[k + 1]
        s0, s1 = self._s_nodes[k], self._s_nodes[k + 1]
        th = t0 + (s - s0) / (s1 - s0) * (t1 - t0)
        for _ in range(50):
            f = s0 + self._gl_length(t0, th) - s
            dth = f / self.speed_theta(th)
            th -= dth
            if abs(dth) < 1e-14 * max(1.0, abs(th)):
                break
        return float(th)

    # ------------------------------------------------------------------
    # arc-length parametrization
    # ------------------------------------------------------------------
    def _eval(self, s: float):
        """(theta, r, T, dT/ds, d2T/ds2) at arc length s (cached for the last s)."""
        if self._cache_s is not None and s == self._cache_s:
            return self._cache_val
        th = self.theta_of_s(s)
        d = self.derivs_theta(th)
        R1, R2, R3 = d[1], d[2], d[3]
        sig = np.linalg.norm(R1)
        sig1 = np.dot(R1, R2) / sig
        sig2 = (np.dot(R2, R2) + np.dot(R1, R3)) / sig - sig1 ** 2 / sig
        T = R1 / sig
        dT = R2 / sig - R1 * sig1 / sig ** 2
        ddT = R3 / sig - 2.0 * R2 * sig1 / sig ** 2 - R1 * sig2 / sig ** 2 + 2.0 * R1 * sig1 ** 2 / sig ** 3
        a = dT / sig
        b = (ddT / sig - dT * sig1 / sig ** 2) / sig
        val = (th, d[0].copy(), T, a, b)
        self._cache_s = s
        self._cache_val = val
        return val

    def theta(self, s: float) -> float:
        return self._eval(s)[0]

    def r(self, s: float) -> np.ndarray:
        return self._eval(s)[1]

    def v(self, s: float) -> np.ndarray:
        return self._eval(s)[2]

    def a(self, s: float) -> np.ndarray:
        return self._eval(s)[3]

    def b(self, s: float) -> np.ndarray:
        return self._eval(s)[4]

    def speed(self, s: float) -> float:
        return 1.0

    def segment_length(self, ta: float, tb: float, nsub: int = 1) -> float:
        return tb - ta

    def length(self) -> float:
        return self.tmax - self.tmin

    def equidistant(self, n: int) -> np.ndarray:
        return np.linspace(self.tmin, self.tmax, n)

    # ------------------------------------------------------------------
    # strip frame
    # ------------------------------------------------------------------
    def frame(self, s: float):
        T = self.v(s)
        b = _EZ
        n = np.cross(b, T)
        return T, n, b

    def curvatures(self, s: float):
        T, n, b = self.frame(s)
        return 0.0, 0.0, float(np.dot(n, self.a(s)))

    def signed_curvature(self, s: float) -> float:
        return self.curvatures(s)[2]

    # ------------------------------------------------------------------
    # diagnostics
    # ------------------------------------------------------------------
    def turn_gap(self, theta: float) -> float:
        """Normal distance between the turn at theta and the next turn (theta + 2 pi).

        Measured along the outward normal of the base curve; equals the pitch up to
        second order in pitch / (2 pi R).
        """
        p0 = self.point_theta(theta)
        p1 = self.point_theta(theta + 2.0 * math.pi)
        nu = self.outward_normal(theta)
        return float(np.dot(p1 - p0, nu))


class LameSpiral(PlanarSpiral):
    """Spiral winding on a Lame curve (superellipse)  |x/a|^n + |y/b|^n = 1.

    Parameters
    ----------
    a, b : semi-axes of the innermost turn (mandrel) along x and y
    n    : Lame exponent (2 = ellipse; 4 ... 6 = rounded rectangle).  Even integers
           give an analytic curve; other exponents are only C^(ceil(n)-1) on the axes.
    pitch : radial growth per turn (tape thickness plus insulation)
    nturns : number of turns (may be fractional)
    theta0 : polar angle where the winding starts
    sense : +1 winds counter-clockwise (seen from +z), -1 clockwise
    """

    def __init__(self, a: float, b: float, n: float, pitch: float, nturns: float,
                 theta0: float = 0.0, sense: int = +1, nodes_per_turn: int = 64):
        if a <= 0.0 or b <= 0.0:
            raise ValueError("semi-axes must be positive")
        if n < 2.0:
            raise ValueError("Lame exponent must be >= 2 (convex curve)")
        self.a_axis = float(a)
        self.b_axis = float(b)
        self.n_exp = float(n)

        theta = sp.Symbol("theta")  # not declared real: keeps (cos^2)^(n/2) free of Abs/sign terms
        if float(n).is_integer() and int(n) % 2 == 0:
            nn = sp.Integer(int(n))
            cn = sp.cos(theta) ** nn
            sn = sp.sin(theta) ** nn
        else:
            nn = sp.Float(n)
            cn = (sp.cos(theta) ** 2) ** (nn / 2)
            sn = (sp.sin(theta) ** 2) ** (nn / 2)
        rho = (cn / sp.Float(a) ** nn + sn / sp.Float(b) ** nn) ** (-1 / nn)
        PlanarSpiral.__init__(self, theta, rho * sp.cos(theta), rho * sp.sin(theta),
                              pitch, nturns, theta0, sense, nodes_per_turn)

    def __repr__(self):
        return ("LameSpiral(a=%g, b=%g, n=%g, pitch=%g, nturns=%g, theta0=%g, sense=%+d)"
                % (self.a_axis, self.b_axis, self.n_exp, self.pitch, self.nturns, self.theta0, self.sense))
