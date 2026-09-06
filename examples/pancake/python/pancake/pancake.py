"""Composite base curve of a pancake coil: inner lead + planar spiral + outer lead.

The whole curve is parametrized by its arc length s in [0, L]:

    [0, L_in]                inner lead (from the inner connector into the winding)
    [L_in, L_in + L_sp]      planar spiral
    [L_in + L_sp, L]         outer lead (from the winding to the outer connector)

Both leads are *described* in the direction that leaves the coil (that is the
natural way to think about a terminal: "go straight, twist by 90 degrees, bend
up, go straight").  For the inner lead this description is integrated
backwards from the start of the spiral with the leaving frame

    T_L = -T(0),  n_L = -n(0),  b_L = b(0)

so that the tape width b still points along +z when the description starts.
The resulting piece is then reversed to obtain the forward parametrization.
Reversal flips the sign of the odd derivatives (v, jerk), of T and n and of
kappa_n, and keeps a, b, kappa_g and tau (see the docstring of `_inner`).
"""

import numpy as np

from .basecurve import Basecurve
from .spiral import PlanarSpiral
from .lead import Lead, Release


def _make_lead(spec, release: float, kappa_n_junction: float) -> Lead:
    if isinstance(spec, Lead):
        segments, kw = list(spec.segments), dict(max_step_angle=spec.max_step_angle, min_steps=spec.min_steps)
    else:
        segments, kw = list(spec), {}
    if release > 0.0:
        segments.insert(0, Release(release, kappa_n=kappa_n_junction))
    return Lead(segments, **kw)


class PancakeCurve(Basecurve):
    """
    Parameters
    ----------
    spiral : PlanarSpiral
    inner_lead, outer_lead : list of Segment (or a Lead), described in the direction
        that leaves the coil; None for a free end
    release : length over which the curvature of the winding is faded to zero
        before the user segments start (0 = tape leaves tangentially with a
        curvature jump; the curve is then only C^1 at the junction)
    """

    def __init__(self, spiral: PlanarSpiral, inner_lead=None, outer_lead=None, release: float = 0.0):
        Basecurve.__init__(self)
        self.spiral = spiral
        self.inner = None
        self.outer = None
        self.release = float(release)

        if inner_lead is not None:
            r0 = spiral.r(spiral.tmin)
            T0, n0, b0 = spiral.frame(spiral.tmin)
            # leaving frame (-T, -n, b): kappa_n changes sign, see `_inner`
            inner = _make_lead(inner_lead, self.release, -spiral.curvatures(spiral.tmin)[2])
            inner.build(r0, -T0, -n0, b0)
            self.inner = inner

        if outer_lead is not None:
            r1 = spiral.r(spiral.tmax)
            T1, n1, b1 = spiral.frame(spiral.tmax)
            outer = _make_lead(outer_lead, self.release, spiral.curvatures(spiral.tmax)[2])
            outer.build(r1, T1, n1, b1)
            self.outer = outer

        self.L_inner = self.inner.length if self.inner is not None else 0.0
        self.L_spiral = spiral.tmax - spiral.tmin
        self.L_outer = self.outer.length if self.outer is not None else 0.0

        self.s_spiral0 = self.L_inner
        self.s_spiral1 = self.L_inner + self.L_spiral
        self.tmin = 0.0
        self.tmax = self.L_inner + self.L_spiral + self.L_outer

        bp = []
        if self.inner is not None:
            bp.extend((self.L_inner - self.inner.breakpoints)[::-1][:-1])
        bp.append(self.s_spiral0)
        bp.append(self.s_spiral1)
        if self.outer is not None:
            bp.extend(self.s_spiral1 + self.outer.breakpoints[1:])
        self.breakpoints = np.array(bp)

    # ------------------------------------------------------------------
    # bookkeeping
    # ------------------------------------------------------------------
    def piece(self, s: float) -> str:
        if s < self.s_spiral0:
            return "inner"
        if s > self.s_spiral1:
            return "outer"
        return "spiral"

    def _inner(self, s: float):
        """Forward quantities of the reversed inner lead at global s.

        With s_L = L_in - s and the leaving-direction lead quantities (r, v_L, a_L, j_L),
        frame (T_L, n_L, b_L) and curvatures (tau_L, kg_L, kn_L):
            r = r,  v = -v_L,  a = a_L,  jerk = -j_L
            T = -T_L,  n = -n_L,  b = b_L        (right handed)
            tau = tau_L,  kappa_g = kg_L,  kappa_n = -kn_L
        """
        sL = self.L_inner - s
        r, vL, aL, jL = self.inner.derivatives(sL)
        _, TL, nL, bL = self.inner.state(sL)
        tauL, kgL, knL = self.inner.omega(sL)
        return r, -vL, aL, -jL, (-TL, -nL, bL), (tauL, kgL, -knL)

    # ------------------------------------------------------------------
    # derivatives
    # ------------------------------------------------------------------
    def r(self, s: float) -> np.ndarray:
        p = self.piece(s)
        if p == "spiral":
            return self.spiral.r(s - self.s_spiral0)
        if p == "outer":
            return self.outer.state(s - self.s_spiral1)[0]
        return self._inner(s)[0]

    def v(self, s: float) -> np.ndarray:
        p = self.piece(s)
        if p == "spiral":
            return self.spiral.v(s - self.s_spiral0)
        if p == "outer":
            return self.outer.derivatives(s - self.s_spiral1)[1]
        return self._inner(s)[1]

    def a(self, s: float) -> np.ndarray:
        p = self.piece(s)
        if p == "spiral":
            return self.spiral.a(s - self.s_spiral0)
        if p == "outer":
            return self.outer.derivatives(s - self.s_spiral1)[2]
        return self._inner(s)[2]

    def b(self, s: float) -> np.ndarray:
        p = self.piece(s)
        if p == "spiral":
            return self.spiral.b(s - self.s_spiral0)
        if p == "outer":
            return self.outer.derivatives(s - self.s_spiral1)[3]
        return self._inner(s)[3]

    # ------------------------------------------------------------------
    # frame
    # ------------------------------------------------------------------
    def frame(self, s: float):
        p = self.piece(s)
        if p == "spiral":
            return self.spiral.frame(s - self.s_spiral0)
        if p == "outer":
            _, T, n, b = self.outer.state(s - self.s_spiral1)
            return T, n, b
        return self._inner(s)[4]

    def curvatures(self, s: float):
        p = self.piece(s)
        if p == "spiral":
            return self.spiral.curvatures(s - self.s_spiral0)
        if p == "outer":
            tau, kg, kn = self.outer.omega(s - self.s_spiral1)
            return float(tau), float(kg), float(kn)
        tau, kg, kn = self._inner(s)[5]
        return float(tau), float(kg), float(kn)

    # ------------------------------------------------------------------
    # arc length is the parameter
    # ------------------------------------------------------------------
    def speed(self, s: float) -> float:
        return 1.0

    def segment_length(self, ta: float, tb: float, nsub: int = 1) -> float:
        return tb - ta

    def length(self) -> float:
        return self.tmax - self.tmin

    def equidistant(self, n: int) -> np.ndarray:
        return np.linspace(self.tmin, self.tmax, n)

    def sample_spacing(self, ds: float, include_breakpoints: bool = True) -> np.ndarray:
        """Sample parameters with spacing ~ds, snapped so that breakpoints are hit exactly."""
        edges = np.concatenate([[self.tmin], self.breakpoints, [self.tmax]]) if include_breakpoints \
            else np.array([self.tmin, self.tmax])
        edges = np.unique(edges)
        pts = [edges[0]]
        for k in range(len(edges) - 1):
            m = max(1, int(np.ceil((edges[k + 1] - edges[k]) / ds)))
            pts.extend(np.linspace(edges[k], edges[k + 1], m + 1)[1:])
        self.t = np.array(pts)
        return self.t

    def connector_states(self):
        """(r, T, n, b) at the two free ends of the curve (inner first)."""
        out = []
        for s in (self.tmin, self.tmax):
            T, n, b = self.frame(s)
            out.append((self.r(s), T, n, b))
        return out
