"""Lead (terminal) sections defined intrinsically by their strip curvatures.

A lead is a chain of segments.  Each segment prescribes the three curvature
parameters of the strip along its arc length x in [0, L]

    omega(x) = ( tau(x), kappa_g(x), kappa_n(x) )

and the space curve together with its frame follows by integrating the
generalized Frenet-Serret equations for strips (Russenschuck Eq. (19.9))

    r' = T,  T' = kappa_n n - kappa_g b,  n' = -kappa_n T + tau b,  b' = kappa_g T - tau n

from the state at the end of the winding.  Because the frame is *integrated*
rather than derived from v x a, it is defined and continuous on straight
sections, through inflection points and across the transition into the bend.
The derivatives v, a, b of the position are then exact expressions in the
frame and the curvature functions:

    v = T
    a = T' = kappa_n n - kappa_g b
    b = T'' = kappa_n' n + kappa_n n' - kappa_g' b - kappa_g b'

Sign conventions (positive amplitudes):
    kappa_n  turns the tangent towards +n   ("easy way" bend about b)
    kappa_g  turns the tangent towards -b   ("hard way" bend about n)
    tau      rotates n towards +b           (right-handed twist about T)

Each amplitude is multiplied by a profile P(x):

    shape = "plateau": rises smoothly from 0 to 1 over the first `ramp` length
                       units and falls back to 0 over the last `ramp`
                       (ramp = 0: constant curvature -> circular arc / helix / plain twist)
    shape = "fade":    1 -> 0 over the whole segment (used to release the winding curvature)
    shape = "rise":    0 -> 1 over the whole segment

The quintic smoothstep makes the curvatures C^2, so frame, curvature and jerk
are continuous everywhere.  On intervals with constant curvature the frame and
position are advanced with the closed-form matrix exponential (exact to
round-off); on ramps a fourth-order Runge-Kutta scheme with a small step is
used.
"""

import math
import numpy as np


def smoothstep(x: float):
    """Quintic smoothstep S(x) = 6x^5 - 15x^4 + 10x^3 on [0, 1] and its derivative."""
    if x <= 0.0:
        return 0.0, 0.0
    if x >= 1.0:
        return 1.0, 0.0
    return x * x * x * (10.0 + x * (-15.0 + 6.0 * x)), 30.0 * x * x * (1.0 - x) ** 2


class Segment:
    """Strip section of given length with curvature amplitudes times a profile."""

    def __init__(self, length: float, tau: float = 0.0, kappa_g: float = 0.0, kappa_n: float = 0.0,
                 ramp: float = 0.0, shape: str = "plateau", name: str = "segment"):
        if length <= 0.0:
            raise ValueError("segment length must be positive")
        if shape not in ("plateau", "fade", "rise"):
            raise ValueError("shape must be 'plateau', 'fade' or 'rise'")
        if shape != "plateau" and ramp != 0.0:
            raise ValueError("ramp is only used with shape='plateau'")
        if ramp < 0.0 or 2.0 * ramp > length * (1.0 + 1e-12):
            raise ValueError("ramp must satisfy 0 <= 2 ramp <= length")
        self.length = float(length)
        self.ramp = float(ramp)
        self.shape = shape
        self.amplitude = np.array([tau, kappa_g, kappa_n], dtype=float)
        self.name = name
        if shape == "plateau":
            knots = [0.0, self.ramp, self.length - self.ramp, self.length]
        else:
            knots = [0.0, self.length]
        self.knots = np.unique(np.clip(np.array(knots), 0.0, self.length))

    def profile(self, x: float):
        """Profile P(x) and dP/dx."""
        if self.shape == "fade":
            p, dp = smoothstep(x / self.length)
            return 1.0 - p, -dp / self.length
        if self.shape == "rise":
            p, dp = smoothstep(x / self.length)
            return p, dp / self.length
        if self.ramp == 0.0:
            return 1.0, 0.0
        if x < self.ramp:
            p, dp = smoothstep(x / self.ramp)
            return p, dp / self.ramp
        if x > self.length - self.ramp:
            p, dp = smoothstep((self.length - x) / self.ramp)
            return p, -dp / self.ramp
        return 1.0, 0.0

    def is_constant(self, x0: float, x1: float) -> bool:
        """True if the profile is constant (== 1) on [x0, x1]."""
        if self.shape != "plateau":
            return False
        if self.ramp == 0.0:
            return True
        tol = 1e-12 * self.length
        return x0 >= self.ramp - tol and x1 <= self.length - self.ramp + tol

    def omega(self, x: float) -> np.ndarray:
        return self.amplitude * self.profile(x)[0]

    def domega(self, x: float) -> np.ndarray:
        return self.amplitude * self.profile(x)[1]

    @property
    def effective_length(self) -> float:
        """Integral of the profile: the length that contributes to the total rotation."""
        if self.shape == "plateau":
            return self.length - self.ramp
        return 0.5 * self.length

    @property
    def total_rotation(self) -> np.ndarray:
        """Integrated (tau, kappa_g, kappa_n) in radians."""
        return self.amplitude * self.effective_length

    def __repr__(self):
        return "%s(L=%g, tau=%g, kappa_g=%g, kappa_n=%g, ramp=%g, shape=%s)" % (
            self.name, self.length, self.amplitude[0], self.amplitude[1], self.amplitude[2], self.ramp, self.shape)


def Straight(length: float) -> Segment:
    return Segment(length, name="straight")


def Twist(angle_deg: float, length: float, ramp: float = 0.0) -> Segment:
    """Twist the tape by angle_deg about its tangent over the given length (base curve straight)."""
    tau = math.radians(angle_deg) / (length - ramp)
    return Segment(length, tau=tau, ramp=ramp, name="twist")


def Bend(angle_deg: float, radius: float, ramp: float = 0.0, hard: bool = False) -> Segment:
    """Bend the tape by angle_deg with the given bending radius.

    hard=False: easy-way bend about b (kappa_n), the tangent turns towards +n for
    positive angles.  hard=True: hard-way bend about n (kappa_g), the tangent turns
    towards -b for positive angles.  A ramp > 0 blends the curvature in and out
    (clothoid-like), which lengthens the segment by `ramp`.
    """
    if radius <= 0.0:
        raise ValueError("bend radius must be positive")
    ang = math.radians(angle_deg)
    length = abs(ang) * radius + ramp
    kappa = math.copysign(1.0 / radius, ang)
    if hard:
        return Segment(length, kappa_g=kappa, ramp=ramp, name="hardbend")
    return Segment(length, kappa_n=kappa, ramp=ramp, name="bend")


def Release(length: float, kappa_n: float = 0.0, kappa_g: float = 0.0, tau: float = 0.0) -> Segment:
    """Fade the given curvatures smoothly to zero over `length` (tape leaving the winding)."""
    return Segment(length, tau=tau, kappa_g=kappa_g, kappa_n=kappa_n, shape="fade", name="release")


# ----------------------------------------------------------------------
# integration of  r' = T,  F' = F [omega]_x   with  F = [T n b]
# ----------------------------------------------------------------------
def _rhs(y: np.ndarray, omega: np.ndarray) -> np.ndarray:
    T = y[3:6]
    n = y[6:9]
    b = y[9:12]
    tau, kg, kn = omega
    dy = np.empty(12)
    dy[0:3] = T
    dy[3:6] = kn * n - kg * b
    dy[6:9] = -kn * T + tau * b
    dy[9:12] = kg * T - tau * n
    return dy


def _orthonormalize(y: np.ndarray) -> np.ndarray:
    T = y[3:6] / np.linalg.norm(y[3:6])
    n = y[6:9] - np.dot(y[6:9], T) * T
    n /= np.linalg.norm(n)
    b = np.cross(T, n)
    y[3:6] = T
    y[6:9] = n
    y[9:12] = b
    return y


def _rk4(y: np.ndarray, seg: Segment, x0: float, h: float) -> np.ndarray:
    k1 = _rhs(y, seg.omega(x0))
    k2 = _rhs(y + 0.5 * h * k1, seg.omega(x0 + 0.5 * h))
    k3 = _rhs(y + 0.5 * h * k2, seg.omega(x0 + 0.5 * h))
    k4 = _rhs(y + h * k3, seg.omega(x0 + h))
    return _orthonormalize(y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4))


def _exact(y: np.ndarray, omega: np.ndarray, h: float) -> np.ndarray:
    """Closed-form step for constant body-frame Darboux vector omega = (tau, kappa_g, kappa_n).

    F(h) = F0 exp(h K),  r(h) = r0 + F0 ( int_0^h exp(x K) dx ) e_1,  K = [omega]_x,
    with exp(x K) = I + sin(w x) K^ + (1 - cos(w x)) K^^2  (Rodrigues), w = |omega|.
    """
    r = y[0:3]
    F = np.column_stack([y[3:6], y[6:9], y[9:12]])
    w = float(np.linalg.norm(omega))
    out = np.empty(12)
    if w * h < 1e-14:
        out[0:3] = r + h * F[:, 0]
        out[3:12] = y[3:12]
        return out
    k = omega / w
    K = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    K2 = K @ K
    sw, cw = math.sin(w * h), math.cos(w * h)
    E = np.eye(3) + sw * K + (1.0 - cw) * K2
    I = h * np.eye(3) + (1.0 - cw) / w * K + (h - sw / w) * K2
    Fn = F @ E
    out[0:3] = r + F @ I[:, 0]
    out[3:6] = Fn[:, 0]
    out[6:9] = Fn[:, 1]
    out[9:12] = Fn[:, 2]
    return _orthonormalize(out)


class Lead:
    """Chain of segments integrated from a start state (r0, T0, n0, b0)."""

    def __init__(self, segments, max_step_angle: float = 0.01, min_steps: int = 32):
        self.segments = list(segments)
        if not self.segments:
            raise ValueError("a lead needs at least one segment")
        self.max_step_angle = float(max_step_angle)
        self.min_steps = int(min_steps)
        self.length = float(sum(seg.length for seg in self.segments))
        self.offsets = np.concatenate([[0.0], np.cumsum([seg.length for seg in self.segments])])
        self._nodes_s = None
        self._nodes_y = None
        self._node_seg = None   # index of the segment owning the interval [s_k, s_{k+1}]
        self.built = False

    # ------------------------------------------------------------------
    def _advance(self, y, seg: Segment, x0: float, h: float, nsub: int = 1):
        if seg.is_constant(x0, x0 + h):
            return _exact(y, seg.amplitude, h)
        hs = h / nsub
        for i in range(nsub):
            y = _rk4(y, seg, x0 + i * hs, hs)
        return y

    def build(self, r0, T0, n0, b0):
        y = np.concatenate([np.asarray(r0, float), np.asarray(T0, float),
                            np.asarray(n0, float), np.asarray(b0, float)])
        y = _orthonormalize(y.copy())
        nodes_s = [0.0]
        nodes_y = [y.copy()]
        node_seg = []
        s = 0.0
        for iseg, seg in enumerate(self.segments):
            amp = float(np.linalg.norm(seg.amplitude))
            for k in range(len(seg.knots) - 1):
                x0, x1 = seg.knots[k], seg.knots[k + 1]
                if x1 <= x0:
                    continue
                if seg.is_constant(x0, x1):
                    nsteps = max(4, int(math.ceil((x1 - x0) * amp / (10.0 * self.max_step_angle))))
                else:
                    nsteps = max(self.min_steps, int(math.ceil((x1 - x0) * amp / self.max_step_angle)))
                h = (x1 - x0) / nsteps
                for i in range(nsteps):
                    y = self._advance(y, seg, x0 + i * h, h)
                    s += h
                    nodes_s.append(s)
                    nodes_y.append(y.copy())
                    node_seg.append(iseg)
            s = float(self.offsets[iseg + 1])
            nodes_s[-1] = s
        self._nodes_s = np.array(nodes_s)
        self._nodes_y = np.array(nodes_y)
        self._node_seg = np.array(node_seg, dtype=int)
        self.built = True
        return self

    # ------------------------------------------------------------------
    def segment_at(self, s: float):
        """(segment, local coordinate x) for global lead coordinate s."""
        i = int(np.searchsorted(self.offsets, s, side="right") - 1)
        i = min(max(i, 0), len(self.segments) - 1)
        return self.segments[i], s - self.offsets[i]

    def omega(self, s: float) -> np.ndarray:
        seg, x = self.segment_at(s)
        return seg.omega(x)

    def domega(self, s: float) -> np.ndarray:
        seg, x = self.segment_at(s)
        return seg.domega(x)

    def state(self, s: float):
        """(r, T, n, b) at lead coordinate s in [0, length]."""
        if not self.built:
            raise RuntimeError("lead has not been built")
        s = float(s)
        k = int(np.searchsorted(self._nodes_s, s, side="right") - 1)
        k = min(max(k, 0), len(self._nodes_s) - 1)
        y = self._nodes_y[k].copy()
        ds = s - self._nodes_s[k]
        if ds != 0.0:
            k = min(k, len(self._node_seg) - 1)
            seg = self.segments[self._node_seg[k]]
            x0 = self._nodes_s[k] - self.offsets[self._node_seg[k]]
            y = self._advance(y, seg, x0, ds, nsub=4)
        return y[0:3], y[3:6], y[6:9], y[9:12]

    def derivatives(self, s: float):
        """(r, v, a, jerk) at lead coordinate s from the frame and the curvature functions."""
        r, T, n, b = self.state(s)
        tau, kg, kn = self.omega(s)
        dtau, dkg, dkn = self.domega(s)
        dT = kn * n - kg * b
        dn = -kn * T + tau * b
        db = kg * T - tau * n
        ddT = dkn * n + kn * dn - dkg * b - kg * db
        return r, T, dT, ddT

    def end_state(self):
        return self.state(self.length)

    @property
    def breakpoints(self) -> np.ndarray:
        return self.offsets.copy()

    def __repr__(self):
        return "Lead(" + ", ".join(repr(s) for s in self.segments) + ")"
