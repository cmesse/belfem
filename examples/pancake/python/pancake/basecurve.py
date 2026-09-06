"""Base class for regular space curves with analytic derivatives and a moving frame.

Notation follows Russenschuck, *Field Computation for Accelerator Magnets* (2010):

* Sec. 3.1  Frenet frame of a space curve with general parameter t,
  Eqs. (3.29)-(3.30):  T = v/|v|,  B = v x a / |v x a|,  N = B x T,
  kappa = |v x a| / |v|^3,  tau = (v x a) . da/dt / |v x a|^2.
* Sec. 19.2 generalized Frenet-Serret equations for strips, Eq. (19.9):

      T' =  kappa_n n - kappa_g b
      n' = -kappa_n T + tau     b
      b' =  kappa_g T - tau     n

  with the strip frame {T, n, b}: T tangent, n normal to the broad face of
  the tape, b along the tape width.  kappa_n is the "easy way" (normal)
  curvature, kappa_g the "hard way" (geodesic) curvature and tau the twist.
* Eqs. (19.30)-(19.35): an additional twist angle theta_T(s) about T
  rotates n, b into n*, b* and changes the curvature parameters to
  tau* = tau + dtheta_T/ds, kappa_g* = cos(theta) kappa_g + sin(theta) kappa_n,
  kappa_n* = cos(theta) kappa_n - sin(theta) kappa_g.

The default `frame` implemented here is the classical Frenet frame, which is
only defined where the curvature does not vanish.  Curves with straight
sections (e.g. the leads of a pancake coil) must override `frame` and
`curvatures` with a frame that stays defined and continuous.
"""

import numpy as np

_GL_X, _GL_W = np.polynomial.legendre.leggauss(7)

_EZ = np.array([0.0, 0.0, 1.0])


class Basecurve:
    """Regular space curve r(t) on [tmin, tmax] with derivatives v, a, b.

    Subclasses implement

        r(t) -> position            (3,)
        v(t) -> dr/dt               (3,)
        a(t) -> d2r/dt2             (3,)
        b(t) -> d3r/dt3 ("jerk")    (3,)

    The parameter t need not be the arc length, but all frame quantities are
    computed such that they are invariant under reparametrization.
    """

    def __init__(self, tmin: float = 0.0, tmax: float = 1.0):
        self.tmin = float(tmin)
        self.tmax = float(tmax)
        # sample parameters, filled by `sample`
        self.t = None

    # ------------------------------------------------------------------
    # derivatives, to be provided by the subclass
    # ------------------------------------------------------------------
    def r(self, t: float) -> np.ndarray:
        raise NotImplementedError()

    def v(self, t: float) -> np.ndarray:
        raise NotImplementedError()

    def a(self, t: float) -> np.ndarray:
        raise NotImplementedError()

    def b(self, t: float) -> np.ndarray:
        raise NotImplementedError()

    # ------------------------------------------------------------------
    # arc length
    # ------------------------------------------------------------------
    def speed(self, t: float) -> float:
        return float(np.linalg.norm(self.v(t)))

    def segment_length(self, ta: float, tb: float, nsub: int = 1) -> float:
        """Arc length between ta and tb by composite 7-point Gauss-Legendre."""
        edges = np.linspace(ta, tb, nsub + 1)
        total = 0.0
        for k in range(nsub):
            t0, t1 = edges[k], edges[k + 1]
            hm = 0.5 * (t1 - t0)
            tm = 0.5 * (t1 + t0)
            acc = 0.0
            for x, w in zip(_GL_X, _GL_W):
                acc += w * self.speed(tm + hm * x)
            total += hm * acc
        return total

    def length(self) -> float:
        return self.segment_length(self.tmin, self.tmax, nsub=max(1, self._nsub_default()))

    def _nsub_default(self) -> int:
        return 200

    def equidistant(self, n: int) -> np.ndarray:
        """Parameters t_0..t_{n-1} such that the arc length between neighbours is constant."""
        m = max(20 * n, 400)
        tt = np.linspace(self.tmin, self.tmax, m + 1)
        s = np.zeros(m + 1)
        for k in range(m):
            s[k + 1] = s[k] + self.segment_length(tt[k], tt[k + 1])
        target = np.linspace(0.0, s[-1], n)
        out = np.empty(n)
        out[0] = self.tmin
        out[-1] = self.tmax
        for i in range(1, n - 1):
            k = int(np.searchsorted(s, target[i], side="right") - 1)
            k = min(max(k, 0), m - 1)
            # Newton on the cumulative length within interval k
            t = tt[k] + (target[i] - s[k]) / (s[k + 1] - s[k]) * (tt[k + 1] - tt[k])
            for _ in range(30):
                f = s[k] + self.segment_length(tt[k], t) - target[i]
                dt = f / self.speed(t)
                t -= dt
                if abs(dt) < 1e-13 * max(1.0, abs(t)):
                    break
            out[i] = t
        return out

    def sample(self, n: int) -> np.ndarray:
        """Store and return n parameters with equal arc-length spacing."""
        self.t = self.equidistant(n)
        return self.t

    # ------------------------------------------------------------------
    # frames and curvatures
    # ------------------------------------------------------------------
    def frame(self, t: float):
        """Strip frame (T, n, b).  Default: Frenet frame (T, N, B), Eq. (3.29)."""
        v = self.v(t)
        a = self.a(t)
        nv = np.linalg.norm(v)
        T = v / nv
        vxa = np.cross(v, a)
        nvxa = np.linalg.norm(vxa)
        if nvxa <= 1e-12 * nv * nv * max(1.0, np.linalg.norm(a)):
            raise ValueError("Frenet frame undefined (vanishing curvature) at t = %g" % t)
        B = vxa / nvxa
        N = np.cross(B, T)
        return T, N, B

    def curvatures(self, t: float):
        """Curvature parameters (tau, kappa_g, kappa_n) of the frame returned by `frame`.

        For the Frenet frame this is (tau, 0, kappa), Eq. (3.30).
        """
        kappa, tau = self.kappa_tau(t)
        return tau, 0.0, kappa

    def kappa_tau(self, t: float):
        """Frenet curvature and torsion, Eq. (3.30)."""
        v = self.v(t)
        a = self.a(t)
        j = self.b(t)
        vxa = np.cross(v, a)
        nvxa = np.linalg.norm(vxa)
        nv = np.linalg.norm(v)
        kappa = nvxa / nv ** 3
        tau = float(np.dot(vxa, j) / nvxa ** 2) if nvxa > 0.0 else 0.0
        return float(kappa), tau

    def twisted_frame(self, t: float, theta_T: float = 0.0):
        """Frame (T, n*, b*) after twisting about T by theta_T, Eqs. (19.31)-(19.33)."""
        T, n, b = self.frame(t)
        if theta_T == 0.0:
            return T, n, b
        c = np.cos(theta_T)
        s = np.sin(theta_T)
        return T, c * n + s * b, c * b - s * n

    def transform(self, t: float, theta_T: float = 0.0) -> np.ndarray:
        """Rotation matrix mapping local cross-section coordinates to space.

        Column 0: b*  (tape width direction)
        Column 1: -n* (tape thickness direction)
        Column 2: T   (tape length direction)

        Local coordinates (x0, y0, 0) of a cross section are mapped by
        r(t) + R @ (x0, y0, 0).  The triple (b, -n, T) is right handed.
        """
        T, n, b = self.twisted_frame(t, theta_T)
        R = np.empty((3, 3))
        R[:, 0] = b
        R[:, 1] = -n
        R[:, 2] = T
        return R

    def strip_curvatures(self, t: float, theta_T: float = 0.0, dtheta_T_ds: float = 0.0):
        """(tau*, kappa_g*, kappa_n*) of the twisted strip, Eqs. (19.30), (19.34), (19.35)."""
        tau, kg, kn = self.curvatures(t)
        c = np.cos(theta_T)
        s = np.sin(theta_T)
        return tau + dtheta_T_ds, c * kg + s * kn, c * kn - s * kg

    def darboux(self, t: float) -> np.ndarray:
        """Darboux vector of the strip frame in space coordinates: tau T + kappa_g n + kappa_n b."""
        T, n, b = self.frame(t)
        tau, kg, kn = self.curvatures(t)
        return tau * T + kg * n + kn * b


class RigidTransform(Basecurve):
    """Curve obtained from another one by a rotation R and a translation d: r* = R r + d."""

    def __init__(self, curve: Basecurve, R=None, d=None):
        Basecurve.__init__(self, curve.tmin, curve.tmax)
        self.curve = curve
        self.R = np.eye(3) if R is None else np.asarray(R, dtype=float)
        self.d = np.zeros(3) if d is None else np.asarray(d, dtype=float)
        if abs(np.linalg.det(self.R) - 1.0) > 1e-10:
            raise ValueError("R must be a proper rotation matrix")
        self.t = curve.t

    def r(self, t):
        return self.R @ self.curve.r(t) + self.d

    def v(self, t):
        return self.R @ self.curve.v(t)

    def a(self, t):
        return self.R @ self.curve.a(t)

    def b(self, t):
        return self.R @ self.curve.b(t)

    def frame(self, t):
        T, n, b = self.curve.frame(t)
        return self.R @ T, self.R @ n, self.R @ b

    def curvatures(self, t):
        return self.curve.curvatures(t)

    def equidistant(self, n):
        return self.curve.equidistant(n)


def rotation_z(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotation_x(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
