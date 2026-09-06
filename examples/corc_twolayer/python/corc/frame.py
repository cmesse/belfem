"""Continuous Frenet-based frame along a Curve.

Algorithm (pinned by the audited plan, refined per the final audit):
- Frenet frame from exact derivatives (t = v/|v|, b = v x a/|v x a|,
  n = b x t) wherever kappa >= KAPPA_TOL.
- Reference-up fallback where kappa < KAPPA_TOL (straight segments), with
  hysteresis so the frame does not chatter at the threshold. The
  Frenet/fallback decision of a query point follows the nearest reference
  sample, so it is piecewise constant in s.
- Sign continuity enforced against a dense precomputed reference sampling:
  n and b flip together so that n stays continuous along s. NOTE: through
  an inflection (kappa -> 0 crossing) this yields a CONTINUOUS frame that
  is deliberately NOT the pure Frenet frame (whose normal would jump);
  what raises FrameError is a rotation between adjacent samples that stays
  above MAX_STEP_ROTATION even after adaptive resampling — i.e. a genuine
  kink, not a smooth inflection.
"""

import math

import numpy as np


class FrameError(Exception):
    pass


class Frame:

    KAPPA_TOL = 1e-8          # below: fallback frame (1/mm)
    MAX_STEP_ROTATION = 0.35  # rad between adjacent reference samples

    def __init__(self, curve, reference_up=(0.0, 0.0, 1.0), nsamples=2048):
        self.curve = curve
        self.reference_up = np.array(reference_up, dtype=float)
        # adaptive density: a smooth admissible curve can legitimately
        # rotate its frame fast (high kappa*L); refine before concluding
        # that a step-rotation violation is a genuine kink
        for attempt in range(3):
            try:
                self._build_reference(nsamples)
                return
            except FrameError as err:
                if "frame rotates" not in str(err) or attempt == 2:
                    raise
                nsamples *= 8

    # ------------------------------------------------------------ raw frame

    def _raw_frame(self, s: float, use_fallback=None):
        """Pointwise frame before continuity fixes: returns (p, t, n, b,
        kappa). The Frenet-vs-fallback decision can be forced so that
        queries follow the hysteresis decision of the reference pass
        instead of chattering around the threshold."""
        tpar = self.curve.t_of_s(s)
        p, v, a, _ = self.curve.derivatives(tpar)
        nv = np.linalg.norm(v)
        if nv < 1e-12:
            raise FrameError("zero velocity at s={:g}".format(s))
        t = v / nv
        va = np.cross(v, a)
        nva = np.linalg.norm(va)
        kappa = nva / nv ** 3
        if use_fallback is None:
            use_fallback = kappa < self.KAPPA_TOL
        if not use_fallback and nva > 0.0:
            b = va / nva
            n = np.cross(b, t)
        else:
            n = self._fallback_normal(t)
            b = np.cross(t, n)
        return p, t, n, b, kappa

    def _fallback_normal(self, t: np.ndarray) -> np.ndarray:
        up = self.reference_up
        n = up - np.dot(up, t) * t
        norm = np.linalg.norm(n)
        if norm < 1e-8:
            up = np.array([1.0, 0.0, 0.0]) if abs(t[0]) < 0.9 \
                else np.array([0.0, 1.0, 0.0])
            n = up - np.dot(up, t) * t
            norm = np.linalg.norm(n)
        return n / norm

    # ------------------------------------------------------- reference pass

    def _build_reference(self, nsamples: int):
        L = self.curve.length
        self._sref = np.linspace(0.0, L, nsamples + 1)
        self._nref = np.zeros((nsamples + 1, 3))
        self._kappa = np.zeros(nsamples + 1)

        self._fallback = np.zeros(nsamples + 1, dtype=bool)

        prev_n = None
        prev_t = None
        fallback = None
        for i, s in enumerate(self._sref):
            # hysteresis: enter the fallback below tol/2, leave above 2*tol
            _, _, _, _, kappa_probe = self._raw_frame(float(s))
            if fallback is None:
                fallback = kappa_probe < self.KAPPA_TOL
            elif fallback and kappa_probe > 2.0 * self.KAPPA_TOL:
                fallback = False
            elif not fallback and kappa_probe < 0.5 * self.KAPPA_TOL:
                fallback = True
            self._fallback[i] = fallback

            _, t, n, b, kappa = self._raw_frame(float(s), use_fallback=fallback)
            if prev_n is not None:
                if np.dot(prev_n, n) < 0.0:
                    n = -n
                # incremental rotation between frames
                cosang = np.clip(np.dot(prev_n, n), -1.0, 1.0)
                cosang_t = np.clip(np.dot(prev_t, t), -1.0, 1.0)
                ang = max(math.acos(cosang), math.acos(cosang_t))
                if ang > self.MAX_STEP_ROTATION:
                    raise FrameError(
                        "frame rotates {:.3f} rad between s={:.6g} and "
                        "s={:.6g} — inflection or kink in the centerline; "
                        "refusing to sweep".format(
                            ang, self._sref[i - 1], s))
            self._nref[i] = n
            self._kappa[i] = kappa
            prev_n = n
            prev_t = t

        self.kappa_max = float(self._kappa.max())

    # ----------------------------------------------------------- public API

    def _nearest_sample(self, s: float) -> int:
        i = int(np.clip(np.searchsorted(self._sref, s), 1, len(self._sref) - 1))
        if abs(self._sref[i - 1] - s) < abs(self._sref[i] - s):
            i -= 1
        return i

    def R(self, s: float) -> np.ndarray:
        """Rotation matrix with columns [n, b, t]: local (x, y, s) ->
        global; the Frenet/fallback mode and the n,b sign follow the
        reference sampling (continuity, no threshold chatter)."""
        i = self._nearest_sample(s)
        _, t, n, b, _ = self._raw_frame(s, use_fallback=bool(self._fallback[i]))
        if np.dot(self._nref[i], n) < 0.0:
            n = -n
            b = -b
        return np.column_stack((n, b, t))

    def origin(self, s: float) -> np.ndarray:
        return self.curve.p(self.curve.t_of_s(s))

    def transform(self, s: float, x: float, y: float) -> np.ndarray:
        """Map local cross-section coordinates (x, y) at arc length s to
        global coordinates."""
        p = self.origin(s)
        R = self.R(s)
        return p + x * R[:, 0] + y * R[:, 1]

    # ---------------------------------------------------------- validations

    def check_domain_admissibility(self, domain_radius: float):
        """Local: the domain tube must not fold (kappa * r < 1). Global:
        distant centerline points must stay 2 r apart."""
        if self.kappa_max * domain_radius >= 1.0:
            raise FrameError(
                "domain tube folds: kappa_max = {:.4g} 1/mm requires "
                "domainRadius < {:.4g} mm (got {:g})".format(
                    self.kappa_max, 1.0 / self.kappa_max, domain_radius))

        # short arc separations are covered by the local check (worst case
        # is a circle at kappa_max: chord >= 2r for arc <= pi*r as long as
        # kappa*r < 1), so only test pairs separated by more than pi*r —
        # otherwise gentle bends trip the check spuriously. For curves
        # shorter than pi*r NO pair is tested, and that is sound: turning
        # back towards oneself needs total curvature >= pi, i.e. arc
        # length >= pi/kappa_max > pi*r under the local condition.
        pts = np.array([self.origin(float(s)) for s in self._sref])
        n = len(pts)
        ds = self._sref[1] - self._sref[0]
        min_sep = int(math.ceil(math.pi * domain_radius / ds)) + 1
        worst = np.inf
        for i in range(n):
            j0 = i + min_sep
            if j0 >= n:
                break
            d = np.linalg.norm(pts[j0:] - pts[i], axis=1)
            worst = min(worst, float(d.min()))
        if worst < 2.0 * domain_radius:
            raise FrameError(
                "domain tube self-intersects: distant centerline points "
                "come within {:.4g} mm < 2 * domainRadius = {:.4g} mm".format(
                    worst, 2.0 * domain_radius))
