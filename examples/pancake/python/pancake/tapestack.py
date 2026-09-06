"""A stack of tapes swept along a base curve.

At every station t of the curve the local frame R = transform(t) has columns
(b, -n, T).  The stack lives in the (b, n) plane: tape width along b (local
x0), the tapes offset along the stack normal (local y0), which follows the
normal n of the base curve.  Every tape is represented like in
tapestack3d.geo by three points across its width: left tip, centre, right tip.
"""

import numpy as np
from scipy.spatial import cKDTree

from .basecurve import Basecurve


class TapeStack:

    def __init__(self, curve: Basecurve, numtapes: int, tapewidth: float, tapedistance: float, t=None):
        if t is None:
            t = curve.t
        if t is None:
            raise ValueError("sample the curve first (curve.sample_spacing(ds))")
        self.curve = curve
        self.t = np.asarray(t, float)
        self.numtapes = int(numtapes)
        self.tapewidth = float(tapewidth)
        self.tapedistance = float(tapedistance)
        self.thickness = self.tapedistance * (self.numtapes - 1)   # distance between outer tapes

        # local cross-section coordinates (x0 along b, y0 along the stack normal)
        self.y = np.linspace(-0.5 * self.thickness, 0.5 * self.thickness, self.numtapes)
        self.x = np.array([-0.5 * self.tapewidth, 0.0, 0.5 * self.tapewidth])

        m = len(self.t)
        self.center = np.zeros((m, 3))
        self.frames = np.zeros((m, 3, 3))
        # points[k, i, j] : tape k, station i, j = 0 left tip, 1 centre, 2 right tip
        self.points = np.zeros((self.numtapes, m, 3, 3))
        for i, ti in enumerate(self.t):
            r = curve.r(ti)
            R = curve.transform(ti)
            self.center[i] = r
            self.frames[i] = R
            for k in range(self.numtapes):
                for j in range(3):
                    self.points[k, i, j] = r + R @ np.array([self.x[j], self.y[k], 0.0])

    # ------------------------------------------------------------------
    def envelope(self):
        """Corners (m, 4, 3) of the rectangle enclosing the stack at every station."""
        local = np.array([[self.x[0], self.y[0], 0.0], [self.x[2], self.y[0], 0.0],
                          [self.x[2], self.y[-1], 0.0], [self.x[0], self.y[-1], 0.0]])
        return self.center[:, None, :] + np.einsum("mab,cb->mca", self.frames, local)

    def clearance(self, min_arc_distance=None):
        """Smallest distance between stack cross sections that are far apart along the curve.

        Returns (distance, i, j) for the closest pair of stations with
        |t_i - t_j| > min_arc_distance (default: twice the stack diagonal).  A
        distance below zero cannot occur with this estimate; a value close to
        zero (below one tape distance, say) means the tapes cross or touch.
        """
        diag = np.hypot(self.tapewidth, self.thickness)
        if min_arc_distance is None:
            min_arc_distance = 2.0 * diag
        corners = self.envelope()
        # boundary samples of every rectangle: corners and edge mid points
        mids = 0.5 * (corners + np.roll(corners, -1, axis=1))
        samples = np.concatenate([corners, mids], axis=1)      # (m, 8, 3)

        tree = cKDTree(self.center)
        pairs = tree.query_pairs(r=diag, output_type="ndarray")
        if len(pairs) == 0:
            return np.inf, -1, -1
        keep = np.abs(self.t[pairs[:, 0]] - self.t[pairs[:, 1]]) > min_arc_distance
        pairs = pairs[keep]
        if len(pairs) == 0:
            return np.inf, -1, -1

        best = (np.inf, -1, -1)
        for i, j in pairs:
            d = min(self._points_to_rectangle(samples[i], corners[j], self.frames[j]),
                    self._points_to_rectangle(samples[j], corners[i], self.frames[i]))
            if d < best[0]:
                best = (float(d), int(i), int(j))
        return best

    @staticmethod
    def _points_to_rectangle(P, corners, R):
        """Minimum distance from points P (p, 3) to the rectangle given by its corners and frame R."""
        e1 = R[:, 0]     # width direction
        e2 = R[:, 1]     # stack normal direction
        nrm = R[:, 2]
        c = corners.mean(axis=0)
        d = P - c
        u = d @ e1
        v = d @ e2
        w = d @ nrm
        hu = 0.5 * np.linalg.norm(corners[1] - corners[0])
        hv = 0.5 * np.linalg.norm(corners[3] - corners[0])
        du = np.maximum(np.abs(u) - hu, 0.0)
        dv = np.maximum(np.abs(v) - hv, 0.0)
        return float(np.sqrt(du ** 2 + dv ** 2 + w ** 2).min())

    # ------------------------------------------------------------------
    def plot(self, ax, tapes=None, every=1, color="k", lw=0.4, sections_every=0, section_color="tab:orange"):
        """Draw the tips of the tapes as lines and, optionally, the stack envelope at some stations."""
        if tapes is None:
            tapes = range(self.numtapes)
        for k in tapes:
            for j in (0, 2):
                p = self.points[k, ::every, j]
                ax.plot(p[:, 0], p[:, 1], p[:, 2], color=color, lw=lw)
        if sections_every > 0:
            corners = self.envelope()
            for i in range(0, len(self.t), sections_every):
                q = np.vstack([corners[i], corners[i][:1]])
                ax.plot(q[:, 0], q[:, 1], q[:, 2], color=section_color, lw=1.0)
