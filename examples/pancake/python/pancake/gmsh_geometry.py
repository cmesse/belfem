"""Gmsh .geo writer: tapes as strip surfaces along a base curve, stacked into blocks.

Ported from the old `frenet` package.  A `CrossSection` holds the local
(width, thickness) coordinates of the left and right edge of every tape;
`Tape` maps them through `Basecurve.transform` to space, and `TapeBlock`
closes the volume between two neighbouring tapes.
"""

import numpy as np

from .basecurve import Basecurve


class Point:
    def __init__(self, x, y, z, res):
        self.id = 0
        self.x, self.y, self.z, self.res = float(x), float(y), float(z), float(res)

    def write(self):
        return "Point(%d) = {%.12f,%.12f,%.12f,%.3f};" % (self.id, self.x, self.y, self.z, self.res)


class Curve:
    def __init__(self, label: str):
        self.id = 0
        self.label = label
        self.points = []

    def write(self):
        ids = [p.id for p in self.points]
        n = len(ids)
        line = "%s(%d) = {" % (self.label, self.id)
        if n > 1 and (all(ids[k] == ids[k - 1] + 1 for k in range(1, n)) or
                      all(ids[k] == ids[k - 1] - 1 for k in range(1, n))):
            line += "%d:%d" % (ids[0], ids[-1])
        else:
            line += ",".join(str(i) for i in ids)
        return line + "};"


class CurveLoop:
    def __init__(self):
        self.id = 0
        self.curves = []
        self.signs = []

    def write(self):
        signs = self.signs if self.signs else [1] * len(self.curves)
        return "Curve Loop(%d) = {" % self.id + ", ".join(
            "%d" % (s * c.id) for s, c in zip(signs, self.curves)) + "};"


class Surface:
    def __init__(self):
        self.id = 0
        self.loops = []
        self.is_plane = False

    def write(self):
        label = "Plane Surface" if self.is_plane else "Surface"
        return "%s(%d) = {" % (label, self.id) + ", ".join("%d" % l.id for l in self.loops) + "};"


class SurfaceLoop:
    def __init__(self):
        self.id = 0
        self.surfaces = []
        self.signs = []

    def write(self):
        signs = self.signs if self.signs else [1] * len(self.surfaces)
        return "Surface Loop(%d) = {" % self.id + ", ".join(
            "%d" % (s * f.id) for s, f in zip(signs, self.surfaces)) + "};"


class Volume:
    def __init__(self):
        self.id = 0
        self.loops = []

    def write(self):
        return "Volume(%d) = {" % self.id + ", ".join("%d" % l.id for l in self.loops) + "};"


class CrossSection:
    """`numtapes` tapes of width `tapewidth` stacked in the thickness direction with `tapedistance`."""

    def __init__(self, numtapes: int, tapewidth: float, tapedistance: float):
        self.numtapes = numtapes
        self.tapewidth = tapewidth
        self.tapedistance = tapedistance
        h = tapedistance * (numtapes - 1)
        m = np.linspace(-0.5 * h, 0.5 * h, numtapes)
        self.leftpoints = np.zeros((numtapes, 2))
        self.rightpoints = np.zeros((numtapes, 2))
        for k in range(numtapes):
            self.leftpoints[k] = (-0.5 * tapewidth, m[k])
            self.rightpoints[k] = (0.5 * tapewidth, m[k])


class Tape:
    def __init__(self, basecurve: Basecurve, cross_section: CrossSection, index: int, resolution: float = 5.0):
        self.index = index
        self.id = index + 1
        self.basecurve = basecurve
        self.cross_section = cross_section
        self.resolution = resolution
        self.points_left = []
        self.points_right = []
        self.curves = []
        self.curveloops = []
        self.surfaces = []
        self._make_points()
        self._make_curves()
        self._make_surface()

    def _make_points(self):
        pl = np.array([*self.cross_section.leftpoints[self.index], 0.0])
        pr = np.array([*self.cross_section.rightpoints[self.index], 0.0])
        for t in self.basecurve.t:
            m = self.basecurve.r(t)
            R = self.basecurve.transform(t)
            p = m + R @ pl
            q = m + R @ pr
            self.points_left.append(Point(p[0], p[1], p[2], self.resolution))
            self.points_right.append(Point(q[0], q[1], q[2], self.resolution))

    def _make_curves(self):
        F = Curve("Line")
        F.points = [self.points_left[0], self.points_right[0]]
        R = Curve("Spline")
        R.points = list(self.points_right)
        B = Curve("Line")
        B.points = [self.points_right[-1], self.points_left[-1]]
        L = Curve("Spline")
        L.points = list(reversed(self.points_left))
        self.curves = [F, R, B, L]

    def _make_surface(self):
        loop = CurveLoop()
        loop.curves = list(self.curves)
        self.curveloops.append(loop)
        S = Surface()
        S.loops.append(loop)
        self.surfaces.append(S)


class TapeBlock:
    def __init__(self, bottom: Tape, top: Tape):
        self.bottom = bottom
        self.top = top
        self.curves = []
        self.curveloops = []
        self.surfaces = []
        self.surfaceloops = []
        self.volumes = []
        self._make_surfaces()

    @staticmethod
    def _line(p, q):
        c = Curve("Line")
        c.points = [p, q]
        return c

    @staticmethod
    def _loop(curves, signs):
        loop = CurveLoop()
        loop.curves = list(curves)
        loop.signs = list(signs)
        return loop

    def _make_surfaces(self):
        bot, top = self.bottom, self.top
        l0 = self._line(bot.points_left[0], top.points_left[0])
        r0 = self._line(bot.points_right[0], top.points_right[0])
        l1 = self._line(bot.points_left[-1], top.points_left[-1])
        r1 = self._line(bot.points_right[-1], top.points_right[-1])
        self.curves = [l0, r0, l1, r1]

        loops = [
            (self._loop([bot.curves[0], r0, top.curves[0], l0], [1, 1, -1, -1]), True),
            (self._loop([bot.curves[2], l1, top.curves[2], r1], [1, 1, -1, -1]), True),
            (self._loop([l0, bot.curves[3], l1, top.curves[3]], [1, 1, -1, -1]), False),
            (self._loop([r0, top.curves[1], r1, bot.curves[1]], [1, 1, -1, -1]), False),
        ]
        for loop, plane in loops:
            self.curveloops.append(loop)
            S = Surface()
            S.is_plane = plane
            S.loops.append(loop)
            self.surfaces.append(S)

        SL = SurfaceLoop()
        SL.surfaces = list(bot.surfaces) + list(top.surfaces) + list(self.surfaces)
        V = Volume()
        V.loops.append(SL)
        self.surfaceloops.append(SL)
        self.volumes.append(V)


class Geometry:
    """Collects tapes and blocks along a sampled base curve and writes a Gmsh .geo file."""

    def __init__(self, basecurve: Basecurve, cross_section: CrossSection, resolution: float = 5.0):
        if basecurve.t is None:
            raise ValueError("sample the base curve first (curve.sample(n) or curve.sample_spacing(ds))")
        self.basecurve = basecurve
        self.cross_section = cross_section
        self.tapes = [Tape(basecurve, cross_section, k, resolution) for k in range(cross_section.numtapes)]
        self.tape_blocks = [TapeBlock(self.tapes[k - 1], self.tapes[k]) for k in range(1, cross_section.numtapes)]

    def save(self, path: str):
        lines = []
        pid = 0
        for tape in self.tapes:
            for p in tape.points_left + tape.points_right:
                pid += 1
                p.id = pid
                lines.append(p.write())
        cid = clid = sid = slid = vid = 0
        curves, loops, surfs, sloops, vols = [], [], [], [], []
        for obj in self.tapes + self.tape_blocks:
            for c in obj.curves:
                cid += 1; c.id = cid; curves.append(c)
            for l in obj.curveloops:
                clid += 1; l.id = clid; loops.append(l)
            for s in obj.surfaces:
                sid += 1; s.id = sid; surfs.append(s)
            for l in getattr(obj, "surfaceloops", []):
                slid += 1; l.id = slid; sloops.append(l)
            for v in getattr(obj, "volumes", []):
                vid += 1; v.id = vid; vols.append(v)
        for group in (curves, loops, surfs, sloops, vols):
            lines.extend(o.write() for o in group)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return path
