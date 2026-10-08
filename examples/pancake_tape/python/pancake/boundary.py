# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Domain boundary: a structured tube on r = R_dom whose slots include
the lead's exit angle and the lead bar's exit angle, and whose rows
include the pack levels, so the lead's end edge and the bar's end
column ARE tube nodes and the bar's end face is the output terminal;
two caps meshed by gmsh in 2D, the one on the lead side with a hole
where the arm bar's end face (the input terminal) sits."""

import math
import os
import subprocess

import numpy as np

from pancake.mshio import read_msh2
from pancake.pack import Pack


class BoundaryError(Exception):
    pass


class Boundary:

    def __init__(self, sheet, pack, domain_resolution, workdir):
        self.S = sheet
        self.P = pack
        self.h = float(domain_resolution)
        self.workdir = workdir
        self.xyz = pack.xyz            # grows as nodes are added
        self.surfaces = []             # (tag, tris)

    def _add_nodes(self, pts):
        start = self.xyz.shape[0] + 1
        self.xyz = np.vstack((self.xyz, pts))
        return np.arange(start, start + pts.shape[0])

    # ---------------------------------------------------------------- tube

    def build_tube(self):
        S = self.S
        ns = max(8, int(math.ceil(2.0 * math.pi * S.Rdom / self.h)))
        step = 2.0 * math.pi / ns
        if S.barExitAngle - S.exitAngle >= step:
            raise BoundaryError(
                "outerBarThickness {:g} mm is not below the tube slot spacing "
                "{:.3g} mm; lower it or domainResolution".format(
                    S.tLead, S.Rdom * step))
        phi = np.concatenate(([S.exitAngle, S.barExitAngle],
                              S.exitAngle + np.arange(1, ns) * step))
        self.phi = phi
        ns1 = len(phi)
        nb = max(1, int(math.ceil((S.H - 0.5 * S.b) / self.h)))
        below = np.linspace(-S.H, -0.5 * S.b, nb + 1)[:-1]
        above = np.linspace(0.5 * S.b, S.H, nb + 1)[1:]
        rows = np.concatenate((below, S.zPack, above))
        self.rows = rows
        nr = len(rows)
        ids = np.zeros((nr, ns1), dtype=int)
        lead = self.P.lead_end_ids()
        bar = self.P.lead_bar_end_ids()
        r0 = nb
        for i, z in enumerate(rows):
            pts = np.column_stack((S.Rdom * np.cos(phi), S.Rdom * np.sin(phi),
                                   np.full(ns1, z)))
            if r0 <= i < r0 + S.nZ + 1:
                ids[i, 0] = lead[i - r0]
                ids[i, 1] = bar[i - r0]
                ids[i, 2:] = self._add_nodes(pts[2:])
            else:
                ids[i] = self._add_nodes(pts)
        for column in (lead, bar):
            worst = np.abs(np.linalg.norm(self.xyz[column - 1][:, :2], axis=1) - S.Rdom).max()
            if worst > 1e-9:
                raise BoundaryError("a lead column is off the domain cylinder")
        self.tubeIds = ids
        a = ids[:-1, :]
        b = np.roll(ids[:-1, :], -1, axis=1)
        c = ids[1:, :]
        d = np.roll(ids[1:, :], -1, axis=1)
        keep = np.ones(a.shape, dtype=bool)
        keep[r0:r0 + S.nZ, 0] = False      # the terminal patch comes from the pack
        tris = np.vstack((np.column_stack((a[keep], b[keep], d[keep])),
                          np.column_stack((a[keep], d[keep], c[keep]))))
        self.surfaces.append((Pack.TUBE, tris))

    # ---------------------------------------------------------------- caps

    def _cap_geo(self, path, hole_xy):
        S = self.S
        lines = []
        pid = [0]

        def point(x, y, h):
            pid[0] += 1
            lines.append("Point({:d}) = {{{:.17g}, {:.17g}, 0, {:.17g}}};".format(
                pid[0], x, y, h))
            return pid[0]

        def loop(pts, h, first_line):
            ids = [point(p[0], p[1], h) for p in pts]
            n = len(ids)
            for i in range(n):
                lines.append("Line({:d}) = {{{:d}, {:d}}};".format(
                    first_line + i, ids[i], ids[(i + 1) % n]))
                lines.append("Transfinite Curve {{{:d}}} = 2;".format(first_line + i))
            return list(range(first_line, first_line + n))

        ring = np.column_stack((S.Rdom * np.cos(self.phi), S.Rdom * np.sin(self.phi)))
        outer = loop(ring, self.h, 1)
        lines.append("Curve Loop(1) = {{{}}};".format(", ".join(map(str, outer))))
        if hole_xy is not None:
            inner = loop(hole_xy, S.leadResolution, 1000)
            lines.append("Curve Loop(2) = {{{}}};".format(", ".join(map(str, inner))))
            lines.append("Plane Surface(1) = {1, 2};")
        else:
            lines.append("Plane Surface(1) = {1};")
        lines.append("Mesh.MeshSizeFromPoints = 1;")
        lines.append("Mesh.MeshSizeExtendFromBoundary = 1;")
        lines.append("Mesh.SaveAll = 1;")
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")

    def build_caps(self):
        S = self.S
        outline = self.P.hole_outline_ids()
        for side in ("top", "bottom"):
            with_hole = (side == "top") == (S.sigma > 0)
            geo = os.path.join(self.workdir, "pancake_cap_{}.geo".format(side))
            out = os.path.join(self.workdir, "pancake_cap_{}.msh".format(side))
            self._cap_geo(geo, self.xyz[outline - 1][:, :2] if with_hole else None)
            r = subprocess.run(["gmsh", "-2", "-format", "msh22", geo, "-o", out],
                               capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                raise BoundaryError("cap meshing failed:\n" + r.stderr[-2000:])
            gids, gxyz, elements = read_msh2(out)
            self._instantiate(side, outline if with_hole else None, gids, gxyz, elements)

    def _instantiate(self, side, outline, gids, gxyz, elements):
        S = self.S
        z = S.H if side == "top" else -S.H
        ring = self.tubeIds[-1] if side == "top" else self.tubeIds[0]
        n = len(gids)
        final = np.zeros(n, dtype=int)
        claimed = np.zeros(n, dtype=bool)

        rad = np.hypot(gxyz[:, 0], gxyz[:, 1])
        on_ring = np.abs(rad - S.Rdom) < 1e-6
        ang = np.arctan2(gxyz[:, 1], gxyz[:, 0])
        # nearest slot by angle (slots are not uniform)
        diff = np.angle(np.exp(1j * (ang[:, None] - self.phi[None, :])))
        slot = np.argmin(np.abs(diff), axis=1)
        exact = np.column_stack((S.Rdom * np.cos(self.phi[slot]), S.Rdom * np.sin(self.phi[slot])))
        off = np.linalg.norm(gxyz[:, :2] - exact, axis=1)
        if np.any(on_ring & (off > 1e-6)) or on_ring.sum() != len(self.phi):
            raise BoundaryError("cap ring does not match the tube slots")
        final[on_ring] = ring[slot[on_ring]]
        claimed |= on_ring

        if outline is not None:
            pts = self.xyz[outline - 1][:, :2]
            for i in range(pts.shape[0]):
                d = np.linalg.norm(gxyz[:, :2] - pts[i], axis=1)
                d[claimed] = np.inf
                k = int(np.argmin(d))
                if d[k] > 1e-6:
                    raise BoundaryError("cap misses hole outline node {:d}".format(i))
                final[k] = outline[i]
                claimed[k] = True

        rest = ~claimed
        pts = np.column_stack((gxyz[rest, 0], gxyz[rest, 1], np.full(rest.sum(), z)))
        final[rest] = self._add_nodes(pts)

        index = {g: i for i, g in enumerate(gids)}
        tags, conn = elements[2]
        tris = np.vectorize(lambda g: final[index[g]])(conn)
        tag = Pack.CAP_TOP if side == "top" else Pack.CAP_BOTTOM
        self.surfaces.append((tag, tris))
