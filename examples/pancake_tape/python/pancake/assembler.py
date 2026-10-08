# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Air volume: the pack's boundary, the domain boundary and the bare
sheet pieces are written as one discrete shell, gmsh meshes the single
air volume around them (embedding the bare sheet), and the result is
read back with our node ids restored by coordinate matching.

The mesh size comes from a background field, not from the boundary
mesh: the pack's top and bottom faces carry 50 um edges that would
otherwise seed millions of air tets."""

import os
import subprocess

import numpy as np

from pancake.mshio import read_msh2, write_msh2
from pancake.pack import Pack


class AssemblerError(Exception):
    pass


class Assembler:

    def __init__(self, sheet, xyz, surfaces, domain_resolution, workdir,
                 pack_air_resolution, air_growth, optimize_netgen=False,
                 algorithm3d=10, size_surfaces=()):
        self.S = sheet
        self.xyz = xyz
        self.surfaces = dict(surfaces)
        self.h = float(domain_resolution)
        self.hFace = float(pack_air_resolution)
        self.grow = float(air_growth)
        self.workdir = workdir
        self.netgen = optimize_netgen
        self.algorithm3d = int(algorithm3d)
        self.sizeSurfaces = [t for t in size_surfaces if self.surfaces[t].shape[0] > 0]

    def _write_geo(self, path):
        S = self.S
        loop = ", ".join(str(t) for t, tris in self.surfaces.items() if tris.shape[0] > 0)
        f = ("max({res:.17g}, min({dmax:.17g}, {res:.17g} + {grow:.17g} * sqrt("
             "(max(0, sqrt(x*x+y*y) - {ro:.17g}))^2 + "
             "(max(0, {ri:.17g} - sqrt(x*x+y*y)))^2 + "
             "(max(0, abs(z) - {hb:.17g}))^2)))").format(
                 res=self.hFace, dmax=self.h, grow=self.grow, ro=S.ro, ri=S.ri,
                 hb=0.5 * S.b)
        lines = [
            'Merge "pancake_shell.msh";',
            "Surface Loop(1) = {{{}}};".format(loop),
            "Volume(1) = {1};",
            "Field[1] = MathEval;",
            'Field[1].F = "{}";'.format(f),
        ]
        fields = "1"
        if self.sizeSurfaces:
            lines += [
                "Field[2] = Distance;",
                "Field[2].SurfacesList = {{{}}};".format(
                    ", ".join(str(t) for t in self.sizeSurfaces)),
                "Field[2].Sampling = 20;",
                "Field[3] = Threshold;",
                "Field[3].InField = 2;",
                "Field[3].SizeMin = {:.17g};".format(S.leadResolution),
                "Field[3].SizeMax = {:.17g};".format(self.h),
                "Field[3].DistMin = {:.17g};".format(S.leadResolution),
                "Field[3].DistMax = {:.17g};".format(4.0 * self.h),
            ]
            fields = "1, 3"
        lines += [
            "Field[4] = Min;",
            "Field[4].FieldsList = {{{}}};".format(fields),
            "Background Field = 4;",
            "Mesh.MeshSizeExtendFromBoundary = 0;",
            "Mesh.MeshSizeFromPoints = 0;",
            "Mesh.MeshSizeFromCurvature = 0;",
            "Mesh.MeshSizeMax = {:.17g};".format(self.h),
            "Mesh.Algorithm3D = {:d};".format(self.algorithm3d),
            "Mesh.OptimizeNetgen = {:d};".format(1 if self.netgen else 0),
            "Mesh.SaveAll = 1;",
        ]
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    def run(self):
        shell = os.path.join(self.workdir, "pancake_shell.msh")
        geo = os.path.join(self.workdir, "pancake_volume.geo")
        out = os.path.join(self.workdir, "pancake_volume.msh")
        groups = [(2, t, t, tris) for t, tris in self.surfaces.items()
                  if tris.shape[0] > 0]
        write_msh2(shell, self.xyz, groups)
        self._write_geo(geo)
        r = subprocess.run(
            ["gmsh", "-3", "-format", "msh22", "pancake_volume.geo", "-o",
             os.path.basename(out)],
            cwd=self.workdir, capture_output=True, text=True, timeout=7200)
        if r.returncode != 0:
            raise AssemblerError("air meshing failed:\n" + r.stderr[-3000:]
                                 + "\n" + r.stdout[-2000:])
        for line in r.stdout.splitlines():
            if "Error" in line or "Warning" in line:
                print("gmsh:", line.strip())
        return self._read_back(out)

    def _read_back(self, path):
        gids, gxyz, elements = read_msh2(path)

        # our shell nodes -> gmsh ids, renumbering-safe, 1e-9 strict
        eps = 1e-9
        keys = {}
        for i in range(len(gids)):
            p = gxyz[i]
            keys.setdefault((round(p[0], 9), round(p[1], 9), round(p[2], 9)), []).append(i)

        # only nodes on the shell surfaces reach gmsh; the pack's
        # interior sheet nodes stay ours untouched
        n = self.xyz.shape[0]
        used = np.unique(np.concatenate([t.ravel() for t in self.surfaces.values()]))
        mine_to_g = np.full(n, -1, dtype=int)
        for a in used - 1:
            p = self.xyz[a]
            best, best_d = -1, eps
            for dx in (0.0, -eps, eps):
                for dy in (0.0, -eps, eps):
                    for dz in (0.0, -eps, eps):
                        key = (round(p[0] + dx, 9), round(p[1] + dy, 9), round(p[2] + dz, 9))
                        for i in keys.get(key, ()):
                            d = np.abs(gxyz[i] - p).max()
                            if d < best_d:
                                best, best_d = i, d
            mine_to_g[a] = best
        matched = mine_to_g[used - 1]
        if np.any(matched < 0):
            raise AssemblerError("{:d} shell nodes not preserved by the air run".format(
                int((matched < 0).sum())))
        if len(np.unique(matched)) != len(used):
            raise AssemblerError("shell node matching is not injective")

        # final ids: ours keep theirs, new air nodes follow
        g_to_final = np.full(len(gids), -1, dtype=int)
        g_to_final[matched] = used
        new = np.where(g_to_final < 0)[0]
        g_to_final[new] = np.arange(n + 1, n + 1 + len(new))
        xyz = np.vstack((self.xyz, gxyz[new]))

        index = np.full(gids.max() + 1, -1, dtype=int)
        index[gids] = np.arange(len(gids))
        tags, conn = elements[4]
        air = conn[tags[:, 1] == 1]
        if air.shape[0] == 0:
            raise AssemblerError("air run produced no tetrahedra")
        air = g_to_final[index[air]]
        if np.any(air <= 0):
            raise AssemblerError("air tetrahedron references an unknown node")

        # surface triangle counts per tag must survive untouched
        stags, _ = elements[2]
        for t, tris in self.surfaces.items():
            got = int((stags[:, 1] == t).sum())
            if got != tris.shape[0]:
                raise AssemblerError(
                    "surface {:d}: {:d} triangles in, {:d} out".format(t, tris.shape[0], got))
        return xyz, air
