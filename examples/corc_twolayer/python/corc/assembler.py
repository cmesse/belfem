"""Volume meshing over the discrete shell boundary: writes the shell msh +
volume .geo (Surface Loops sharing the strip surfaces, Spike A/B recipe),
runs gmsh -3, reads the result back, maps our node ids through gmsh's
renumbering and restores the exact analytic shell coordinates.

Volume regions (same decomposition as the old generator):
  1        inner        cap disks + layer-0 strips
  2..nL    interlayer   cap annuli + adjacent layer strips
  nL+1     outer        cap outer pieces + last layer strips + domain tube
"""

import os
import subprocess

import numpy as np

from corc.mshio import read_msh2


class AssemblerError(Exception):
    pass


class Assembler:

    def __init__(self, builder, workdir: str, domain_resolution: float):
        self.b = builder
        self.workdir = workdir
        self.domain_resolution = domain_resolution

    # ------------------------------------------------------------- volumes

    def volume_loops(self):
        b = self.b
        surf = b.surfaces

        def tags(role, **match):
            out = []
            for t, m in surf.items():
                if m["role"] != role:
                    continue
                if all(m.get(k) == v for k, v in match.items()):
                    out.append(t)
            if not out:
                raise AssemblerError(
                    "no surfaces of role '{}' ({}) — Surface Loop would be "
                    "open".format(role, match))
            return out

        strip_layers = [m["layer"] for m in surf.values()
                        if m["role"] in ("tape", "gap")]
        if not strip_layers:
            raise AssemblerError("builder holds no strip surfaces")
        nlayers = max(strip_layers) + 1

        def layer_strips(l):
            return tags("tape", layer=l) + tags("gap", layer=l)

        loops = []
        loops.append(tags("cap_disk", side="front")
                     + tags("cap_disk", side="back") + layer_strips(0))
        for l in range(nlayers - 1):
            loops.append(tags("cap_annulus", side="front", annulus=l)
                         + tags("cap_annulus", side="back", annulus=l)
                         + layer_strips(l) + layer_strips(l + 1))
        loops.append(tags("cap_outer", side="front")
                     + tags("cap_outer", side="back")
                     + layer_strips(nlayers - 1) + tags("tube"))
        return loops

    # ----------------------------------------------------------------- run

    def run(self):
        b = self.b
        shell = os.path.join(self.workdir, "corc_shell.msh")
        geo = os.path.join(self.workdir, "corc_volume.geo")
        out = os.path.join(self.workdir, "corc_volume.msh")
        b.write_shell_msh(shell)
        loops = self.volume_loops()

        with open(geo, "w") as f:
            f.write('Merge "corc_shell.msh";\n')
            for i, loop in enumerate(loops):
                f.write("Surface Loop({:d}) = {{{:s}}};\n".format(
                    i + 1, ", ".join(str(t) for t in loop)))
                f.write("Volume({:d}) = {{{:d}}};\n".format(i + 1, i + 1))
            f.write("Mesh.MeshSizeExtendFromBoundary = 1;\n")
            f.write("Mesh.MeshSizeMax = {:.17g};\n".format(self.domain_resolution))
            # the Netgen pass removes most front-collision slivers the
            # default optimizer leaves behind (p5 quality ~10x better)
            f.write("Mesh.OptimizeNetgen = 1;\n")
            f.write("Mesh.SaveAll = 1;\n")

        r = subprocess.run(
            ["gmsh", "-3", "-format", "msh22", "corc_volume.geo", "-o",
             os.path.basename(out)],
            cwd=self.workdir, text=True, timeout=3600)

        if r.returncode != 0:
            raise AssemblerError("volume meshing failed:\n" + r.stderr[-3000:]
                                 + "\n" + r.stdout[-2000:])

        return self._read_back(out, len(loops))

    # ------------------------------------------------------------ read-back

    def _read_back(self, path, nvolumes):
        """Returns (nodes, tets, shell_tris, old_to_gmsh):
        nodes {gmsh_id: [x,y,z]} with shell coordinates restored exactly,
        tets [(region, [n1..n4])], shell_tris [(surface_tag, [n1..n3])]
        mapped to gmsh node ids, old_to_gmsh mapping for our builder ids."""
        b = self.b
        nodes, elements = read_msh2(path)

        # match our shell nodes to gmsh's (renumbering-safe): hash into
        # 1e-9 bins but probe neighboring bins too, so nodes straddling a
        # bin boundary still match; verify with a true distance check
        lookup = {}
        for gid, p in nodes.items():
            key = (round(p[0], 9), round(p[1], 9), round(p[2], 9))
            lookup.setdefault(key, []).append(gid)

        eps = 1e-9

        def find(p):
            best = None
            best_d = eps   # strict 1e-9 preservation policy
            for dx in (0.0, -eps, eps):
                for dy in (0.0, -eps, eps):
                    for dz in (0.0, -eps, eps):
                        key = (round(p[0] + dx, 9), round(p[1] + dy, 9),
                               round(p[2] + dz, 9))
                        for gid in lookup.get(key, ()):  # noqa: B007
                            q = nodes[gid]
                            d = max(abs(p[0] - q[0]), abs(p[1] - q[1]),
                                    abs(p[2] - q[2]))
                            if d < best_d:
                                best = gid
                                best_d = d
            return best

        old_to_gmsh = {}
        missing = 0
        for i, p in enumerate(b.nodes):
            gid = find(p)
            if gid is None:
                missing += 1
            else:
                old_to_gmsh[i + 1] = gid
        if missing:
            raise AssemblerError(
                "{:d} shell nodes were not preserved by the volume run "
                "(boundary must stay fixed)".format(missing))
        if len(set(old_to_gmsh.values())) != len(old_to_gmsh):
            raise AssemblerError("shell node matching is not injective")

        # restore exact analytic coordinates on the shell nodes
        nodes = {gid: np.array(p) for gid, p in nodes.items()}
        for old, gid in old_to_gmsh.items():
            nodes[gid] = b.nodes[old - 1].copy()

        tets = [(g, nn) for (t, p, g, nn) in elements if t == 4]
        regions = sorted(set(g for g, _ in tets))
        if regions != list(range(1, nvolumes + 1)):
            raise AssemblerError(
                "expected volume regions 1..{:d}, got {}".format(
                    nvolumes, regions))

        # per-surface triangle counts must be untouched (a same-total
        # reclassification would silently corrupt the physical tags)
        shell_tris = [(g, nn) for (t, p, g, nn) in elements if t == 2]
        got = {}
        for g, _ in shell_tris:
            got[g] = got.get(g, 0) + 1
        want = {}
        for g, *_ in b.tris:
            want[g] = want.get(g, 0) + 1
        if got != want:
            diff = {g: (want.get(g, 0), got.get(g, 0))
                    for g in set(want) | set(got)
                    if want.get(g, 0) != got.get(g, 0)}
            raise AssemblerError(
                "per-surface triangle counts changed in the volume run "
                "(surface: expected, got): {}".format(diff))

        # full connectivity identity per surface: the multiset of node-id
        # triples must survive unchanged (counts alone would let a
        # same-count retessellation slip through)
        want_tris = {}
        for g, a, bb, c in b.tris:
            want_tris.setdefault(g, set()).add(frozenset(
                (old_to_gmsh[a], old_to_gmsh[bb], old_to_gmsh[c])))
        got_tris = {}
        for g, nn in shell_tris:
            got_tris.setdefault(g, set()).add(frozenset(nn))
        if got_tris != want_tris:
            bad = [g for g in want_tris if got_tris.get(g) != want_tris[g]]
            raise AssemblerError(
                "shell connectivity changed on surfaces {}".format(bad))

        return nodes, tets, shell_tris, old_to_gmsh
