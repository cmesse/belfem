# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""End caps: meshed ONCE in the plane, then instantiated at z = 0 and,
translated by the cable length, at z = length. Cap periodicity holds
node-for-node by construction.

The layer rings are polylines through the exact slot positions of the
z = 0 row (one gmsh Point per slot, one two-node Line per interval), so
the cap boundary nodes ARE the strip ring nodes, whatever the slot
spacing does across the flat/round junctions. The domain circle keeps
the transfinite quarter arcs. Boundary nodes are snapped to their exact
positions after the gmsh 2D run and REUSE the ring node ids of the
strips, so the shell is watertight with a single global node id space.
"""

import math
import os
import subprocess

import numpy as np

from rutherford.mshio import read_msh2


class CapError(Exception):
    pass


class CapGenerator:

    def __init__(self, shellgen, inner_resolution: float, workdir: str):
        self.sg = shellgen
        self.w = shellgen.w
        self.inner_resolution = inner_resolution
        self.workdir = workdir

    # ------------------------------------------------------------ geometry

    def build_cap_geo(self, path):
        """Write the planar cap .geo. Returns metadata needed to interpret
        the mesh: ring slot positions and surface ids."""
        w = self.w
        sg = self.sg
        meta = {"rings": [], "surfaces": {}}
        next_ids = {"point": 1, "curve": 1, "loop": 1, "surface": 1}

        lines = []

        def point(x, y, res):
            pid = next_ids["point"]
            next_ids["point"] += 1
            lines.append("Point({:d}) = {{{:.17g}, {:.17g}, 0, {:.17g}}};".format(
                pid, x, y, res))
            return pid

        center = point(0.0, 0.0, self.inner_resolution)

        ring_loops = []
        for l, layer in enumerate(w.layers):
            slots, features = sg.layer_slots(l)
            nslots = len(slots)
            x, y = layer.stadium.points(layer.sigma * slots)
            slot_xy = np.column_stack((x, y))
            meta["rings"].append({"stadium": layer.stadium, "xy": slot_xy})

            sizes = np.full(nslots, w.tapeResolution)
            for kind, k, start, count in features:
                if kind == "gap":
                    sizes[start:start + count] = w.gapResolution
            pids = [point(slot_xy[j][0], slot_xy[j][1], float(sizes[j]))
                    for j in range(nslots)]

            curves = []
            for j in range(nslots):
                cid = next_ids["curve"]
                next_ids["curve"] += 1
                lines.append("Line({:d}) = {{{:d}, {:d}}};".format(
                    cid, pids[j], pids[(j + 1) % nslots]))
                lines.append("Transfinite Curve {{{:d}}} = 2;".format(cid))
                curves.append(cid)
            # loops run counter-clockwise; the arc length runs clockwise,
            # so the slot walk of a positive-sense layer is reversed
            if layer.sigma > 0:
                curves = [-c for c in reversed(curves)]
            lid = next_ids["loop"]
            next_ids["loop"] += 1
            lines.append("Curve Loop({:d}) = {{{:s}}};".format(
                lid, ", ".join(str(c) for c in curves)))
            ring_loops.append(lid)

        # domain circle: quarter arcs, nslots divisible by 4
        ndom = sg.tube_slot_count()
        R = sg.domain_radius
        dom_angles = np.arange(ndom) * (2.0 * math.pi / ndom)
        dom_xy = np.column_stack((R * np.cos(dom_angles), R * np.sin(dom_angles)))
        meta["dom"] = {"radius": R, "xy": dom_xy}
        qids = []
        for q in range(4):
            k = q * ndom // 4
            qids.append(point(dom_xy[k][0], dom_xy[k][1],
                              self.sg.domain_resolution))
        dom_curves = []
        for q in range(4):
            cid = next_ids["curve"]
            next_ids["curve"] += 1
            lines.append("Circle({:d}) = {{{:d}, {:d}, {:d}}};".format(
                cid, qids[q], center, qids[(q + 1) % 4]))
            lines.append("Transfinite Curve {{{:d}}} = {:d};".format(
                cid, ndom // 4 + 1))
            dom_curves.append(cid)
        lid = next_ids["loop"]
        next_ids["loop"] += 1
        lines.append("Curve Loop({:d}) = {{{:s}}};".format(
            lid, ", ".join(str(c) for c in dom_curves)))
        dom_loop = lid

        # plane surfaces: core section, annuli, outer
        sid = next_ids["surface"]
        lines.append("Plane Surface({:d}) = {{{:d}}};".format(sid, ring_loops[0]))
        lines.append("Point{{{:d}}} In Surface{{{:d}}};".format(center, sid))
        meta["surfaces"]["disk"] = sid
        next_ids["surface"] += 1

        annuli = []
        for l in range(len(w.layers) - 1):
            sid = next_ids["surface"]
            next_ids["surface"] += 1
            lines.append("Plane Surface({:d}) = {{{:d}, {:d}}};".format(
                sid, ring_loops[l + 1], ring_loops[l]))
            annuli.append(sid)
        meta["surfaces"]["annuli"] = annuli

        sid = next_ids["surface"]
        next_ids["surface"] += 1
        lines.append("Plane Surface({:d}) = {{{:d}, {:d}}};".format(
            sid, dom_loop, ring_loops[-1]))
        meta["surfaces"]["outer"] = sid

        lines.append("Mesh.SaveAll = 1;")
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return meta

    # ------------------------------------------------------------- meshing

    def mesh_cap(self):
        geo = os.path.join(self.workdir, "rutherford_cap.geo")
        out = os.path.join(self.workdir, "rutherford_cap.msh")
        meta = self.build_cap_geo(geo)
        r = subprocess.run(
            ["gmsh", "-2", "-format", "msh22", geo, "-o", out],
            capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            raise CapError("cap 2D meshing failed:\n" + r.stderr[-2000:])
        nodes, elements = read_msh2(out)
        return self._interpret(meta, nodes, elements)

    def _interpret(self, meta, nodes, elements):
        """Snap boundary nodes to exact slot positions and classify every
        cap node: ('ring', layer, slot) / ('dom', slot) / ('center',) /
        ('interior',). Returns dict with node arrays and tris per piece."""
        ids = sorted(nodes.keys())
        xy = np.array([[nodes[i][0], nodes[i][1]] for i in ids])
        id_index = {nid: k for k, nid in enumerate(ids)}

        kind = [("interior",)] * len(ids)
        claimed = {}

        def claim(target_xy, on_curve, tag_fn, ring_label):
            for slot, exy in enumerate(target_xy):
                # only nodes on this ring's curve are candidates, so a
                # close interior node can never steal a boundary slot
                d = np.where(on_curve, np.linalg.norm(xy - exy, axis=1),
                             np.inf)
                k = int(np.argmin(d))
                if d[k] > 1e-6:
                    raise CapError(
                        "no cap node within 1e-6 of {} slot {:d} "
                        "(nearest {:.3e})".format(ring_label, slot, d[k]))
                if k in claimed:
                    raise CapError(
                        "cap node {:d} claimed twice ({} slot {:d} and "
                        "{})".format(ids[k], ring_label, slot, claimed[k]))
                claimed[k] = "{} {:d}".format(ring_label, slot)
                kind[k] = tag_fn(slot)
                xy[k] = exy   # snap to exact position

        for l, ring in enumerate(meta["rings"]):
            on_curve = np.abs(ring["stadium"].distance(xy[:, 0], xy[:, 1])) < 1e-6
            claim(ring["xy"], on_curve, lambda slot, l=l: ("ring", l, slot),
                  "ring {:d}".format(l))
        radii = np.linalg.norm(xy, axis=1)
        claim(meta["dom"]["xy"], np.abs(radii - meta["dom"]["radius"]) < 1e-6,
              lambda slot: ("dom", slot), "dom")

        d0 = np.linalg.norm(xy, axis=1)
        kc = int(np.argmin(d0))
        if d0[kc] > 1e-6:
            raise CapError("no cap node at the center")
        xy[kc] = (0.0, 0.0)
        kind[kc] = ("center",)

        tris = {"disk": [], "annuli": [[] for _ in meta["surfaces"]["annuli"]],
                "outer": []}
        for (etype, phys, geo, nn) in elements:
            if etype != 2:
                continue
            tri = tuple(id_index[i] for i in nn)
            if geo == meta["surfaces"]["disk"]:
                tris["disk"].append(tri)
            elif geo == meta["surfaces"]["outer"]:
                tris["outer"].append(tri)
            else:
                for a, sid in enumerate(meta["surfaces"]["annuli"]):
                    if geo == sid:
                        tris["annuli"][a].append(tri)
                        break

        if not tris["disk"] or not tris["outer"]:
            raise CapError("cap mesh misses disk or outer piece")
        return {"xy": xy, "kind": kind, "tris": tris}

    # -------------------------------------------------------- instantiation

    def instantiate(self, builder, cap):
        """Create front and back cap instances in the builder, wiring
        boundary nodes to the strip ring ids. Fills builder.cap_pairs,
        the cap surface tags and builder.vertex_nodes."""
        sg = self.sg
        L = self.w.length

        xy = cap["xy"]
        kind = cap["kind"]
        n = len(xy)

        front_id = np.zeros(n, dtype=int)
        back_id = np.zeros(n, dtype=int)

        worst_reuse = 0.0
        for k in range(n):
            kd = kind[k]
            x, y = xy[k]
            pf = np.array([x, y, 0.0])
            pb = np.array([x, y, L])
            if kd[0] == "ring":
                _, l, slot = kd
                front_id[k] = builder.rings[l][0][slot]
                # the last row is the first row advanced by a whole number
                # of tape positions, so the node under planar slot k is
                # ring slot k - shift (the reuse check below verifies it)
                nring = len(builder.rings[l][-1])
                back_id[k] = builder.rings[l][-1][(slot - sg.end_slot_shift(l)) % nring]
            elif kd[0] == "dom":
                slot = kd[1]
                front_id[k] = builder.rings["dom"][0][slot]
                back_id[k] = builder.rings["dom"][-1][slot]
            else:
                front_id[k] = builder.add_node(*pf)
                back_id[k] = builder.add_node(*pb)
                continue
            # reused ring ids must BE the planar positions: the geometric
            # closure of the whole cap construction
            worst_reuse = max(
                worst_reuse,
                float(np.linalg.norm(builder.nodes[front_id[k] - 1] - pf)),
                float(np.linalg.norm(builder.nodes[back_id[k] - 1] - pb)))
        if worst_reuse > 1e-9:
            raise CapError(
                "cap boundary reuse mismatch: ring nodes deviate {:.3e} "
                "from the planar positions (closure inconsistency)".format(
                    worst_reuse))

        builder.cap_pairs = [(int(front_id[k]), int(back_id[k])) for k in range(n)]

        # surfaces: front pieces face -z, back pieces +z; every triangle's
        # normal is asserted against the cap plane normal (the 2D mesher
        # does not guarantee uniform winding per piece)
        ez = np.array([0.0, 0.0, 1.0])

        def emit(piece_tris, role, side, meta_extra):
            tag = builder.new_surface(role, side=side, **meta_extra)
            ids = front_id if side == "front" else back_id
            outward = -ez if side == "front" else ez
            for (a, b, c) in piece_tris:
                p1 = builder.nodes[ids[a] - 1]
                p2 = builder.nodes[ids[b] - 1]
                p3 = builder.nodes[ids[c] - 1]
                proj = float(np.dot(np.cross(p2 - p1, p3 - p1), outward))
                area2 = float(np.linalg.norm(np.cross(p2 - p1, p3 - p1)))
                if area2 <= 0.0 or abs(proj) < 1e-6 * area2:
                    raise CapError(
                        "degenerate or tilted cap triangle on {} {}".format(
                            role, side))
                if proj > 0.0:
                    builder.add_tri(tag, int(ids[a]), int(ids[b]), int(ids[c]))
                else:
                    builder.add_tri(tag, int(ids[a]), int(ids[c]), int(ids[b]))
            return tag

        for side in ("front", "back"):
            emit(cap["tris"]["disk"], "cap_disk", side, {})
            for aidx, atris in enumerate(cap["tris"]["annuli"]):
                emit(atris, "cap_annulus", side, {"annulus": aidx})
            emit(cap["tris"]["outer"], "cap_outer", side, {})

        # the ten vertices: front center, 4 front quarter points, back
        # center, 4 back quarter points (quarter points live on the tube)
        kc = next(k for k in range(n) if kind[k][0] == "center")
        ndom = sg.tube_slot_count()
        quarters = [0, ndom // 4, ndom // 2, 3 * ndom // 4]
        builder.vertex_nodes = (
            [int(front_id[kc])]
            + [builder.rings["dom"][0][q] for q in quarters]
            + [int(back_id[kc])]
            + [builder.rings["dom"][-1][q] for q in quarters])
