# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Structured shell strips: tapes, gaps and the domain tube, generated in
the unrolled (s, z) plane and mapped through the stadium. Every node sits
exactly on its layer's stadium.

Perimeter layout per layer ("slots"): tape k occupies n_t uniform
intervals, the gap after it n_g intervals; slot positions are arc
lengths at z = 0. The node at (row i, slot j) has arc length
sigma * (slot_j + shift_i) with shift_i = numPitches * P * i / nz, so
tape and gap strips share their seam nodes by construction, and the last
row is the first row advanced by a whole number of tape positions.

Rows sit at constant z. A strip quad therefore has one side along the
perimeter and one along the tape path, and the end planes z = 0 and
z = length are rows of the mesh: the periodic cut never slices an
element.
"""

import math

import numpy as np


class ShellGenerator:

    def __init__(self, winding, domain_radius: float,
                 domain_resolution: float):
        self.w = winding
        self.domain_radius = domain_radius
        self.domain_resolution = domain_resolution

    # ------------------------------------------------------------ helpers

    def layer_slots(self, l: int):
        """Slot arc lengths (z = 0, before the sense is applied) for layer
        l and the per-feature slot ranges: returns (slots, features) where
        features is a list of ("tape"|"gap", tape_index, slot_start,
        slot_count)."""
        w = self.w
        L = w.layers[l]
        nt = w.n_s(l, "tape")
        ng = w.n_s(l, "gap")
        slots = []
        features = []
        for k, T in enumerate(L.tapes):
            start = len(slots)
            for j in range(nt):
                slots.append(T.s0 + j * L.tape_span / nt)
            features.append(("tape", k, start, nt))
            start = len(slots)
            for j in range(ng):
                slots.append(T.s1 + j * L.gap_span / ng)
            features.append(("gap", k, start, ng))
        return np.array(slots), features

    def end_slot_shift(self, l: int) -> int:
        """Slot index shift between a layer's last row and its first: ring
        slot j at z = length sits at the z = 0 position of slot
        j + q * (slots per sector), q = winding.tape_shift(), the same for
        both winding senses."""
        q = self.w.tape_shift()
        return q * (self.w.n_s(l, "tape") + self.w.n_s(l, "gap"))

    def z_rows(self):
        nz = self.w.z_partition_count()
        return np.linspace(0.0, self.w.length, nz + 1)

    def tube_rows(self):
        ns = max(2, int(math.ceil(self.w.length / self.domain_resolution)))
        return np.linspace(0.0, self.w.length, ns + 1)

    def tube_slot_count(self):
        n = int(math.ceil(2.0 * math.pi * self.domain_radius
                          / (4.0 * self.domain_resolution)))
        return 4 * max(1, n)

    # -------------------------------------------------------------- build

    def build(self, builder):
        w = self.w
        zrows = self.z_rows()
        nz = len(zrows) - 1

        for l, layer in enumerate(w.layers):
            slots, features = self.layer_slots(l)
            nslots = len(slots)
            sigma = layer.sigma
            P = layer.perimeter
            st = layer.stadium

            rows = []
            for i, z in enumerate(zrows):
                shift = w.numPitches * P * (float(i) / nz)
                x, y = st.points(sigma * (slots + shift))
                ids = [builder.add_node(float(x[j]), float(y[j]), float(z))
                       for j in range(nslots)]
                rows.append(ids)
            builder.rings[l] = rows

            for kind, k, start, count in features:
                tag = builder.new_surface(kind, layer=l, tape=k)
                for i in range(nz):
                    a_row, b_row = rows[i], rows[i + 1]
                    axis = np.array([0.0, 0.0, 0.5 * (zrows[i] + zrows[i + 1])])
                    for j in range(start, start + count):
                        j1 = (j + 1) % nslots
                        a, b = a_row[j], a_row[j1]
                        c, d = b_row[j], b_row[j1]
                        # split along the short diagonal b-c: the quad is
                        # sheared by the tape advance per row, so a-d is
                        # the long one. Outward orientation: the arc
                        # length runs clockwise, so a positive-sense slot
                        # walk crossed with +z points inward and the
                        # triangles are wound the other way; sigma flips
                        # the handedness
                        if sigma > 0:
                            builder.add_tri(tag, a, c, b)
                            builder.add_tri(tag, b, c, d)
                        else:
                            builder.add_tri(tag, a, b, c)
                            builder.add_tri(tag, b, d, c)
                        self._assert_outward(builder, axis, -2)
                        self._assert_outward(builder, axis, -1)

        # ---------------------------------------------------- domain tube
        rows_z = self.tube_rows()
        nslots = self.tube_slot_count()
        phi = np.arange(nslots) * (2.0 * math.pi / nslots)
        R_dom = self.domain_radius
        xs = R_dom * np.cos(phi)
        ys = R_dom * np.sin(phi)

        tag = builder.new_surface("tube")
        tube_rows = []
        for z in rows_z:
            ids = [builder.add_node(float(xs[j]), float(ys[j]), float(z))
                   for j in range(nslots)]
            tube_rows.append(ids)
        builder.rings["dom"] = tube_rows

        for i in range(len(rows_z) - 1):
            a_row, b_row = tube_rows[i], tube_rows[i + 1]
            axis = np.array([0.0, 0.0, 0.5 * (rows_z[i] + rows_z[i + 1])])
            for j in range(nslots):
                j1 = (j + 1) % nslots
                a, b = a_row[j], a_row[j1]
                c, d = b_row[j], b_row[j1]
                builder.add_tri(tag, a, b, d)
                builder.add_tri(tag, a, d, c)
                self._assert_outward(builder, axis, -2)
                self._assert_outward(builder, axis, -1)

    def _assert_outward(self, builder, axis, which):
        """The triangle added at builder.tris[which] must face away from
        the axis point of its row pair. The stadium is convex, so this is
        a valid outwardness test for the strips as well as the tube."""
        tag, a, b, c = builder.tris[which]
        p1 = builder.nodes[a - 1]
        p2 = builder.nodes[b - 1]
        p3 = builder.nodes[c - 1]
        n = np.cross(p2 - p1, p3 - p1)
        radial = (p1 + p2 + p3) / 3.0 - axis
        if float(np.dot(n, radial)) <= 0.0:
            raise AssertionError(
                "inward-facing shell triangle on surface {:d}".format(tag))

    # -------------------------------------------------------- validations

    def check_surface(self, builder, tol=1e-9):
        """Every ring node must sit exactly on its layer's stadium."""
        worst = 0.0
        for l, layer in enumerate(self.w.layers):
            for row in builder.rings[l]:
                p = np.array([builder.nodes[nid - 1] for nid in row])
                d = layer.stadium.distance(p[:, 0], p[:, 1])
                worst = max(worst, float(np.abs(d).max()))
        if worst > tol:
            raise AssertionError(
                "shell surface error {:.3e} exceeds {:.1e}".format(worst, tol))
        return worst
