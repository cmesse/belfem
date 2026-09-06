"""Structured shell strips: tapes, gaps and the domain tube, generated
analytically in (s, phi) space and mapped through the frame. Every node
sits exactly at its prescribed radius from the centerline.

Angular layout per layer ("slots"): tape k occupies n_t uniform intervals,
the gap after it n_g intervals; slot angles are body angles at s = 0. The
node at (row i, slot j) has body angle sigma * slot_j + sigma * omega * s_i,
so tape and gap strips share their seam nodes by construction.

The domain tube twists by the common cap rotation theta over the length so
its end ring matches the rotated back cap (theta = 0 for integer turns).
"""

import math

import numpy as np


class ShellGenerator:

    def __init__(self, winding, frame, domain_radius: float,
                 domain_resolution: float):
        self.w = winding
        self.frame = frame
        self.domain_radius = domain_radius
        self.domain_resolution = domain_resolution
        self.theta = winding.cap_rotation()

    # ------------------------------------------------------------ helpers

    def layer_slots(self, l: int):
        """Slot boundary angles (body frame, s=0) for layer l and the
        per-feature slot ranges: returns (angles, features) where features
        is a list of ("tape"|"gap", tape_index, slot_start, slot_count)."""
        w = self.w
        L = w.layers[l]
        nt = w.n_phi(l, "tape")
        ng = w.n_phi(l, "gap")
        angles = []
        features = []
        for k, T in enumerate(L.tapes):
            start = len(angles)
            for j in range(nt):
                angles.append(T.phi0 + j * L.tape_angle / nt)
            features.append(("tape", k, start, nt))
            start = len(angles)
            for j in range(ng):
                angles.append(T.phi1 + j * L.gap_angle / ng)
            features.append(("gap", k, start, ng))
        return np.array(angles), features

    def end_slot_shift(self, l: int) -> int:
        """Ring-slot index shift between a layer's last row and its first:
        with the n-fold closure (winding.tape_shift() = q) ring slot j at
        s=L sits at the s=0 position of slot j + q * (slots per sector),
        the same for both winding senses. 0 when the caps close by rotation
        instead."""
        q = self.w.tape_shift()
        if q is None or self.w.numLayers == 1:
            return 0
        return q * (self.w.n_phi(l, "tape") + self.w.n_phi(l, "gap"))

    def s_rows(self):
        ns = self.w.s_partition_count()
        return np.linspace(0.0, self.w.length, ns + 1)

    def tube_rows(self):
        ns = max(2, int(math.ceil(self.w.length / self.domain_resolution)))
        return np.linspace(0.0, self.w.length, ns + 1)

    def tube_slot_count(self):
        n = int(math.ceil(2.0 * math.pi * self.domain_radius
                          / (4.0 * self.domain_resolution)))
        return 4 * max(1, n)

    # -------------------------------------------------------------- build

    def build(self, builder):
        frame = self.frame
        w = self.w

        # frame samples are shared by all layers (same s rows)
        srows = self.s_rows()
        pR = [self._pR(float(s)) for s in srows]

        for l, layer in enumerate(w.layers):
            angles, features = self.layer_slots(l)
            nslots = len(angles)
            sigma = layer.sigma
            r = layer.radius

            rows = []
            for i, s in enumerate(srows):
                p, R = pR[i]
                phi = sigma * (angles + w.omega * float(s))
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                ids = [builder.add_node(*(p + R @ np.array([x[j], y[j], 0.0])))
                       for j in range(nslots)]
                rows.append(ids)
            builder.rings[l] = rows

            for kind, k, start, count in features:
                tag = builder.new_surface(kind, layer=l, tape=k)
                for i in range(len(srows) - 1):
                    a_row, b_row = rows[i], rows[i + 1]
                    axis = 0.5 * (pR[i][0] + pR[i + 1][0])
                    for j in range(start, start + count):
                        j1 = (j + 1) % nslots
                        a, b = a_row[j], a_row[j1]
                        c, d = b_row[j], b_row[j1]
                        # split along the short diagonal b-c: the quad is
                        # sheared by omega*ds along the helix, so a-d is the
                        # long one and yields ~128 deg triangles that gmsh
                        # cannot fill the inter-layer slab over. Outward
                        # (radial) orientation; sigma flips the handedness
                        if sigma > 0:
                            builder.add_tri(tag, a, b, c)
                            builder.add_tri(tag, b, d, c)
                        else:
                            builder.add_tri(tag, a, c, b)
                            builder.add_tri(tag, b, c, d)
                        self._assert_outward(builder, axis, -2)
                        self._assert_outward(builder, axis, -1)

        # ---------------------------------------------------- domain tube
        rows_s = self.tube_rows()
        nslots = self.tube_slot_count()
        slot = np.arange(nslots) * (2.0 * math.pi / nslots)
        R_dom = self.domain_radius
        L = w.length

        tag = builder.new_surface("tube")
        tube_rows = []
        for s in rows_s:
            p, R = self._pR(float(s))
            phi = slot + self.theta * (float(s) / L)
            x = R_dom * np.cos(phi)
            y = R_dom * np.sin(phi)
            ids = [builder.add_node(*(p + R @ np.array([x[j], y[j], 0.0])))
                   for j in range(nslots)]
            tube_rows.append(ids)
        builder.rings["dom"] = tube_rows

        for i in range(len(rows_s) - 1):
            a_row, b_row = tube_rows[i], tube_rows[i + 1]
            axis = 0.5 * (self._pR(float(rows_s[i]))[0]
                          + self._pR(float(rows_s[i + 1]))[0])
            for j in range(nslots):
                j1 = (j + 1) % nslots
                a, b = a_row[j], a_row[j1]
                c, d = b_row[j], b_row[j1]
                builder.add_tri(tag, a, b, d)
                builder.add_tri(tag, a, d, c)
                self._assert_outward(builder, axis, -2)
                self._assert_outward(builder, axis, -1)

    def _pR(self, s: float):
        R = self.frame.R(s)
        p = self.frame.origin(s)
        return p, R

    def _assert_outward(self, builder, axis, which):
        """The triangle added at builder.tris[which] must face away from
        the local centerline point (radially outward)."""
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

    def check_radius(self, builder, tol=1e-9):
        """Every ring node must sit exactly at its layer radius from the
        centerline (measured in the cross-section plane it was built in)."""
        worst = 0.0
        srows = self.s_rows()
        for l, layer in enumerate(self.w.layers):
            rows = builder.rings[l]
            for i, s in enumerate(srows):
                p, R = self._pR(float(s))
                for nid in rows[i]:
                    d = builder.nodes[nid - 1] - p
                    loc = R.T @ d
                    rr = math.hypot(loc[0], loc[1])
                    worst = max(worst, abs(rr - layer.radius))
        if worst > tol:
            raise AssertionError(
                "shell radius error {:.3e} exceeds {:.1e}".format(worst, tol))
        return worst
