"""Final mesh assembly and validation.

Output follows the established conventions exactly:
- elements 1..10 are the vertex elements "k 15 2 0 k k" on nodes 1..10
  (front center, 4 front quarter points, back center, 4 back quarter pts)
- physical tags: Volume_1..nV, Tape_1..nT, DomainBoundary, DomainFront,
  DomainBack (in this order); gap strip elements are dropped
- one $Periodic dim-2 link DomainBack -> DomainFront with all cap node
  pairs (BELFEM derives the affine map from the vertex triples; the pairs
  are kept for checking)
"""

import math

import numpy as np

import mesh


class ValidationError(Exception):
    pass


class Postprocessor:

    def __init__(self, builder, winding, frame, shellgen):
        self.b = builder
        self.w = winding
        self.frame = frame
        self.sg = shellgen

    # ------------------------------------------------------------ assembly

    def build_final(self, nodes, tets, shell_tris, old_to_gmsh):
        b = self.b

        # mesh/ Node objects per gmsh id
        node_obj = {}
        for gid, p in nodes.items():
            N = mesh.Node(p[0], p[1], p[2])
            N.id = gid
            node_obj[gid] = N

        def make(cls, gmsh_ids, geo, phys):
            E = cls()
            E.nodes = [node_obj[i] for i in gmsh_ids]
            E.geometryTag = geo
            E.physicalTag = phys
            return E

        M = mesh.Mesh()
        elements = []

        # ---- ten vertices first
        vertex_gids = [old_to_gmsh[v] for v in b.vertex_nodes]
        for k, gid in enumerate(vertex_gids):
            V = mesh.Vertex()
            V.nodes = [node_obj[gid]]
            V.geometryTag = k + 1
            V.physicalTag = 0
            elements.append(V)

        # ---- volumes
        nV = max(g for g, _ in tets)
        for region in range(1, nV + 1):
            M.add_physical(3, region, "Volume_{:d}".format(region))
        for g, nn in tets:
            elements.append(make(mesh.Tet4, nn, g, g))

        # ---- tape sidesets (layer-major, tape order)
        tape_tags = []
        for l in range(self.w.numLayers):
            for m_tag, meta in sorted(b.surfaces.items()):
                if meta["role"] == "tape" and meta["layer"] == l:
                    tape_tags.append(m_tag)
        tris_by_surface = {}
        for g, nn in shell_tris:
            tris_by_surface.setdefault(g, []).append(nn)

        tag = nV
        self.tape_physicals = []
        for t_index, m_tag in enumerate(tape_tags):
            tag += 1
            M.add_physical(2, tag, "Tape_{:d}".format(t_index + 1))
            self.tape_physicals.append(tag)
            for nn in tris_by_surface[m_tag]:
                elements.append(make(mesh.Tri3, nn, tag, tag))

        # ---- domain boundary (tube), front, back
        def emit_group(tags_list, name):
            nonlocal tag
            tag += 1
            M.add_physical(2, tag, name)
            for m_tag in tags_list:
                for nn in tris_by_surface[m_tag]:
                    elements.append(make(mesh.Tri3, nn, tag, tag))
            return tag

        # cap pieces by role: cap_disk (bore) and cap_outer (outer air) stay
        # in DomainFront / DomainBack, so the tape ids, the curves
        # ( cap @ tape ) and the vertex-triple periodicity are unchanged;
        # the solder annuli ( cap_annulus, one per interlayer bed ) become
        # their own physical surfaces AFTER them. BELFEM's current terminal
        # for a soldered stack is the solder END FACE as a sideset: the bulk
        # branch of suggest_Homology takes its boundary, the envelope of the
        # conductor cross-section ( tapestack3d.geo, "The current terminal" ).
        # A cap that also carries bore and outer-air faces has the domain
        # boundary circle as its boundary instead, and the imposed current
        # then constrains the periodic axial loop, not the loop around the
        # cable ( corc_solder, 2026-09-02 ).
        def cap_tags(side, roles):
            return sorted(t for t, m in b.surfaces.items()
                          if m["role"] in roles and m["side"] == side)

        self.boundary_physical = emit_group(b.surface_tags("tube"), "DomainBoundary")
        self.front_physical = emit_group(cap_tags("front", ("cap_disk", "cap_outer")), "DomainFront")
        self.back_physical = emit_group(cap_tags("back", ("cap_disk", "cap_outer")), "DomainBack")
        self.solder_front_physicals = []
        self.solder_back_physicals = []
        for aidx in range(self.w.numLayers - 1):
            f = [t for t, m in b.surfaces.items()
                 if m["role"] == "cap_annulus" and m["side"] == "front" and m["annulus"] == aidx]
            bk = [t for t, m in b.surfaces.items()
                  if m["role"] == "cap_annulus" and m["side"] == "back" and m["annulus"] == aidx]
            self.solder_front_physicals.append(emit_group(sorted(f), "SolderFront_{:d}".format(aidx + 1)))
            self.solder_back_physicals.append(emit_group(sorted(bk), "SolderBack_{:d}".format(aidx + 1)))
        self.volume_count = nV

        M.elements = elements

        # ---- node ordering: ten vertex nodes first, then the rest
        used = set()
        for E in elements:
            for N in E.nodes:
                used.add(N.id)
        first = [node_obj[g] for g in vertex_gids]
        first_set = set(vertex_gids)
        rest = [node_obj[g] for g in sorted(used) if g not in first_set]
        M.nodes = first + rest

        # final ids (mesh.save re-derives these, same order)
        for i, N in enumerate(M.nodes):
            N.id = i + 1

        # ---- periodic link with final ids
        pairs = []
        for f_old, b_old in self.b.cap_pairs:
            fg = old_to_gmsh[f_old]
            bg = old_to_gmsh[b_old]
            pairs.append((node_obj[bg].id, node_obj[fg].id))  # slave, master
        M.periodic = [mesh.PeriodicLink(2, self.back_physical,
                                        self.front_physical, None, pairs)]
        return M

    # ---------------------------------------------------------- validation

    def validate(self, M):
        self._check_vertices(M)
        T = self._vertex_transform(M)
        self._check_cap_pairs(M, T)
        self._check_tets(M)
        self._check_end_ring_planes(M)
        print("Validation OK: {:d} nodes, {:d} elements, {:d} cap pairs".format(
            len(M.nodes), len(M.elements), len(M.periodic[0].node_pairs)))

    def _check_vertices(self, M):
        for k in range(10):
            E = M.elements[k]
            if E.GMSH != 15 or E.geometryTag != k + 1 or E.physicalTag != 0 \
                    or E.nodes[0].id != k + 1:
                raise ValidationError(
                    "vertex element {:d} violates 'k 15 2 0 k k'".format(k + 1))

    def _vertex_transform(self, M):
        """Affine map from the front vertex triple (1,2,3) to the back
        triple (6,7,8), the way BELFEM will derive it."""
        pos = {N.id: np.array([N.x, N.y, N.z]) for N in M.nodes[:10]}
        P, A, B = pos[1], pos[2], pos[3]
        Q, E, F = pos[6], pos[7], pos[8]
        u1, u2 = A - P, B - P
        u3 = np.cross(u1, u2)
        v1, v2 = E - Q, F - Q
        v3 = np.cross(v1, v2)
        U = np.column_stack((u1, u2, u3))
        V = np.column_stack((v1, v2, v3))
        R = V @ np.linalg.inv(U)
        # rigid map check
        if np.abs(R.T @ R - np.eye(3)).max() > 1e-6:
            raise ValidationError("vertex-defined map is not a rotation")
        return lambda x: R @ (x - P) + Q

    def _check_cap_pairs(self, M, T, tol=1e-6):
        pos = {N.id: np.array([N.x, N.y, N.z]) for N in M.nodes}
        link = M.periodic[0]
        slaves = set()
        masters = set()
        worst = 0.0
        for s, m in link.node_pairs:
            if s in slaves or m in masters:
                raise ValidationError("cap pairing is not bijective")
            slaves.add(s)
            masters.add(m)
            worst = max(worst, float(np.linalg.norm(T(pos[m]) - pos[s])))
        if worst >= tol:
            raise ValidationError(
                "cap periodicity residual {:.3e} under the vertex-defined "
                "transform".format(worst))
        # every DomainBack / DomainFront node must be covered
        front_nodes = set()
        back_nodes = set()
        for E in M.elements:
            if E.GMSH != 2:
                continue
            if E.geometryTag == self.front_physical:
                front_nodes.update(N.id for N in E.nodes)
            elif E.geometryTag == self.back_physical:
                back_nodes.update(N.id for N in E.nodes)
        if front_nodes - masters or back_nodes - slaves:
            raise ValidationError(
                "cap nodes missing from the periodic link: {:d} front, "
                "{:d} back".format(len(front_nodes - masters),
                                   len(back_nodes - slaves)))
        print("cap periodicity: {:d} pairs, max residual {:.3e}".format(
            len(link.node_pairs), worst))

    def _check_tets(self, M):
        worst_gamma = float("inf")
        inverted = 0
        count = 0
        gammas = []
        for E in M.elements:
            if E.GMSH != 4:
                continue
            count += 1
            p = np.array([[N.x, N.y, N.z] for N in E.nodes])
            a, b, c, d = p
            v = float(np.dot(np.cross(b - a, c - a), d - a)) / 6.0
            if v <= 0:
                inverted += 1
                continue
            faces = [(a, b, c), (a, b, d), (a, c, d), (b, c, d)]
            area = sum(0.5 * np.linalg.norm(np.cross(f[1] - f[0], f[2] - f[0]))
                       for f in faces)
            rin = 3 * v / area
            aa = np.linalg.norm(b - a) * np.linalg.norm(c - d)
            bb = np.linalg.norm(c - a) * np.linalg.norm(b - d)
            cc = np.linalg.norm(d - a) * np.linalg.norm(b - c)
            s = 0.5 * (aa + bb + cc)
            rc = math.sqrt(max(s * (s - aa) * (s - bb) * (s - cc), 0.0)) / (6 * v)
            g = 3 * rin / rc if rc > 0 else 0.0
            gammas.append(g)
            worst_gamma = min(worst_gamma, g)
        if inverted:
            raise ValidationError("{:d} inverted tets".format(inverted))
        gammas = np.array(gammas)
        print("tets: {:d}, gamma min {:.3f} / p5 {:.3f} / median {:.3f}".format(
            count, worst_gamma, float(np.percentile(gammas, 5)),
            float(np.median(gammas))))
        if worst_gamma < 1e-3:
            raise ValidationError("degenerate tets (gamma < 1e-3)")

    def _check_end_ring_planes(self, M, tol=1e-9):
        for s, label in ((0.0, "front"), (self.w.length, "back")):
            p0 = self.frame.origin(s)
            t = self.frame.R(s)[:, 2]
            worst = 0.0
            for l in range(self.w.numLayers):
                row = self.b.rings[l][0 if s == 0.0 else -1]
                for nid in row:
                    worst = max(worst, abs(float(
                        np.dot(self.b.nodes[nid - 1] - p0, t))))
            if worst > tol:
                raise ValidationError(
                    "{:s} end ring off its cap plane by {:.3e}".format(
                        label, worst))
