# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Structured mesh of the solder pack, the two lead bars and the tape
sheet.

Node (k, m) sits on column k at pack level m; k < 0 are the sliver's
inner boundary columns, k >= 0 the sheet columns. A pack cell k is the
hexahedron between columns (k, k+1) inside and (k + nTheta, k + nTheta + 1)
outside, one level tall; the sliver cells are k in [-nTheta, 0). The
inner lead is a strip in the unrolled plane of the small cylinder, node
(q, j) at forward index q and width index j; its last nWi + 1 columns
ARE the sliver's inner boundary columns -nWi..0, so the sliver cells
pair with the departing tape. A two-cell-wide bar (inner boundary,
middle line one pitch inside the tape, tape) rides inside the strip;
at the departure point its middle column is the sliver's column
-nTheta, so the bar's start face shares the sliver's start face. The
outer lead bar is one cell wide between its inner boundary and the
sheet over the pack levels.

Every hexahedron is split into the six Kuhn tetrahedra sharing its
main diagonal (0,0,0) -> (1,1,1) in its local (along, outward, up) frame;
each quad face is then cut from its (0,0) to its (1,1) corner. Every
cell uses the same frame directions (along = forward along the tape,
outward = away from the winding's or the small cylinder's axis, up =
width index), so neighbors agree on the diagonal and the surface
triangles below are faces of the tetrahedra.

BELFEM opens a thin shell by block: the slave block of each sideset
gets the duplicate nodes. That is why the pack is cut into half-turn
blocks with the tape sidesets split on the same columns, why a solder
sliver sits inside the first turn, and why both leads carry a bar: a
tape facet must always separate two different blocks, a block must
never be inside one piece of a sideset and outside another, and where
air touches the tape it must be on the outer side, so that "master =
inner" holds along the whole sheet.
"""

import numpy as np


KUHN = (((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)),
        ((0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)),
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)),
        ((0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)),
        ((0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)),
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)))

# parity of the six paths above: (i,j,l) +, (i,l,j) -, (j,i,l) -,
# (j,l,i) +, (l,i,j) +, (l,j,i) -
KUHN_PARITY = np.array([1, -1, -1, 1, 1, -1])


def kuhn(v):
    """Six tets per hexahedron; v maps (i, j, l) to an id array."""
    return np.stack([np.column_stack([v[c] for c in p]) for p in KUHN],
                    axis=1).reshape(-1, 4)


def face_tris(v00, v10, v01, v11):
    """The two triangles of a quad face, cut along (0,0)-(1,1) in the
    face's two local coordinates."""
    return np.vstack((np.column_stack((v00, v10, v11)),
                      np.column_stack((v00, v11, v01))))


class Pack:

    # surfaces handed to gmsh as the air boundary, plus the terminals
    (TOP, BOTTOM, SLIVER_INNER, TAPE_OUTSIDE,
     LEAD_IN_OUTSIDE, LEAD_IN_INNER, LEAD_IN_SIDES, LEAD_IN_START, TERMINAL_IN,
     LEAD_BAR_INNER, LEAD_BAR_TOP, LEAD_BAR_BOTTOM, TERMINAL_OUT,
     TUBE, CAP_TOP, CAP_BOTTOM) = range(1, 17)
    BAR_SURFACES = (LEAD_IN_OUTSIDE, LEAD_IN_INNER, LEAD_IN_SIDES, LEAD_IN_START,
                    LEAD_BAR_INNER, LEAD_BAR_TOP, LEAD_BAR_BOTTOM)
    # the pack faces where solder or a bar touches the air; the tape
    # facets facing the air (TAPE_OUTSIDE, LEAD_IN_OUTSIDE) belong to the
    # thin shell and the terminals are domain faces
    SOLDER_AIR_SURFACES = (TOP, BOTTOM, SLIVER_INNER,
                           LEAD_IN_INNER, LEAD_IN_SIDES, LEAD_IN_START,
                           LEAD_BAR_INNER, LEAD_BAR_TOP, LEAD_BAR_BOTTOM)

    def __init__(self, sheet):
        self.S = S = sheet
        nT = S.nTheta
        self.nz1 = nz1 = S.nZ + 1
        self.nCol = S.M + 1 + nT

        blocks = []
        allcols = np.vstack((S.vcols, S.cols))
        blocks.append(np.column_stack((np.repeat(allcols, nz1, axis=0),
                                       np.tile(S.zPack, self.nCol))))
        self.nPack = self.nCol * nz1
        base = self.nPack

        # inner lead: tape nodes for q < nU - nWi, middle nodes for
        # q < nU, inner nodes for all q; each (q, j) row-major in (q, j)
        nQt = S.nU - S.nWi
        self.offTape = base
        base += nQt * nz1
        tape = S.lead_point(S.leadSZ[:nQt], S.r1)
        blocks.append(tape.reshape(-1, 3))
        self.offMid = base
        base += S.nU * nz1
        mid = S.lead_point(S.leadSZ[:S.nU], S.r1 - S.pitch)
        blocks.append(mid.reshape(-1, 3))
        self.offInner = base
        base += (S.nU + 1) * nz1
        inner = S.lead_point(S.leadSZ, S.r1 - S.tIn)
        blocks.append(inner.reshape(-1, 3))

        # outer lead bar
        self.offLead = base
        nlb = S.M - S.kLead
        base += nlb * nz1
        blocks.append(np.column_stack((np.repeat(S.leadBar[1:], nz1, axis=0),
                                       np.tile(S.zPack, nlb))))
        self.xyz = np.vstack(blocks)

        self.blockNames = (["Solder_{:d}".format(p + 1) for p in range(S.nPieces)]
                           + ["LeadIn", "LeadOut"])
        self.LEAD_IN = S.nPieces + 1
        self.LEAD_OUT = S.nPieces + 2
        self.innerPiece = S.nSpiralTapePieces       # 0-based tape piece of the inner lead
        self.nTapePieces = S.nSpiralTapePieces + 1

        self._build_tets()
        self._build_sheet()
        self._build_surfaces()

    # ------------------------------------------------------------ node ids

    def pid(self, k, m):
        return (np.asarray(k) + self.S.nTheta) * self.nz1 + np.asarray(m) + 1

    def tid(self, q, j):
        """Inner lead tape node: the sliver boundary column q - nU for
        the wedge part, own nodes before it."""
        q, j = np.asarray(q), np.asarray(j)
        S = self.S
        return np.where(q >= S.nU - S.nWi, self.pid(q - S.nU, j),
                        self.offTape + q * self.nz1 + j + 1)

    def mid(self, q, j):
        q, j = np.asarray(q), np.asarray(j)
        S = self.S
        return np.where(q == S.nU, self.pid(-S.nTheta, j),
                        self.offMid + q * self.nz1 + j + 1)

    def iid(self, q, j):
        return self.offInner + np.asarray(q) * self.nz1 + np.asarray(j) + 1

    def lid(self, k, m):
        k, m = np.asarray(k), np.asarray(m)
        S = self.S
        return np.where(k == S.kLead, self.pid(S.kLead - S.nTheta, m),
                        self.offLead + (k - S.kLead - 1) * self.nz1 + m + 1)

    def hole_outline_ids(self):
        """Closed outline of the inner bar's face on the cap (q = 0)."""
        j = np.arange(self.nz1)
        nZ = self.S.nZ
        return np.concatenate((
            self.iid(0, j),
            [self.mid(0, nZ)],
            self.tid(np.zeros(self.nz1, dtype=int), j[::-1]),
            [self.mid(0, 0)]))

    def lead_end_ids(self):
        return self.pid(self.S.M, np.arange(self.nz1))

    def lead_bar_end_ids(self):
        return self.lid(np.full(self.nz1, self.S.M), np.arange(self.nz1))

    # ---------------------------------------------------------------- tets

    def _build_tets(self):
        S = self.S
        nT = S.nTheta
        parts = []
        blocks = []

        kk, mm = np.meshgrid(np.arange(-nT, S.nCells), np.arange(S.nZ), indexing="ij")
        kk, mm = kk.ravel(), mm.ravel()
        v = {(i, j, l): self.pid(kk + i + j * nT, mm + l)
             for i in (0, 1) for j in (0, 1) for l in (0, 1)}
        parts.append(kuhn(v))
        blocks.append(np.repeat(S.piece(kk) + 1, 6))

        qq, jj = np.meshgrid(np.arange(S.nU), np.arange(S.nZ), indexing="ij")
        qq, jj = qq.ravel(), jj.ravel()
        v = {(i, j, l): (self.iid(qq + i, jj + l) if j == 0 else self.mid(qq + i, jj + l))
             for i in (0, 1) for j in (0, 1) for l in (0, 1)}
        parts.append(kuhn(v))
        v = {(i, j, l): (self.mid(qq + i, jj + l) if j == 0 else self.tid(qq + i, jj + l))
             for i in (0, 1) for j in (0, 1) for l in (0, 1)}
        parts.append(kuhn(v))
        blocks.append(np.full(2 * 6 * len(qq), self.LEAD_IN))

        kk, mm = np.meshgrid(np.arange(S.kLead, S.M), np.arange(S.nZ), indexing="ij")
        kk, mm = kk.ravel(), mm.ravel()
        v = {(i, j, l): (self.lid(kk + i, mm + l) if j == 0 else self.pid(kk + i, mm + l))
             for i in (0, 1) for j in (0, 1) for l in (0, 1)}
        parts.append(kuhn(v))
        blocks.append(np.full(6 * len(kk), self.LEAD_OUT))

        tets = np.vstack(parts)
        self.tetBlock = np.concatenate(blocks)
        P = self.xyz[tets - 1]
        vol = np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                        P[:, 3] - P[:, 0])
        if np.any(vol == 0.0):
            raise AssertionError("degenerate structured tetrahedra")
        # the six Kuhn tets of a sound cell carry the parity of their
        # path times the handedness of the cell's frame; a cell breaking
        # that pattern is twisted or inverted and must not be repaired by
        # flipping
        sign = np.sign(vol).reshape(-1, 6) * KUHN_PARITY[None, :]
        mixed = np.any(sign != sign[:, :1], axis=1)
        if np.any(mixed):
            raise AssertionError(
                "{:d} structured cells are twisted or inverted (first at tet {:d}); "
                "check the bar thicknesses against the geometry".format(
                    int(mixed.sum()), int(np.argmax(mixed) * 6)))
        flip = vol < 0.0
        tets[flip, 1], tets[flip, 2] = tets[flip, 2], tets[flip, 1].copy()
        self.tets = tets

    # --------------------------------------------------------------- sheet

    def _oriented(self, tris, ref):
        """Flip triangles whose normal does not point along ref."""
        P = self.xyz[tris - 1]
        n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
        flip = np.einsum("ij,ij->i", n, ref) < 0.0
        tris[flip, 1], tris[flip, 2] = tris[flip, 2], tris[flip, 1].copy()
        return tris

    def _travel_ref(self, k):
        """Right of the direction of travel: radially outward on the
        spiral, away from the coil on the outer lead."""
        d = self.S.col(k + 1) - self.S.col(k)
        return np.column_stack((d[:, 1], -d[:, 0], np.zeros(len(k))))

    def _segments(self, k0, k1):
        S = self.S
        kk, mm = np.meshgrid(np.arange(k0, k1), np.arange(S.nZ), indexing="ij")
        kk, mm = kk.ravel(), mm.ravel()
        tris = face_tris(self.pid(kk, mm), self.pid(kk + 1, mm),
                         self.pid(kk, mm + 1), self.pid(kk + 1, mm + 1))
        seg = np.concatenate((kk, kk))
        return self._oriented(tris, self._travel_ref(seg)), seg

    def _build_sheet(self):
        S = self.S
        tris, seg = self._segments(0, S.M)
        piece = S.piece(seg) - 2
        # inner lead: normal away from the small cylinder's axis
        qq, jj = np.meshgrid(np.arange(S.nU), np.arange(S.nZ), indexing="ij")
        qq, jj = qq.ravel(), jj.ravel()
        lead = face_tris(self.tid(qq, jj), self.tid(qq + 1, jj),
                         self.tid(qq, jj + 1), self.tid(qq + 1, jj + 1))
        c = self.xyz[lead - 1].mean(axis=1)
        ref = np.column_stack((c[:, 0] - S.C1[0], c[:, 1] - S.C1[1], np.zeros(len(c))))
        lead = self._oriented(lead, ref)
        self.sheetPack, self.sheetPackSegment = tris, seg
        self.sheetLead = lead
        self.sheetLeadQ = np.concatenate((qq, qq))
        self.sheetTris = np.vstack((tris, lead))
        self.sheetPiece = np.concatenate((piece, np.full(lead.shape[0], self.innerPiece)))
        self.sliverInner, _ = self._segments(-S.nTheta, -S.nWi)

    # ------------------------------------------------------------ surfaces

    def _build_surfaces(self):
        S = self.S
        nT, nZ, nU = S.nTheta, S.nZ, S.nU
        kk = np.arange(-nT, S.nCells)

        def cap(m):
            return face_tris(self.pid(kk, m), self.pid(kk + 1, m),
                             self.pid(kk + nT, m), self.pid(kk + nT + 1, m))

        top = cap(nZ)
        bottom = cap(0)
        outside = self.sheetPack[self.sheetPackSegment >= S.nCells]

        # inner lead
        lead_outside = self.sheetLead[self.sheetLeadQ < nU - S.nWi]
        qq, jj = np.meshgrid(np.arange(nU), np.arange(nZ), indexing="ij")
        qq, jj = qq.ravel(), jj.ravel()
        lead_inner = face_tris(self.iid(qq, jj), self.iid(qq + 1, jj),
                               self.iid(qq, jj + 1), self.iid(qq + 1, jj + 1))
        q = np.arange(nU)
        sides = []
        for j in (0, nZ):
            sides.append(face_tris(self.iid(q, j), self.iid(q + 1, j),
                                   self.mid(q, j), self.mid(q + 1, j)))
            sides.append(face_tris(self.mid(q, j), self.mid(q + 1, j),
                                   self.tid(q, j), self.tid(q + 1, j)))
        lead_sides = np.vstack(sides)
        j = np.arange(nZ)
        lead_start = face_tris(self.iid(nU, j), self.mid(nU, j),
                               self.iid(nU, j + 1), self.mid(nU, j + 1))
        terminal_in = np.vstack((
            face_tris(self.iid(0, j), self.mid(0, j), self.iid(0, j + 1), self.mid(0, j + 1)),
            face_tris(self.mid(0, j), self.tid(0, j), self.mid(0, j + 1), self.tid(0, j + 1))))

        # outer lead bar
        kk3, mm3 = np.meshgrid(np.arange(S.kLead, S.M), np.arange(nZ), indexing="ij")
        kk3, mm3 = kk3.ravel(), mm3.ravel()
        bar_inner = face_tris(self.lid(kk3, mm3), self.lid(kk3 + 1, mm3),
                              self.lid(kk3, mm3 + 1), self.lid(kk3 + 1, mm3 + 1))
        k = np.arange(S.kLead, S.M)
        bar_top = face_tris(self.lid(k, nZ), self.lid(k + 1, nZ),
                            self.pid(k, nZ), self.pid(k + 1, nZ))
        bar_bottom = face_tris(self.lid(k, 0), self.lid(k + 1, 0),
                               self.pid(k, 0), self.pid(k + 1, 0))
        mm = np.arange(nZ)
        terminal_out = face_tris(self.lid(S.M, mm), self.pid(S.M, mm),
                                 self.lid(S.M, mm + 1), self.pid(S.M, mm + 1))

        self.surfaces = [
            (self.TOP, top), (self.BOTTOM, bottom),
            (self.SLIVER_INNER, self.sliverInner), (self.TAPE_OUTSIDE, outside),
            (self.LEAD_IN_OUTSIDE, lead_outside), (self.LEAD_IN_INNER, lead_inner),
            (self.LEAD_IN_SIDES, lead_sides), (self.LEAD_IN_START, lead_start),
            (self.TERMINAL_IN, terminal_in),
            (self.LEAD_BAR_INNER, bar_inner), (self.LEAD_BAR_TOP, bar_top),
            (self.LEAD_BAR_BOTTOM, bar_bottom), (self.TERMINAL_OUT, terminal_out)]

    def shell_surfaces(self):
        return list(self.surfaces)

    def solder_air_tris(self):
        """The SOLDER_AIR_SURFACES as one triangle array, each normal
        pointing out of the pack tet that owns the face."""
        S = dict(self.surfaces)
        tris = np.vstack([S[t] for t in self.SOLDER_AIR_SURFACES]).copy()
        base = self.xyz.shape[0] + 1

        def keys(t):
            t = np.sort(t, axis=1).astype(np.int64)
            return (t[:, 0] * base + t[:, 1]) * base + t[:, 2]

        tets = self.tets
        faces = np.vstack((tets[:, [0, 1, 2]], tets[:, [0, 1, 3]],
                           tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]))
        owner = np.tile(np.arange(tets.shape[0]), 4)
        fk = keys(faces)
        order = np.argsort(fk, kind="stable")
        fk, owner = fk[order], owner[order]
        k = keys(tris)
        pos = np.minimum(np.searchsorted(fk, k), len(fk) - 1)
        if np.any(fk[pos] != k):
            raise AssertionError("a solder/air triangle is not a face of a pack tet")
        nxt = np.minimum(pos + 1, len(fk) - 1)
        if np.any((fk[nxt] == k) & (nxt != pos)):
            raise AssertionError("a solder/air triangle is shared by two pack tets")
        out = self.xyz[tris - 1].mean(axis=1) - self.xyz[tets[owner[pos]] - 1].mean(axis=1)
        return self._oriented(tris, out)
