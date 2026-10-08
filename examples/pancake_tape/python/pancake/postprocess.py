# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Validation of the assembled mesh: tetrahedron orientation and
quality, node valence, closedness of the boundary, the terminals, the
solder/air contact sideset, and the three rules BELFEM's block-wise thin-shell opening imposes on the
tape sidesets."""

import numpy as np


class ValidationError(Exception):
    pass


def tet_quality(xyz, tets):
    """Signed volumes and gamma = 3 r_in / r_circ."""
    P = xyz[tets - 1]
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    vol = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0
    area = 0.0
    for f in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
        area = area + 0.5 * np.linalg.norm(np.cross(f[1] - f[0], f[2] - f[0]), axis=1)
    rin = 3.0 * np.abs(vol) / area
    aa = np.linalg.norm(b - a, axis=1) * np.linalg.norm(c - d, axis=1)
    bb = np.linalg.norm(c - a, axis=1) * np.linalg.norm(b - d, axis=1)
    cc = np.linalg.norm(d - a, axis=1) * np.linalg.norm(b - c, axis=1)
    s = 0.5 * (aa + bb + cc)
    rc = np.sqrt(np.maximum(s * (s - aa) * (s - bb) * (s - cc), 0.0)) / (6.0 * np.abs(vol))
    gamma = np.where(rc > 0.0, 3.0 * rin / np.maximum(rc, 1e-300), 0.0)
    return vol, gamma


def _keys(tris, base):
    t = np.sort(tris, axis=1).astype(np.int64)
    return (t[:, 0] * base + t[:, 1]) * base + t[:, 2]


def _edge_keys(tris, base):
    t = tris.astype(np.int64)
    e = np.vstack((t[:, [0, 1]], t[:, [1, 2]], t[:, [0, 2]]))
    e = np.sort(e, axis=1)
    return e[:, 0] * base + e[:, 1]


def validate(sheet, xyz, tets, tet_block, block_names, air_block,
             sheet_tris, sheet_piece, sidesets, inner_piece, outer_piece):
    """sidesets: dict name -> tris for TerminalIn, TerminalOut,
    DomainBoundary, DomainTop, DomainBottom and SolderAir."""
    S = sheet
    base = xyz.shape[0] + 1

    vol, gamma = tet_quality(xyz, tets)
    if np.any(vol <= 0.0):
        raise ValidationError("{:d} inverted tets".format(int((vol <= 0).sum())))
    if gamma.min() < 1e-4:
        raise ValidationError("degenerate tets (gamma < 1e-4)")
    flat = int((gamma < 1e-3).sum())
    if flat:
        print("WARNING: {:d} air tets with gamma below 1e-3 (min {:.1e})".format(
            flat, gamma.min()))
    groups = [("Solder", (tet_block < air_block - 2)),
              ("LeadIn", tet_block == air_block - 2),
              ("LeadOut", tet_block == air_block - 1),
              ("Air", tet_block == air_block)]
    for name, mask in groups:
        g = gamma[mask]
        print("{:>8s}: {:8d} tets, gamma min {:.3f} / p5 {:.3f} / median {:.3f}".format(
            name, int(mask.sum()), g.min(), float(np.percentile(g, 5)), float(np.median(g))))

    # BELFEM counts a node's elements in a byte
    valence = np.bincount(tets.ravel())
    if valence.max() >= 256:
        raise ValidationError(
            "{:d} nodes carry 256 or more tets (max {:d}); lower "
            "packAirResolution".format(int((valence >= 256).sum()), int(valence.max())))
    print("node valence: max {:d}, {:d} nodes above 128".format(
        int(valence.max()), int((valence >= 128).sum())))

    # faces of the tets, with the owning tet of each face occurrence
    faces = np.vstack((tets[:, [0, 1, 2]], tets[:, [0, 1, 3]],
                       tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]))
    owner = np.tile(np.arange(tets.shape[0]), 4)
    fk = _keys(faces, base)
    order = np.argsort(fk, kind="stable")
    fk_sorted = fk[order]
    owner_sorted = owner[order]
    uniq, first, count = np.unique(fk_sorted, return_index=True, return_counts=True)
    if np.any(count > 2):
        raise ValidationError("a face is shared by more than two tets")

    def lookup(tris):
        k = _keys(tris, base)
        pos = np.searchsorted(uniq, k)
        pos = np.minimum(pos, len(uniq) - 1)
        found = uniq[pos] == k
        return pos, found

    # every boundary face belongs to a domain sideset, no more, no less
    boundary = uniq[count == 1]
    dom = np.concatenate([_keys(sidesets[n], base) for n in
                          ("DomainBoundary", "DomainTop", "DomainBottom",
                           "TerminalIn", "TerminalOut")])
    if len(np.setdiff1d(boundary, dom)) or len(np.setdiff1d(dom, boundary)):
        raise ValidationError(
            "domain sidesets and the mesh boundary differ ({:d} boundary faces "
            "untagged, {:d} tagged faces not on the boundary)".format(
                len(np.setdiff1d(boundary, dom)), len(np.setdiff1d(dom, boundary))))

    # tape facets: two tets, the outer one on the normal side
    pos, found = lookup(sheet_tris)
    if not found.all() or np.any(count[pos] != 2):
        raise ValidationError("a tape facet is not shared by exactly two tets")
    t1 = owner_sorted[first[pos]]
    t2 = owner_sorted[first[pos] + 1]
    P = xyz[sheet_tris - 1]
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    c = P.mean(axis=1)
    c1 = xyz[tets[t1] - 1].mean(axis=1)
    side1 = np.einsum("ij,ij->i", c1 - c, n) > 0.0
    outer = np.where(side1, t1, t2)
    inner = np.where(side1, t2, t1)
    bo = tet_block[outer]
    bi = tet_block[inner]
    if np.any(bo == bi):
        raise ValidationError(
            "{:d} tape facets have the same block on both sides".format(int((bo == bi).sum())))
    air_in = bi == air_block
    if np.any(air_in):
        raise ValidationError(
            "{:d} tape facets have air on their inner side; the master side "
            "would flip along the sheet".format(int(air_in.sum())))
    for p in range(sheet_piece.max() + 1):
        m = sheet_piece == p
        clash = np.intersect1d(np.unique(bo[m]), np.unique(bi[m]))
        if len(clash):
            raise ValidationError(
                "tape piece {:d}: block(s) {} lie inside and outside the same "
                "sideset; the block-wise node relinking would tear the sheet".format(
                    p + 1, [block_names[b - 1] for b in clash]))
    print("tape: {:d} facets in {:d} pieces; every facet separates two blocks, air is "
          "always outside, no block is on both sides of a piece".format(
              sheet_tris.shape[0], int(sheet_piece.max()) + 1))

    # SolderAir: every triangle has one solder or bar tet and one air tet,
    # and every such face that is not a tape facet is in it, so the cut
    # pruning finds the conductor/air contact closed
    contact = sidesets["SolderAir"]
    pos, found = lookup(contact)
    if not found.all() or np.any(count[pos] != 2):
        raise ValidationError("a SolderAir triangle is not shared by exactly two tets")
    b1 = tet_block[owner_sorted[first[pos]]]
    b2 = tet_block[owner_sorted[first[pos] + 1]]
    bad = (b1 == air_block) == (b2 == air_block)
    if np.any(bad):
        raise ValidationError("{:d} SolderAir triangles do not have exactly one air "
                              "and one solder tet".format(int(bad.sum())))
    interior = uniq[count == 2]
    i1 = tet_block[owner_sorted[first[count == 2]]]
    i2 = tet_block[owner_sorted[first[count == 2] + 1]]
    mixed = interior[(i1 == air_block) != (i2 == air_block)]
    missing = np.setdiff1d(np.setdiff1d(mixed, _keys(sheet_tris, base)),
                           _keys(contact, base))
    if len(missing):
        raise ValidationError("{:d} solder/air faces are neither tape facets nor in "
                              "SolderAir".format(len(missing)))
    print("SolderAir: {:d} facets, each between one solder and one air tet; together "
          "with the {:d} tape facets facing air they cover the conductor/air "
          "contact".format(contact.shape[0], len(np.intersect1d(mixed, _keys(sheet_tris, base)))))

    # terminals: their edges shared with the inner lead piece and the
    # last spiral piece
    for name, piece, expect, where in (
            ("TerminalIn", inner_piece, S.nZ, lambda p: np.abs(p[:, 2] - S.sigma * S.H) < 1e-9),
            ("TerminalOut", outer_piece, S.nZ,
             lambda p: np.abs(np.hypot(p[:, 0], p[:, 1]) - S.Rdom) < 1e-9)):
        tri = sidesets[name]
        nodes = np.unique(tri)
        if not where(xyz[nodes - 1]).all():
            raise ValidationError("{} is not on the domain boundary".format(name))
        shared = np.intersect1d(_edge_keys(tri, base),
                                _edge_keys(sheet_tris[sheet_piece == piece], base))
        if len(shared) != expect:
            raise ValidationError("{} shares {:d} edges with tape piece {:d}, expected "
                                  "{:d}".format(name, len(shared), piece + 1, expect))
        print("{}: {:d} facets, {:d} edges on tape piece {:d}".format(
            name, tri.shape[0], len(shared), piece + 1))
