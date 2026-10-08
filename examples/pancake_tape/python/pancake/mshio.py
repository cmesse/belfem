# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""msh2 ASCII reader and writer on numpy arrays. Meshes of a million
elements pass through here, so nothing is stored per element as a
Python object."""

import numpy as np


def read_msh2(path: str):
    """Returns (node_ids, xyz, elements): node_ids (n,), xyz (n, 3),
    elements a dict etype -> (tags (m, 2) [physical, geometry],
    connectivity (m, nn))."""
    with open(path) as f:
        text = f.read()
    lines = text.split("\n")

    i = lines.index("$Nodes")
    n = int(lines[i + 1])
    block = np.fromstring(" ".join(lines[i + 2:i + 2 + n]), sep=" ")
    block = block.reshape(n, 4)
    node_ids = block[:, 0].astype(int)
    xyz = block[:, 1:4].copy()

    i = lines.index("$Elements")
    m = int(lines[i + 1])
    per_type = {}
    for line in lines[i + 2:i + 2 + m]:
        d = line.split()
        etype = int(d[1])
        ntags = int(d[2])
        per_type.setdefault(etype, []).append(
            (int(d[3]) if ntags > 0 else 0, int(d[4]) if ntags > 1 else 0,
             *[int(x) for x in d[3 + ntags:]]))
    elements = {}
    for etype, rows in per_type.items():
        arr = np.array(rows, dtype=int)
        elements[etype] = (arr[:, :2], arr[:, 2:])
    return node_ids, xyz, elements


def write_msh2(path: str, xyz, groups, physicals=None):
    """xyz (n, 3) with implicit ids 1..n; groups is a list of
    (etype, physical, geometry, connectivity (m, nn)) written in order;
    physicals is a list of (dimension, tag, name)."""
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        if physicals:
            f.write("$PhysicalNames\n{:d}\n".format(len(physicals)))
            for dim, tag, name in physicals:
                f.write('{:d} {:d} "{:s}"\n'.format(dim, tag, name))
            f.write("$EndPhysicalNames\n")
        n = xyz.shape[0]
        f.write("$Nodes\n{:d}\n".format(n))
        ids = np.arange(1, n + 1)
        np.savetxt(f, np.column_stack((ids, xyz)),
                   fmt=["%d", "%.17g", "%.17g", "%.17g"])
        f.write("$EndNodes\n")
        total = sum(g[3].shape[0] for g in groups)
        f.write("$Elements\n{:d}\n".format(total))
        eid = 1
        for etype, phys, geo, conn in groups:
            m = conn.shape[0]
            if m == 0:
                continue
            head = np.column_stack((
                np.arange(eid, eid + m), np.full(m, etype), np.full(m, 2),
                np.full(m, phys), np.full(m, geo)))
            np.savetxt(f, np.hstack((head, conn)), fmt="%d")
            eid += m
        f.write("$EndElements\n")
