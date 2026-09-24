# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Minimal msh2 ASCII reader used by the cap and assembler stages (the
final export still goes through the mesh/ package)."""


def read_msh2(path: str):
    """Returns (nodes, elements): nodes is {id: (x, y, z)}, elements is a
    list of (etype, phys, geo, [node ids])."""
    with open(path) as f:
        lines = f.read().splitlines()

    i = lines.index("$Nodes")
    n = int(lines[i + 1])
    nodes = {}
    for line in lines[i + 2:i + 2 + n]:
        d = line.split()
        nodes[int(d[0])] = (float(d[1]), float(d[2]), float(d[3]))

    i = lines.index("$Elements")
    n = int(lines[i + 1])
    elements = []
    for line in lines[i + 2:i + 2 + n]:
        d = [int(x) for x in line.split()]
        ntags = d[2]
        elements.append((d[1], d[3] if ntags > 0 else 0,
                         d[4] if ntags > 1 else 0, d[3 + ntags:]))
    return nodes, elements
