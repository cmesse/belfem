# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

class PeriodicLink:

    def __init__(self, dim: int, slave_tag: int, master_tag: int,
                 affine=None, node_pairs=None):
        self.dim = dim
        self.slave_tag = slave_tag
        self.master_tag = master_tag
        self.affine = affine
        self.node_pairs = node_pairs if node_pairs is not None else []

    def with_node_pairs(self, node_pairs):
        return PeriodicLink(self.dim, self.slave_tag, self.master_tag,
                            self.affine, node_pairs)
