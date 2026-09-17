# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

class Physical :

    def __init__(self, dimension: int, phystag: int, label: str ):
        self.Dimension = dimension
        self.physicalTag = phystag
        self.Label = label

    def gmsh(self):
        return f'{self.Dimension:d} {self.physicalTag:d} "{self.Label:s}"'