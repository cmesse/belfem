# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

from mesh.Element import Element

## 2-Node Line Element
class Line2( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 1, 3, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ] * 2

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        return '{:d} {:d} 2 {:d} {:d} {:d} {:d}'.format(
            self.id,
            self.GMSH,
            self.physicalTag,
            self.geometryTag,
            self.nodes[ 0 ].id,
            self.nodes[ 1 ].id )
