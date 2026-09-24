# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

from mesh.Element import Element

## 9-Node Quadrilateral Element
class Quad9( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 10, 28, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ] * 9

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        return '{:d} {:d} 2 {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d}'.format(
            self.id,
            self.GMSH,
            self.physicalTag,
            self.geometryTag,
            self.nodes[ 0 ].id,
            self.nodes[ 1 ].id,
            self.nodes[ 2 ].id,
            self.nodes[ 3 ].id,
            self.nodes[ 4 ].id,
            self.nodes[ 5 ].id,
            self.nodes[ 6 ].id,
            self.nodes[ 7 ].id,
            self.nodes[ 8 ].id )