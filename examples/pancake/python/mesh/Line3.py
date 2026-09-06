from mesh.Element import Element

## 3-Node Line Element
class Line3( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 8, 21, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ] * 3

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        return '{:d} {:d} 2 {:d} {:d} {:d} {:d} {:d}'.format(
            self.id,
            self.GMSH,
            self.physicalTag,
            self.geometryTag,
            self.nodes[ 0 ].id,
            self.nodes[ 1 ].id,
            self.nodes[ 2 ].id )
