from mesh.Element import Element

## 1-Node Element
class Vertex( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 15, 1, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ]

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        return '{:d} {:d} 2 {:d} {:d} {:d}'.format(
            self.id,           \
            self.GMSH,         \
            self.physicalTag,  \
            self.geometryTag,  \
            self.nodes[ 0 ].id )

# -------------------------------------------------------------------------------