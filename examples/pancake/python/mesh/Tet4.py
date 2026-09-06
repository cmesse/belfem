from mesh.Element import Element

## 4-Node Tet Element
class Tet4( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 4, 10, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ] * 4

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        return '{:d} {:d} 2 {:d} {:d} {:d} {:d} {:d} {:d}'.format(
            self.id,
            self.GMSH,
            self.physicalTag,
            self.geometryTag,
            self.nodes[0].id,
            self.nodes[1].id,
            self.nodes[2].id,
            self.nodes[3].id)

# -------------------------------------------------------------------------------