from mesh.Element import Element

## 18-Node Pentaedron Element
class Penta15( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 18, 26, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ] * 15

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        return '{:d} {:d} 2 {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d}'.format(
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
            self.nodes[ 8 ].id,
            self.nodes[ 9].id,
            self.nodes[ 10 ].id,
            self.nodes[ 11 ].id,
            self.nodes[ 12 ].id,
            self.nodes[ 13 ].id,
            self.nodes[ 14 ].id,
            self.nodes[ 15 ].id )