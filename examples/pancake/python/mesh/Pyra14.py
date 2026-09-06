from mesh.Element import Element

## 14-Node Pyramid Element
class Pyra14( Element ):

    def __init__( self, Geo=0, Phys=0 ):

        # call parent constructor
        Element.__init__( self, 14, 14, Geo, Phys )

        # Allocate Node Array
        self.nodes = [ None ] * 14

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
        s = '{:d} {:d} 2 {:d} {:d}'.format( self.id,
            self.GMSH,
            self.physicalTag,
            self.geometryTag )

        for k in range(14):
            s = s + ' ' + str( self.nodes[k].id )

        return s
