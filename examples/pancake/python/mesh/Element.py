#===============================================================================
## An abstract element class
class Element:
#-------------------------------------------------------------------------------

    ## Element constructor#
    # @param: GMSH: Element Type in GMSH
    # @param: VTK:  Element Type in VTK
    # @param: GEO: Geometry Tag in GMSH
    # @param: PHYS:  Physisical Tag in GMSH

    def __init__( self, GMSH=0, VTK=0, Geo=0, Phys=0 ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#       Member Variables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ## ID of this element
        self.id = -1

        ## Array of nodes connected to this element
        ## Ordered according to the GMSH order
        self.nodes = []

        ## Geometry Tag of this element ( Block Number )
        self.geometryTag = Geo

        ## Physical Tag of this element ( Material ID )
        self.physicalTag = Phys

        ## GMSH Tag of Element
        self.GMSH = GMSH

        ## VTK Tag of Element
        self.VTK = VTK

# -------------------------------------------------------------------------------

    ## create a gmsh string of this element
    def gmsh(self):
            pass


# -------------------------------------------------------------------------------

    ## flag all nodes that belong to thos element
    def flag_nodes(self):
        for n in self.nodes :
            n.flag = True