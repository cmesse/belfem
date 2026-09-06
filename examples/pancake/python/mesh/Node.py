#===============================================================================
## The node class for this element
class Node:
#-------------------------------------------------------------------------------

    def __init__( self, x=0.0, y=0.0, z=0.0 ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#       Member Variables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # ID
        self.id = -1

        # index in memory ( needed for orientation vector computation )
        self.index = -1

        # Coordinates in m
        self.x = x
        self.y = y
        self.z = z

        # Flag for a node
        self.flag = False

# -------------------------------------------------------------------------------

    def gmsh(self):
            return "{:g} {:.15f} {:.15f} {:.15f}".format(self.id, self.x, self.y, self.z)
