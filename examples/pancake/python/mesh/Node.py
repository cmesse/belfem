# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

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
