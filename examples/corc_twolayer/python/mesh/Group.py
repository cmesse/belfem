
# SideSet for the Mesh definition
class Group:
    def __init__(self, label: str, dimension: int ):
        self.label = label
        self.dimension = dimension

        self.MaterialID = 0

        # container for elements
        self.elements = []

        # ids for lines
        self._lines = 0

        # ids for triangles
        self._tris = 0

        # ids for quad
        self._quads = 0

        # ids for tets
        self._tets = 0

        # ids for pentas
        self._pentas = 0

        # ids for hexs
        self._hexs = 0

        # Element Counters
        self._num_lines  = 0
        self._num_tris   = 0
        self._num_quads  = 0
        self._num_tets   = 0
        self._num_pentas = 0
        self._num_hexs   = 0

    # -------------------------------------------------------------------------------

    def _count_elements(self):

        # reset all placeholders
        self._num_lines  = 0
        self._num_tris   = 0
        self._num_quads  = 0
        self._num_tets   = 0
        self._num_pentas = 0
        self._num_hexs   = 0

        # loop over all elements in group
        for E in self.elements:
            # Check GMSH number of element
            if E.GMSH == 1 or E.GMSH == 8 :
                self._num_lines += 1
            elif E.GMSH == 2 or E.GMSH == 9 :
                self._num_tris += 1
            elif E.GMSH == 3 or E.GMSH == 10 or E.GMSH == 16 :
                self._num_quads += 1
            elif E.GMSH == 4 or E.GMSH == 11 :
               self._num_tets +=1
            elif E.GMSH == 6 or E.GMSH == 12 or E.GMSH == 18 :
                self._num_pentas +=1
            elif E.GMSH == 5 or E.GMSH == 12 or E.GMSH == 17 :
                self._num_hexs +=1

    # -------------------------------------------------------------------------------

    def _compute_block_numbers(self, BlockCount: int ):

        if self._num_lines > 0 :
            BlockCount += 1
            self._lines = BlockCount
        if self._num_tris > 0 :
            BlockCount += 1
            self._tris = BlockCount
        if self._num_quads > 0 :
            BlockCount += 1
            self._quads = BlockCount
        if self._num_tets > 0 :
            BlockCount += 1
            self._tets = BlockCount
        if self._num_pentas > 0 :
            BlockCount += 1
            self._pentas = BlockCount
        if self._num_hexs > 0 :
            BlockCount += 1
            self._hexs = BlockCount

        return BlockCount

    # -------------------------------------------------------------------------------

    def assign_materials(self):
        for E in self.elements :
            E.physicalTag = self.MaterialID

    # -------------------------------------------------------------------------------

    def assign_block_numbers(self, BlockCount: int ):

        # count elements per type
        self._count_elements()

        # count block IDs
        BlockCount = self._compute_block_numbers( BlockCount )

        # loop over all elements in group and assign values
        for E in self.elements:
            # Check GMSH number of element
            if E.GMSH == 1 or E.GMSH == 8 :
                E.geometryTag = self._lines
            elif E.GMSH == 2 or E.GMSH == 9 :
                E.geometryTag = self._tris
            elif E.GMSH == 3 or E.GMSH == 10 or E.GMSH == 16 :
                E.geometryTag = self._quads
            elif E.GMSH == 4 or E.GMSH == 11 :
                E.geometryTag = self._tets
            elif E.GMSH == 6 or E.GMSH == 12 or E.GMSH == 18 :
                E.geometryTag = self._pentas
            elif E.GMSH == 5 or E.GMSH == 12 or E.GMSH == 17 :
                E.geometryTag = self._hexs

        # return the counter
        return BlockCount