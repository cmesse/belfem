class Physical :

    def __init__(self, dimension: int, phystag: int, label: str ):
        self.Dimension = dimension
        self.physicalTag = phystag
        self.Label = label

    def gmsh(self):
        return f'{self.Dimension:d} {self.physicalTag:d} "{self.Label:s}"'