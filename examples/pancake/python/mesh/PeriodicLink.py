class PeriodicLink:

    def __init__(self, dim: int, slave_tag: int, master_tag: int,
                 affine=None, node_pairs=None):
        self.dim = dim
        self.slave_tag = slave_tag
        self.master_tag = master_tag
        self.affine = affine
        self.node_pairs = node_pairs if node_pairs is not None else []

    def with_node_pairs(self, node_pairs):
        return PeriodicLink(self.dim, self.slave_tag, self.master_tag,
                            self.affine, node_pairs)
