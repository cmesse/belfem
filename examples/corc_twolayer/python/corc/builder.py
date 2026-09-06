"""Shared mesh container for the discrete shell construction: nodes with a
single global id space, tri3 elements classified on surface tags, and the
registries the cap/assembler stages need. Writes the shell msh2 file the
volume run consumes."""

import numpy as np


class Builder:

    def __init__(self):
        self.nodes = []          # list of np.array([x, y, z]); id = index+1
        self.tris = []           # (surface_tag, n1, n2, n3)
        self.surfaces = {}       # tag -> {"role": str, ...}
        self._next_surface = 0

        # rings[layer_index][row] -> list of node ids per angular slot;
        # "dom" key for the domain tube. Filled by shells, used by caps.
        self.rings = {}

        # cap registries (filled by caps.py)
        self.cap_pairs = []      # (front_node_id, back_node_id)
        self.vertex_nodes = []   # ten node ids: front C,A,B,C,D back C,E,F,G,H

    # ------------------------------------------------------------- entities

    def add_node(self, x: float, y: float, z: float) -> int:
        self.nodes.append(np.array([x, y, z], dtype=float))
        return len(self.nodes)

    def new_surface(self, role: str, **meta) -> int:
        self._next_surface += 1
        self.surfaces[self._next_surface] = dict(role=role, **meta)
        return self._next_surface

    def add_tri(self, tag: int, a: int, b: int, c: int):
        self.tris.append((tag, a, b, c))

    def surface_tags(self, role: str):
        return [t for t, m in self.surfaces.items() if m["role"] == role]

    # ---------------------------------------------------------------- files

    def write_shell_msh(self, path: str):
        with open(path, "w") as f:
            f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            f.write("$Nodes\n{:d}\n".format(len(self.nodes)))
            for i, p in enumerate(self.nodes):
                f.write("{:d} {:.17g} {:.17g} {:.17g}\n".format(
                    i + 1, p[0], p[1], p[2]))
            f.write("$EndNodes\n$Elements\n{:d}\n".format(len(self.tris)))
            for i, (g, a, b, c) in enumerate(self.tris):
                f.write("{:d} 2 2 {:d} {:d} {:d} {:d} {:d}\n".format(
                    i + 1, g, g, a, b, c))
            f.write("$EndElements\n")
