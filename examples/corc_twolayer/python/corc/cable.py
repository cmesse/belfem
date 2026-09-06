"""Cable facade: the one object main.py talks to.

    C = Cable()
    C.tapeWidth = 4
    ...
    C.centerline = Lame(76.1, 25, 1.55)   # or StraightLine() (default)
    C.build()                              # shells + caps + volumes + checks
    C.save("corc.msh")                     # relative: the run directory, where
    C.write_topology("corc_topology.conf") # the deck looks for the mesh

Show-off mode (closed curves like Trefoil): build() generates the bent
tape shells only — no caps, volumes, periodicity or topology.
"""

import math
import os

from corc.curve import StraightLine
from corc.frame import Frame
from corc.winding import Winding
from corc.builder import Builder
from corc.shells import ShellGenerator
from corc.caps import CapGenerator
from corc.assembler import Assembler
from corc.postprocess import Postprocessor
from corc import belfem


class Cable:

    def __init__(self):
        # parameters and units exactly as the old generator
        self.tapeWidth = 12          # mm
        self.pitch = 10              # mm z-advance per radian
        self.numTurns = 1
        self.gap = 1                 # mm
        self.tapeThickness = 100     # µm (radial slot reserved per layer)
        self.solderThickness = 10    # µm
        self.stackThickness = None   # µm physical tape stack (BELFEM layers
                                     # sum) for the solder density correction;
                                     # None falls back to tapeThickness
        self.numTapesPerLayer = 3
        self.numLayers = 1           # >1 needs tapeResolution ~ interlayer gap
        self.delta = 30              # degrees of winding per element row
        self.domainRadius = 100      # mm

        self.tapeResolution = 0.5    # mm
        self.innerResolution = 1.0   # mm
        self.domainResolution = 10.0  # mm

        self.centerline = StraightLine()
        self.workdir = "/tmp"

        self._mesh = None
        self._post = None

    # ------------------------------------------------------------- building

    def build(self):
        curve = self.centerline

        if isinstance(curve, StraightLine):
            # straight cable: the curve length follows from the winding
            curve.set_length(self.pitch * 2.0 * math.pi * self.numTurns)
            numTurns = self.numTurns
        else:
            # bent cable: the centerline sets the length
            numTurns = curve.length / (self.pitch * 2.0 * math.pi)

        self.winding = Winding(
            tapeWidth=self.tapeWidth, gap=self.gap,
            tapeThickness=self.tapeThickness,
            solderThickness=self.solderThickness,
            numTapesPerLayer=self.numTapesPerLayer,
            numLayers=self.numLayers, pitch=self.pitch,
            numTurns=numTurns, delta=self.delta,
            tapeResolution=self.tapeResolution)

        outer = self.winding.layers[-1].radius
        if not curve.closed and self.domainRadius <= outer + self.gap:
            raise ValueError(
                "domainRadius {:g} mm must clear the outermost layer "
                "radius {:.4g} mm".format(self.domainRadius, outer))

        self.frame = Frame(curve)

        self.shellgen = ShellGenerator(
            self.winding, self.frame, self.domainRadius,
            self.domainResolution)

        self.builder = Builder()

        if curve.closed:
            # show-off mode: shells only
            self.shellgen.build(self.builder)
            self.shellgen.check_radius(self.builder)
            print("show-off shells: {:d} nodes, {:d} triangles".format(
                len(self.builder.nodes), len(self.builder.tris)))
            return

        self.frame.check_domain_admissibility(self.domainRadius)

        self.shellgen.build(self.builder)
        worst = self.shellgen.check_radius(self.builder)
        print("shell radius exactness: {:.3e} mm".format(worst))

        capgen = CapGenerator(self.shellgen, self.innerResolution, self.workdir)
        cap = capgen.mesh_cap()
        capgen.instantiate(self.builder, cap)

        assembler = Assembler(self.builder, self.workdir, self.domainResolution)
        nodes, tets, shell_tris, old_to_gmsh = assembler.run()

        self._post = Postprocessor(self.builder, self.winding, self.frame,
                                   self.shellgen)
        self._mesh = self._post.build_final(nodes, tets, shell_tris, old_to_gmsh)
        self._post.validate(self._mesh)

    def density_corrections(self):
        """BELFEM 'density correction' per interlayer volume (2..nL), if
        those volumes are meshed as solder. See Winding.solder_fractions."""
        t = self.stackThickness if self.stackThickness is not None \
            else self.tapeThickness
        return self.winding.solder_fractions(t)

    # --------------------------------------------------------------- output

    def save(self, path: str):
        if self._mesh is None:
            raise RuntimeError("call build() first")
        self._mesh.save(path)

    def write_topology(self, path: str):
        if self._post is None:
            raise RuntimeError("call build() first (FEM mode)")
        belfem.write_topology(
            path,
            volume_count=self._post.volume_count,
            tape_physicals=self._post.tape_physicals,
            boundary_physical=self._post.boundary_physical,
            front_physical=self._post.front_physical,
            back_physical=self._post.back_physical,
            solder_corrections=self.density_corrections()
                if self.winding.numLayers > 1 else None,
            solder_front_physicals=getattr(self._post, "solder_front_physicals", None),
            solder_back_physicals=getattr(self._post, "solder_back_physicals", None))

    def save_stl(self, path: str):
        """Write the tape shells as ASCII STL (show-off / inspection)."""
        b = self.builder
        tape_tags = set(b.surface_tags("tape"))
        with open(path, "w") as f:
            f.write("solid tapes\n")
            for tag, a, bb, c in b.tris:
                if tag not in tape_tags:
                    continue
                import numpy as np
                p1 = b.nodes[a - 1]
                p2 = b.nodes[bb - 1]
                p3 = b.nodes[c - 1]
                n = np.cross(p2 - p1, p3 - p1)
                norm = float(np.linalg.norm(n))
                if norm > 0:
                    n = n / norm
                f.write("  facet normal {:.6e} {:.6e} {:.6e}\n".format(*n))
                f.write("    outer loop\n")
                for p in (p1, p2, p3):
                    f.write("      vertex {:.6e} {:.6e} {:.6e}\n".format(*p))
                f.write("    endloop\n  endfacet\n")
            f.write("endsolid tapes\n")

    def print(self):
        if self._post is None:
            return
        p = self._post
        print("Volumes: 1..{:d}".format(p.volume_count))
        print("Tapes:", p.tape_physicals)
        print("DomainBoundary:", p.boundary_physical)
        print("DomainFront:", p.front_physical)
        print("DomainBack:", p.back_physical)
        if self.winding.numLayers > 1:
            for l, f in enumerate(self.density_corrections()):
                print("density correction (solder in Volume_{:d}): "
                      "{:.4f}".format(l + 2, f))
