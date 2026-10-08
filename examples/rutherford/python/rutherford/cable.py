# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Cable facade: the one object main.py talks to.

    C = Cable()
    C.tapeWidth = 4
    ...
    C.build()                                  # shells + caps + volume + checks
    C.save("rutherford.msh")                   # relative: the run directory,
    C.write_topology("rutherford_topology.conf")  # where the deck looks

The cable axis is z. The tape is wound around a flat core with a stadium
cross-section (width cableWidth along x, thickness 2 edgeRadius along y)
at the lay angle layAngle measured from the axis. The core itself is not
meshed as a conductor: its section is the air block 1.
"""

import math
import os
import tempfile

from rutherford.winding import Winding
from rutherford.builder import Builder
from rutherford.shells import ShellGenerator
from rutherford.caps import CapGenerator
from rutherford.assembler import Assembler
from rutherford.postprocess import Postprocessor
from rutherford import belfem


class Cable:

    def __init__(self):
        # geometry in mm and degrees, the layer gap in um
        self.tapeWidth = 4           # 2, 4 or 6 in the reference study
        self.layAngle = 30           # from the cable axis; 30 or 45
        self.edgeRadius = 0.9        # half the core thickness
        self.cableWidth = 12         # placeholder: the core width is not
                                     # documented in the reference study
        self.numTapes = 5            # tapes per layer, equally spaced
        self.numLayers = 1           # odd layers are counter-wound
        self.layerGap = 200          # um, radial step between layers
        self.numPitches = 0.2        # cable length in pitches; times
                                     # numTapes must be whole
        self.sense = 1               # +1 or -1: winding sense of layer 0
        self.domainRadius = 50       # mm, air cylinder around the axis

        # resolutions in mm
        self.tapeResolution = 0.4    # along and across the tape
        self.gapResolution = 0.8     # across the bare core between tapes
        self.innerResolution = 0.5   # cap mesh inside the core section
        self.domainResolution = 4.0  # air cylinder and cap outer piece

        self.tapeCurrent = 30.0      # A per tape; the written amplitude is
                                     # the total through one conductor group

        self.workdir = None          # gmsh scratch; None = a fresh temp dir

        self._mesh = None
        self._post = None

    # ------------------------------------------------------------- building

    def build(self):
        self.winding = Winding(
            tapeWidth=self.tapeWidth, layAngle=self.layAngle,
            edgeRadius=self.edgeRadius, cableWidth=self.cableWidth,
            numTapes=self.numTapes, numLayers=self.numLayers,
            layerGap=self.layerGap, numPitches=self.numPitches,
            tapeResolution=self.tapeResolution,
            gapResolution=self.gapResolution, sense=self.sense)

        outer = self.winding.layers[-1].stadium
        if self.domainRadius <= 0.5 * outer.width + self.gapResolution:
            raise ValueError(
                "domainRadius {:g} mm must clear the outermost layer's "
                "half width {:.4g} mm".format(self.domainRadius,
                                               0.5 * outer.width))

        if self.workdir is None:
            self.workdir = tempfile.mkdtemp(prefix="rutherford_")
        os.makedirs(self.workdir, exist_ok=True)

        self.shellgen = ShellGenerator(self.winding, self.domainRadius,
                                       self.domainResolution)
        self.builder = Builder()
        self.shellgen.build(self.builder)
        worst = self.shellgen.check_surface(self.builder)
        print("shell surface exactness: {:.3e} mm".format(worst))

        capgen = CapGenerator(self.shellgen, self.innerResolution, self.workdir)
        cap = capgen.mesh_cap()
        capgen.instantiate(self.builder, cap)

        assembler = Assembler(self.builder, self.workdir, self.domainResolution)
        nodes, tets, shell_tris, old_to_gmsh = assembler.run()

        self._post = Postprocessor(self.builder, self.winding, self.shellgen)
        self._mesh = self._post.build_final(nodes, tets, shell_tris, old_to_gmsh)
        self._post.validate(self._mesh)

    # --------------------------------------------------------------- output

    def save(self, path: str):
        if self._mesh is None:
            raise RuntimeError("call build() first")
        self._mesh.save(path)

    def write_topology(self, path: str):
        if self._post is None:
            raise RuntimeError("call build() first")
        belfem.write_topology(
            path,
            volume_count=self._post.volume_count,
            tape_physicals=self._post.tape_physicals,
            boundary_physical=self._post.boundary_physical,
            front_physical=self._post.front_physical,
            back_physical=self._post.back_physical,
            conductor_groups=self.winding.conductor_groups(),
            tape_current=self.tapeCurrent)

    def save_stl(self, path: str):
        """Write the tape shells as ASCII STL (inspection)."""
        import numpy as np
        b = self.builder
        tape_tags = set(b.surface_tags("tape"))
        with open(path, "w") as f:
            f.write("solid tapes\n")
            for tag, a, bb, c in b.tris:
                if tag not in tape_tags:
                    continue
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
        w = self.winding
        print("core: width {:g} mm, thickness {:g} mm, flat length {:.4g} mm, "
              "perimeter {:.4g} mm".format(
                  w.cableWidth, 2.0 * w.edgeRadius, w.ws,
                  w.layers[0].perimeter))
        print("pitch: {:.4g} mm ({:.4g} mm by the face-only rule "
              "2 w / tan theta); cable length {:.4g} mm, {:d} rows".format(
                  w.pitch, 2.0 * w.cableWidth / math.tan(math.radians(w.layAngle)),
                  w.length, w.z_partition_count()))
        for L in w.layers:
            print("layer {:d}: edge radius {:.4g} mm, lay angle {:.4g} deg, "
                  "sense {:+d}, {:d} tape(s) covering {:.1f} % of the surface; "
                  "tape bending radius on the edges {:.3g} mm, twist rate "
                  "{:.3g} rad/mm".format(
                      L.index, L.radius, math.degrees(L.theta), L.sigma,
                      w.numTapes, 100.0 * L.coverage(), L.bending_radius(),
                      L.geodesic_torsion()))
        if self._post is None:
            return
        p = self._post
        print("Volumes: 1..{:d}".format(p.volume_count))
        print("Tapes:", p.tape_physicals)
        print("DomainBoundary:", p.boundary_physical)
        print("DomainFront:", p.front_physical)
        print("DomainBack:", p.back_physical)
