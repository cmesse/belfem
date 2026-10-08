# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Pancake facade: the one object main.py talks to.

    C = Pancake()
    C.numTurns = 200
    ...
    C.build()                                    # sheet + pack + air + checks
    C.save("pancake.msh")                        # relative: the run directory
    C.write_topology("pancake_topology.conf")

A single tape of width tapeWidth is wound numTurns times around the z
axis from innerRadius to outerRadius; the solder between the turns is
a solid that touches itself and the contact surface is the tape, a thin
shell. A solder sliver sits inside the first turn. The inner lead climbs
the bore wall to the domain cap on a solder bar, the outer lead bends
easy-way in the plane and leaves radially to the domain cylinder on a
solder bar; the bars' end faces are the terminals. Lengths in mm.
"""

import os
import tempfile

import numpy as np

from pancake.sheet import Sheet
from pancake.pack import Pack
from pancake.boundary import Boundary
from pancake.assembler import Assembler
from pancake.postprocess import validate, ValidationError
from pancake.mshio import write_msh2
from pancake import belfem


class Pancake:

    def __init__(self):
        # geometry in mm
        self.tapeWidth = 4.0
        self.innerRadius = 7.0
        self.outerRadius = 17.0
        self.numTurns = 200
        self.innerLeadDirection = 1     # +1: the inner lead climbs to z = +H
        self.leadRadius = 3.0           # departure circle inside the bore,
                                        # tangent to the spiral at its inner end
        self.leadWrap = 90.0            # degrees wrapped on that circle before
                                        # the riser bends the tape toward z
        self.riserRadius = None         # hard-way bend radius in the unrolled
                                        # plane of the small cylinder;
                                        # None = 2 x tapeWidth
        self.innerBarThickness = 1.0    # solder bar inside the inner lead
        self.bendRadius = 5.0           # easy-way bend of the outer lead
        self.outerBarThickness = None   # solder bar inside the outer lead;
                                        # None = wedgeGap, so the wedge's end
                                        # face is the bar's start face
        self.domainRadius = 80.0
        self.domainHalfHeight = 50.0

        # resolutions in mm; the pack is always one cell thick
        self.tapeResolution = 1.0       # along the spiral and across the width
        self.leadResolution = 1.0       # along the two leads
        self.domainResolution = 8.0     # far field
        self.packAirResolution = 0.2    # air size on the pack faces. A node's
                                        # element count is a byte in BELFEM, and
                                        # a coarse air tet fanning over the 50 um
                                        # strips overflows it
        self.airGrowth = 0.5            # size gradient away from the pack: a
                                        # node at distance d wants size
                                        # packAirResolution + airGrowth * d.
                                        # 0.3 breeds slivers, 1.0 lets the first
                                        # air layer fan over hundreds of strips
        self.wedgeGap = None            # solder wedge closes at this gap;
                                        # None = tapeResolution

        self.tapeCurrent = 100.0        # A, written into the deck
        self.optimizeNetgen = False     # segfaulted on embedded surfaces once
        self.algorithm3d = 1            # gmsh 3-D algorithm for the air:
                                        # Delaunay (1) or HXT (10); if the one
                                        # chosen leaves slivers, the other is
                                        # tried once
        self.workdir = None             # gmsh scratch; None = fresh temp dir

        self._xyz = None

    # ------------------------------------------------------------- building

    def build(self):
        wedge = self.tapeResolution if self.wedgeGap is None else self.wedgeGap
        self.sheet = Sheet(
            tapeWidth=self.tapeWidth, innerRadius=self.innerRadius,
            outerRadius=self.outerRadius, numTurns=self.numTurns,
            bendRadius=self.bendRadius, domainRadius=self.domainRadius,
            domainHalfHeight=self.domainHalfHeight,
            tapeResolution=self.tapeResolution,
            leadResolution=self.leadResolution, wedgeGap=wedge,
            innerLeadDirection=self.innerLeadDirection,
            leadRadius=self.leadRadius, leadWrap=self.leadWrap,
            riserRadius=2.0 * self.tapeWidth if self.riserRadius is None else self.riserRadius,
            innerBarThickness=self.innerBarThickness,
            outerBarThickness=wedge if self.outerBarThickness is None else self.outerBarThickness)
        self.print_sheet()

        if self.workdir is None:
            self.workdir = tempfile.mkdtemp(prefix="pancake_")
        os.makedirs(self.workdir, exist_ok=True)

        self.pack = Pack(self.sheet)
        boundary = Boundary(self.sheet, self.pack, self.domainResolution, self.workdir)
        boundary.build_tube()
        boundary.build_caps()

        surfaces = self.pack.shell_surfaces() + boundary.surfaces
        S = dict(surfaces)
        # the terminals are faces of the bars, not of the air
        air_surfaces = [(t, tris) for t, tris in surfaces
                        if t not in (Pack.TERMINAL_IN, Pack.TERMINAL_OUT)]
        self.blockNames = self.pack.blockNames + ["Air"]
        self.airBlock = len(self.blockNames)
        self.nTape = self.pack.nTapePieces
        self._sidesets = [("Tape_{:d}".format(p + 1),
                           self.pack.sheetTris[self.pack.sheetPiece == p])
                          for p in range(self.nTape)]
        self._sidesets += [("TerminalIn", S[Pack.TERMINAL_IN]),
                           ("TerminalOut", S[Pack.TERMINAL_OUT]),
                           ("DomainBoundary", S[Pack.TUBE]),
                           ("DomainTop", S[Pack.CAP_TOP]),
                           ("DomainBottom", S[Pack.CAP_BOTTOM]),
                           ("SolderAir", self.pack.solder_air_tris())]

        # the air is meshed with the chosen algorithm; slivers are a
        # matter of luck with the near-coplanar pack faces, so the other
        # algorithm gets one try before giving up
        other = 10 if self.algorithm3d == 1 else 1
        for attempt, algorithm in enumerate((self.algorithm3d, other)):
            assembler = Assembler(
                self.sheet, boundary.xyz, air_surfaces, self.domainResolution,
                self.workdir, self.packAirResolution, self.airGrowth,
                self.optimizeNetgen, algorithm, Pack.BAR_SURFACES)
            xyz, air = assembler.run()
            self._xyz = xyz
            self._tets = np.vstack((self.pack.tets, air))
            self._tetBlock = np.concatenate((self.pack.tetBlock,
                                             np.full(air.shape[0], self.airBlock)))
            try:
                validate(self.sheet, xyz, self._tets, self._tetBlock, self.blockNames,
                         self.airBlock, self.pack.sheetTris, self.pack.sheetPiece,
                         dict(self._sidesets), self.pack.innerPiece,
                         self.sheet.nSpiralTapePieces - 1)
                break
            except ValidationError as e:
                if attempt == 1 or "degenerate" not in str(e):
                    raise
                print("air mesh with algorithm {:d} rejected ({}); retrying with "
                      "algorithm {:d}".format(algorithm, e, other))
        print("mesh: {:d} nodes, {:d} tets, {:d} tape facets, {:d} blocks, {:d} sidesets".format(
            xyz.shape[0], self._tets.shape[0], self.pack.sheetTris.shape[0],
            len(self.blockNames), len(self._sidesets)))

    # --------------------------------------------------------------- output

    def physical_ids(self):
        names = self.blockNames + [n for n, _ in self._sidesets]
        return {name: i + 1 for i, name in enumerate(names)}

    def save(self, path):
        if self._xyz is None:
            raise RuntimeError("call build() first")
        ids = self.physical_ids()
        groups = []
        physicals = []
        for b, name in enumerate(self.blockNames, start=1):
            groups.append((4, ids[name], ids[name], self._tets[self._tetBlock == b]))
            physicals.append((3, ids[name], name))
        for name, tris in self._sidesets:
            groups.append((2, ids[name], ids[name], tris))
            physicals.append((2, ids[name], name))
        write_msh2(path, self._xyz, groups, physicals)

    def write_topology(self, path):
        if self._xyz is None:
            raise RuntimeError("call build() first")
        belfem.write_topology(path, ids=self.physical_ids(),
                              n_solder=self.sheet.nPieces, n_tape=self.nTape,
                              inner_piece=self.pack.innerPiece,
                              outer_piece=self.sheet.nSpiralTapePieces - 1,
                              tape_current=self.tapeCurrent)

    def print_sheet(self):
        S = self.sheet
        tets, facets = S.estimate()
        print("winding: {:d} turns from r = {:g} to {:g} mm, pitch {:.4g} um, tape "
              "{:g} mm wide; spiral length {:.4g} m".format(
                  S.N, S.ri, S.ro, 1e3 * S.pitch, S.b,
                  1e-3 * np.pi * S.N * (S.ri + S.ro)))
        print("columns: {:d} per turn, step {:.3g} mm at r = {:g} and {:.3g} mm at "
              "r = {:g}; {:d} levels across the width; {:d} solder pieces, {:d} tape "
              "pieces".format(S.nTheta, S.ri * S.dtheta, S.ri, S.ro * S.dtheta, S.ro,
                              S.nZ, S.nPieces, S.nSpiralTapePieces + 1))
        print("inner lead: departure circle r1 = {:g} mm about ({:g}, 0), wrap {:g} deg, "
              "riser R_h = {:g} mm, {:.3g} rad on the small cylinder in all, {:d} columns "
              "({:d} in the solder wedge), up to z = {:+g} mm, on a {:g} mm solder bar".format(
                  S.r1, S.C1[0], np.degrees(S.psi1), S.Rh, S.leadWrap, S.nU, S.nWi,
                  S.sigma * S.H, S.tIn))
        print("outer lead: bend {:.4g} deg at R = {:g} mm, then {:.4g} mm radial to "
              "r = {:g} at {:.4g} deg; solder wedge over {:d} columns ({:.3g} mm), "
              "then a {:g} mm solder bar".format(
                  np.degrees(S.bendAngle), S.Rb, S.straightLength, S.Rdom,
                  np.degrees(S.exitAngle), S.nW, S.wedgeLength, S.tLead))
        print("estimate: {:d} structured tets, {:d} tape facets".format(tets, facets))

    def print(self):
        ids = self.physical_ids()
        print("blocks: Solder_1..Solder_{:d} = 1..{:d}, LeadIn {:d}, LeadOut {:d}, "
              "Air {:d}".format(self.sheet.nPieces, self.sheet.nPieces,
                                ids["LeadIn"], ids["LeadOut"], ids["Air"]))
        print("sidesets: Tape_1..Tape_{:d} = {:d}..{:d}".format(
            self.nTape, ids["Tape_1"], ids["Tape_{:d}".format(self.nTape)]))
        for name in ("TerminalIn", "TerminalOut", "DomainBoundary", "DomainTop", "DomainBottom",
                     "SolderAir"):
            print("{:>15s}: {:d}".format(name, ids[name]))
