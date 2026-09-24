# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Winding schedule for the tape-wrapped flat core: layer cross-sections,
lay angles, tape/gap partitions of the perimeter, discretization counts
and the periodic closure rule. Pure parameters, no geometry entities.

Conventions:
- lengths in mm, the layer gap in um, the lay angle in degrees
- the cable axis is z; the lay angle theta is measured from it, so a
  tape advances ds/dz = tan theta along the perimeter and one turn takes
  the pitch L = P / tan theta in z
- the pitch is common to all layers. Layer l sits on a stadium with the
  same flat length and the edge radius r + l * gap, so its perimeter and
  hence its lay angle differ: tan theta_l = P_l / L. layAngle is the
  layer-0 value
- odd layers are counter-wound (sense sigma = -1 relative to layer 0)
- tape k of a layer is centered at s = sigma * (k * P / n + z * tan theta)
  with s = 0 at the top-face center
- the tape width wt is measured across the tape; at fixed z the band
  spans wt / cos theta of perimeter and wt / sin theta of z
"""

import math


class WindingError(Exception):
    pass


class Tape:
    def __init__(self, layer, index: int, center: float):
        self.layer = layer
        self.index = index
        self.center = center            # center arc length at z = 0
        half = 0.5 * layer.tape_span
        self.s0 = center - half         # leading edge at z = 0
        self.s1 = center + half         # trailing edge at z = 0


class Layer:
    def __init__(self, winding, index: int, radius: float):
        from rutherford.stadium import Stadium

        self.winding = winding
        self.index = index
        self.radius = radius
        self.stadium = Stadium(winding.ws, radius)
        self.perimeter = self.stadium.perimeter
        self.sigma = winding.sense * (1 if index % 2 == 0 else -1)

        # common pitch, layer-specific lay angle
        self.tan_theta = self.perimeter / winding.pitch
        self.theta = math.atan(self.tan_theta)
        self.cos_theta = math.cos(self.theta)
        self.sin_theta = math.sin(self.theta)

        self.tape_span = winding.tapeWidth / self.cos_theta
        self.sector = self.perimeter / winding.numTapes
        self.gap_span = self.sector - self.tape_span
        if self.gap_span <= 0.0:
            raise WindingError(
                "layer {:d}: tapes overlap ({:d} x {:.4g} mm of perimeter "
                "per tape >= {:.4g} mm); reduce tapeWidth or numTapes, or "
                "raise the lay angle".format(
                    index, winding.numTapes, self.tape_span,
                    self.perimeter))

        self.tapes = [Tape(self, k, k * self.sector)
                      for k in range(winding.numTapes)]

    def coverage(self) -> float:
        """Fraction of the layer surface under tape."""
        return self.winding.numTapes * self.tape_span / self.perimeter

    def bending_radius(self) -> float:
        """Radius of the tape centerline's normal curvature on the round
        edges: the surface curvature 1/r taken in the direction of the
        oblique tape path (Euler's theorem, kappa_n = sin^2 theta / r)."""
        return self.radius / (self.sin_theta * self.sin_theta)

    def geodesic_torsion(self) -> float:
        """Twist rate of the tape on the round edges, sin cos / r."""
        return self.sin_theta * self.cos_theta / self.radius


class Winding:

    def __init__(self, tapeWidth: float, layAngle: float, edgeRadius: float,
                 cableWidth: float, numTapes: int, numLayers: int,
                 layerGap: float, numPitches: float, tapeResolution: float,
                 gapResolution: float, sense: int = 1):
        if numTapes < 1 or numLayers < 1:
            raise WindingError("need at least one tape and one layer")
        if tapeWidth <= 0.0:
            raise WindingError("tapeWidth must be positive")
        if not 0.0 < layAngle < 90.0:
            raise WindingError(
                "layAngle {:g} deg must lie strictly between 0 and 90 "
                "(measured from the cable axis)".format(layAngle))
        if numPitches <= 0.0:
            raise WindingError("numPitches must be positive")
        if tapeResolution <= 0.0 or gapResolution <= 0.0:
            raise WindingError("resolutions must be positive")
        if sense not in (1, -1):
            raise WindingError("sense must be +1 or -1")
        if numLayers > 1 and layerGap <= 0.0:
            raise WindingError("layerGap must be positive for stacked layers")

        self.tapeWidth = tapeWidth
        self.layAngle = layAngle
        self.edgeRadius = edgeRadius
        self.cableWidth = cableWidth
        self.numTapes = numTapes
        self.numLayers = numLayers
        self.numPitches = numPitches
        self.tapeResolution = tapeResolution
        self.gapResolution = gapResolution
        self.sense = sense

        self.ws = cableWidth - 2.0 * edgeRadius
        if edgeRadius <= 0.0 or self.ws < 0.0:
            raise WindingError(
                "edgeRadius {:g} mm must be positive and at most half the "
                "cableWidth {:g} mm".format(edgeRadius, cableWidth))

        # layer-0 perimeter fixes the pitch; the layers stack outward
        theta0 = math.radians(layAngle)
        P0 = 2.0 * self.ws + 2.0 * math.pi * edgeRadius
        self.pitch = P0 / math.tan(theta0)
        self.length = numPitches * self.pitch
        self.interlayer = layerGap * 0.001   # um -> mm

        self.layers = []
        for l in range(numLayers):
            self.layers.append(Layer(self, l, edgeRadius + l * self.interlayer))

        if self.tape_shift() is None:
            raise WindingError(
                "numPitches * numTapes = {:.12g} is not a whole number: the "
                "end planes are not translated images of each other. Use a "
                "whole number of pitches, or a fraction whose product with "
                "numTapes is whole (the half-pitch closure by a 180 degree "
                "rotation is not implemented)".format(numPitches * numTapes))

        for L in self.layers:
            if L.gap_span < tapeResolution:
                print("WARNING: layer {:d}: the gap between tapes ({:.3g} mm "
                      "of perimeter) is narrower than tapeResolution {:g} mm; "
                      "expect sliver triangles in the gap strip".format(
                          L.index, L.gap_span, tapeResolution))

        # inter-layer thin-gap sizing (the CORC generator's Spike D rule):
        # boundary element size on layer interfaces must stay within a few
        # gaps, or gmsh's boundary recovery fails on long anisotropic
        # triangles threading the slab
        if numLayers > 1:
            if tapeResolution > 2.0 * self.interlayer:
                raise WindingError(
                    "tapeResolution {:g} mm exceeds 2x the inter-layer gap "
                    "{:g} mm; tet quality in the gap collapses. Reduce the "
                    "tape resolution to about the gap size".format(
                        tapeResolution, self.interlayer))
            if gapResolution > 5.0 * self.interlayer:
                raise WindingError(
                    "gapResolution {:g} mm exceeds 5x the inter-layer gap "
                    "{:g} mm; gmsh cannot recover the boundary in the thin "
                    "slab".format(gapResolution, self.interlayer))
            if tapeResolution > self.interlayer:
                print("WARNING: tapeResolution {:g} mm > inter-layer gap "
                      "{:g} mm; expect flat (valid but low-quality) tets "
                      "between layers".format(tapeResolution, self.interlayer))
            if tapeResolution < 0.5 * self.interlayer:
                print("WARNING: tapeResolution {:g} mm is finer than half "
                      "the inter-layer gap {:g} mm; expect sliver tets in "
                      "this regime".format(tapeResolution, self.interlayer))

    # -------------------------------------------------------- FEM closure

    def tape_shift(self):
        """Number of tape positions every tape advances over the cable
        length, when that is a whole number: q = numPitches * numTapes.
        Then the partition of EVERY layer at z = length coincides with
        its partition at z = 0 (counter-wound layers advance by -q, the
        same set), the end planes are pure translations of each other,
        and tape k's far end is tape (k + q mod n)'s near end. Returns
        None when q is not whole."""
        q = self.numPitches * self.numTapes
        if abs(q - round(q)) > 1e-9 * max(1.0, abs(q)):
            return None
        return int(round(q)) % self.numTapes

    def conductor_groups(self):
        """Galvanically connected tape sets, as lists of tape indices in
        the layer-major, tape-order numbering of the sidesets. Across the
        periodic planes tape k continues into tape k + q (q = numPitches *
        numTapes), so a layer splits into gcd(q, n) conductors, each one
        tape of n / gcd(q, n) times the length. Every group is ONE current
        condition and consumes one cohomology generator; listing its
        curves separately asks for more conditions than the topology has
        generators."""
        n = self.numTapes
        q = int(round(self.numPitches * n))
        g = math.gcd(q, n)
        groups = []
        for l in range(self.numLayers):
            for c in range(g):
                groups.append([l * n + k for k in range(n) if k % g == c])
        return groups

    # ---------------------------------------------------- discretizations

    def z_partition_count(self) -> int:
        """Number of z-intervals shared by all strips: the step along the
        layer-0 tape path is tapeResolution, so dz = tapeResolution *
        cos theta_0. With the slot width across the tape also equal to
        tapeResolution, a strip quad is a parallelogram with the angle
        90 - theta between its sides, and the short-diagonal split gives
        near-equilateral triangles at 30 degrees."""
        dz = self.tapeResolution * self.layers[0].cos_theta
        return max(1, int(math.ceil(self.length / dz)))

    def n_s(self, layer_index: int, kind: str) -> int:
        """Slot count of one tape ('tape') or one gap ('gap') strip of a
        layer, at least 2 so a strip is never a single sliver."""
        L = self.layers[layer_index]
        if kind == "tape":
            return max(2, int(math.ceil(L.tape_span / self.tapeResolution)))
        return max(2, int(math.ceil(L.gap_span / self.gapResolution)))
