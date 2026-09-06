"""Cross-section schedule for the CORC cable: layer radii, winding sense,
tape/gap angular partitions, discretization counts and the FEM closure
rule. Pure parameters — no geometry entities.

Conventions (identical to the old generator):
- lengths in mm, tape/solder thickness in µm
- pitch = z-advance per radian, so ω = dφ/ds = 1/pitch and the cable
  length is L = pitch * 2π * numTurns
- layer 0 radius from the packing perimeter: 2πr0 = n * (width + gap)
- odd layers are counter-wound (sense σ = -1)
- tape k of layer l is centered at φ = σ_l (l·π/n + k·2π/n + ω s)
"""

import math


class WindingError(Exception):
    pass


class Tape:
    def __init__(self, layer, index: int, center: float):
        self.layer = layer
        self.index = index
        self.center = center          # center angle at s = 0 (body frame)
        half = 0.5 * layer.tape_angle
        self.phi0 = center - half     # leading edge at s = 0
        self.phi1 = center + half     # trailing edge at s = 0


class Layer:
    def __init__(self, winding, index: int, radius: float):
        self.winding = winding
        self.index = index
        self.radius = radius
        self.sigma = 1 if index % 2 == 0 else -1
        self.tape_angle = winding.tapeWidth / radius
        self.gap_angle = 2.0 * math.pi / winding.numTapesPerLayer - self.tape_angle
        if self.gap_angle <= 1e-12:
            raise WindingError(
                "layer {:d}: tapes overlap (tape angle {:.4g} rad >= "
                "pitch angle {:.4g} rad); reduce tapeWidth or "
                "numTapesPerLayer".format(
                    index, self.tape_angle,
                    2.0 * math.pi / winding.numTapesPerLayer))

        dphi = 2.0 * math.pi / winding.numTapesPerLayer
        self.tapes = [Tape(self, k, index * 0.5 * dphi + k * dphi)
                      for k in range(winding.numTapesPerLayer)]

    def phi(self, tape_or_edge_angle: float, s: float) -> float:
        """Body angle of a cross-section feature at arc length s."""
        return self.sigma * (tape_or_edge_angle + self.winding.omega * s)


class Winding:

    def __init__(self, tapeWidth: float, gap: float, tapeThickness: float,
                 solderThickness: float, numTapesPerLayer: int,
                 numLayers: int, pitch: float, numTurns: float,
                 delta: float, tapeResolution: float):
        if numTapesPerLayer < 1 or numLayers < 1:
            raise WindingError("need at least one tape and one layer")
        if gap <= 0.0:
            raise WindingError("gap must be positive (got {:g})".format(gap))
        if numTurns <= 0.0 or pitch <= 0.0:
            raise WindingError("pitch and numTurns must be positive")

        self.tapeWidth = tapeWidth
        self.gap = gap
        self.numTapesPerLayer = numTapesPerLayer
        self.numLayers = numLayers
        self.pitch = pitch
        self.numTurns = numTurns
        self.delta = delta
        self.tapeResolution = tapeResolution

        # winding rate and cable length (straight-equivalent)
        self.omega = 1.0 / pitch
        self.length = pitch * 2.0 * math.pi * numTurns

        # layer radii: packing perimeter for layer 0, then radial stacking
        r = numTapesPerLayer * (tapeWidth + gap) / (2.0 * math.pi)
        self.interlayer = (tapeThickness + solderThickness) * 0.001  # mm
        self.layers = []
        for l in range(numLayers):
            self.layers.append(Layer(self, l, r))
            r += self.interlayer

        # inter-layer thin-gap sizing law (Spike D): boundary element size
        # on layer interfaces must stay within a few gaps IN BOTH
        # DIRECTIONS — gmsh's boundary recovery fails on long anisotropic
        # triangles threading a thin slab (observed at ds = 16 * gap)
        if numLayers > 1:
            if tapeResolution > 2.0 * self.interlayer:
                raise WindingError(
                    "tapeResolution {:g} mm exceeds 2x the inter-layer gap "
                    "{:g} mm; tet quality in the gap collapses (Spike D). "
                    "Reduce the tape resolution to about the gap size.".format(
                        tapeResolution, self.interlayer))
            ds = self.length / self.s_partition_count()
            if ds > 5.0 * self.interlayer:
                # suggest a delta that survives the round() in
                # s_partition_count: derive it from the integer row count
                rows = int(math.ceil(self.length / (5.0 * self.interlayer)))
                good = numTurns * 360.0 / (rows + 0.5)
                raise WindingError(
                    "s-spacing {:.3g} mm exceeds 5x the inter-layer gap "
                    "{:g} mm; gmsh cannot recover the boundary in the thin "
                    "slab (use delta <= {:.3g} degrees)".format(
                        ds, self.interlayer, good))
            if tapeResolution > self.interlayer:
                print("WARNING: tapeResolution {:g} mm > inter-layer gap "
                      "{:g} mm; expect flat (valid but low-quality) tets "
                      "between layers".format(tapeResolution, self.interlayer))
            if tapeResolution < 0.5 * self.interlayer:
                print("WARNING: tapeResolution {:g} mm is finer than half "
                      "the inter-layer gap {:g} mm; Spike D saw sliver tets "
                      "in this regime".format(tapeResolution, self.interlayer))

    # --------------------------------------------------- solder fraction

    def solder_fractions(self, stackThickness: float):
        """Real-solder fraction of each inter-layer annulus: the BELFEM
        'density correction' for a deck that meshes those volumes as
        solder while the tapes ride on the layer surfaces as thin shells.

        stackThickness is the physical tape stack in um (the sum of the
        BELFEM layers block), as opposed to tapeThickness, which is the
        radial slot the winding schedule reserves per layer.

        Convention: each layer's stack lies radially outward of its shell
        surface -- the same one-slot-per-layer rule the layer radii above
        are stacked by -- so the annulus between layers l and l+1 carries
        exactly layer l's tapes, and the outermost stack extends into the
        outer air volume, where there is no solder mass to correct.

        A helical strip of arc width w covers w * L of cylinder surface
        regardless of pitch, so per unit length the annulus holds
        pi (r1^2 - r0^2) of meshed solder and n * w * t of real tape and
        the cable length drops out. For bent centerlines this is the
        straight-equivalent value: bending stretches tape and annulus by
        the same local factor to first order.

        Returns one fraction per annulus, i.e. per BELFEM volume 2..nL
        (assembler ordering: 1 = inner, 2..nL = interlayer, nL+1 = outer).
        """
        t = stackThickness * 0.001  # um -> mm
        fractions = []
        for l in range(self.numLayers - 1):
            r0 = self.layers[l].radius
            r1 = self.layers[l + 1].radius
            annulus = math.pi * (r1 * r1 - r0 * r0)
            tape = self.numTapesPerLayer * self.tapeWidth * t
            if tape >= annulus:
                raise WindingError(
                    "annulus {:d}: tape stack volume {:.4g} mm^2 per unit "
                    "length exceeds the annulus cross-section {:.4g} mm^2; "
                    "stackThickness {:g} um cannot be physical".format(
                        l, tape, annulus, stackThickness))
            fractions.append(1.0 - tape / annulus)
        return fractions

    # -------------------------------------------------------- FEM closure

    def tape_shift(self):
        """Number of tape positions every tape advances over the cable
        length, when that is a whole number: q = numTurns * numTapesPerLayer.
        With q whole the partition of EVERY layer at s=L coincides with its
        partition at s=0 (the pattern is n-fold symmetric, and counter-wound
        layers advance by -q, which is the same set), so the caps close with
        zero rotation and tape k's far end is tape (k+q mod n)'s near end --
        the periodic twist: with n = 3 and numTurns = 1/3 (or 4/3, ...)
        sideset 4 continues into 5, 5 into 6, 6 into 4, and the inner layer
        likewise. Returns None when q is not whole."""
        q = self.numTurns * self.numTapesPerLayer
        if abs(q - round(q)) > 1e-9 * max(1.0, abs(q)):
            return None
        return int(round(q)) % self.numTapesPerLayer

    def cap_rotation(self):
        """The single in-plane rotation that maps every layer's partition
        at s=0 to its partition at s=L. Raises if none exists (then the
        caps cannot be periodic images and FEM mode must refuse)."""
        wl = self.omega * self.length  # = 2 pi numTurns
        if self.numLayers == 1:
            return wl % (2.0 * math.pi)
        # n-fold closure: whole tape positions advanced -> zero rotation,
        # tapes permute (see tape_shift)
        if self.tape_shift() is not None:
            return 0.0
        # counter-wound layers: need wl == -wl (mod 2 pi), i.e. integer or
        # half-integer turns (tolerance relative to the turn count)
        half_turns = 2.0 * self.numTurns
        if abs(half_turns - round(half_turns)) > 1e-9 * max(1.0, half_turns):
            raise WindingError(
                "counter-wound layers close only for integer or "
                "half-integer numTurns (got {:.12g}); adjust numTurns or "
                "pitch".format(self.numTurns))
        return math.pi * (round(half_turns) % 2)

    # ---------------------------------------------------- discretizations

    def s_partition_count(self) -> int:
        """Number of s-intervals shared by all strips of one layer (and the
        domain tube): from delta (degrees of winding per row, as today)."""
        return max(1, int(round(self.numTurns * 360.0 / self.delta)))

    def n_phi(self, layer_index: int, kind: str) -> int:
        """Single source of truth for angular interval counts (Grok #13).
        kind is 'tape' or 'gap'. The count also guarantees each angular
        slot spans <= pi/2 so cap circle arcs can always be chunked below
        the built-in kernel's pi limit (matters for numTapesPerLayer=1,
        where a single gap can span nearly 2 pi)."""
        L = self.layers[layer_index]
        angle = L.tape_angle if kind == "tape" else L.gap_angle
        arclen = angle * L.radius
        return max(2, int(math.ceil(arclen / self.tapeResolution)),
                   int(math.ceil(angle / (0.5 * math.pi))))
