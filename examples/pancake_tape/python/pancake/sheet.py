# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Footprint of the tape sheet and the column indexing the structured
pack mesh is built on.

Frame: t is the tape direction, n the face normal, b = t x n the
binormal across the width. In the winding n is radial and b = e_z.

Column k of the pack runs, in this order:

    sliver     k in [-nTheta, 0)          the spiral continued one turn
                                          inward: inner boundary of the
                                          solder sliver (not tape)
    spiral     k in [0, kEnd]             r = ri + p theta / 2 pi,
                                          theta = k dtheta
    wedge      k in [kEnd, kEnd + nW]     the easy-way bend, sampled at
                                          the polar angles j dtheta so the
                                          last turn pairs with it radially
    lead       k in (kEnd + nW, M]        rest of the bend, then the
                                          radial straight to r = R_dom

with p = (ro - ri) / N the turn pitch and dtheta = 2 pi / nTheta a
uniform step, so column k and column k + nTheta sit at the same polar
angle and one turn apart. The solder pack is the set of cells between
column pairs (k, k+1) and (k + nTheta, k + nTheta + 1) for k in
[-nTheta, nCells).

The inner lead is a separate branch that leaves the spiral at the
departure point D = column 0 (r = ri, theta = 0): the tape follows a
circle of radius r1 tangent to the spiral at D about the center
C1 = (ri - r1, 0), wraps by leadWrap into the bore, then bends in the
unrolled plane of that small cylinder with radius R_h until t runs along
z, and climbs to the cap. The first columns of the branch are sampled at
the polar angles of the sliver's inner boundary, so the sliver cells
pair with them radially and the thin gap between the departing tape and
the first turn is solder by construction until it reaches wedgeGap. A
two-cell-wide solder bar rides inside the branch on the small cylinder.
"""

import math

import numpy as np


class SheetError(Exception):
    pass


class Sheet:

    def __init__(self, tapeWidth, innerRadius, outerRadius, numTurns,
                 bendRadius, domainRadius, domainHalfHeight,
                 tapeResolution, leadResolution, wedgeGap,
                 innerLeadDirection, leadRadius, leadWrap, riserRadius,
                 innerBarThickness, outerBarThickness):
        if tapeWidth <= 0.0:
            raise SheetError("tapeWidth must be positive")
        if not 0.0 < innerRadius < outerRadius:
            raise SheetError("need 0 < innerRadius < outerRadius")
        if numTurns < 1:
            raise SheetError("numTurns must be at least 1")
        if bendRadius <= 0.5 * tapeWidth:
            raise SheetError(
                "bendRadius {:g} mm must exceed half the tape width "
                "{:g} mm".format(bendRadius, 0.5 * tapeWidth))
        if domainRadius <= outerRadius + bendRadius:
            raise SheetError("domainRadius must clear outerRadius + bendRadius")
        if domainHalfHeight <= 0.5 * tapeWidth:
            raise SheetError("domainHalfHeight must exceed half the tape width")
        if tapeResolution <= 0.0 or leadResolution <= 0.0 or wedgeGap <= 0.0:
            raise SheetError("resolutions and wedgeGap must be positive")
        if innerLeadDirection not in (1, -1):
            raise SheetError("innerLeadDirection must be +1 or -1")
        if not 0.0 < leadRadius < innerRadius:
            raise SheetError("leadRadius must lie strictly between 0 and innerRadius")
        if riserRadius <= 0.5 * tapeWidth:
            raise SheetError("riserRadius must exceed half the tape width")
        if riserRadius >= domainHalfHeight:
            raise SheetError("riserRadius must be below domainHalfHeight")
        if innerBarThickness <= 0.0 or outerBarThickness <= 0.0:
            raise SheetError("bar thicknesses must be positive")
        if innerBarThickness >= leadRadius:
            raise SheetError("innerBarThickness must be below leadRadius")
        if leadWrap <= 0.0:
            raise SheetError("leadWrap must be positive")

        self.b = float(tapeWidth)
        self.ri = float(innerRadius)
        self.ro = float(outerRadius)
        self.N = int(numTurns)
        self.Rb = float(bendRadius)
        self.Rdom = float(domainRadius)
        self.H = float(domainHalfHeight)
        self.tapeResolution = float(tapeResolution)
        self.leadResolution = float(leadResolution)
        self.wedgeGap = float(wedgeGap)
        self.sigma = int(innerLeadDirection)
        self.r1 = float(leadRadius)
        self.psi1 = math.radians(float(leadWrap))
        self.Rh = float(riserRadius)
        self.tIn = float(innerBarThickness)
        self.tLead = float(outerBarThickness)

        self.pitch = (self.ro - self.ri) / self.N
        if self.tIn <= 1.2 * self.pitch:
            raise SheetError(
                "innerBarThickness {:g} mm must exceed the solder sliver, about "
                "1.1 x the turn pitch {:.4g} mm".format(self.tIn, self.pitch))
        self.nTheta = int(math.ceil(2.0 * math.pi * self.ro / self.tapeResolution))
        self.dtheta = 2.0 * math.pi / self.nTheta
        self.half = self.nTheta // 2

        self.kEnd = self.N * self.nTheta
        self._build_columns()

        # z levels of the pack (tape width)
        self.nZ = max(1, int(math.ceil(self.b / self.tapeResolution)))
        self.zPack = np.linspace(-0.5 * self.b, 0.5 * self.b, self.nZ + 1)

        self.nWinding = (self.N - 1) * self.nTheta
        self.nCells = self.nWinding + self.nW
        self.M = self.cols.shape[0] - 1        # last column index
        self._build_sliver_and_outer_bar()
        self._build_inner_lead()

    # ------------------------------------------------------------ curves

    def spiral(self, theta):
        r = self.ri + self.pitch * theta / (2.0 * math.pi)
        return np.column_stack((r * np.cos(theta), r * np.sin(theta)))

    def _bend_point(self, phi):
        C, n, t, R = self._C, self._n, self._t, self.Rb
        return C + R * (-n * math.cos(phi) + t * math.sin(phi))

    def _bend_dir(self, phi):
        return self._n * math.sin(phi) + self._t * math.cos(phi)

    @staticmethod
    def _polar(P):
        return math.atan2(P[1], P[0])

    def _build_columns(self):
        nT, N, dth = self.nTheta, self.N, self.dtheta
        theta = np.arange(0, self.kEnd + 1) * dth
        cols = [self.spiral(theta)]
        E = cols[0][-1]
        E[:] = (self.ro, 0.0)   # exact outer end

        # spiral tangent at the outer end, right-hand normal, bend center
        thE = 2.0 * math.pi * N
        t = np.array([self.pitch / (2.0 * math.pi) * math.cos(thE) - self.ro * math.sin(thE),
                      self.pitch / (2.0 * math.pi) * math.sin(thE) + self.ro * math.cos(thE)])
        t /= np.linalg.norm(t)
        n = np.array([t[1], -t[0]])
        self._t, self._n = t, n
        self._C = E + self.Rb * n

        # bend angle: the lead leaves radially, i.e. its direction is
        # parallel to the position vector at the end of the bend
        def cross(beta):
            P = self._bend_point(beta)
            d = self._bend_dir(beta)
            return P[0] * d[1] - P[1] * d[0]
        lo, hi = 0.0, 0.5 * math.pi
        if cross(lo) <= 0.0 or cross(hi) >= 0.0:
            raise SheetError("no radial exit within a quarter bend")
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if cross(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        self.bendAngle = 0.5 * (lo + hi)
        self.bendEnd = self._bend_point(self.bendAngle)
        self.exitAngle = self._polar(self.bendEnd)
        self.leadEnd = self.bendEnd / np.linalg.norm(self.bendEnd) * self.Rdom
        self.straightLength = self.Rdom - float(np.linalg.norm(self.bendEnd))

        # wedge: sample the bend at the polar angles of the last turn's
        # columns until the gap to that turn reaches wedgeGap
        polar_max = self.exitAngle
        wedge = []
        phis = []
        j = 0
        while True:
            j += 1
            target = j * dth
            if target >= polar_max:
                raise SheetError(
                    "the wedge does not close before the lead turns radial "
                    "(gap below wedgeGap {:g} mm); raise bendRadius or lower "
                    "wedgeGap".format(self.wedgeGap))
            lo, hi = phis[-1] if phis else 0.0, self.bendAngle
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if self._polar(self._bend_point(mid)) < target:
                    lo = mid
                else:
                    hi = mid
            phi = 0.5 * (lo + hi)
            phis.append(phi)
            wedge.append(self._bend_point(phi))
            if self._gap(j, phis) >= self.wedgeGap:
                break
        self.nW = j
        self.wedgePhis = phis
        cols.append(np.array(wedge))

        # rest of the lead: the remaining bend, then the radial straight
        arc = self.Rb * (self.bendAngle - phis[-1])
        L = arc + self.straightLength
        nR = max(1, int(math.ceil(L / self.leadResolution)))
        rest = []
        for q in range(1, nR + 1):
            s = L * q / nR
            if s <= arc:
                rest.append(self._bend_point(phis[-1] + s / self.Rb))
            else:
                rest.append(self.bendEnd + (s - arc) * self._bend_dir(self.bendAngle))
        rest[-1] = self.leadEnd.copy()
        self.nRest = nR
        cols.append(np.array(rest))
        self.cols = np.vstack(cols)
        self.wedgeLength = self.Rb * phis[-1]
        self.leadLength = self.Rb * self.bendAngle + self.straightLength

    def _gap(self, j, phis):
        P = self._bend_point(phis[j - 1])
        S = self.spiral(np.array([2.0 * math.pi * (self.N - 1) + j * self.dtheta]))[0]
        return float(np.linalg.norm(P - S))

    # ------------------------------------------------- sliver, outer bar

    def _build_sliver_and_outer_bar(self):
        nT, dth = self.nTheta, self.dtheta
        # sliver: the spiral continued one turn inward; column -nTheta is
        # at (ri - p, 0), one pitch inside the departure point
        k = np.arange(-nT, 0)
        theta = k * dth
        r = self.ri + self.pitch * theta / (2.0 * math.pi)
        self.vcols = np.column_stack((r * np.cos(theta), r * np.sin(theta)))

        # outer lead bar: inner (coil) side of the lead beyond the wedge.
        # Its first column is the last turn's point at the wedge end, so
        # the wedge's end face becomes the bar's start face; its last
        # column sits on the domain cylinder one bar thickness to the left
        self.kLead = self.kEnd + self.nW
        ks = np.arange(self.kLead, self.M + 1)
        pts = np.empty((len(ks), 2))
        pts[0] = self.cols[self.kLead - nT]
        for i, kk in enumerate(ks[1:-1], start=1):
            d = self.cols[kk + 1] - self.cols[kk - 1]
            d /= np.linalg.norm(d)
            left = np.array([-d[1], d[0]])
            pts[i] = self.cols[kk] + self.tLead * left
        self.barExitAngle = self.exitAngle + self.tLead / self.Rdom
        pts[-1] = self.Rdom * np.array([math.cos(self.barExitAngle),
                                        math.sin(self.barExitAngle)])
        self.leadBar = pts
        self.nLeadBar = len(ks) - 1            # cells along the outer bar

    # ---------------------------------------------------------- inner lead

    def _small_circle(self, psi):
        """Point on the departure circle at angle psi about C1 (psi = 0
        at D, negative going away from the spiral)."""
        return self.C1 + self.r1 * np.array([math.cos(psi), math.sin(psi)])

    def _unrolled(self, lam):
        """Centerline of the inner lead in the unrolled (s, z) plane of
        the small cylinder at backward arc length lam from D, with the
        in-plane width direction nu (nu = e_z on the departure arc)."""
        La = self.r1 * self.psi1
        Lb = La + 0.5 * math.pi * self.Rh
        sg = self.sigma
        if lam <= La:
            return np.array([-lam, 0.0]), np.array([0.0, 1.0])
        if lam <= Lb:
            phi = (lam - La) / self.Rh
            s = -La - self.Rh * math.sin(phi)
            z = sg * self.Rh * (1.0 - math.cos(phi))
            return np.array([s, z]), np.array([sg * math.sin(phi), math.cos(phi)])
        z = sg * (self.Rh + (lam - Lb))
        return np.array([-La - self.Rh, z]), np.array([sg, 0.0])

    def lead_point(self, sz, radius):
        """Map unrolled (s, z) on the small cylinder to 3-D at the given
        radius about C1 (s is arc length on the tape's cylinder r1)."""
        psi = sz[..., 0] / self.r1
        return np.stack((self.C1[0] + radius * np.cos(psi),
                         self.C1[1] + radius * np.sin(psi), sz[..., 1]), axis=-1)

    def _build_inner_lead(self):
        dth = self.dtheta
        self.C1 = np.array([self.ri - self.r1, 0.0])
        La = self.r1 * self.psi1
        Lend = La + 0.5 * math.pi * self.Rh + (self.H - self.Rh)
        wrap = (La + self.Rh + 0.5 * self.b) / self.r1
        if wrap >= 2.0 * math.pi - 0.05:
            raise SheetError(
                "the inner lead wraps {:.3g} rad around the small cylinder; lower "
                "leadWrap or riserRadius, or raise leadRadius".format(wrap))

        # wedge columns: departure points at the polar angles -j dtheta of
        # the sliver's inner boundary, until the gap to that boundary
        # reaches wedgeGap
        lams = [0.0]
        j = 0
        psi_prev = 0.0
        while True:
            j += 1
            target = -j * dth
            lo, hi = -self.psi1, psi_prev
            if self._polar(self._small_circle(lo)) > target:
                raise SheetError(
                    "the inner wedge does not close within leadWrap; raise "
                    "leadWrap or lower wedgeGap")
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if self._polar(self._small_circle(mid)) > target:
                    hi = mid
                else:
                    lo = mid
            psi = 0.5 * (lo + hi)
            psi_prev = psi
            lams.append(-self.r1 * psi)
            gap = float(np.linalg.norm(self._small_circle(psi) - self.vcols[-j]))
            if gap >= self.wedgeGap:
                break
        self.nWi = j
        # the rest of the branch at leadResolution, ending exactly on the cap
        L0 = lams[-1]
        n = max(1, int(math.ceil((Lend - L0) / self.leadResolution)))
        lams += [L0 + (Lend - L0) * q / n for q in range(1, n + 1)]
        lams = np.array(lams)
        self.nU = len(lams) - 1                 # cells along the branch

        # nodes of the branch in the unrolled plane, forward order: index
        # q = nU - j runs from the cap (q = 0) to the departure point
        sz = np.empty((self.nU + 1, self.nZ + 1, 2))
        for j, lam in enumerate(lams):
            c, nu = self._unrolled(lam)
            q = self.nU - j
            sz[q] = c[None, :] + self.zPack[:, None] * nu[None, :]
        sz[self.nU, :, 0] = 0.0                  # D exactly
        sz[self.nU, :, 1] = self.zPack
        self.leadSZ = sz
        self.leadWrap = wrap
        self.leadLengthIn = float(Lend)

    def col(self, k):
        """Column k of the sheet (k >= 0) or of the sliver's inner
        boundary (k < 0)."""
        k = np.asarray(k)
        return np.where((k < 0)[:, None], self.vcols[np.clip(k, -self.nTheta, -1) + self.nTheta],
                        self.cols[np.clip(k, 0, self.M)])

    def piece(self, k):
        """Half-turn piece of cell or segment k (k >= -nTheta): the
        sliver is pieces 0 and 1, the first turn 2 and 3, ..."""
        k = np.asarray(k) + self.nTheta
        return 2 * (k // self.nTheta) + (k % self.nTheta >= self.half)

    @property
    def nPieces(self):
        return int(self.piece(self.nCells - 1)) + 1

    @property
    def nSpiralTapePieces(self):
        return int(self.piece(self.M - 1)) - 2 + 1

    def estimate(self):
        tets = 6 * ((self.nCells + self.nTheta) * self.nZ + 2 * self.nU * self.nZ
                    + self.nLeadBar * self.nZ)
        facets = 2 * (self.M + self.nU) * self.nZ
        return tets, facets
