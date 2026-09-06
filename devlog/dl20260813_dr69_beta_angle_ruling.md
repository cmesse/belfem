# DR-69 Beta-Angle Convention: Christian's Ruling

**Date:** 2026-08-13
**Purpose:** Record the DR-69 ruling on the tape-normal jc angle convention, plus the
database probe that widened the finding to sst-1. Discussion session; no source edits.
**Module:** `src/fem/kernel` (Calculator angle helpers), `src/physics/materials`

## Background

DR-69 (found earlier on 2026-08-13, three-AI audit in
`tmp/ai_exchange/bn_angle_signed_audit.md`) established that the solve folds the
tape-normal angle to [0,90°] via `abs(dot(n,b))` in `Calculator::bn_angle`, while the
thin-shell postprocessor uses the signed dot ([0,180°]) and the measured sp-ap.hdf5
table is asymmetric about 90° — so the solve samples only one lobe of measured
asymmetric pinning. The row was parked "needs Christian's ruling on the convention."

## The ruling (Christian, 2026-08-13)

There are **two essential definitions of beta** and they are settled separately:

1. **ReBCO (tape): angle between tape normal and B.** For **table-based Jc**
   (`JcFunctionDatabase`) the solve switches to the signed dot —
   `acos( dot(n,b)/|b| )`, the [0,180°] range the databases natively span. Nothing
   beyond 180° is meaningful (it is the angle between two vectors); the one genuine
   sign choice is the **orientation of n itself**, which maps θ↔180°−θ and therefore
   becomes load-bearing material/deck input that must match the measurement
   convention (film side vs substrate side). **Kim-type analytic laws stay folded** —
   they are even in θ, so the material must advertise whether its Jc(θ) is even; one
   helper cannot serve both.

2. **Metals: angle between current and field driving the Kohler law.** A different
   definition (`bj_angle`, dependency `angleBxJ`), and **no action**: Pippard's
   angular interpolation `Along·cos²β + Atrans·sin²β` (`cl_Material_Copper.cpp`,
   `Copper::kohler`) is even about 90°, so the fold in `bj_angle_3d` is provably
   lossless. The two definitions keep their separate helpers.

3. **Bulk HTS: no angle.** There is no tape normal, so `beta_dummy()` = π/2 stays.
   Academic case (real bulk superconductors are rare in our decks); the only residue
   is a doc footnote that an angle-dependent table read by a bulk material silently
   uses its 90° column.

Implementation remains **gated** on the two pre-registered gates in the DR-69 row:
(1) the constant-jc/n control run in which the J/Jc asymmetry must vanish, and
(2) provenance of the θ↔180°−θ asymmetry in the raw Robinson (SuperCurrent) sweep,
ruling out a smoothing/export artifact.

## Database probe run this session

Both `share/material/sp-ap.hdf5` and `share/material/sst-1.hdf5` were read directly
(h5py):

- **Identical grids:** T = 4–92 K in 4 K steps (23 points), log₁₀B = −2…1 in 0.1
  steps (31 points), θ = **0–180° in exactly 1° steps** (181 points) — confirming
  the tables natively carry the unfolded range the ruling adopts.
- **sst-1 is also asymmetric about 90°**, which the original DR-69 audit had not
  established (it probed sp-ap only). At the raw-coefficient level ~95% of grid
  points differ by >20% between θ and 180°−θ in both files. Caveat, stated in the
  register too: the stored `values` are **order-2 spline control points**, not nodal
  samples, so these percentages are qualitative until evaluated through BELFEM's
  spline; sp-ap's spline-evaluated ratios (1.42 Claude / 1.461 Codex) remain the
  quantitative reference. Confidence: high on the grid ranges, medium-high on the
  sst-1 asymmetry magnitude.

Also confirmed in code this session: every metal Kohler consumer routes through the
`angleBxJ` dependency, and the 2D `bj_angle_2d` returns π/2 unconditionally (current
out of plane) — consistent with the "no action on metals" leg of the ruling.

## Register changes

`todo/debt_register.md` DR-69: "Ruling needed" → **RULED** with the three-part
convention above; status column now reads "ruled 2026-08-13, implementation NOT
started, gated"; ruling-session addendum records the sst-1 finding. The row stays
open and unstruck — the code change, its β-tangent consequence (signed ∂θ/∂q chain,
never the `bj_angle` one), and both gates are all outstanding.
