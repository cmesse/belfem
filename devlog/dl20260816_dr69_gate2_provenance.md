# DR-69 Gate 2 Closed: the Raw Robinson Sweep Carries the Asymmetry

**Date:** 2026-08-16
**Purpose:** Record the provenance check that unblocks half of the DR-69
signed-β implementation
**Module:** maxwell, physics/materials

## What was checked

The 2026-08-13 ruling (signed dot for table-based Jc) was gated on two
pre-registered checks. Gate 2 asked whether the θ↔180°−θ asymmetry in
`sp-ap.hdf5` exists in the RAW Robinson (SuperCurrent) angular Ic data,
or was introduced by the smoothing/export — in which case the current
fold would be accidentally right.

The raw dataset is public: figshare 10.6084/m9.figshare.4256624.v3
(Wimbush/Strickland, Robinson Research Institute), SuperPower Advanced
Pinning wire SCS12050-AP M3-1386-3, sample SP066, measured Jan–Feb 2021,
CC-BY. Its angle column is defined as the sample **normal** with respect
to the applied field — the exact `bn_angle` quantity. Downloaded the
0.3 T and 0.5 T angle-dependence workbooks and evaluated the 80 K and
85 K sheets (the register's 82–86 K / 0.3–0.5 T reference window).

## Result: PASSED, verified

- **The asymmetry is in the raw measurement.** Ic(θ)/Ic(180°−θ) at
  0.3/0.5 T, 80/85 K runs 1.27–1.39, peaking near θ ≈ 80°. BELFEM's
  spline-evaluated reference numbers (1.42 Claude / 1.461 Codex at 84 K,
  interpolated between these sheets) sit right on top of the raw trend.
- **Artifact control:** the raw sweeps extend to 240°, which permits a
  decisive physics check — Ic must be 180°-periodic under field
  reversal. It is, to ≤2.6% over all 13 available (θ, θ+180°) pairs per
  sheet, while the 90°-mirror deviates 22–39% at the same angles.
  Field-reversal symmetry holds and the mirror symmetry does not: the
  asymmetry is real asymmetric pinning, not rig drift and not the
  smoothed export.

"Verified" is used deliberately: this is an executable check that ran,
on the raw public data, with a built-in control.

Residue, noted in the register: `sst-1.hdf5`'s raw counterpart was not
matched — figshare hosts three Shanghai Superconductor characterisations
and the product designation behind `sst-1` needs Christian's say before
the same check runs there. The gate's quantitative reference was sp-ap
and the ruling stands on it.

## What remains before implementation

- **Gate 1** (Christian's run): constant-jc/n control deck in which the
  J/Jc asymmetry must vanish; doubles as the bit-identical regression
  check for the DR-07 tangent patch.
- Then the implementation campaign per the ruling: signed
  `acos(dot(n,b)/|b|)` for table-based Jc only; the material must
  advertise whether its Jc(θ) is even (Kim stays folded); tape-normal
  sign graduates to load-bearing deck input; a future β-tangent needs
  its own signed ∂θ/∂q chain.

Analysis artifacts (downloaded workbooks + scripts) live in the session
scratchpad only; the numbers above and the dataset DOI are the durable
record.
