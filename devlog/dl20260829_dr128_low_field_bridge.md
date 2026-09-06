# DR-128: Reader-Side Self-Field Bridge for the Jc Tables

**Date:** 2026-08-29
**Purpose:** Remove the 10 mT clamp floor of the HTS jc/n lookup tables — the last open row of
the jjc-noise cluster — with a measured-anchor transition in the table reader.
**Module:** physics/materials, io

## The defect

The jc tables' B axis starts at log10 B = −2 because log10(0) is undefined, and
`JcFunctionDatabase` clamped below it: no field dependence at all across the entire self-field
operating range of mT-scale decks (clamp error up to 23 % near Tc), and a ∂/∂B Newton-tangent
kink exactly at 10 mT, crossed in production at t ≈ 2.7 s.

## The fix (Christian's Kohler instinct, literature-grounded)

No clamp — a transition function, as `Copper::create_kohler` does for magnetoresistance, in
LINEAR B where B = 0 is an ordinary point:

    jc(B) = J0(T) + c(θ,T)·B + b(θ,T)·B²      below Bmin, C1 at the join

- **J0 is measured, not extrapolated.** Two designs that extrapolated it from the table's edge
  failed the data: measured jc at B = 0 is angle-independent to 0.25–0.73 % while any
  edge-derived level keeps ~2.6–4 % angular shape (ill-conditioned by ~4 % value vs ~200 %
  slope spreads). The rebuilt tables embed their raw rows; the loader bins the 946 B = 0 rows
  by temperature, averages over stage angle, and **calibrates** against the table's own
  `meta/Icw_77p5K_sf_A_per_m / t_eff_m` — no vendor unit assumption (the raw column is A/cm;
  assuming A/m had left the bridge silently inert behind the monotonicity guard).
- **The linear term is kept.** Kim-family route (Riva 2023 Eq. 2, Denis 2026 Eq. 11,
  Messe 2023 §2.6 naming Kim 1962): d jc/d|B| at 0 is finite; the anisotropy lives inside
  |B_eff| and vanishes with the field — exactly what the measured angular collapse shows.
- **Monotone limiter** (code audit): the exact-match quadratic is non-monotone at 60 % of
  (T,θ) — worst hump 6.1 % near Tc — so where c > 0 the bridge falls to the pure parabola:
  value match kept, monotone by construction, tangents still differentiate the value.
- **Consistent degradation:** inconsistent (T,θ) slices (measured level below table edge; low
  T, tiny lift) keep the historical clamp WITH its tangents; induced jc(0) spread ≤ 0.14 %,
  inside measurement scatter. Tables without `/source` keep the exact historical clamp — old
  files bit-unchanged.
- θ/T tangents below Bmin are endpoint-exact with a declared mixed-partial deferral (bound
  |Sθ|·Bmin/4, Grok-verified); `n` reads its own column, unscaled.

Also fixed en route: `load_strings_from_file` read vlen strings as ASCII while h5py writes
UTF-8 — HDF5 refuses that conversion, so any Python-written string dataset was unreadable
(the runtime crash that surfaced it). The memtype cset is now copied from the dataset.

## Ceremony and evidence

Two plan rounds (axis-extension design rejected by both voices → Christian's transition
ruling), one code round (both voices request-changes; all findings actioned same evening —
the tangent-consistency guard, the monotone limiter, schema bounds, eval-routed edge value
for lift-factor compatibility, self_field assert). Full trail:
`tmp/ai_exchange/review_dr128_spap_low_field.md`; design history
`todo/dr128_spap_low_field_axis.md`.

**Executed gates** (standalone probes against the shipped sp-ap table, MPI initialized):
C0/C1 join ~3e-7 / ~3e-5; jc(B→0) angle-free to 1.3e-12; d jc/dB(0) = −3.0e10 finite
negative; monotone at 181 angles × 3 temperatures; FD-vs-analytic dB tangent to 1e-9;
jc(0)/jc(Bmin) = 1.0278 at 77 K — the measured lift.

**Production gate:** tapestack3d (fully coupled, warm restart, DR-128 binary verified by its
loader string) crossed 10 mT at t = 2.725 s with frame ratios ×1.0465 → ×1.0463 through the
crossing — the Newton kink is gone by construction and nothing kicked: zero cuts, zero
resets, 2 Picard iterates per step, run continuing past t = 3.575 s.

DR-128 struck and archived. This closes the jjc-noise cluster (DR-126, DR-127, DR-128).

## Files

- `src/physics/materials/cl_JcFunction_Database.hpp` — the bridge, loader, calibration
- `src/io/hdf5_tools.hpp` — vlen-UTF-8 string reads
- `todo/debt_register.md` / `debt_register_closed.md`, `todo/dr128_spap_low_field_axis.md`
