# C1 continuation of the tabulated jc and n below 10 mT

**Date:** 2026-09-23
**Purpose:** Session record of the DR-189 fix: `JcFunctionDatabase` continues n and REBCO jc from
Bmin = 10 mT down to B = 0, C1 and monotone, with exact derivatives, without touching the tables
**Module:** `src/physics/materials`, `src/physics/database`, `tests/physics`, `share/material`,
`share/python/database`
**AIs involved:** Claude Opus 5.5 (plan v1–v4, implementation, measurement, mutations), Codex
(`gpt-6-astra`/high: plan rounds 1–3, code jury, confirmation round), Grok `grok-4.7`/high (plan
rounds 1–4, code jury)

## Summary

DR-189 was an evaluator defect, not a table defect: below ~55 K the reader's low-field bridge fell
back to a flat continuation under 10 mT, with a slope kink at 10 mT. The cause was a one-point
calibration of the measured self-field curve at 77.5 K plus a guard that returned the clamp
whenever the measured level sat below the table's edge value.

Christian's rulings, taken through the plan rounds: keep the measured, angle-free anchor with the
Kim finite slope (DR-128); C1 and monotone everywhere wins over the measured anchor near Tc; n
continues from the table alone as a parabola with zero slope at B = 0; the HDF5 databases are not
changed; the derivatives must be exact; below Tmin the clamp stays for now.

## What changed

| File | Change |
|---|---|
| `cl_Database.{hpp,cpp}` | tensor index map of the element nodes; exact mixed derivatives `evaluate_derivxy` / `evaluate_derivyz`; certified Bernstein bounds of df/dy and of f on the face y = min(1), full tensor net, orders 1–3; `dimension()` |
| `cl_JcFunction_Database.hpp` | load refusals in `open_checked` ( 3 axes, angle coverage, face slope bound ≤ 1e-10, n floor 1.01 ) under a `unique_ptr`; measured curve projected to its running maximum from high T, Fritsch & Butland 1984 knot slopes, C1 Hermite with linear continuation; edge struct with V, S and their T and angle tangents; jc Hermite to J = smax( a, V − S Bmin/3 ) with the origin slope 2Δ − S for β ≤ 2 and 0 above; n parabola; all four evaluators exact on every branch |
| `test_Database.cpp` | mixed derivatives at orders 1–3 with distinct steps; the face bounds on Codex's HEX27 counterexample and two controls |
| `test_JcDatabaseTable.cpp` | rewritten: refusals ( folded, half, 2-D, rising inside an element, n floor ), closed forms on synthetic tables ( n parabola, jc quadratic to its anchor, projected curve with Fritsch–Butland slopes ), and contracts over all six shipped tables x { jc, n } ( monotone in B and T, C1 at 10 mT by exact slopes, angle-free anchor where the curve dominates, tangents vs central differences, n(0) > 1.01, clamp below Tmin ). The two `DISABLED_` pins are gone |
| docs | `share/material/README.md` window bullets; `howtoread.py` docstrings ( it clamps below Bmin, the solver does not ); two materials docs; Fritsch & Carlson 1980 and Fritsch & Butland 1984 in `doc/literature_references.md` |

## Measurements

| Quantity | Value |
|---|---|
| rise with \|B\| ( six tables, jc and n, B from 0 ) | ≤ 2.1e-14 relative |
| rise of jc with T ( 2 K to Tmax + 1, seven fields ) | ≤ 1.8e-14 |
| exact vs central-difference tangents, below and above Bmin | ≤ 2e-8 of the scale |
| n(0) minimum | 1.29 ( superox near Tc ) |
| anchor lift J/a − 1, 20–84 K | ≤ 2.4 % ( fysc ) |
| anchor lift, 3 K ( below the lowest measured knot, 12.5–17.5 K ) | up to 54 % ( fesc ) |
| anchor lift near Tc | up to 0.92 ( sp-ap ) and 5.2 ( fesc, where the curve's continuation reaches ~0 ) |
| share of sampled points with β > 2 ( origin slope 0 ) | 16–63 % per table |

The earlier estimate of the lift ( ≤ 1.3 % over 20–84 K ) held for four tables; fysc reaches 2.4 %.
The β > 2 branch is routine, not the rare exception the plan called it, and the comment says so.

The secant test at 10 mT is useless where the table is flat in B: at 3.5 K the slope of sp-ap is
~1e-11 of the scale, below the rounding noise of any difference quotient, so the old relative
secant jump read 2.0 from noise alone. The C1 test compares the exact slopes instead.

The isotonic projection never changes a shipped table inside its T window; the only rising knot is
superox at 95 K, above Tmax = 91 K. A synthetic table exercises it.

## Code jury

Blind, on a restricted diff ( `tmp/ai_exchange/review_jc_c1_code.diff` ), with nine pre-registered
claims. Grok: all nine confirmed, no P0/P1. Codex: C1, C3–C8 confirmed, C2 and C9 refuted in part,
verdict "revise". Reconciled and fixed:

| Finding | From | Fix |
|---|---|---|
| angle snapped onto a window end still returned the endpoint's ∂/∂θ ( value is flat there ) | Codex P1, Grok P2 | `angle_is_snapped`; `deval_dbeta` returns 0 in the 1e-6 band. Pre-existing, inherited by the continuation |
| NaN table values passed both face certificates ( `std::max( bound, NaN )` keeps the bound ) | Codex P1 | the reductions return NaN on a non-finite coefficient, which every caller comparison refuses |
| smooth-maximum width `ε( \|a\| + \|b\| )/2` kinks where the continued curve crosses 0, and J then rises by ~1e-7 in T | Codex P2 | width `ε b`; b ≥ V > 0, and J_T ≤ 0 whenever a_T, b_T ≤ 0 |
| `release()` precedes the rest of construction | Codex P2 | `mDatabase` is a `unique_ptr` member |
| 2-axis refusal survived its mutant ( a debug-only assert refused too ) | Grok P2 | the test pins the refusal message |
| Fritsch & Carlson locator | Codex ( from a transcription of the paper ) | §4 ( region S₁, the square ); Fritsch & Butland Eq. 5 |
| wording: n-floor message, a test comment's −0.075 ( is −0.15 ), resistivity page "continuous slope", Python docstring ( "B ≥ Bmin", "none of these are read by the solver" ), two review-history remarks | both | corrected |

Confirmation round ( Codex astra/high, on the fix delta ): all five fixes confirmed, including the
proof that J_a, J_b > 0 with the new width so J_T ≤ 0, and the Fritsch & Carlson §4 and Fritsch &
Butland Eq. 5 locators against copies of the papers. Two follow-ups, both applied: the Python
reader's `deval_dbeta` gets the same snap mask ( C++/Python parity ), and a direct Database test
pins the value bound's NaN guard, which the constructor test cannot reach ( the slope refusal fires
first ). The user-facing passages had a Codex language sweep ( luna/medium ); it narrowed "C1" to
"C1 in B".

Not changed: the β = 2 seam ( ruled ); the Bernstein bound's conservatism ( both: the safe direction;
subdivide only if a real table trips it ); `pow( 0, 0 )` in the Bernstein matrix ( 1 under C Annex F );
calibration before the projection ( Grok open risk: a future file rising through 77.5 K would leave
its knot off the declared level ).

## Mutations

Twenty-seven mutants, each caught by the test that claims it: edge slope dropped from the jc and the
n field derivative; the Hermite S term dropped from the value; the f_Tu and f_uθ terms; e_q; s0_q;
a_T; n floor loosened; clamp zero removed; T-window zero removed; n parabola made linear; face
slope refusal off; Bernstein matrix set to identity; mixed derivative on the wrong axes ( xy, yz );
face bound in reference units; projection removed; Fritsch–Butland weights swapped; the measured
curve's T axis reversed; the field axis reversed; the anchor on the requirement alone; the 3-axis
refusal removed; the NaN refusal removed; the |a|-based width; the angle-snap zero removed; the value bound's NaN guard removed.

## Open

- The Bernstein bound does not subdivide, so it refuses a face slope that is negative everywhere
  but concave ( −0.1 − ξ² ). The shipped tables pass; subdivision is the fix if a table ever trips it.
- `howtoread.py` does not reproduce the continuation; its docstring and the README say so. The
  copies embedded in the HDF5 files cannot be updated without changing the files.
- "Fritsch & Carlson 1980, §4" rests on Codex's reading of a transcription; the literature library
  does not carry the paper.
- Below Tmin the value is still held ( ruled for now ).

## Gates

`make check` 17/17 in both trees ( Armadillo `build/`, Blaze `cmake-build-debug/` ) before the jury
and again after the jury fixes; the physics suite 145 passed, 1 skipped ( magnesia has no ρ, by
design ).
