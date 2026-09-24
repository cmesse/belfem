# Piecewise power law: one clamped Bézier root for the value and both tangents

**Date:** 2026-09-24
**Purpose:** Nightly pipeline 41 (job 150, RELEASE/Armadillo) failed `HtsLaws.PiecewiseDegenerateBlendStaysFinite`: `rho_piecewise` and `drho_piecewise_dT` return NaN at J = j3 for n = 4. Close it with one helper for the Bézier root that every path shares.
**Module:** `src/physics/materials` (`powerlaws.hpp`), `tests/physics/test_MaterialDerivatives.cpp`
**AIs involved:** Claude (diagnosis relay from the CI session, plan, implementation); Codex + Grok (plan jury, code jury)
**Status:** CLOSED 2026-09-24 — landed in `aa468c3` behind both juries; R1–R5 done, A5–A6 applied, code-round fixes F1–F4 applied and re-gated (test_physics 151/151 release and debug Armadillo, Blaze syntax-clean). Owed: R6's shared-tree `check-fast` and the manual pipeline (Christian).

> **Scope guards:**
> - No change to the knots, the control points, the `|a| > eps` assert, or the root's formula. The change is the clamp: a discriminant below zero by roundoff is projected to zero, where the old code returned NaN.
> - `drho_piecewise_dJ` keeps its early return on the discriminant (landed 2026-09-23); it only switches to the helper for the root it computes after that guard.
> - No new `input.conf` key. No doc other than the header comments that describe the changed lines.

## 1. Defect

With `p = n − 1` and `q = mNff = 3`, the discriminant of the Bézier root at the upper knot is
`b² + ac = R² (1/p − 1/q)²` (comment at `powerlaws.hpp:2300`), which is exactly 0 for n = 4. In
floating point it lands a few ulps either side of 0. On the negative side `std::pow( disc, 0.5 )` is
NaN.

| path | root site | guarded? |
|---|---|---|
| 8 × `rho_piecewise` | `powerlaws.hpp:469,538,613,688,756,823,889,954`: `(b+std::pow(b*b+a*c,0.5))/a` | **no** → NaN value |
| 8 × `drho_piecewise_dJ` | `:1392..1953`, two `pow(…,0.5)` each | yes: early return on `!( tDisc > eps (b²+1) )` at `:1389` etc., before the root |
| 2 × `drho_piecewise_dT` | `:2262`, `:2403`: `tSqrt = pow( tDisc, 0.5 )` | **no**: the guard at `:2333` comes *after* `tSqrt`, `tParam`, `rhoFF` and `tFrozen` are computed from it, and returns `tFrozen`, which is already NaN. The comment at `:2326-2332` believes the negated comparison suffices; it protects the division, not the value. |

Only the RELEASE/Armadillo job failed; the debug jobs and RELEASE/Blaze passed. Which side of 0
the roundoff lands on depends on the compiled evaluation order, so this is expected to vary by
configuration. The C++ release flags are `-O2` without `-march=native`; the probe in R0 measures the sign.

## 2. Mechanism

One private inline helper on `Material`, declared in `cl_Material.hpp` next to `mNff` and defined
at the top of `powerlaws.hpp`:

```cpp
// root of the Bézier discriminant b² + ac, which is exactly 0 at the upper
// knot when n − 1 == mNff and lands a few ulps below it in floating point;
// the clamp returns 0 there, which gives t = 1 and rho = rho_n
static inline real
piecewise_bezier_sqrt( const real a, const real b, const real c )
{
    return std::sqrt( std::max( b * b + a * c, 0.0 ) ) ;
}
```

- **value (8 sites):** `real tParam = ( b + piecewise_bezier_sqrt( a, b, c ) ) / a ;`
  At disc = 0 this is `b/a`; at the degenerate knot `a == b` (comment `:2313`), so t = 1 and
  `rhoFF = rho3 = rho_n`, continuous with the `J > j3` branch that returns `rhon`.
- **dJ (8 sites):** keep the guard; below it `tSqrt = piecewise_bezier_sqrt( a, b, c )` replaces the two
  roots (`tParam` and `dt_dc`): two real `pow(…,0.5)` in seven overloads, two complex `std::sqrt( d )` in the
  `( J, T, B, β )` overload at `:1867-1871`, whose `cplx` intermediates become `real`. Behavior unchanged
  where it runs (disc > eps there, so the complex root was real).
- **dT (2 sites):** `tSqrt = piecewise_bezier_sqrt( a, b, c )`; the guard stays on `tDisc` and its
  `tFrozen` return is now finite, which is what its own comment already claims. The comment block
  `:2326-2332` is rewritten: the negated comparison still takes NaN but no NaN arises any more; the
  guard's job is the `1/tSqrt` in `dSdT`.

`std::sqrt` replaces `std::pow( x, 0.5 )`. `sqrt` is correctly rounded; `pow` need not be, so a
last-ulp drift on the positive branch is possible. The overload-pair test at 1e-14 relative
(`test_MaterialDerivatives.cpp:745`) and the finite-difference tests are the gate for that (R6).

## 3. Gap table

| # | class | gap | step |
|---|---|---|---|
| G1 | (a) defect | value overloads take the root of a possibly negative discriminant | R2 |
| G2 | (a) defect | dT overloads compute the root before the guard | R4 |
| G3 | (b) duplication | 26 copies of the root expression across 18 functions | R1–R4 |
| G4 | (c) evidence | the failing configuration is only on the CI runner | R6 |

## 4. Ordered steps

- [x] **R0** Probe (scratchpad `disc/disc.cpp`, no repo edit) — done 2026-09-24. With the test's `ec`, `jc`, n = 4 and
  rho_n from 1e-7 to 1e-4 Ω·m, `b² + ac` at J = j3 is exactly 0 for rho_n ≤ 3e-6 and **−3.6e-15 to −4.4e-15 for rho_n ≥ 1e-5** under `-O2`, `-Og`, and `-O2 -ffp-contract=off` alike; under `-march=native` it is negative from 3e-7 up (−2e-16). `a − b` is −1.8e-15 at the negative cases. `pow( disc, 0.5 )` is NaN in every negative case. The sign is a property of the log10 arithmetic at this rho_n, not of the optimizer alone; which side the runner's rho_n( 77 K ) lands on is the only configuration-dependent part.
- [x] **R1** Declare and define `piecewise_bezier_sqrt`.
- [x] **R2** Eight `rho_piecewise` sites → helper.
- [x] **R3** Eight `drho_piecewise_dJ` sites → helper for the post-guard root.
- [x] **R4** Two `drho_piecewise_dT` sites → helper; rewrite the guard comment.
- [x] **R5** Test: a second test `PiecewiseDegenerateBlendFiniteAcrossJc` (§7): n = 4, jc log-spaced
  1e9…1e15 × T = 20…90 K; at J = j3 the value (`( J, T, B, β )` and `( J )` overloads) and the dT tangent
  are finite, and `rho_piecewise( j3, T, B, β )` equals `rho( T )` to 1e-9 relative. dT is finiteness
  only: the two sides of the knot have different tangents (Grok's algebra: blend side 0 for constant jc, n;
  normal side `drhodT`). Red on pristine HEAD in the release Armadillo build (63 of 488 points).
- [ ] **R6** Gate (scratch part done by Claude: release + debug Armadillo `test_physics` 151/151, Blaze syntax; shared trees and pipeline open): `test_physics` in a **release Armadillo** build (the failing configuration) —
  scratch worktree, own build dir; then `make check-fast` in both shared trees (Christian); then a
  manual pipeline on devel before tonight's schedule (Christian).

## 5. Open questions

- [x] **O1** R5 asserts the endpoint value at 1e-9 relative and dT finiteness only (jury, both voices).
- [x] **O2** Private static member, declared in the class and defined in `powerlaws.hpp` as
  `inline real Material::piecewise_bezier_sqrt( … )` (no `static` at the definition).

## 6. Definition of done

- [x] All 18 functions use the helper: no `pow( …, 0.5 )` of `b*b + a*c` or `tDisc`, and no `std::sqrt( d )`
  of a `cplx` discriminant, outside the helper (`grep -nE 'pow\( *(b\*b|tDisc)|sqrt\( *d *\)' powerlaws.hpp` is empty).
- [x] `test_physics` green in a release Armadillo build (and debug); [ ] `make check-fast` green in both trees.
- [ ] Pipeline green on all four jobs.
- [ ] Devlog; this plan closed.

## 7. Addendum after dispatch (2026-09-24, 14:20) — the local reproducer

Not in the version the plan jury read.

- A release Armadillo build of pristine HEAD `ad8064f` (scratch worktree `fix/`, the runner's
  configuration) **passes** `PiecewiseDegenerateBlendStaysFinite`. The runner's failure does not
  reproduce at the test's own point on this machine.
- Reason, by binade arithmetic: with `logj2 == log10( j3 )` bitwise (n − 1 == mNff, same
  `pow` argument and exponent), `a = ( logj1 − 2·logj2 ) + log10( j3 )` is exact when
  `|logj1 − 2·logj2| < 16` (the partial sum stays in the binade of `logj1`, then Sterbenz), so
  `a == b` and `c == −b` bitwise and `b² + a·c` is exactly 0. At the test's `jc = 3.1e10` that
  holds. It breaks when the partial sum crosses 16 (`j3` above ~1e13), and independently
  whenever the compiler contracts `b*b + a*c` into a fused multiply-add (then the result is the
  rounding error of `b²`, of either sign). The runner's side of zero is therefore a property of its
  compiler or libm, not of the source revision; the CI session's "only RELEASE/Armadillo" is
  consistent with either.
- **Reproducer in the same build (probe test, scratch only):** n = 4, `ec = 1e-4`, jc log-spaced
  from 1e9 to 1e15 (61 values) × T = 20…90 K (8 values): **63 of 488 points return NaN** from
  `rho_piecewise( J, T, B, β )`, `rho_piecewise( J )` and `drho_piecewise_dT`, the first at
  jc = 2.5e11 A/m², a physical value for REBCO at 20 K. Every NaN point has all three NaN.
- **R5 becomes:** replace the single-point finiteness loop's role with this sweep as a second
  test, `PiecewiseDegenerateBlendFiniteAcrossJc`, expecting finiteness and, at the knot,
  `rho_piecewise( j3 ) == rho( T )` within 1e-9 relative. Red on pristine HEAD in this build (63
  points), green after the helper. The original test stays as it is.

## 8. Amendments from the plan jury (2026-09-24, `tmp/ai_exchange/review_piecewise_bezier_root.md`)

- [x] **A1** (Codex P0, Grok 2) Helper declared in the class, defined as a qualified `inline` member. O2.
- [x] **A2** (Codex P1, Grok 1) The eighth dJ overload's complex roots join R3; §6's check names all
  26 root evaluations.
- [x] **A3** (Codex P1, Grok 3) R5 uses the value endpoint only; `drhodT` is the tangent's real name. O1.
- [x] **A4** (Codex P2, Grok 5) The rounding claim is withdrawn; R6 gates the `pow`→`sqrt` drift.
- [x] **A5** (Grok 6) The helper takes `std::max( disc, 0.0 )` in that order (a NaN stays NaN) and
  carries `BELFEM_ASSERT( disc > -BELFEM_EPSILON * ( b*b + 1.0 ) )`: a few ulps below zero is the
  legal state, an order-one negative is a logic error and aborts in debug.
- [x] **A6** (Grok 6) Both dT comment copies are rewritten, not one.
- [x] **A7** (Grok 4, 6) Scope guard reworded: the clamp is the change; §7 records that the tree's
  own fixture point is exactly zero here and the runner's sign is unexplained.

## 9. Code round (2026-09-24, same exchange file)

Codex terra/high and Grok 4.7/high: no blocking defect. Applied as a separate step and re-gated:
- [x] **F1** helper: `d < 0.0 ? 0.0 : d` in place of `std::max` (specified for NaN); comment names the band, not "order one".
- [x] **F2** class comment: t = 1 and rho = rho_n "to roundoff".
- [x] **F3** sweep: one material per T, jc via `set_constant` (1.5 s → 24 ms); index `i`; dJ finiteness added.
- [x] **F4** `PiecewiseBezierRootClampsRoundoff`: the clamp's contract through a derived-type seam; assert throw under debug.
- Accepted residue, not changed: a discriminant below the band aborts in debug on the value and dT paths and is the silent flux-flow limit on the dJ paths (their guard returns first). Unreachable for legal knots; release agrees across all three.
