# Lobatto Two-Point Rule Returns Wrong Weights

**Date:** 2026-11-15
**Purpose:** `intpoints_lobatto` returns the weights {1, 4/3} for n = 2 instead of {1, 1}, so every
`IntegrationScheme::LOBATTO` line rule requested at order 0 or 1 gives wrong integrals and raises
no error. Fix the Fortran branch, make the rule size safe to compute for every `uint` order, and
lock the LOBATTO and GAUSSCLASSIC dispatch paths under the existing monomial-exactness test.
**Module:** `src/numerics/integration` (+ `tests/fem`)
**AIs involved:** Claude Opus 5.5 (plan), Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high
(plan jury, round 1)
**Status:** ✅ COMPLETE (2026-11-15). R1-R5 and R7-R9 landed. Gates: `make check-fast` 17/17
(Christian, pasted output: 8 fast and 4 mpi tests; the tree was not named), and the reporting
session's Fortran driver and moment check (HEAD fails exactly orders 0 and 1 and Gauss-quadruple
n = 3; the fix passes orders 0…20 and n = 1…6). Not run: the in-tree falsifier, i.e. the new gtests
against the unfixed line 233. The Lobatto tests' power to detect the bug rests on the reporter's
Fortran falsifier and on source trace. Residue: no test reaches `intpoints_gauss_quadruple`
(DR-174 residue, closed register).

> **Scope guards:**
> - No change to the Gauss tables or the Gauss-Legendre Fortran routines.
> - The `quadruple` kind (`intpoints.f90:16`, `selected_real_kind( r = 60 )`, which resolves to
>   kind 8 = binary64) is a known residue of the closed DR-174
>   (`todo/closed/materials_tangent_cluster.md:69-74`) and is not touched here. Widening it would
>   change both Gauss routines and the Debye model (`debye.f90:461`).
> - The uncapped Newton loops (`intpoints.f90:121, 173, 219`) are pre-existing and stay as they are.

## Source

An external bug report from the `intpoints` repository, `tmp/belfem_lobatto_n2/BUG_REPORT.md`
(revised version; reproducer `driver.f90`, exactness probe `moments.f90`, `fix.diff`, kind probe
`kind.f90`). The reporting session compiled the unmodified `intpoints.f90` with SCLS gfortran 16.1
and found that the original fails exactly orders 0 and 1, and that `fix.diff` passes orders 0…20.
That run happened outside this checkout. Here the claim rests on source trace.

## Defects

- **D1 (HIGH): n = 2 falls through the n = 3 branch.** `intpoints.f90:207-208` set
  `w(1) = w(n) = 2/(n(n-1))`, which is the correct value 1 for n = 2. The `else` at
  `intpoints.f90:231-232` then sets `w(2) = 4/3`, Simpson's middle weight. For n = 2 that slot is
  the right endpoint, so the weights sum to 7/3 instead of 2, the length of [−1, 1]. Confidence:
  high (source trace, three independent readers).
  - Reference values: Bronshtein §19.3.3.2, Eq. 19.84a (3 points: x = −1, 0, 1; c = 1/3, 4/3,
    1/3, the Simpson rule) and Eq. 19.84b (4 points: x = ±1, ±1/√5; c = 1/6, 5/6, 5/6, 1/6). The
    2-point rule is the one-panel trapezoid formula, Bronshtein §19.3.2.2, Eq. 19.75 with
    h = 2: weights {1, 1}. Bronshtein's index is one less than the point count, so N points are
    exact to degree 2N−3, which matches `fn_intpoints.cpp:546`.
  - Reach: `intpoints( LOBATTO, LINE, p )` for p ∈ {0, 1}. In the FEM layer only p = 1 can occur,
    because order 0 is rewritten first (`cl_IF_IntegrationData.cpp:78-79`,
    `cl_FEM_SideSet.cpp:506-508`). **No deck can select LOBATTO today:**
    `string_to_integration_scheme` and `KernelParameters::set_integration_scheme` have no caller,
    and `doc/input_schema.yaml` has no such key (Grok). The defect is latent and reachable only
    from code. The bug report's "an input deck can select it" does not hold.
- **D2 (MEDIUM, latent): the rule size can be 0 or 1, and n = 1 corrupts memory.**
  `tN = ceil( 0.5 * ( aOrder + 3 ) )` (`fn_intpoints.cpp:546`) adds in `unsigned int`, so it wraps
  for aOrder ≥ UINT_MAX − 2 and gives tN ∈ {0, 1}. Near UINT_MAX / 2 the `double` → `int`
  conversion overflows (Codex). For n = 1 the Fortran routine divides by zero (line 207) and then
  writes `w(2)` past the end of a length-1 array (line 232) (Grok, and the reporter). R1 removes
  the out-of-bounds store for every n ≠ 3; the division remains. No present caller passes such an
  order.
- **D3 (LOW, comments only): wrong labels in `legendre`.** Lines 51-52 label P1 and P2 as P0 and
  P1. Lines 54-55 repeat the off-by-one (`! dP1dx` on P2′ = 3x, `! d2P1dx2` on `h0 = 2`, which
  is neither P1″ = 0 nor P2″ = 3). `case(1)` and `case(2)` (lines 42-48) return 0 and 1. All of
  these initial values are dead for every call this file makes: Gauss enters with n ≥ 4 and
  Lobatto with n−1 ≥ 3, and lines 67-72 overwrite both derivatives.
- **D4 (test gap):** `test_IntegrationExactness.cpp:157` hard-codes `IntegrationScheme::GAUSS`, so
  neither LOBATTO nor GAUSSCLASSIC is checked for exactness. GAUSSCLASSIC is **not** the kernel
  default (`cl_FEM_KernelParameters.hpp:53` is GAUSS, and `cl_FEM_Block.cpp:129-131` takes the
  parent's scheme). Its real caller is `cl_DatabaseProjector.cpp:183-188`. On QUAD orders 4-13 and
  HEX orders 3-7, GAUSS uses dedicated tables, so the tensor-product path `gauss_quad` /
  `gauss_hex` is not checked at those orders today.

- **D5 (MEDIUM, latent; found during implementation):** `intpoints_gauss_quadruple` n = 3
  (`intpoints.f90:164-166`) wrote `5_QP/9_QP` and `8_QP/9_QP`. These are integer literals of kind
  QP, so the division truncates and all three weights were 0. gfortran `-Wall` flags this as
  `-Winteger-division`. The only caller, `debye.f90:461`, requests n = 10
  (`cl_Material_Metal.cpp:318`) and n = 7 (`debye.f90:212`), so the branch is not reached today.
  The scratch probe on HEAD showed `gaussQ n=3 sum=0`.

## Plan

- [x] **R1 — Fortran fix (D1).** `intpoints.f90:231`: `else` → `else if( n .eq. 3 ) then`. The
  single-line `if` at line 230 does not take an `else`, so the new branch binds to line 210.
  Among n ≥ 2 only n = 2 changes.
- [x] **R2 — safe rule size (D2), after R1.** Form per O3. Recommended: compute the size in integer
  arithmetic, `const uint tN = aOrder / 2 + 2 ;`. This equals `ceil( (p+3)/2 )` for every p (p = 2m
  gives m+2; p = 2m+1 gives m+2), cannot wrap, and is ≥ 2 by construction. Then
  `BELFEM_ERROR( tN <= INT_MAX, … )` before the `int` conversion for the Fortran call. The check
  stays active in release: `intpoints()` runs on the setup path (coding philosophy, error tiers),
  and this rejects a caller's bad input; it is not a debug invariant.
- [x] **R3 — comment fix (D3).** Labels at lines 51-52 and 54-55 corrected to P1, P2, P2′, and
  "initial value, overwritten below for n ≥ 3". The dead `case(1)` / `case(2)` values stay as code
  (outside the scope of this fix); one comment states that `legendre` is called only for n ≥ 3.
- [x] **R4 — exactness test parameterized by scheme (D4).** In
  `tests/fem/test_IntegrationExactness.cpp`, `check_exactness` and `check_all_orders` take an
  `IntegrationScheme` as a **trailing** argument defaulting to `GAUSS`, so the seven existing calls
  stay as they are. Add `LINE_LOBATTO` (orders 0…20, n = 2…12), `QUAD_GAUSSCLASSIC` and
  `HEX_GAUSSCLASSIC` (orders 0…20; HEX 20 is 11³ = 1331 points, the size the tolerance comment
  already names).
- [x] **R5 — literature-anchored reference rules (D1).** `TEST( IntegrationExactness,
  LOBATTO_reference_rules )`: orders 1, 3 and 5 through the public dispatch must give n = 2, 3 and
  4 with the nodes and weights of Bronshtein Eq. 19.75 (h = 2), 19.84a and 19.84b, to a tolerance
  of 1e-14. **One permutation** that sorts the nodes ascending is applied to nodes and weights
  together, and each pair is then compared by index, so a rule only counts as correct if every
  node carries its own weight. This is what the moment test cannot show: a same-size rule that is
  exact to the required degree but has other nodes passes R4 and fails here.
- [x] **R7 — D5 fix (plan amendment, approved 2026-11-15).** `5._QP/9._QP`, `8._QP/9._QP`. The
  routine is not exported to C++, so no gtest covers it; the probe and the reporter's driver do.
- [x] **R8 — `n*(n-1)` overflow (code jury, Codex + Grok; approved 2026-11-15).** The product was
  formed in `c_int` and overflowed for n ≥ 46342 (order ≥ 92680), inside the R2 guard. Both sites
  (`intpoints.f90:209, 229`) now convert to real before the multiply. Pre-existing and far outside
  any FEM order.
- [x] **R9 — changelog and register.** `CHANGELOG.md` `[Unreleased]` → Fixed (three bullets;
  Codex `gpt-6-luna`/medium language sweep applied). The missing test for
  `intpoints_gauss_quadruple` goes onto the DR-174 residue in `todo/debt_register_closed.md`
  (Christian: widen, not a new row).
- [x] **R6 — gate.** Christian builds and runs `make check-fast` (`tests/fem` is labeled `fast`,
  `tests/fem/CMakeLists.txt:2`). Expected: before R1, `LINE_LOBATTO` fails at orders 0 and 1, and
  `LOBATTO_reference_rules` fails at n = 2. After R1, all pass. Second gate: the reporting session
  re-runs its `driver.f90` / `moments.f90` against the edited `intpoints.f90`, including the
  falsifier (the original must still fail exactly at orders 0 and 1).
  Probe already run (Claude, scratchpad, gfortran 16.1, `-O0 -fcheck=bounds`): fixed file gives
  Lobatto n = 2…6 weight sums of 2 and the reference values, and Gauss-quadruple n = 1…5 sums
  of 2. HEAD gives Lobatto n = 2 sum 7/3 and Gauss-quadruple n = 3 sum 0. This is the probe rung
  only; the gtest suite has not run.
  **Second gate ran 2026-11-15 (reporting session intpoints-b5, SCLS gfortran 16.1, read-only
  copies of HEAD `b81b309` and the working tree):** HEAD gives Lobatto n = 2 {1, 4/3}, fails exactly
  three moment checks (order 0 x⁰, order 1 x⁰, x¹) under both sizings, and Gauss-quadruple n = 3
  {0, 0, 0}. The fix gives Lobatto n = 2…6 sums of 2 with n = 3…6 identical to HEAD, 0 moment
  failures for orders 0…20 under both sizings, and Gauss-quadruple n = 1…6 all exact to degree
  2n−1. `aOrder/2+2 == ceil((aOrder+3)/2)` checked for orders 0…1000. Not run there: the C++
  dispatch and the n ≥ 46342 overflow (source trace only).

## What could break

- R1 changes only n = 2 among valid sizes. Confidence: high.
- R2 changes the size computation. It gives identical results for every order a caller can
  meaningfully pass; this is checked by the R4 sweep, which goes through the dispatch.
- R4: n = 4…12 relies on the uncapped Newton loop meeting the 1e-12 tolerance. The reporter's
  Fortran probe says it does; if it does not, that is a pre-existing finding to report, not a
  regression. A red GAUSSCLASSIC test would likewise be a new finding, reported rather than patched
  inside this plan.

## Open questions

- **O1:** ~~the maximum order for `HEX_GAUSSCLASSIC`~~ **RESOLVED 2026-11-15 → 20.** `gauss_hex`
  (`fn_intpoints.cpp:640-650`) is generic: a tensor product of `tNi = ceil( p/2 + 1 )` points,
  exact to degree ≥ p+1.
- **O2:** ~~`pi` precision~~ **WITHDRAWN 2026-11-15**: the kind is binary64, so nothing narrows
  (Grok, reporter). The misleading name of the kind remains a DR-174 residue (scope guard).
- **O3 (Christian):** the R2 guard. (a) Recommended: integer sizing plus `BELFEM_ERROR` in C++ only.
  (b) In addition, a Fortran-side `n < 2` guard in the exported routine. Its options are
  `error stop`, which bypasses the `BELFEM_ERROR` and syslog path, or an early return, which is
  silent. (c) `BELFEM_ASSERT` only, as in v1. Grok preferred this, but on the "uint implies
  tN ≥ 2" premise that Codex refuted.

## Audit trail

- `tmp/ai_exchange/review_lobatto_two_point.md` — plan jury round 1 (verification and
  reconciliation table at its end); the code jury will append to the same file.
