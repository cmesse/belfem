# Lobatto two-point weights, Gauss-quadruple three-point weights

**Date:** 2026-11-15
**Purpose:** An external bug report (reporting session intpoints-b5, `intpoints` repository,
`tmp/belfem_lobatto_n2/BUG_REPORT.md`) said that `intpoints_lobatto` returns the weights {1, 4/3}
for n = 2. Confirm it, fix it, and add tests so that the integration rules outside GAUSS cannot
drift unnoticed again.
**Module:** `src/numerics/integration` (+ `tests/fem`)
**AIs involved:** Claude Opus 5.5 (primary). Plan jury: Codex `gpt-5.6-terra`/high and Grok
`grok-4.7`/high. Code jury: Codex `gpt-6.1-terra`/high and Grok `grok-4.7`/high. Changelog language
sweep: Codex `gpt-6-luna`/medium.
**Status:** COMPLETE — plan `todo/closed/lobatto_two_point_weights.md`
**Literature References:** Bronshtein §19.3.2.2 Eq. 19.75 (trapezoid, h = 2), §19.3.3.2
Eq. 19.84a/b (3- and 4-point Lobatto)
**Verification:** `make check-fast` 17/17 (Christian, pasted output; tree not named). Reporting
session, SCLS gfortran 16.1 on read-only copies: HEAD `b81b309` fails exactly orders 0 and 1 of
the moment check (three monomials) and gives Gauss-quadruple n = 3 weights {0, 0, 0}; the fix
passes orders 0…20 and n = 1…6. Claude scratch probe, same result. **Not run:** the new gtests
against the unfixed line 233 (in-tree falsifier).
**Exchange:** `tmp/ai_exchange/review_lobatto_two_point.md`

## Summary

`intpoints.f90:207-208` set both endpoint weights to `2/(n(n-1))`, which is 1 for n = 2. The
`else` of `if(n .gt. 3)`, written for Simpson's middle weight, also ran for n = 2 and overwrote
the right endpoint with 4/3, so the weights summed to 7/3 instead of 2, the length of [−1, 1]. Now
`else if(n .eq. 3)`. The defect was latent. Grok showed that no deck can select LOBATTO:
`string_to_integration_scheme` and `set_integration_scheme` have no caller, and no input key
exists. This corrects the bug report's claim that an input deck reaches it. In the FEM layer only
order 1 would hit n = 2, because order 0 is rewritten first.

## What landed

- **Lobatto n = 2** (R1), as above.
- **Rule size** (R2, Christian's ruling (a)). `ceil( 0.5*( aOrder + 3 ) )` added in `unsigned
  int`: it wrapped to n ∈ {0, 1} for aOrder ≥ UINT_MAX − 2, and its `int` conversion overflowed near
  UINT_MAX/2 (Codex). The size is now `aOrder/2 + 2`, which is equal for every order, ≥ 2 and
  cannot wrap, and a release-active `BELFEM_ERROR` checks that it fits in `int`. For n = 1 the old
  Fortran would have divided by zero and written `w(2)` out of bounds (Grok, reporter). No Fortran
  guard was added: Christian rejected `error stop`, which bypasses the syslog path.
- **Gauss-quadruple n = 3** (D5, found by gfortran `-Wall` during the syntax gate).
  `5_QP/9_QP` was integer division of kind-8 integer literals, so all three weights were 0. Now
  `5._QP/9._QP`. Its only caller, `debye.f90`, requests n = 10 and n = 7, so the bug was latent.
- **`n*(n-1)` in `c_int`** (R8, code jury, both seats). The product overflowed for n ≥ 46342,
  inside the R2 bound. It is now converted to real before the multiply (`intpoints.f90:209, 229`).
- **Comments** (R3): the `legendre` seeds are labeled P1, P2, P2′. The dead `case(1)` / `case(2)`
  values are unchanged and unreachable, since every caller passes n ≥ 3.
- **Tests** (R4, R5): the scheme is a trailing default argument of `check_all_orders`, so the
  existing seven sweeps are unchanged. New: `LINE_LOBATTO` (orders 0…20), `QUAD_GAUSSCLASSIC`
  and `HEX_GAUSSCLASSIC` (orders 0…20; these reach the tensor-product path that GAUSS bypasses
  with dedicated tables at QUAD 4-13 and HEX 3-7), and `LOBATTO_reference_rules`, which compares
  nodes and weights through one shared permutation against Bronshtein for n = 2, 3, 4.
- `CHANGELOG.md` `[Unreleased]` → Fixed; the version is 0.9.5 in `CMakeLists.txt` (Christian),
  `CITATION.cff` stays at the published 0.9.4.

## Juries

The plan jury raised 10 findings, all CONFIRMED and folded into plan v2. Codex showed that the
v1 claim "`uint` implies n ≥ 2" was false (the wrap above) and that the O2 rationale named the
wrong routine. Grok showed the n = 1 out-of-bounds store, that decks cannot reach LOBATTO, that
the GAUSSCLASSIC justification was wrong (the kernel default is GAUSS,
`cl_FEM_KernelParameters.hpp:53`; the real caller is `cl_DatabaseProjector.cpp:183-188`), and
that the `pi` item was moot because `quadruple` is binary64 (already a DR-174 residue). Both
caught that R5 had to sort node–weight pairs together. In the code jury neither seat refuted a
functional claim. Both raised the `n*(n-1)` overflow (→ R8) and the missing test for
Gauss-quadruple.

## Open

- No gtest reaches `intpoints_gauss_quadruple` (not exported to C++). This is added to the DR-174
  residue in `todo/debt_register_closed.md` (Christian: widen, no new row).
- In-tree falsifier not run (see Verification).
- Not done: the reporter's side note that `cl_TensorMeshFactory.cpp:1141-1144` cites
  "Zienkiewicz et al. 2013, Eq. 17.4" has not been checked against the book.
