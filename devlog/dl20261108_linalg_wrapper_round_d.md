# linalg wrapper hardening, round D: the consumer escape hatches come out

**Date:** 2026-11-08
**Purpose:** Round D of `todo/linalg_wrapper_hardening_plan.md` (v3): after rounds A–C gave the
`Matrix`/`Vector` wrappers complete expression overload sets on both backends, retire the
`matrix_data()`/`vector_data()` lines the consumers had written to get around the gaps. This entry
covers the inventory (R4.0) and the two sub-rounds done in this session, R4.1 (`cl_IWG_Timestep.cpp`)
and R4.2 (`src/fem/kernel`, `src/fem/postproc`). R4.3 was handed to the second session.
**Module:** `src/fem/iwg`, `src/fem/kernel`, `src/fem/postproc`, `src/linalg` (one comment), `tests/linalg`
**AIs involved:** Claude Fable 5.1 (inventory, edits, gates, reconciliation); code juries per sub-round,
Codex `gpt-6-astra`/high + Grok `grok-4.7`/high, parallel blind (`scripts/cross_review.sh --jury`)
**Status:** R4.1 and R4.2 VERIFIED — Christian's `make check` green on Blaze and on Armadillo (2026-11-08); committed
**Verification:** compile gate (`g++ -c` with each TU's `flags.make`, `-Werror`) on every touched TU, Blaze
and Armadillo, release and debug: clean. Blaze allocation-parity probe on every R4.1 expression shape:
wrapper form == hatch form. Christian's `make check` on Blaze with R4.1 + R4.2 in the tree: every suite
green (before the two post-jury edits D21/D22, a comment and a new test); then on Armadillo with D21/D22
in the tree: green. Afterwards, executed by Claude at Christian's request (he was away): the full suite
in the shared tree, Blaze **debug** (17/17, 09:58, the first debug-mode run since round A) and Armadillo
**debug** after a reconfigure (`USE_MATRIX_ARMADILLO=ON`, 17/17, 12:07, log `Testing/Temporary/LastTest.log`);
the D22 test and the `CheckPolicyDebug` cases ran and passed on both. Round D is gated on Blaze release,
Blaze debug and Armadillo debug.
**Exchange:** `tmp/ai_exchange/review_linalg_wrapper_hardening.md` (R4.1 and R4.2 pre-registrations,
verifications, reconciliations; diffs `…_roundD1.diff`, `…_roundD2.diff`)

## Summary

Of the 110 `matrix_data()`/`vector_data()` lines outside `src/linalg`, 41 are not linalg hatches
(`Cell::vector_data()` moves, the HDF5 string accessor, `Spline::matrix_data()`) and stay. The
Timestep IWG held 79; 70 of them are live and are now wrapper expressions (`aJ += M()`,
`aRHS += aJ * tSwap`, `aSwap = b0 * qold(0)`, …), the nine in Crank–Nicolson and Galerkin stay because
both schemes are refused at construction and the refusal is tested. Twelve more lines in the kernel
Calculator, the Gradient post-processor, the edge-length routine and the eigenvalue power iteration are
replaced. One of those is not mechanical: the eigen residual on Blaze.

## What the juries found

**R4.1 (Timestep).** No finding from either seat. Codex traced the history-buffer constraint —
`qold()` overwrites and returns the shared `mq0`, so every history term must be consumed in its own
statement — and confirmed the rewrite keeps it; Grok re-derived the sign pattern of `collect_qhist`,
the `=`-versus-`+=` split, and the BDF order (history product, then `*= alpha`, then `+= h K`). Grok's two
open items: Crank–Nicolson/Galerkin still on raw accessors (intended, dead code), and the
`explicit_euler` product `( aJ − h K ) * qold()` passing through the `operator*( const A&, const Vector& )`
catch-all that carries no debug size assert — pre-existing, carried as a round-E candidate.

**R4.2 (kernel, postproc, eigen).** No P0/P1. Both seats confirmed from the SCLS headers that
`norm( mY − λ mZ )` on a complex vector is the Euclidean norm on both backends (Blaze: `SqrAbs` then
`Sqrt`, result complex with zero imaginary part; Armadillo: real `pod_type`), that the Blaze expression
does not allocate (`CompositeType` is a reference), and that the `fn_norm.hpp` comment states the overload
rule correctly: an rvalue cannot bind `norm( ET & )`, so argument-dependent lookup reaches the backend's
`norm`; a `const ET &` overload would be an exact match against both backends' own `norm` and ambiguous
(it was tried in this round and was). Grok then made the one substantive point: **the old Blaze branch
computed the square root of the 1-norm**, `sqrt( Σ|rᵢ| )`, so with `mEpsilon = 1e-4` Blaze stopped at
`|r| < 1e-8` while Armadillo, already on the 2-norm, stopped at `|r| < 1e-4`. Both backends now use the
2-norm. This is the right definition and it is a change of the Blaze stop test; no test drives the power
iteration, so it is recorded here rather than claimed covered.

Two findings applied after the jury: **D21** (Grok) — my comment at the eigen site carried history
("used to take the square root of the 1-norm"), which `doc/commenting_guidelines.md` §8 forbids; trimmed to
the return-type contract. **D22** (Codex) — the `norm` tests exercised real vectors only; added
`Vector.ComplexExpressionNormIsTheTwoNorm`: complex `y`, `z`, `λ` with nonzero imaginary parts, modulus
against the hand-summed 2-norm, imaginary part zero, and a `static_assert` that the expression `norm` is a
scalar, not a materialized `Vector` — the latter pins Grok's open risk that a later overload change could
make the implicit `Vector` constructor the better match and allocate an `mN` temporary per iteration.

Predictions: R4.1 P1 (Armadillo parity measured, not inferred) and P2 (extra shape compare in
`collect_qhist`) were not raised; R4.2 P1 (behaviour change needs recording) was raised, with the exact
factor; P2 (`inv` materialization on Blaze) was not. P3 (no P0/P1) held both times.

## What is open

- R4.3 (`cl_Gas.cpp` ×5, `cl_GM_EoS_Cubic.cpp:115`, `fn_ODE_RK45.cpp:166`, `fn_ODE_DOP853.cpp:286`): handed
  to the second session (belfem-5e) with the gate and jury instructions; it is holding until Christian
  approves source edits in that session. It also corrected me: the tooling fixes (`cross_review.sh` mktemp
  template, `ask_copilot.sh` under bash 3.2) owe a round-1 jury under the 2026-09-22 ruling, not an exemption.
- R4.4 (`nonfree/cl_CN_Reaction.cpp:101`, a nonfree session) and R4.5 (`cl_Cohomology.cpp:1143`, report to
  Gregory).
- No numerical test for `fn_Mesh_compute_edge_lengths` or the eigenvalue power iteration (Appendix D gap).
- The unchecked `operator*( const A&, const Vector& )` catch-all (round E).

## Files Updated

`src/fem/iwg/cl_IWG_Timestep.cpp`, `src/fem/kernel/cl_FEM_Calculator.hpp`,
`src/fem/kernel/cl_FEM_DofMgr_EigenValues.cpp`, `src/fem/postproc/cl_Gradient.cpp`,
`src/fem/postproc/fn_Mesh_compute_edge_lengths.cpp`, `src/linalg/fn_norm.hpp` (comment),
`tests/linalg/test_Vector.cpp` (D22), `todo/linalg_wrapper_hardening_plan.md` (Status, R4.1/R4.2 ticks,
round-D paragraph, Appendix D outcome), this devlog, `devlog/README.md`.
