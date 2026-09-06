# Devlog 2026-08-22 — Coulomb Gauge Task, Steps B and C.1/C.2 (G-operator)

**Date:** 2026-08-22 (overnight session, continuation of dl20260822_coulomb_gauge_stepA.md)
**Topic:** EdgeFunction::G() interface landed (Step B, audited); G-operator math
and implementation plan for TET4+TRI3 audited through two blind rounds (C.1,
C.2); stopped at the plan sign-off gate.
**AIs involved:** Claude (implementation/derivation), Codex + Grok (blind audits,
four rounds this session: B-code, C.1-math, C.2-plan ×1 each vendor per round)
**Claude Confidence:** high
**Codex Audit Confidence:** high (all rounds)
**Grok Audit Confidence:** high ~85–90% (all rounds)
**Literature References:** layout/identities tie back to the Step A note
(Monk §5.5.1, §7.2.1; theory note §7); no new literature this session
**Verification:** Step B build gate VERIFIED — `make libbelfem_interpolation.a`
and `make libbelfem_kernel.a` green in cmake-build-debug (run twice, incl.
after the post-round comment fixes). Everything else this session is reviewed,
not verified; the class-level executable gate is C.4 (`make check`), not yet
reached.

## Summary

Christian signed off Step A and clarified the goal (G-operator + a penalty
factor that actually stabilizes; the theory note's answer stands: on simplex
meshes the stabilizing operator is the weak divergence, the pointwise G serves
the hex-mesh variant, diagnostics, and the test chain). Step B then landed the
interface: pure virtual `G( const uint aIndex = 0 )` in `EdgeFunction`,
protected member `mGrad` carrying the one-layout contract
(`G(i+d*j,e) = ∂(w_e)_j/∂x_i`, 9×n 3D / 4×n 2D, divergence = trace rows,
explicit curl-tie rows), hard-fail `BELFEM_ERROR` stubs in all ten subclasses,
ctor allocations in nine (LINE3 deferred). Deviation from the brief, flagged
and audit-endorsed: the member is `mGrad`, not `mG` — seven subclasses already
own a private `mG` that would shadow it.

Step C then ran its first two stages for the affine Whitney pair TET4+TRI3:
C.1 math (pair tables re-read from `E()`, gradient = s·(A⊗B − B⊗A) constant
per element, no Hessian term on affine maps, curl-tie and trace identities)
and the C.2 implementation plan (fill in `link()` beside `mC`, header-inline
`G()` like `C()`, zero per-Gauss-point cost, C.4 battery inside
`tests/fem/test_EdgeFunctions.cpp` reusing its FD stencil). Both C rounds
confirmed the substance and produced material corrections (below). Work is
STOPPED at the brief's mandated sign-off on the first element's plan.

## Key Findings

- Step B audits: zero functional defects; Grok caught the "and cyclic"
  curl-tie comment as a footgun (cycling row integers gives wrong pairs —
  fixed to three explicit rows) and the unnamed LINE3 sizing exception
  (fixed). Codex verified C++17 default-argument-through-base-pointer
  semantics for the HEX-family overrides.
- C.1 audits: pair tables and identities confirmed independently (Grok
  re-derived from `mG`/`mH` packing, expanded all six TET4 curl columns;
  INC-243's edge-2 site clean in the current tree). The FD test recipe in the
  math note was the weak point — E() aliases mE, precompute() clobbers the
  quadrature tables, mJ is transposed — retired in favor of the existing
  test_curl stencil.
- C.2 audits: both "needs revision", all paragraph-sized, folded into plan
  v3: TRI3 loop written out explicitly (pair-table silence IS the INC-243
  class), delete-the-stub + header-inline G() (link trap otherwise), a
  flipped-edge gradient test (the only mS falsifier — all default fixtures
  have s=+1), pinned non-constant φ for the gradient-mode test (L-01: v2's
  unspecified φ admitted a vacuously passing q=0), Calculator::G explicitly
  deferred as out-of-scope.
- Process: one dispatch failed on a relative script path after a cwd drift
  (redispatched, no side effects); one earlier redirect error gave Codex an
  empty prompt (failed instantly on usage, no pollution).

## Changes Made / Proposed

- `src/fem/interpolation/nedelec/cl_EF_EdgeFunction.hpp` — G() pure virtual,
  mGrad + layout contract (edit approved by Christian's go-ahead).
- All ten `cl_EF_*.hpp/.cpp` — G() declarations, BELFEM_ERROR stubs, ctor
  allocations (nine).
- Exchange: `coulomb_gauge_stepB.md` (+ per-vendor files),
  `coulomb_gauge_stepC1.md` (C.1 math, C.2 plan v1→v3, reconciliations,
  + per-vendor files ×2 rounds).
- No edits outside `src/fem/interpolation/nedelec/`.

## Open Questions

- Christian's sign-off on plan v3 (TET4+TRI3 implementation + C.4 battery),
  explicitly covering the edit to `tests/fem/test_EdgeFunctions.cpp`
  (outside the nedelec scope rule).
- LINE3 G contract (§12.1 of the theory note) — stub stands until decided.
- Later: Calculator::G() passthrough (outside scope, needs clearance);
  HEX8's C.1 must derive the DJ term (grad does NOT inherit the curl's
  Piola cancellation — Grok, C.1 round).

## Files Updated

- src/fem/interpolation/nedelec/ (21 files, see git diff)
- devlog/dl20260822_coulomb_gauge_stepBC.md (this file)
- devlog/README.md (index line)
