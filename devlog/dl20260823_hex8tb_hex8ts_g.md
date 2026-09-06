# Devlog 2026-08-23 — HEX8TB + HEX8TS G-operators complete (campaign element work done)

**Date:** 2026-08-23
**Topic:** The final element pair's G-operators through the full C.1–C.4
pipeline; the HEX8TS term1 gradient contract; HEX8TB's first compiled
curl-sign gate; 9 of 10 elements done — the campaign's element work is
finished (LINE3 remains a contract decision).
**AIs involved:** Claude (derivation, implementation), Codex + Grok (blind
audits — two rounds: C.1 math with folded C.2 plan, C.3 implementation).
The C.3 Grok audit was dispatched MANUALLY by Christian from a
self-contained embedded-source prompt; the report was pasted back and
recorded verbatim in the exchange.
**Claude Confidence:** high
**Codex Audit Confidence:** high (both rounds; numeric probes for the
fixture misclassification and break (c))
**Grok Audit Confidence:** high (both rounds; full hand re-derivation of
the wall frame, both dF tables, and the Stokes numbers)
**Literature References:** none new — both elements transcribe their own
class conventions; the thin-shell approximation context is the class's
existing design, not a new formulation.
**Verification:** executed gates — `test_fem` 145/145, `make check` 14/14
(100%) after every edit; four deliberate breaks red with pre-registered
signatures, restores grep-verified.

## Summary

HEX8TB (4-dof side-connector wall) is the campaign's simplest 3D case
once its frame is understood: link() imposes an exact orthonormal cuboid
frame (t, b, n) from block data — deliberately not from the node metric —
so the map is affine, G = s(∇F)⊗(∇ξ) with dF rows 1–2 only, and
div e ≡ 0 identically by orthonormality. Both auditors re-derived the
handedness (b = n×t, t×b = n) and confirmed the implemented C is +curl —
no DR-98 sibling. The battery gave the element its **first compiled
curl-sign gate**: a per-dof Stokes check on the bottom face (E-loop vs
face integrals of both C·n and antisym(G)·n). Per-dof matters: the
obvious coefficient vector (1,1,0,0) telescopes to a vacuous 0 = 0 — a
Grok catch that killed my planned version of the gate.

HEX8TS (8-dof thin-shell machinery, factory- and Calculator-live) carried
the round's one genuine formulation decision. Its per-point nablas are a
thin-shell **convention** (pseudo-inverse of the 2×3 mid-surface Jacobian
plus a normal column), not the gradients of an exact inverse map. The
audited contract: **G = term1** — the reference derivatives chained
through those nablas held fixed, exactly the ingredients C() uses — so
antisym(G) reproduces the implemented C on every geometry. Scope, stated
in the header docstring in three cases: exact gradient on affine
(parallelogram) mid-surfaces; on planar non-parallelograms C is the exact
curl but G omits the symmetric inverse-map Hessian; on warped quads E, C
and G all share the thin-shell approximation. Nobody demanded the
pseudo-inverse-derivative machinery.

The round's headline catch (three independent finders: both auditors and
my own in-round prep): the claim document called `TS_TestPrism`'s
mid-surface "rectangular" — it is a flat **general** quad, on which term1
is not the full gradient (Codex quantified: |FD − term1| = 0.116, trace
= 0.042). The FD and trace gates therefore run on a **rotated rectangle**
(θ = 0.35, sides 1.3 × 0.8) via a new optional corner override on the
fixture; Grok sharpened the requirement — a parallelogram buys
FD-exactness, but trace = 0 additionally needs orthogonality. The default
general quad stays in the battery for the curl tie, which the contract
makes geometry-independent.

## Gate record (all executed)

- 6 new tests, first run green (the top-flip variant added post-C.3);
  max |G − FD(E)| = 7.1e-12 (HEX8TS rectangle), 4.1e-12 (HEX8TB) — exact
  identities at roundoff on the affine fixtures.
- test_fem 145/145; make check 14/14, 100%.
- Breaks: (a) global mS drop → only the bottom-flip test red; (a′)
  top-columns-only mS drop → only the NEW top-flip test red (proving the
  post-audit test closes exactly the hole both vendors named); (b) HEX8TB
  nabla-column swap → FD + y/z ties + Stokes-vs-antisym(G) red while
  trace stayed green (the trace weakness demonstrated) and Stokes-vs-C
  stayed green; (c) HEX8TS fill transpose → all three TS gradient tests
  red via the curl tie (a rank-1 gradF⊗nabla is never antisymmetric —
  Grok's parity argument, Codex's numeric check).

## Changes Made

- src/fem/interpolation/nedelec/cl_EF_HEX8TB.{hpp,cpp} — G() body
  replacing the stub; scope docstring.
- src/fem/interpolation/nedelec/cl_EF_HEX8TS.{hpp,cpp} — G() body
  (update_nabla + full mNx = mInvJ·mFxi product; axis table); three-case
  scope docstring.
- tests/fem/support/cl_TS_TestStack.hpp — TS_TestPrism appended optional
  aCorners mid-surface override (z = 0 plane); header prose notes the
  default quad is deliberately not a parallelogram.
- tests/fem/test_InterfaceOrientation.cpp — hex8ts_gradient_battery,
  hex8tb_gradient_battery + 6 TESTs (Hex8TsGradient, …GeneralQuad,
  …FlippedBottom, …FlippedTop, Hex8TbGradient, Hex8TbGradientFlipped).
- Exchange: coulomb_gauge_stepC_tb.md (claim doc + two reconciliations),
  per-vendor files (Grok's C.3 entry recorded from the manual dispatch).

## Open / Remaining

- LINE3: the only element without G — contract decision with Christian;
  standing recommendation is a permanent hard-fail stub (theory note
  §12.1: the ambient gradient is not defined from a 1D manifold element
  alone).
- Recorded low-severity residuals (both vendors concur): the
  skew-parallelogram FD case is claimed by the docstring but gated only
  on the rectangle subclass; warped quads unexercised (shared thin-shell
  convention); the HEX8TB fixture is axis-aligned; the HEX8TS kernel node
  pairs are hardcoded (match the element table today).
- Calculator::G() passthrough still deferred (out of approved scope).
- Final campaign report: todo/coulomb_gauge_penalty_report.md — next.

## Files Updated

- see Changes Made; devlog/README.md index line added.
