# Devlog 2026-08-23 — G-operator Step C: TET4+TRI3 and PENTA6TS+QUAD4TS complete

**Date:** 2026-08-23
**Topic:** First four edge-function G-operators implemented, tested, and closed
through the full C.1–C.4 pipeline with blind two-vendor audits at every stage.
**AIs involved:** Claude (implementation), Codex + Grok (blind audits — six
rounds this session: TET4/TRI3 C.3, TS C.1, TS C.2, TS C.3, each ×2 vendors
counted per stage)
**Claude Confidence:** high
**Codex Audit Confidence:** high (all rounds; re-ran the volume probes where
sandbox allowed)
**Grok Audit Confidence:** high ~85–92% (all rounds; independent re-derivation
in read-only sandboxes)
**Literature References:** none new — implements the already-audited C.1 math
(Whitney product rule; layout contract in `cl_EF_EdgeFunction.hpp`)
**Verification:** executed gates, all named and quoted below — `make check`
14/14 (100%) and `test_fem` 122/122 at session close; two deliberate-break
(L-01) gates went red as designed and green after restore.

## Summary

Christian signed off plan v3 (mGrad name confirmed — also avoids the mGram
adjacency) and later authorized the thin-shell pair. Four G-operators landed:

- **TET4 / TRI3** (affine Whitney): constant `mGrad` filled in `link()`
  beside `mC` from the post-division nablas; header-inline `G()`;
  cpp stubs deleted. Battery in `tests/fem/test_EdgeFunctions.cpp`
  (6 TESTs): FD vs compiled E (max dev 1.4e-14 TRI3 / 1.7e-13 TET4),
  curl tie ≤ 3.5e-17, Frobenius ‖G‖² = ½‖C‖², trace, gradient-mode
  blindness, shape ASSERTs, flipped-edge variants. Break gate: dropping
  `mS` → flipped variant red (curl-tie dev 0.457), unflipped green — the
  battery discriminates exactly as designed.
- **PENTA6TS / QUAD4TS** (thin shells): QUAD4TS G is constant
  (s·F′·∇η⊗∇ξ, link-filled, header-inline); PENTA6TS is the first
  point-dependent G (s·[F·W ∓ ½∇τ⊗w], new member `real mGradW[9]` — ONE
  tensor, since W is identical for all three pairs by ∇ζ = −∇ξ−∇η, an
  auditor-supplied identity). Battery in
  `tests/fem/test_InterfaceOrientation.cpp` (4 TESTs): FD through a
  NON-transposed test-side prism Jacobian (max dev 3.0e-12), curl tie at
  τ = 0.4 AND at a second cluster point (aIndex-dependence hardening),
  trace, corrected kernel test (φ = x+2y+3z circulations, Σc = 0, plus the
  all-s Ampère mode as a negative control), flip variants incl. a BOTTOM
  edge. Break gate: F₀↔F₁ swap → FD dev 0.299 / curl-tie 0.597, red as
  designed; restored green. NO Frobenius test — the identity is
  simplex-only (PENTA6TS G has a symmetric part; QUAD4TS is rank-1 with
  ‖G‖² = ‖C‖²).

## Key audit findings that shaped the code

- Both auditors independently REFUTED the first TS kernel-test sketch:
  equal signed dofs on the three triangle edges are the constant-curl
  (Ampère) mode, not a constant field — the correct curl-null vector needs
  c₀+c₁+c₂ = 0. The wrong vector became the negative control.
- Grok caught the FD "Jacobian rows" wording as the transpose defect class
  the TET4 round had already closed (class mJ is the TRANSPOSED geometry
  Jacobian) — plan rewritten before implementation.
- Convention caveat documented: the class nabla vectors are not always the
  geometric parameter gradients (QUAD4TS mNablaEta follows the STACKING
  direction). FD gates run on fixtures where class convention == geometry;
  the G-vs-C tie is the geometry-independent gate.
- Theory-note leftover fixed: §7 no longer claims QUAD4TS divergence is
  identically zero (it is stacking-dependent; kernel stays blind).
- TET4 fill loop renamed tEdge/tRow/tCol — link()'s scope carries live
  Jacobian reference aliases named e and i.

## Changes Made

- src/fem/interpolation/nedelec/: cl_EF_TET4.{hpp,cpp},
  cl_EF_TRI3.{hpp,cpp}, cl_EF_PENTA6TS.{hpp,cpp}, cl_EF_QUAD4TS.{hpp,cpp}
- tests/fem/test_EdgeFunctions.cpp (+6 TESTs), test_InterfaceOrientation.cpp
  (+4 TESTs, +2 helpers, +FD constants)
- todo/coulomb_gauge_penalty_theory.md §7 correction
- Exchange threads: coulomb_gauge_stepC1.md (TET4/TRI3, closed),
  coulomb_gauge_stepC_ts.md (TS pair, closed), + 8 per-vendor audit files

## Open / Remaining

- Elements without a real G (stubs hard-fail): HEX8 (its C.1 must derive
  the DJ term — grad does NOT inherit the curl's Piola cancellation, both
  auditors), TRI6, TET10 (straight + curved paths), HEX8TB, HEX8TS,
  LINE3 (contract decision with Christian, theory note §12.1).
- Calculator::G() passthrough deferred (outside approved scope).
- Final report (todo/coulomb_gauge_penalty_report.md) owed when Step C
  coverage is settled.

## Files Updated

- see Changes Made; devlog/README.md index line added.
