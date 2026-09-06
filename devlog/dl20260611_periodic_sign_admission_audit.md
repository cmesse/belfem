# Devlog 2026-06-11 - Periodic Sign Admission Audit

**Date:** 2026-06-11
**Topic:** Periodic thin-cut seam faces and sign-blind cut-face admission
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high
**Literature References:** N/A

## Summary

Audited the proposal to admit thin-cut faces by `abs(tCase)` in
`CutProcessor::collect_facets()` after the corc reproducer showed parity-random
positive and negative seam-adjacent face cases.

## Key Findings

- `determine_cut_case_3d()` selects the face from the unsigned case table first,
  then applies the cochain/orientation sign. The sign is orientation information,
  not an independent face selector.
- The positive-only admission path is at least the intended single-imposition
  mechanism for clean interior faces, but the code does not assert the invariant.
- `DomainType::Cut` sidesets do not appear to carry default HPhi cut dofs in the
  current Maxwell constructor, so a simple "two facets means two surface
  integrals" argument is too strong.
- `mThinCutFaces` is still algebraically important because it feeds node
  flagging, duplicate creation, and element relinking before cut sidesets are
  emitted.
- Target periodic faces become slave-only faces after `fix_face_slaves()`;
  admitting them into `add_thin_cut_sidesets_to_mesh()` is unsafe because
  `Facet::set_master(..., true)` dereferences `tFace->master()`.

## Changes Made / Proposed

- Appended the detailed counter-check to `todo/ai_exchange.md`.
- No source code was modified.

## Open Questions

- Runtime instrumentation should classify missing seam faces by `tCase`,
  `abs(tCase)`, whether the selected face is the seam face, and whether the face
  is source-side master-owned or target-side slave-only.
- A one-sided seam-duplicate policy would be a design change; it needs a proof
  before replacing the current symmetric periodic duplicate invariant.

## Files Updated

- `todo/ai_exchange.md`
- `devlog/dl20260611_periodic_sign_admission_audit.md`
