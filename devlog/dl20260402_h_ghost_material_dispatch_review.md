# Devlog 2026-04-02 — h_ghost Material Dispatch Review

**Date:** 2026-04-02
**Topic:** Read-only review of how difficult exact resistivity recovery would be inside `h_ghost()`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Checked the `h_hts*` family and the generic `h()` / `h_ts()` kernels to determine whether `h_ghost()` can cheaply recover the same effective resistivity logic.

Conclusion: the user is right. Exact recovery is not trivial because the coefficient depends on the material model family, thermal coupling, defect support, shell-vs-bulk handling, and in the thin-shell case the reconstructed total field and tape-normal angle.

## Key Findings

- `IWG_Maxwell` selects many separate HTS kernels depending on thermal coupling, defects, and piecewise-vs-powerlaw mode ([src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L328), [src/fem/maxwell/cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L458)).
- The generic `h()` kernel already contains a nontrivial material/model dispatch tree for coefficient evaluation, including special handling for HTS, user-defined materials, defects, and thermal coupling ([src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2029)).
- The generic `h_ts()` kernel adds thin-shell-specific field reconstruction and angle handling through `bt + bn` and the tape normal, so exact shell-style resistivity recovery needs more than just `q` and `T` ([src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2175), [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2252)).
- Because of that complexity, a mesh field is not a good primary source of truth. A shared helper extracted from the generic kernels is a better long-term path if exact side recovery is required.

## Changes Made / Proposed

- No source changes made.
- Logged the review in `todo/ai_exchange.md`.

## Open Questions

- Should `h_ghost()` deliberately use only an interface proxy coefficient, or is the goal to match the side kernels as closely as possible?
- If exact matching is required, should the shared helper return only `rho`, or both `rho` and `drho_dJ` so the volume kernels can migrate to it too?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260402_h_ghost_material_dispatch_review.md
