# Devlog 2026-04-01 — Ghost Thin-Shell Review

**Date:** 2026-04-01
**Topic:** Read-only audit of `tmp/ghost_wip_analysis.md` for thin-shell ghost stabilization
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** `todo/thinshell_selective_nitsche_coupling.md`, Alves et al. 2022 (paper6) as routed in local design notes

## Summary

Reviewed the current ghost-stabilization work-in-progress note against the live Maxwell and mesh code. No source files were modified.

## Key Findings

- The current ghost coefficient logic zeroes the interface operator for any HTS-adjacent interface below `T_crit`, not only HTS/HTS pairs. See `src/fem/maxwell/matrices/mt_maxwell_h.cpp:1829-1840`.
- The proposed `aCalc->normal()`-based `Cm/Cs` fix has a hidden prerequisite: ghost sidesets carry `PENTA6TS` master/slave types, but `Calculator` has no `GeometryType::PENTA` normal callback. See `src/fem/kernel/cl_FEM_SideSet.cpp:41-44` and `src/fem/kernel/cl_FEM_Calculator.cpp:708-767`.
- Selective ghost creation cannot currently be implemented in `ThinShellFactory::create_ghost_facets()` based on `mMaterialBlockAssignment`, because materials are attached after shell creation. See `src/mesh/cl_ThinShellFactory.cpp:109-113`, `src/mesh/cl_ThinShellFactory.cpp:262`, and `src/fem/maxwell/cl_MaxwellFactory.cpp:883-888`.
- Filtering ghost interfaces later in `MaxwellFactory` would not restore CG behavior, because the mesh already duplicates layer edges before ghost registration. See `src/mesh/cl_ThinShellFactory.cpp:1459-1510` and `src/mesh/cl_ThinShellFactory.cpp:1595-1625`.
- The abandoned facet-`E` path was also dimensionally inconsistent: `EF_TRI3::E()` is `2x3`, while the draft `Dm/Ds` assembly assumed a `3x3` face operator. See `src/fem/interpolation/nedelec/cl_EF_TRI3.cpp:17-25` and `src/fem/maxwell/matrices/mt_maxwell_h.cpp:1878-1891`.

## Changes Made / Proposed

- Proposed only: fix interface resistivity scaling before continuing the ghost assembly.
- Proposed only: make selective DG decisions before edge duplication if CG must be preserved on same-material interfaces.
- Proposed only: add `PENTA6TS` normal support before using `aCalc->normal()` in `h_ghost`.
- Proposed only: verify `Cm/Cs` versus through-thickness `E`-difference equivalence on an affine prism before changing formulation.

## Open Questions

- Whether the intended ghost formulation should stay with the local 1D through-thickness derivative from `todo/thinshell_selective_nitsche_coupling.md` or move to the full `curl x n` form used in the WIP note.
- What positive lagged interface resistivity should be used for HTS layers in the first ghost prototype.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260401_ghost_thinshell_review.md
