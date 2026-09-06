# Devlog 2026-04-07 — Thin-Shell H-Phi Formulation Note Review

**Date:** 2026-04-07
**Topic:** Read-only review of `todo/thinshell_hphi_formulation.md`
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Alves et al. 2022a, Alves et al. 2022b, Schnaubelt et al. 2023

## Summary

Reviewed Claude's thin-shell H-phi proposal against the current BELFEM Maxwell code and the local literature summaries.

Verdict: the motivation is directionally sound, but the implementation sketch understates the amount of structural work required and contains one formulation-level mismatch with the current BELFEM H-phi path. The strongest issue is that the proposed `phi_ts_insulator` kernel is not actually analogous to the existing air-block `phi` kernels: BELFEM's present Maxwell `phi` kernels assemble into `M`, not `K`.

## Key Findings

- The literature-side motivation is broadly credible:
  - [literature/papers/fem/schnaubelt2023.md](/home/christian/codes/belfem/literature/papers/fem/schnaubelt2023.md)
  - [literature/papers/fem/alves2022a.md](/home/christian/codes/belfem/literature/papers/fem/alves2022a.md)
  - [literature/papers/fem/alves2022b.md](/home/christian/codes/belfem/literature/papers/fem/alves2022b.md)
- The proposed internal H-phi reuse of the existing hanging-edge path is not established by the current code. The existing `hang_thinshell_edges_on_nodes_*` logic is built around external shell-to-volume facets in [src/fem/maxwell/cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1113) and [src/fem/maxwell/cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1220), not around internal layer-to-layer interfaces.
- The proposed `phi_ts_insulator` kernel in the note uses `aMatrices->K()`, but current Maxwell `phi` kernels assemble `trans(B)*B` into `aMatrices->M()` in [src/fem/maxwell/matrices/mt_maxwell_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi.cpp#L24) and [src/fem/maxwell/matrices/mt_maxwell_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi.cpp#L113). So the note's "same operator as existing air blocks" claim is not currently true.
- The current field-list/block-type machinery is more rigid than the note suggests. `DomainType::ThinShell` blocks are hard-mapped to the conductor DOF set in [src/fem/maxwell/cl_Maxwell_FieldList.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_Maxwell_FieldList.cpp#L255), and thin-shell sidesets are separately mapped through the `ThinShell` list in [src/fem/maxwell/cl_Maxwell_FieldList.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_Maxwell_FieldList.cpp#L387). This is not a "small" per-block tweak.
- The recommendation to cap bulk `rho` for `rint` but still let postprocessing "see" the real `rho` is internally inconsistent. Once the solve uses capped material data, postprocessing with uncapped `rho` no longer corresponds to the solved state.

## Open Questions

- Whether an internal H-phi interface inside a thin-shell stack should really be modeled through static condensation / hanging constraints, or instead through a new dedicated interface formulation.
- Whether the correct scalar-phi thin-shell operator should follow the existing BELFEM Maxwell `phi` mass-like assembly or a different reduced model.

## Files Updated

- devlog/dl20260407_thinshell_hphi_formulation_review.md
