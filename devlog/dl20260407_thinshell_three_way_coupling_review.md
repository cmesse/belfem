# Devlog 2026-04-07 — Thin-Shell Three-Way Coupling Review

**Date:** 2026-04-07
**Topic:** Read-only review of the revised thin-shell design direction: condensation for continuity, `phi` for insulators, possible Nitsche/Robin path for contact impedance
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high
**Literature References:** Schnaubelt et al. 2023 (paper `schnaubelt2023`) §III

## Summary

Reviewed the revised design discussion in [todo/thinshell_hphi_formulation.md](/home/christian/codes/belfem/todo/thinshell_hphi_formulation.md) against the current BELFEM implementation. The external hanging-edge / T-matrix path does support the "condensation for continuity" direction. The `H-phi` internal case still needs new topology plumbing, but the algebraic core is already present. The contact-impedance idea is plausible and literature-backed at a high level, but it is not yet equivalent to BELFEM's current continuity-enforcing `h_ghost()` kernel.

## Key Findings

- The external `H-phi` condensation mechanism is real: shell edges are hung on volume nodes in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1364) and [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1420), and the edge-to-node T-matrix resolves to `[1,-1]` in [cl_FEM_DofMgr_DofData.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofMgr_DofData.cpp#L3712).
- The internal `H-phi` case is still not a literal reuse of the external helper path. The existing hang functions hardcode shell/volume facet roles via `aFacet->master()` / `aFacet->slave()` in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1379), [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1429), [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1497), and [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1584). New internal dispatch/helpers are still required.
- The scalar-`phi` operator correction remains valid: current Maxwell `phi` kernels assemble the elliptic operator into `M()`, not `K()`, in [mt_maxwell_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi.cpp#L25), [mt_maxwell_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi.cpp#L107), and [mt_maxwell_phi.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_phi.cpp#L113).
- The field-list scope warning also remains valid: `DomainType::ThinShell` is still hard-wired to conductor-like block DOFs in [cl_Maxwell_FieldList.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_Maxwell_FieldList.cpp#L255) and to the thin-shell sideset table in [cl_Maxwell_FieldList.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_Maxwell_FieldList.cpp#L387).
- Schnaubelt et al. 2023 §III does support collapsing the T2TCL to a surface and introducing two tangential `H` traces plus an additional surface contribution; see the local paper text around [schnaubelt2023.txt](/home/christian/codes/belfem/literature/papers/fem/schnaubelt2023.txt#L220). But this is a dedicated thin-shell interface formulation, not just a relabeling of BELFEM's current continuity penalty.
- BELFEM's current `h_ghost()` remains a continuity-oriented Nitsche kernel: it uses a penalty term weighted by `alpha` plus consistency terms weighted by `rho_harm` in [mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1838) and [mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1901). A contact-impedance / Robin law would need its own derivation and should not be described as "just reinterpreting alpha" without that derivation.

## Changes Made / Proposed

- No source edits.
- Wrote this devlog to capture the distinction between continuity coupling and contact-impedance coupling for the next design pass.

## Open Questions

- What is the exact weak form BELFEM should use for a collapsed contact-impedance interface? Schnaubelt 2023 motivates the surface-collapse idea, but the mapping onto BELFEM's `h_ghost()` matrices is not yet derived.
- If `H-H` internal interfaces default to condensation, should the existing `h_ghost()` path remain a global fallback only, or a per-layer / per-interface option?

## Files Updated

- /home/christian/codes/belfem/devlog/dl20260407_thinshell_three_way_coupling_review.md
