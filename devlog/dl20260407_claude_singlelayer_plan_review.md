# Devlog 2026-04-07 — Review of Claude's Single-Layer Regression Plan

**Date:** 2026-04-07
**Topic:** Read-only audit of Claude's consolidated plan for isolating the single-layer thin-shell regression
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed Claude's suspect list and isolation plan against the current code. Most elimination steps are sound: the thin-shell T-matrix path, `BearingData`, and hanging-DOF resolution are unchanged versus `devel`, and the hanging `BELFEM_UINT_MAX` indices are deliberate. The two strongest next diagnostics remain printing one resolved shell-edge hanging relation and dumping the first assembled `K`/`rhs` for direct `devel` vs `ghost` comparison.

Two hypotheses should be downgraded or removed:
- the `element_rho` field-index-shift theory is unsupported by `create_field_map()`
- the DEBUG node/edge index wipe is not a new ghost-branch change

The bearing/gauge check remains a cheap sanity check, but not a strong regression lead because the bearing-imposition code is byte-identical to `devel`.

## Key Findings

- `DofData::create_field_map()` uses only `aIwg->dof_fields()` and skips `EntityType::ELEMENT`, so `mFields.NonDof.push("element_rho")` does not renumber existing DOF field mappings.
- The bearing-imposition block in `src/fem/maxwell/cl_MaxwellFactory.cpp` is unchanged relative to `devel`, as are `cl_FEM_DofMgr_BearingData.{cpp,hpp}`.
- The DEBUG wipe of node/edge indices in `create_hanging_edges_and_facets()` already existed in `devel`; it is not a new branch-specific suspect.
- The new `PENTA6TS::is_thinshell()` / `slave_integration_penta()` path is real, but it appears more relevant to ghost-facet groups than to the current single-layer air/shell/air path.
- The connectivity rewrite is real, but with `connect_edges_to_ghost_facets()` still commented out, the new ghost-facet traversal in `connect_edges_to_edges()` is likely inert for the present case.

## Changes Made / Proposed

- No source-code changes.
- Added this devlog and a matching audit summary to `todo/ai_exchange.md`.

## Open Questions

- Does the first resolved shell `edge_h` hanging relation on `ghost` match `devel` exactly?
- If not, is the mismatch already present before `DofData`, in the shell-edge source assignment from `MaxwellFactory`?
- If the hanging relation matches, does the first assembled `K`/`rhs` still differ between branches?

## Files Updated

- `devlog/dl20260407_claude_singlelayer_plan_review.md`
- `todo/ai_exchange.md`
