# Devlog 2026-04-07 — Thin-Shell `devel` Comparison for Single-Layer Patch Failure

**Date:** 2026-04-07
**Topic:** Read-only comparison of current thin-shell creation / hanging assembly against `devel` for the single-layer patch-test failure
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Compared the current branch against local `devel` with focus on the single-layer thin-shell patch-test failure. The thin-shell T-matrix and hanging-weight construction in `DofData` are unchanged relative to `devel`. The relevant branch deltas are upstream: `MaxwellFactory` hanging-edge processing order and `ThinShellFactory` front-end changes for duplicates/ghosts plus the recent `collect_nodes()` sort.

The user’s `h_ts_metal` diagnostics indicate the single-layer debug case is using node-based hanging shell edges tied to adjacent `phi` DOFs, and that the shell/master/slave integration points coincide geometrically. That makes the T-matrix path and point matching less likely leads.

## Key Findings

- No diff vs `devel` in:
  - `src/fem/maxwell/cl_Maxwell_TMatrix.cpp`
  - `src/fem/maxwell/cl_Maxwell_TMatrix.hpp`
  - `src/fem/kernel/cl_FEM_Tmatrix.cpp`
  - `src/fem/kernel/cl_FEM_Tmatrix.hpp`
  - `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp`
- `src/fem/maxwell/cl_MaxwellFactory.cpp` differs from `devel` in `create_hanging_edges_and_facets()`:
  - current branch processes air/node coupling before conductor/edge coupling,
  - `devel` processed conductor/edge coupling first.
- `src/mesh/cl_ThinShellFactory.cpp` differs substantially from `devel`, but most additions are duplicate-layer / ghost-facet infrastructure that should be inactive for a true single-layer case.
- In the user’s debug print, the shell local `edge_h` DOFs have `BELFEM_UINT_MAX` indices, which is consistent with deliberately hanging DOFs and not suspicious by itself.
- The printed adjacent master/slave volume elements show only `phi` DOFs, consistent with an air/shell/air single-layer setup using node-based hanging of shell edges.
- The `point k 0 0` diagnostics from `h_ts_metal` show the normal calculator maps master/slave/shell integration points to the same physical coordinates for the inspected facet.

## Changes Made / Proposed

- No source-code changes.
- Added this devlog and a corresponding summary to `todo/ai_exchange.md`.

## Open Questions

- If the failing single-layer case is air/shell/air, what front-end change outside the unchanged T-matrix path causes the node-based hanging constraints to fail the patch test?
- Is the recent `collect_nodes()` sort truly neutral in practice, or does it subtly affect shell-edge canonical orientation for this path?
- If `devel` is known to pass this exact single-layer case, can the regression be reproduced by bisecting the `MaxwellFactory` / `ThinShellFactory` deltas only?

## Files Updated

- `devlog/dl20260407_thinshell_devel_compare_singlelayer.md`
- `todo/ai_exchange.md`
