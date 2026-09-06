# Devlog 2026-04-29 — Connector Weight Physics Audit

**Date:** 2026-04-29
**Topic:** Read-only physical audit of `ThinShellFactory::preprocess_binomial_edges()` side-connector hanging weights
**AIs involved:** Codex
**Codex Audit Confidence:** medium-high
**Literature References:** Monk 2003 §5.5.1; Messe et al. 2023 (paper1) §3; Alves et al. 2022b (paper6) §3.1

## Summary

Audited the current horizontal side-connector hanging weights as a physical/static-condensation constraint. The static-condensation direction is consistent with BELFEM's native hanging-DOF machinery and with the literature preference for eliminating interface DOFs over Lagrange multipliers. The current weight formula is only partially physically justified.

## Key Findings

- The raw factor `mConnectorWidth / dot(b, r)` in `src/mesh/cl_ThinShellFactory.cpp:2961` is defensible for the off-curve Whitney-edge contribution at a boundary vertex. For a curve-edge triangle, the off-curve Nedelec basis contribution to `H · b` scales like `1 / dot(b, r)`, and the connector horizontal DOF adds the wrap thickness factor.
- The final normalization `tWeights /= sum(tWeights)` in `src/mesh/cl_ThinShellFactory.cpp:2966` is physically suspect. It converts geometric point-evaluation coefficients into a partition-of-unity average and can erase the thickness/projection magnitude that the raw edge-circulation relation was trying to encode.
- The target horizontal edge is created from inner node to outer node in `create_binomial_edges()`, while the outer node is displaced by `aSign * b * mConnectorWidth` in `create_outer_nodes()`. Since `preprocess_binomial_edges()` receives only `b`, not `aSign * b`, one side of the tape appears to miss an orientation sign unless `b` is already outward by construction for that side. Current evidence points to an omitted `aSign` factor.
- Excluding curve-edge sources is exact only for locally orthogonal boundary triangles or if the actual basis cleanly separates curve-tangent and binormal components. The current `EF_PENTA6TS` implementation uses the standard Whitney form, so on skew triangles the curve-edge basis can contribute to `H · b`; omitting it leaves tangential-H contamination in the off-curve edge term.
- `preprocess_binomial_edges()` selects every incident edge from a curve node to a flagged interior node. That is stronger than "one off-curve edge per adjacent boundary facet" and may include extra interior spokes on an unstructured surface mesh.

## Changes Made / Proposed

- No source-code changes made.
- Proposed next check: test the formula on a simple skew TRI3 patch with a pure along-curve field. The physically correct horizontal binormal DOF should be zero; the current off-curve-only formula generally produces a nonzero value unless the mesh is orthogonal or symmetric cancellation occurs.

## Open Questions

- Should the production mesh generator enforce that off-curve edges at side curves are normal to the curve? If yes, document that as a precondition. If no, the weight computation should use facet-aware Nedelec point-evaluation weights and include curve-edge terms where needed.
- Should multiple adjacent facet contributions be averaged arithmetically, area/angle weighted, or treated as separate constraints? The current `sum()` normalization is not clearly derived from the Nedelec interpolation.

## Follow-Up: Current Gram-Solve Implementation

- The revised Gram-matrix direction is physically better than the earlier per-edge inverse projection, but the current implementation is not correct yet.
- `preprocess_binomial_edges()` forms `G = trans(r) * r` but calls `posv(G, c)` with `c` still holding the 3-vector target `aSign * mConnectorWidth * b`. The correct RHS is `trans(r) * c_target`, with length equal to the number of source edges.
- `posv()` requires a positive-definite square system matching the RHS length. The current call asserts in debug for the common two-source case and risks invalid LAPACK dimensions in release.
- The current source selection admits any number of off-curve incident edges. If more than two are selected, `r^T r` is singular because the vectors live in a local tangent plane. The implementation either needs to select exactly two independent facet-aware source edges or use a solver/constraint rule appropriate for an underdetermined multi-source system.

## Files Updated

- `devlog/dl20260429_connector_weight_physics_audit.md`
- `devlog/README.md`
