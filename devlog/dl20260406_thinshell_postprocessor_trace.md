# Devlog 2026-04-06 — Thin-Shell Postprocessor Trace

**Date:** 2026-04-06
**Topic:** Read-only trace of why `J/Jc` appears on the Hastelloy side of a thin-shell interface
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** Internal routing consulted; see `literature/papers/fem/index.md` (thin-shell entries) and `src/fem/maxwell/doc/postprocessor_recovery_theory.md`

## Summary

Traced thin-shell construction, ghost/Nitsche coupling, and Maxwell postprocessing. The current evidence says the thin-shell layers are decoupled correctly at the H-based DOF level, but the Maxwell postprocessor reconnects them for visualization because it writes global nodal recovery fields on interfaces where the shell layers still share node objects.

## Key Findings

- Maxwell conductor/thin-shell unknowns are edge-based (`edge_h`, plus `face_h` for higher order), not nodal conductor DOFs (`src/fem/maxwell/cl_IWG_Maxwell.cpp:51`).
- Thin-shell inter-layer decoupling is driven by `hasDuplicates` when adjacent materials differ (`src/mesh/cl_ThinShellFactory.cpp:186`).
- For such interfaces, `ThinShellFactory` duplicates edges/faces and creates ghost facets between neighboring shell blocks (`src/mesh/cl_ThinShellFactory.cpp:1563`, `src/mesh/cl_ThinShellFactory.cpp:1655`).
- However, adjacent shell blocks still use the same geometric node layers at their common interface (`src/mesh/cl_ThinShellFactory.cpp:1254`).
- Maxwell postprocessors for thin-shell conductors and superconductors are created separately, but they still write into shared global mesh node fields (`src/fem/maxwell/cl_MaxwellFactory.cpp:1917`, `src/fem/maxwell/cl_MaxwellPostprocessor.cpp:111`).
- The FEM postprocessor is a nodal SPR recovery method and scatters element contributions onto nodes (`src/fem/maxwell/doc/postprocessor_recovery_theory.md:39`, `src/fem/kernel/cl_FEM_Postprocessor.cpp:910`).
- Because shell layers share node objects across material interfaces, a superconducting layer can write `J/Jc` to a node that is also used by the neighboring Hastelloy layer, making `J/Jc` visible on the Hastelloy side even though the solver-side H DOFs are decoupled.

## Changes Made / Proposed

- No source changes made.
- Appended the investigation result to `todo/ai_exchange.md`.
- Proposed likely remedy directions:
  - duplicate thin-shell interface nodes for postprocessing as well, if discontinuous nodal output is required across Nitsche interfaces, or
  - switch `J/Jc` visualization for thin shells to element/discontinuous fields rather than shared nodal recovery fields.

## Open Questions

- Whether the reported “Hastelloy block” is the thin-shell Hastelloy layer block or an adjacent volume block. The tracing strongly explains the thin-shell layer case.
- Whether the current visualization pipeline can consume element fields cleanly enough to make that the preferred fix over node duplication.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_thinshell_postprocessor_trace.md

## Update — same day, after Christian's review

Christian confirmed the diagnosis is consistent and explicitly chose **not**
to implement the duplicate-node fix at this time. Reasons recorded:

- The math (assembled K, solver behaviour) is correct; this is a
  visualization-only artifact.
- A duplicate-node fix would require a corresponding treatment in the
  *thermal* module (which currently relies on the shared layer-interface
  nodes for temperature continuity). Adding a Nitsche-style coupling there is
  not worth the cost since we plan to migrate the thermal solver to the
  finite-volume-element method.

Action taken: documented the issue, the decision, and the future fix sketch
in `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md`, and
registered the new doc in `src/fem/maxwell/doc/README.md`. No source changes.
