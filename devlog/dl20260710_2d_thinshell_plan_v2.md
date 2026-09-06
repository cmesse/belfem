# Devlog 2026-07-10 — 2D Thin-Shell Plan, Second Iteration (v2)

**Date:** 2026-07-10
**Topic:** Second iteration over `todo/2d_thinshell_todo.md` — verify-first items
resolved, plan re-baselined to HEAD, two-week first-order sprint scoped.
**AIs involved:** Claude (orchestration + 4 verification subagents), Codex (technical
audit + prose pass), Grok (third voice, two passes — first run truncated, retry full)
**Claude Confidence:** high (verification results), medium (schedule estimate)
**Codex Audit Confidence:** high
**Grok Audit Confidence:** high (spot-checks), medium ~70% (feasibility)
**Literature References:** Messe et al. 2023 (paper1) — EF_QUAD4TS curl-coefficient
derivation still owed under E5(b); no other literature consulted this session
(planning/verification pass).

## Summary

Planning-only session (no source edits). All five verify-first items from the
2026-07-01 gap analysis were resolved by parallel verification passes, the todo was
rewritten as v2 with a two-week sprint plan targeting first-order (linear) 2D thin
shells, and all line references were re-baselined to HEAD. The three AIs converged on
the headline re-scoping: **milestone M1 (serial, linear, single-layer) does not need
the E1/E2 h_ghost fixes at all**, and **B1 was already fixed** in the 2026-07-02
CW-mesh session.

## Key Findings

**Resolved verify-first items:**

- **B1 — already FIXED on HEAD** (found by Grok, verified by Claude): the 2D branch of
  `compute_element_adjacencies` (cl_CutFactory.cpp:2761-2773) resets element
  containers and calls the new `Mesh::reset_connectivity(ElementToElement)`
  (cl_Mesh.hpp:841) instead of stamping an empty adjacency. M3 (parallel) becomes a
  smoke test, not a fix.
- **E3 — answered:** single-layer 2D shells create NO Ghost sideset (`nb==1` early
  return, cl_ThinShellFactory.cpp:1837; `g>0` gate :289). E1/E2 gate multi-layer (M2)
  only. The v1 milestone M1 exercise list was wrong.
- **D2 — CONFIRMED cross-wired** (two independent traces): QUAD4TS layer insertion
  order n2=top(F1), n3=top(F0) (cl_ThinShellFactory.cpp:1218-1221) +
  `get_top_nodes(QUAD)` facet-CCW order + `to_master_orientation(LINE2)` flip pairs
  top(F0) with the volume node at F1. Fix in sprint days 3–5.
- **A4 — CONFIRMED mirrored:** `orient_terminal_curves_2D` computes
  `cross(shell, tape)` where 3D computes `cross(boundary, shell)` after the a/b swap —
  anticommuted operands, so the 2D terminal-curve direction is flipped relative to the
  validated 3D convention. One numeric sanity run, then flip.
- **B3 — LARGELY REFUTED:** the BC parser doubles 2D terminal lists
  (`tDomainGroupsOut = tDomainGroupsIn`, cl_MaxwellBoundaryConditionFactory.cpp:125,
  concatenated :144-151), so `size()/2 ≥ 1` for every parsed deck. Downgraded to LOW
  defensive validation.
- **E4 — CONFIRMED, worse than suspected:** silent wrong-physics in debug AND release.
  Armadillo expression assignment silently resizes the target
  (cl_AR_Vector.hpp:375-381); BELFEM's asserting `op_VectorPlus` is bypassed by
  expression templates. `b = bt + bn` (mt_maxwell_h.cpp:137/:2271) grows `b` from 2 to
  3 with no diagnostic.
- **E5 — ESCALATED (medium-high):** no live writer of layer-element edge directions was
  found; the legacy tape-roller path (cl_FEM_Element.cpp:90-111 → :1415-1431) is dead,
  so `EF_QUAD4TS::E()`'s assumed `mS[1] = -mS[0]` invariant is violated by the default
  bitset ({-1,-1}). Open caveat → new task T1: 3D PENTA6TS works in production, so
  either a writer was missed or 3D tolerates the defaults; resolve before fixing
  EF_QUAD4TS.

**Auditor deltas:**

- Codex: cut sidesets are excluded from assembly by
  `set_block_types_in_magnetic_equation` (cl_MaxwellFactory.cpp:1737), so the B5
  `mCuts` no-op may be deliberate and C3 (master-only 2D cut facets, confirmed defect)
  is M1-contingent rather than day-1. E6 should be done in the same pass as E4.
- Grok: D1 (`to_master_orientation` 1-edge case, confirmed at
  fn_to_master_orientation.cpp:594-681) gates conductor-above geometries only — the
  air-above path uses the node overload which already handles LINE2. Top schedule
  risk: no 2D thin-shell reference deck exists in the tree → new day-1 task T0.
- Conflict resolved against Codex: Codex called A3 fixed, but it checked the arclength
  LINE3 switch (cl_CutFactory.cpp:1962-1987, which has its `break`); the
  node-container LINE3 switch (:2015-2030) still falls through into the `default:`
  error (Claude direct read; Grok concurs). A3 stays open (order-2, off-sprint).

## Changes Made / Proposed

- `todo/2d_thinshell_todo.md` — rewritten as v2: sprint scope section, day-sequenced
  two-week plan with S/M/L effort tags, new items T0 (author 2D reference deck) and
  T1 (PENTA6TS edge-direction cross-check), B1/E3 ticked, A4/D2/B3/E4/E5 statuses
  updated with evidence, milestone exercise lists corrected, all line refs
  re-baselined to HEAD of 2026-07-10.
- `todo/2d_thinshell_gap_analysis.md` — second-iteration addendum prepended (the
  2026-07-01 body kept as the snapshot record).
- No source-code edits (planning session; user approval pending for the sprint).

## Open Questions

- T1: how does 3D PENTA6TS survive the apparently missing edge-direction
  initialization? (Writer missed vs. defaults tolerated.)
- E5(b): derivation of the EF_QUAD4TS curl coefficients `B(0,0)=B(0,1)=+0.5` against
  the PENTA6TS reduction (Messe et al. 2023, paper1) — owed before M1 physics is
  trusted.
- C3/C4: whether any live 2D consumer dereferences `slave()` on Cut-sideset facets
  (Codex's :1737 exclusion suggests not, but no consumer audit yet).
- D2 numerical impact: code-level cross-wiring is confirmed; the two-element run that
  shows wrong hanging-DOF weights is still pending (sprint days 3–5).
- Two-week feasibility is medium confidence; the deciding factor is T0 (reference
  deck) landing on day 1 (Grok: integration discovery dominates otherwise).

## Files Updated

- todo/2d_thinshell_todo.md
- todo/2d_thinshell_gap_analysis.md
- devlog/dl20260710_2d_thinshell_plan_v2.md (this file)
- devlog/README.md
