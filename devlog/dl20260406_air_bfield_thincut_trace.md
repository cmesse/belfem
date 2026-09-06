# Devlog 2026-04-06 — Air B-Field Thin-Cut Trace

**Date:** 2026-04-06
**Topic:** Read-only trace of how air-domain `B` is computed and why thin-cut postprocessing is currently inconsistent
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high on the execution path, medium-high on the likely recent regression point
**Literature References:** Internal routing consulted; see `src/homology/doc/cohomology_theory_and_implementation.md` and `src/fem/maxwell/doc/postprocessor_recovery_theory.md`

## Summary

Traced the full air-domain postprocessing path. In the current code, air `B` is not recovered directly; the Maxwell air postprocessor first recovers `H = -grad(phi)` from the solved nodal `phi` field and then sets `B = mu0 * H`. Thin-cut healing is incomplete in that path: the old `create_postprocessors()` hook rewires mesh node source containers, but the air postprocessor does not consume those source containers when it builds the element `phi` vectors. The only recent cut-specific postprocessing change is the March 9, 2026 duplicate-aware patch in `cl_FEM_Postprocessor.cpp`, which merges original/duplicate patches but still operates on raw cut-discontinuous element data.

## Key Findings

- Air `B` is computed through the air Maxwell postprocessor, not through a separate air-specific `B` projector:
  - `MaxwellFactory::create_postprocessors()` creates an `Air` `MaxwellPostprocessor` (`src/fem/maxwell/cl_MaxwellFactory.cpp:1863`).
  - `MaxwellPostprocessor::compute_air()` computes `mH = mCalculator->B(aK) * mDOFs` and then multiplies by `-1` (`src/fem/maxwell/cl_MaxwellPostprocessor.cpp:433`).
  - `MaxwellPostprocessor::run()` then writes `Bx = mu0 * Hx`, `By = mu0 * Hy`, and `Bz = mu0 * Hz` (`src/fem/maxwell/cl_MaxwellPostprocessor.cpp:355`).

- The air postprocessor reads raw solved nodal `phi` values directly from the mesh field:
  - `update_dofs_lagrange()` uses `mCalculator->node_data("phi")` (`src/fem/maxwell/cl_MaxwellPostprocessor.cpp:417`).
  - `Calculator::node_data()` forwards to `IWG::collect_node_data()` (`src/fem/kernel/cl_FEM_Calculator.cpp:1978`).
  - `IWG::collect_node_data()` just copies `mMesh->field_data("phi")` at the element node indices (`src/fem/iwg/cl_IWG.cpp:808`).

- Those `phi` field values are written from the magnetic DOFs, including hanging cut DOFs, before postprocessing:
  - The main magnetic `DofManager` and its DOFs are created in `Kernel::create_field()` / `DofManager::set_equation()` before `MaxwellFactory::create_postprocessors()` is called (`src/fem/kernel/cl_FEM_Kernel.cpp:580`, `src/fem/kernel/cl_FEM_DofManager.cpp:91`, `src/fem/maxwell/cl_MaxwellFactory.cpp:475`, `src/fem/maxwell/cl_MaxwellFactory.cpp:631`).
  - After solve, `SolverData::compute_hanging_dofs()` computes hanging-node field values from the DOF source relations and writes them back into the mesh fields (`src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp:2525`).

- Therefore the old `create_postprocessors()` cut-healing hook does not modify the actual `phi` values used by the air postprocessor:
  - That hook rewires mesh node source containers by removing abstract-node sources from hanging nodes (`src/fem/maxwell/cl_MaxwellFactory.cpp:1760`).
  - But the air postprocessor does not read mesh node source containers; it reads the already-written `phi` mesh field through `collect_node_data()` (`src/fem/iwg/cl_IWG.cpp:808`).
  - Since the DOF source relations were created earlier, the mesh-source rewrite is solver-neutral here but also ineffective for the current SPR-based air recovery.

- Thin-cut handling in the current air recovery therefore relies on the recent duplicate-aware patch in `cl_FEM_Postprocessor.cpp`:
  - On March 9, 2026 (`71d521d`), the postprocessor was changed to merge original nodes and their duplicates when selecting owned nodes/elements and when assembling recovery patches.
  - The new logic is in `select_elements_and_owned_nodes()`, `compute_node_matrices()`, and `recover_fields()` (`src/fem/kernel/cl_FEM_Postprocessor.cpp:223`, `src/fem/kernel/cl_FEM_Postprocessor.cpp:606`, `src/fem/kernel/cl_FEM_Postprocessor.cpp:868`).
  - `CutFactory::link_node_duplicates_and_originals()` does populate the `original()/duplicate()` relation for thin-cut duplicates, so this new path is active for cohomology cuts (`src/homology/cl_CutFactory.cpp:2620`).

- But that duplicate-aware patch is not a true cut-healing step:
  - It merges both sides of the cut into the same SPR patch, but it still computes raw element fields from `B * phi` using cut-discontinuous nodal `phi`.
  - For cut-adjacent elements, the abstract jump contribution is carried by the duplicated cut-side nodal values, so the raw per-element `H` is still contaminated by the thin-cut representation before SPR sees it.
  - The current recovery therefore blends raw data from both cut sides rather than first reconstructing a cut-healed `phi` or `H`.

## Changes Made / Proposed

- No source changes made.
- Appended the investigation result to `todo/ai_exchange.md`.
- Proposed likely next fix directions:
  - add a real thin-cut healing stage for air/ferro postprocessing before `compute_air()` evaluates `B * phi`, or
  - build a dedicated postprocessing `phi` field where duplicate cut nodes have the abstract jump removed, and recover `H` from that healed field instead of the raw solved one.

## Open Questions

- Whether the observed regression is exactly the March 9 duplicate-aware patch, or whether that patch only exposed a pre-existing mismatch between the old mesh-source rewrite and the current SPR pipeline.
- Whether the intended long-term solution should heal `phi` first or reconstruct `H` directly from a cut-aware projection.

## Update — after serial rerun

Christian reran the case in serial and did **not** observe the bad `B` field on
the thin cut. That materially changes the priority of the hypotheses.

Revised assessment:

- High confidence the visible bug is now primarily **parallel-only**.
- High confidence the strongest concrete candidate is in the March 9 duplicate-aware
  postprocessor ownership logic, not in the old serial/mesh-source “healing” hook.

Reason:

- In `Postprocessor::select_elements_and_owned_nodes()`, once an owned original
  node is found, **all** of its duplicates are inserted into the
  would-be-owned node bitset without checking duplicate ownership
  (`src/fem/kernel/cl_FEM_Postprocessor.cpp:246-254`,
  `src/fem/kernel/cl_FEM_Postprocessor.cpp:276-284`).
- Those indices become `mMyOwnedNodeIndices` and are sent to rank 0 as if this
  rank owned them (`src/fem/kernel/cl_FEM_Postprocessor.cpp:293`,
  `src/fem/kernel/cl_FEM_Postprocessor.cpp:368-408`).
- But in `recover_fields()`, element contributions are accumulated only for
  duplicates whose `owner() == mCommRank`
  (`src/fem/kernel/cl_FEM_Postprocessor.cpp:982-995`).
- Therefore a rank can advertise a duplicate node as “owned” while never
  computing a coefficient matrix for it. The output row in `tMyData` remains
  unset/default for that slot (`src/fem/kernel/cl_FEM_Postprocessor.cpp:885`,
  `src/fem/kernel/cl_FEM_Postprocessor.cpp:1000-1026`), yet
  `synch_target_fields()` still writes it back on rank 0 using the claimed
  owned-node list (`src/fem/kernel/cl_FEM_Postprocessor.cpp:146-183`).

Implication:

- In parallel, the same cut duplicate can be claimed by multiple ranks.
- The real owner may send the correct recovered value, while another rank sends
  an empty/default row for the same node.
- Rank-0 gather order then decides which value survives, producing a
  partition-dependent postprocessing artifact.

This behavior is fully consistent with “bad in parallel, clean in serial.”

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_air_bfield_thincut_trace.md

## Update — review of Claude's `cl_FEM_Postprocessor.cpp` fix

I reviewed Claude's current patch in
`src/fem/kernel/cl_FEM_Postprocessor.cpp`.

Verdict:

- High confidence: I agree with the fix.

Reason:

- The earlier diagnosis was that
  `select_elements_and_owned_nodes()` advertised duplicate cut nodes as
  locally owned without checking `duplicate->owner()`.
- `recover_fields()` never accumulated coefficients for those off-rank
  duplicates because it already gates accumulation by local ownership.
- `synch_target_fields()` then wrote back rows using the inconsistent
  owned-node list, allowing an off-rank/default row to overwrite the real
  owner's recovered value on rank 0.

Claude's change fixes exactly that mismatch:

- duplicate nodes are now inserted into `mNodeBitset` only when
  `tDup->owner() == mCommRank`
  (`src/fem/kernel/cl_FEM_Postprocessor.cpp:251-263`,
  `src/fem/kernel/cl_FEM_Postprocessor.cpp:291-299`).

Assessment:

- This is the right minimal fix for the currently visible parallel thin-cut
  artifact.
- It makes the ownership definition used during node selection consistent with
  the ownership check already used during coefficient accumulation.
- Residual conceptual questions about thin-cut healing in the serial air
  postprocessor still exist, but Christian's serial rerun strongly suggests
  those are secondary to this MPI ownership bug for the present symptom.
