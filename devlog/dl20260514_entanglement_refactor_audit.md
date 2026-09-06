# Devlog 2026-05-14 - Periodic DOF Entanglement Refactor Audit

**Date:** 2026-05-14
**Topic:** Audit of refactored `DofData::entangle()` periodic DOF handling
**AIs involved:** Codex
**Codex Audit Confidence:** high for representative-direction finding; medium for both-hanging algebra risk

## Summary

Performed a read-only source audit of the latest Dof/DofData changes and built `hphirun`.
The refactor compiles, and `Dof::reset_sources()` is a reasonable helper for controlled source replacement, but the new hanging-case logic changes the intended smaller-ID representative semantics.

## Key Findings

- `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:3568` enters a symmetric merge whenever either periodic DOF already hangs, then resets both DOFs at `:3599-3602`. This can make the smaller-ID/free representative hang on the other DOF's old sources instead of making the larger-ID DOF hang on the smaller-ID representative.
- `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:3588-3597` combines both source lists and normalizes by `sum(aWeights)`. If both DOFs already have distinct source constraints, this creates an averaged constraint that neither DOF originally represented; if the summed weights cancel, it can divide by zero.
- `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:3565` checks only entity type, not DOF type id. The caller loops matching local DOF index today, but a `type_id()` assertion would document and defend the invariant.

## Changes Made / Proposed

- No source changes made.
- Added this devlog entry and updated `devlog/README.md`.
- Verification run: `make -C cmake-build-debug hphirun -j2` passed.

## Open Questions

- The intended behavior should be restated in code as: choose the smaller-ID DOF as representative; only the larger-ID DOF receives/replaces sources, inheriting the representative's unfolded non-hanging sources when the representative is already hanging.

## Follow-Up

- Revised policy: hanging constraints should override raw smaller-ID direction when one or both periodic DOFs already hang, because h-phi interface edge DOFs can legitimately be periodic and hanging on both sides.
- Latest implementation fixes the temporary-container compaction bug and changes the zero guard to a norm check, but still normalizes source coefficients by `norm(aWeights)`. That changes valid signed edge/Nedelec constraint scaling and remains a correctness blocker.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed.

## Second Follow-Up

- Latest counter-based version removes norm scaling, but has a loop bound bug: the `aB` accumulation loop is bounded by `aA->number_of_sources()`. This breaks one-hanging/free cases and can either skip `aB` sources or read past `aB` sources.
- Counter division is only applied inside the zero-pruning branch. If no coefficient is pruned, duplicate-source contributions remain summed rather than averaged.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed; these are runtime/semantic issues, not build failures.

## Third Follow-Up

- Latest version fixes the asymmetric loop-bound issue and models free DOFs as unit rows in the merge, but this creates a self-source in exactly-one-hanging cases when the merged row is assigned back to the free DOF.
- Counter division is still only applied in the zero-pruning branch; identical duplicate rows can remain doubled if no coefficient is pruned.
- `aWork.reserve(...)` is used without clearing the reused work cell, so behavior depends on later zero-pruning to remove stale entries and can leave temporary flag bit 1 set when no pruning happens.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed; remaining findings are semantic/runtime issues.

## Fourth Follow-Up

- Latest version clears the reused work cell and adds a guard rejecting exactly-one-hanging periodic entanglement, which prevents the immediate self-source bug.
- The guard is only acceptable if periodic pairing after `create_dofwise_t_matrices_master()` is guaranteed to produce both-free or both-hanging pairs. If asymmetric hanging can occur for valid periodic pairs, this becomes a hard runtime failure.
- Counter division remains inside only the zero-pruning branch; both-hanging rows with overlapping nonzero sources are still summed rather than averaged when no zero entry is pruned.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed.

## Fifth Follow-Up

- Latest version applies counter division in both the zero-pruning and no-pruning paths, resolving the doubled-row issue for overlapping source DOFs.
- Remaining minor risks are mostly clarity/robustness: the exactly-one-hanging guard should be documented as a geometric invariant, and the temporary use of DOF flag bit 1 should leave all retained entries unflagged in both branches.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed.

## Sixth Follow-Up

- Cross-review with Claude's notes: the exactly-one-hanging guard remains defensible for core Maxwell h-phi/thin-shell edge geometry.
- The both-hanging union/average rule has an unresolved semantic risk. Periodic node DOFs are entangled before periodic edge DOFs, so an edge's node-source DOFs can already be hanging by the time the edge merge reuses them. `Dof::set_sources()` rejects hanging source DOFs, so the edge merge likely needs source-row unfolding/canonicalization rather than simply unioning raw source pointers.
- The current counter scheme is per-source averaging, not row averaging. If the two rows have disjoint periodic-image source sets, it preserves each coefficient at full strength instead of forming `(rowA + rowB) / 2`; this needs an explicit physics decision or a regression that exercises the both-hanging periodic thin-shell/interface case.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed.

## Seventh Follow-Up

- Claude's recursive unfolding proposal addresses the concrete hanging-source assertion path, but needs two implementation corrections.
- `Vector<real>` does not expose `push()`, so dynamic accumulation either needs `Cell<real>` work buffers or row-local maps followed by sizing/copying into `Vector<real>`.
- Counters should not increment for every recursive leaf occurrence before row consolidation. A single row can unfold to the same terminal source through multiple branches; those weights must first sum within that row, then the two rows can be averaged or merged by policy.
- The core Maxwell edge ordering appears favorable: periodic nodes are processed before edge rows that source nodes, and conductor/source edges are matched before thin-shell edge modes. This reduces but does not eliminate the need for canonical unfolding in the merge.

## Eighth Follow-Up

- Claude's implemented patch uses row-local `Cell<real>` buffers and row-local maps, so the previous `Vector::push()` and leaf-counter concerns are resolved.
- Remaining issue: after per-row averaging, near-zero merged coefficients are not pruned and a fully cancelled row is not rejected. `mMerged.size() > 0` does not catch rows whose coefficients all cancel to zero.
- Low-level robustness note: `unfold_row()` uses `BELFEM_ASSERT` for recursion depth, which is disabled in release builds. A corrupted/cyclic hanging graph would stack-overflow in release instead of raising a controlled error.
- Follow-up verification run: `make -C cmake-build-debug hphirun -j2` passed.
