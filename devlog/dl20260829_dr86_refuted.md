# Devlog 2026-08-29 -- DR-86 Refuted and Struck

**Date:** 2026-08-29
**Topic:** the aura-phi lag hypothesis has no write path; row struck as refuted
**AIs involved:** Claude (source trace); the refuted claim was Grok's, from the 2026-08-16 DR-84 round
**Verification:** static source trace against the tree at HEAD; no probe run, no solver run

## Summary

DR-86 held that on partition boundaries the aura copies of the `phi` node field feeding
`compute_hn` could lag one step, because `phi` was written per accepted iterate for *owned*
nodes but distributed to auras only at step finalize. Christian challenged the premise
directly -- "but we do distribute the fields at the end of each iteration" -- and the trace
confirms him. Both halves of the mechanism are false and both of the row's line anchors were
misread. The row is struck as refuted and its probe retired.

## The trace

**The distribution is per iterate, not per step.** `DofManager::solve()` ends with
`mFieldData->distribute( mIWG->all_fields() )` followed by `comm_barrier()`
(`cl_FEM_DofManager.cpp:1039-1054`), and so does `solve_from_residual()` (`:1074-1088`).
The latter is the h-phi path: `IwgMode::Iterative` with one RHS column is intercepted at the
top of `solve()` and routed through the two-phase `compute_residual` / `solve_from_residual`
pair (`cl_FEM_DofMgr_SolverData.cpp:2071-2078`). `phi` is a dof field, so it is in
`mAllFields` (`cl_IWG.cpp:1343-1368`).

**There is no owned/aura split to desynchronize.** This is the half that kills the row
outright. Every field write in the solver sits behind `if ( mKernel->is_master() )` --
`SolverData::solve()` at `:2092`, `SolverData::solve_from_residual()` at `:2416`, the latter
covering both the NewtonRaphson update loop (`:2456-2459`) and the Picard one (`:2508-2511`).
Workers never write `phi` at all. On the receiving side, `FieldData::distribute` does
`receive( mMesh->field_data( ... ) )` (`cl_FEM_DofMgr_FieldData.cpp:1206-1215`), replacing the
whole vector; rank 0 builds each worker's payload from `comm_table( aTarget )->nodes()`
(`:1229-1232`), which is that proc's complete node list, aura nodes included. So a worker's
`phi` is wholesale overwritten from the master's freshly-solved global field every iterate --
there is no partial write that could leave auras behind.

Closing the loop, `compute_hn` reads exactly that vector: `get_normal_calculator` pulls
`this->group()->parent()->mesh()->field_data( "phi" )` (`cl_FEM_Calculator.cpp:2009`) and fills
`phi_m` / `phi_s` by node index.

## Both anchors were misread

- `cl_FEM_DofMgr_SolverData.cpp:2236` -- cited as the owned-nodes-only write. It is inside the
  `else // right hand side is matrix` branch: the multi-column RHS path, asserted
  `IwgMode::Direct`, which is not the h-phi Newton/Picard path. It is rank-0-only besides, and
  it loops `mParent->mesh()->nodes()`, which on the master is every node.
- `cl_FEM_Controller.cpp:4100` -- cited as the step-finalize aura distribution. That line is
  deck timestep parsing. The controller's three `distribute_fields` calls (`:4540`, `:4879`,
  `:4892`) are all in the memdump-restore path. There is no step-finalize aura synch anywhere;
  the synch was always in `DofManager::solve`.

## How the row survived four sessions

The row recorded its own weakness accurately and nobody acted on it: "plausibility-checked not
traced-to-ground by the other voices." It was registered on 2026-08-16 as by-catch of the
DR-84 round, carried through the DR-84 strike on 2026-08-27 as a deliberately separate row, and
sat on the watch list with no symptom and no probe. Its `[W]` tag was doing exactly what a
watch tag is for -- keeping an untraced hypothesis cheap -- but a `[W]` row is still a claim in
the register, and this one cost a reader the assumption that MPI thin-shell assembly had a
known open hole. The cheap thing was never the probe; it was the ten minutes of reading that
the row's own confidence note said had not happened.

## Surviving caveat, deliberately kept

The refutation rests on the solve being master-centric. `cl_FEM_Controller.cpp:4534-4539`
already carries the same warning for the history-restore twin: "a future fully-distributed
solver would break this and must revisit." If the solve ever goes fully distributed, this
concern becomes real for the first time -- and wants a new row, not a revival of this one.

## Evidence posture

Static trace only. No probe was run and no solver was launched, so this is **reviewed, not
verified** in the ladder's sense. What makes it strong enough to strike on is that the argument
is structural rather than statistical: a rank guard and a whole-vector replace, both readable in
full. The registered probe (`norm(phi_aura - phi_owned_peer)` mid-step) would have been
measuring a quantity that no code path can make nonzero.

## Files touched

- `todo/debt_register.md` -- DR-86 row removed from the live table; `[W]` count 7 -> 6 with the
  refutation recorded in the preamble.
- `todo/debt_register_closed.md` -- struck DR-86 row appended, full trace in the status cell.
- `doc/lessons_learned_evidence.md` -- INC-331's residue clause updated: the aura-phi hypothesis
  is refuted rather than open. The same clause's second half was also stale and is corrected --
  the normal-component fix's executable gate passed on 2026-08-23, per DR-84's closure record.

No source files were touched.
