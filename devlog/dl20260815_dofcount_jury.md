# Dof-Count Jury: the Matrix Size Depends on the Rank Count — RETRACTED

**RETRACTION (same day, Christian):** the premise was false. The 8-proc
figure of 3,445,520 was a combined count including the temperature problem's
dofs: 2,465,712 + 979,808 = 3,445,520 exactly. The two numbers compared in
the brief were never the same quantity. There are no phantom dofs, no
rank-dependent matrix size, no defect. DR-73 is struck.

**The lesson, which is the reason this record stays:** the false number
entered through the jury brief's "established facts — do not re-litigate"
section, placed there by Claude from a secondhand reading of a log that had
already been overwritten. Blind jurors cannot catch contamination in the
shared brief — that section must contain only facts the brief's author has
verified firsthand. The tell was available: all three tracers independently
concluded that rank-0 dof creation SHOULD be partition-independent, and that
unanimous static conclusion deserved more weight than an unverifiable
observation. When every juror says the code cannot produce the phenomenon,
re-examine the phenomenon before hunting the mechanism.

**What survives (verified, still useful):** the printed dof count IS the
allocated SpMatrix dimension; the 4-proc magnetic count is exactly right
(memdump closure); the by-catch list at the end (uninitialized ctor member,
dead PerProc tables, no MPI dof-invariance test, and — newly relevant — the
dof report not naming which kernel it belongs to, which is what enabled the
misreading).

The original record follows as written.

---


**Date:** 2026-08-15
**Topic:** Three-juror blind sweep of the magnetic dof count differing between
4 and 8 MPI ranks on the same mesh — verdict: functional, not cosmetic
**AIs involved:** Claude, Codex, Grok — parallel blind (jury mode), reconciled
**Verification:** read-only static traces + memdump field forensics on the
running job; "reviewed" not "verified"; the discriminating probe is pre-planned
but NOT run (machine occupied by the weekend quench study)
**Register:** DR-73 (P1, blocking)

## The observation

Same remeshed tapestack3d mesh, same deck:

| | 4 procs | 8 procs |
|---|---|---|
| magnetic free dofs | 2,465,712 | 3,445,520 (+979,808) |
| condensed (hanging) | 1,225,062 | 1,225,062 |
| thermal free dofs | 904,949 | 904,949 |

## Ground truth from the memdump

Fields live on mesh entities, so they are partition-independent; at t = 1.7 s
every live dof carries a nonzero value:

    nonzero(phi) + nonzero(edge_h) = 319,750 + 3,371,025 = 3,690,775
                                   = 2,465,712 + 1,225,062 + 1   exactly

The 4-proc count is exactly right — closing to the single fixed dof — and the
8-proc count carried 979,808 phantom entries.

## Unanimous verdict: functional, not cosmetic

All three jurors traced the same chain independently: `reorder_dofs` counts
the free-dof graph on rank 0 (`cl_FEM_DofMgr_DofData.cpp:3279-3296`, the only
nonzero write), `allocate_matrices` builds `SpMatrix(n,n)` from that value
(`cl_FEM_DofMgr_SolverData.cpp:351-406`), and the report prints the same live
reference (`:552`). There is no stale intermediate. **The printed number IS
the matrix dimension passed to the solver** — which answers the original
question ("are we printing something else?") with: no, and that is exactly why
this matters. The 8-proc system was genuinely allocated at 3,445,520².

## Mechanism: NOT pinned — and the honest tension is recorded

All three tracers independently predict rank-0 dof creation to be
partition-independent (rank 0 counts from the full master mesh; the
partitioner only rewrites ownerships; workers' aura appends never reach
rank 0's container). The observation contradicts the prediction, and no juror
found the mutation statically. Two candidate natures survive:

- **(a) Duplicate-ID Dof objects** (Grok's lead): `count_edge_dofs` emits per
  edge *object*; twin edges on enriched meshes share node pairs (established
  house knowledge). Duplicates lose in `mDofMap` (last-wins), keep empty
  adjacency, and become **empty CSR rows — the 8-proc matrix would be
  singular** and fail its first factorization. Never tested: that run was
  OOM-killed before factorizing.
- **(b) Genuinely extra unique edge dofs**, assembled but with their
  condensation constraints dropped — a partition-dependent discretization
  (DR-08's defect class at the mesh level).

## Claims killed during reconciliation

- *"The phantoms are exactly the hanging edge dofs"* (Codex, relayed by
  Claude as confirmed) — **withdrawn**: the memdump closure holds for ANY
  split of the hanging count into node/edge parts, so it carries no
  information about the entity class. Numerology.
- *"The parallel path leaves hanging dofs in mDOFs"* (the brief's landmark) —
  refuted by both external jurors: `remove_hanging_dofs_from_container()` runs
  on all ranks at `:2839` (in the tree since 2024/2025 per git blame).
- *"The old 8-proc campaign converged, so the phantoms can't be empty rows"*
  — over-claim: that was the previous mesh; the remeshed 8-proc system has
  never been factored.
- Claude's leftover-edge-flag hypothesis — dead: the T-matrix pass pre-clears
  (`unflag_all_edges()` at `:3490`); the flag is an intra-pass visit marker.
- The brief's phantom count 979,807 — arithmetic error (Claude); it is
  979,808.

## The discriminating probe (unanimous, pre-planned, not yet run)

One rank-0 debug print before hanging removal — `mDOFs.size()`,
`mDofMap.size()`, `count_edge_dofs` return, `mMesh->number_of_edges()` — then
an 8-proc SETUP-ONLY run killed after the print (no factorization, modest
memory). `mDOFs.size() != mDofMap.size()` names (a); equality names (b).
Machine is occupied by the weekend quench run; the probe fits any pause.

## Fix direction (Grok J5, to be re-audited when implemented)

Authoritative global count = unique IDs in rank-0 `mDofMap` after hanging
removal, excluding fixed; emit-once per `(id, type)` in `count_edge_dofs`
step 5; debug asserts `mDOFs.size() == mDofMap.size()` after creation and
after removal. Regression gate: a small thin-shell fixture through DofManager
at 1/2/4 ranks asserting identical free/hanging/fixed counts — **no existing
test covers dof-count invariance under MPI** — plus an unchanged
`make check-fast`.

## By-catch (recorded, not acted on)

- `mNumberOfFreeDofs` is uninitialized in the DofData constructor — harmless
  today (zeroed before use), a landmine if anything reads the SolverData
  reference earlier.
- `mNumberOfFreeDofsPerProc` / `mNumberOfFixedDofsPerProc`
  (SolverData.hpp:93-94) are populated and never read — dead weight.
- The memory-sizing chain explains part of why the 8-proc attempts were so
  much heavier than rank-scaling predicted: they were building a 40 % larger
  system on top of the per-rank mesh replication.

---

## Addendum: by-catch cleaned up (2026-08-15, Christian's go-ahead)

Three of the surviving items were fixed the same evening (syntax-gated, not
yet compiled into a run):

1. **Dead PerProc tables removed** — `mNumberOfFreeDofsPerProc` /
   `mNumberOfFixedDofsPerProc` (SolverData.hpp) plus their collect/send
   exchange in `compute_rhs_sizes`. Verified written-never-read before
   removal; both the master and the worker side of the exchange went
   together, so collective consistency is preserved. A one-line comment
   marks the removal site.
2. **DofData counters now have in-class initializers** ( = 0 ) — SolverData
   binds const references to them at DofManager construction, before
   `create_dofs()` runs, so an early read must see zero, not garbage. The
   comment states that reason.
3. **The dof report now names its system**: `Number of Degrees of Freedom
   ( edge_h, phi, ... ):` — built from the IWG's dof field labels, capped at
   four plus ellipsis. This is the ambiguity that enabled the retracted
   DR-73 comparison; two kernels print this block back to back and the
   counts were compared as if they were the same quantity.

Not acted on: the MPI dof-invariance test gap (a real fixture, needs its own
session).