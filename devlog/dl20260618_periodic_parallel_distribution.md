# Devlog 2026-06-18 — Periodic CORC in Parallel: Mesh-Distribution Fixes

**Date:** 2026-06-18
**Topic:** Two fixes that get the z-periodic CORC reproducer through MPI mesh distribution on 4 procs (it solved in serial on 2026-06-18, see `dl20260618_periodic_corc_solves.md`). Parallel now reaches the linear solve.
**AIs involved:** Claude (diagnosis + fixes), Codex + Grok (audited the facet-ghosting fix).
**Claude Confidence:** high on both root causes (diagnostic-confirmed) and the fixes (compile clean, run advances past both crashes).
**Codex/Grok Audit Confidence:** high (facet-ghosting fix); the volume-exchange fix was diagnostic-driven, not yet independently audited.
**Literature References:** N/A (MPI distribution mechanics, not a formulation).

## Summary

Parallel CORC hit two sequential aborts in `Kernel::distribute_mesh` / Kernel construction. Both are now fixed; the 4-proc run proceeds through distribution, assembly, and into `MUMPS::solve()`, where it fails with a *separate* solver error (next frontier, below).

## Fix 1 — periodic facet source not ghosted (`cl_Mesh_Distributor.cpp`)

**Crash:** `Tried to access invalid facet id: 87025` on a worker, in `ProtoMesh::create_t_matrices`.

**Root cause:** `Periodicity::set_entity_dependencies` gives a non-hanging periodic *slave* facet B its *master* facet A as a source (`B->add_source(A,…)`, a FACET-typed source). The Distributor's owned-facet loop selected each facet's *node* sources but never called `select_sources(tFacet)` on the facet itself, so master facet A was never ghosted onto the rank owning B → `facet(87025)` missing during t-matrix rebuild. `select_sources` already handles the FACET case (sets `mFacetBitset`); it was simply never invoked for facets.

**Fix:** added `this->select_sources( tFacet )` in the owned-facet loop. The later aura loop then pulls A's nodes and master/slave elements, so A becomes fully usable. **Audited by Codex + Grok (both VERIFIED high).** They also refuted a related claim I made (that the facet/face *positional pairing* in `set_entity_dependencies` was buggy like the edge case) — `map_facets` constructs `mMasterFacets`/`mSlaveFacets` already index-aligned, and `match_nodes` relies on that alignment, so no change was made there. Latent note: element- and controlpoint-typed sources have the same un-ghosted gap, but no current path attaches such sources, so it doesn't bite.

## Fix 2 — volume-exchange array over-sizing with thin shells (`cl_FEM_Kernel.cpp`)

**Crash:** `Tried to access invalid element id: 0` on proc 0, in `compute_element_volumes`.

**Root cause (diagnostic-confirmed):** the per-proc exchange arrays `tMyIDs`/`tNotMyIDs` were **sized** from pass-1 counters (`tNumMyElements`/`tNumNotMyElements`) accumulated over `tMesh->blocks()` **plus** `tShell->blocks()` (thin-shell elements), but **filled** in pass 2 over `tMesh->elements()`, which does **not** contain thin-shell-block elements. So the arrays were over-sized by the owned thin-shell count, leaving trailing **zero** ids. Those zeros were sent to the master, which did `element(0)` → abort. Instrumented proof: worker rank 1 allocated 34392, filled 32518, and the master's zero-id lookups began exactly at position 32518. Triggers only in parallel *with thin shells* — CORC is the first such case to reach distribution; the earlier facet crash had masked it.

**Fix:** size the send arrays from the same `tMesh->elements()` enumeration that fills them (new local `tNumMine`/`tNumOther`); fix the master's loop-2 to size `tProcVolumes` to the received `tProcIDs.length()` rather than the master's own `tNumNotMyElements`; remove the now-dead pass-1 counters. **Behavior-preserving:** thin-shell volumes were never actually exchanged on this path before (the fill never saw them — the tail was just zeros), so this only stops the bogus zeros, it does not change thin-shell volume handling. Compiles clean under `-Werror`.

## Current state — reaches the solve, fails in MUMPS (next frontier)

4-proc CORC now runs the full pipeline up to the linear solve and aborts in `MUMPS::solve()` (`cl_SolverMUMPS.cpp:305`):

```
MUMPS has thrown the error: -1
An error occurred on proc 3.
```

MUMPS `INFO(1) = -1` means "an error occurred on processor `INFO(2)`" (= 3). The wrapper only inspects **rank 0's** `mInfo(0)`, so **proc 3's actual error code is never surfaced** — we don't yet know whether it's a singular/structurally-singular matrix (`-10`/`-6`), insufficient workspace (`-9`), allocation failure (`-13`), etc. This is a *separate* problem from mesh distribution and likely involves the distributed periodic/hanging-DOF matrix structure.

**Recommended next step (with the user):** temporarily print `mInfo(0)`/`mInfo(1)` on *all* ranks (not just rank 0) in `cl_SolverMUMPS.cpp` to recover proc 3's real `INFO(1)`; that single number will route the diagnosis (numeric vs structural vs resource). Not pursued autonomously — it's a new investigation class.

## Open Questions

- Proc 3's true MUMPS `INFO(1)` (blocked on surfacing it).
- Whether the distributed periodic constraint / hanging-DOF assembly is consistent across ranks (prime suspect if the matrix is structurally singular).
- Latent: element/controlpoint source ghosting gap in the Distributor (no current trigger).

## Files Updated

- src/mesh/cl_Mesh_Distributor.cpp (Fix 1)
- src/fem/kernel/cl_FEM_Kernel.cpp (Fix 2)
