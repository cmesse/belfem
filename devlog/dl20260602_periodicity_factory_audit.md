# Devlog 2026-06-02 — PeriodicityFactory audit and fixes

**Date:** 2026-06-02
**Purpose:** Record the audit of `PeriodicityFactory::create_periodicity()` (now node→edge→facet matching) and the correctness/robustness fixes applied during the session.
**Module:** `src/mesh/` (`cl_Mesh_PeriodicityFactory.cpp`)
**AIs involved:** Claude (primary audit + edits), Codex (independent cross-audit)
**Confidence:** high for the applied fixes; medium for the deferred thin-shell item until tested on a mesh with in-plane coincident duplicates.

## Context

`PeriodicityFactory::create_periodicity()` is the pre-cohomology geometry matcher invoked from
`CutFactory::run()` (`cl_CutFactory.cpp:134`). It sets `set_periodic()` on nodes/edges/facets
only — it never touches master/slave. The session started from a `hphirun` diagnostic showing
only 313 of 536 source nodes matched, and progressed through several rewrites of the factory.
It remains WIP: the function still ends in a debug print loop + `exit(0)`, and edge/face
derivation is inline scaffolding that per the plan (C1) must eventually move into
`Periodicity::update()` (it has to run twice — pre- and post-cohomology). See
`todo/periodic_bc_fix_plan.md`.

## Root cause of the 313/536 diagnostic

`match_nodes` (then `map_nodes`) bounded its facet-pair loop by the **node** count
(`n = aSourceNodes.size()`) while indexing the **facet** arrays. With ~1000+ facet pairs and
only 536 nodes, the loop processed the first 536 facet pairs; the 313 distinct nodes living in
those facets matched (a contiguous prefix), and the remaining ~223 nodes — first appearing on
later facets — were never visited. Fixed by bounding the loop with `aSourceFacets.size()` and
adding a `tCount == n` completeness assert.

## Fixes applied this session

- **Loop bound + completeness** — facet loop now bounded by `aSourceFacets.size()`; assert that
  all `n` nodes matched.
- **`tElement` typo** in `create_facet_map` debug block → `tFacet` (was a debug-build compile
  error; release skipped it under `-DNDEBUG`).
- **`create_facet_map` self-indexes** its node argument (was relying implicitly on indices set
  by `match_edges`).
- **3D facet key** uses `number_of_corner_nodes()` (3 nodes give a unique key; quad orientation
  is acknowledged as unknowable).
- **Assert message args** in `match_nodes` un-swapped (`need n, have tCount`).
- **Dead `tWork`** removed from `match_edges`; **`flag(2)` leak** in `create_facet_map` cleared
  via `aNodes` at the end.
- **MED-1 — facet-centroid match tolerance.** `find_closest_node` accumulates a **squared**
  distance; the acceptance test compared it against the linear `BELFEM_MESH_EPSILON` (effective
  ~32 µm vs the 1 nm used for node matching). Fixed by comparing against
  `BELFEM_MESH_EPSILON * BELFEM_MESH_EPSILON`. An intermediate fix had instead `sqrt`-ed inside
  `find_closest_node`, which broke the kd-tree pruning test (`tDelta*tDelta < aBestDist` mixed
  squared vs linear — latent over-pruning for meshes with distances > 1 m); reverted to a fully
  squared formulation so accumulate/compare/prune are all consistent.
- **MED-3 — positional-correspondence assert.** Added a debug check that
  `aSourceNodes(k)->index() == k` after matching, documenting and guarding the unstated invariant
  that `match_edges`/`match_facets` rely on (match order == `collect_nodes` order). The
  `tCount == n` assert catches a count mismatch but not a reordering.
- **`set_size(n, nullptr)` no-op** in `match_nodes` documented (it must not be "simplified" to
  clear the array — the following loop dereferences the collected target nodes still stored
  there).

## N10 / MED-2 — clobbered global indices (addressed separately by the user)

`create_periodicity` calls `update_node/edge/facet_indices()` up front, then overwrites
periodic-plane **node** and boundary **facet** indices with local `0..n-1` match orderings.
Verified that downstream consumers in the cut window need canonical node indices —
`SimplicialComplex` keys chains/cochains directly on `node->index()`
(`cl_SimplicialComplex.cpp:139,193-194,…`), and there is no re-index between
`create_periodicity()` and `compute_cohomologies()`. The user added
`update_node_indices()/update_edge_indices()/update_face_indices()` at the top of
`compute_cohomologies()`. Notes:

- Node re-index is the load-bearing one; placement there is fine because the only earlier
  consumers (RCM/Poisson ordering) are performance-only and self-heal index uniqueness
  (`symrcm`, `rearrange_nodes`).
- `update_face_indices()` (interior `Face`) and `update_edge_indices()` are harmless no-ops here
  — the factory never clobbers Face or Edge indices. The genuinely-dirtied **boundary
  `Facet`** indices come solely from `A->set_index(k)/B->set_index(k)` in `map_facets`, which is
  **unread** anywhere in the cut window (facet code uses `index_on_master()`, not `index()`; the
  only `Facet::index()` reader, `SideSetFactory.cpp:107`, is outside `CutFactory::run()`). Those
  two lines were left in place by choice (potential facet→pair-index hook); consequence is that
  the facet-index restore stays necessary.

## Decisions / deferred

- **N2** (unchecked `Map` lookups in `match_edges`/`match_facets`) and **N7** (validity guards as
  `BELFEM_ASSERT` rather than `BELFEM_ERROR`): closed by decision — assume a well-formed periodic
  mesh; garbage twins should not be constructible.
- **Double-match in `map_facets`** (two source centroids resolving to one target): not added —
  same class as N7, only bites on a degenerate mesh.
- **N5 — thin-shell coincident-duplicate node mispairing:** deferred and logged in
  `todo/periodic_bc_fix_plan.md`. `match_nodes` takes the first unflagged in-plane coordinate
  match, so where a tape pierces a periodic plane `P_top` can bind to `Q_bottom`. The
  completeness assert guarantees the pair count, not topological correctness; a wrong pairing
  propagates silently into the edge/facet keys. Likely fix: side-aware / global-kdtree
  disambiguation. Confirm whether `corc` actually contains coincident in-plane duplicates first.

## Still WIP (not bugs)

- Debug print loop + `exit(0)` at the end of `create_periodicity()` must be removed before the
  path goes live; at that point the `compute_cohomologies` re-index becomes load-bearing.
- Inline `match_edges`/`match_facets` to move into `Periodicity::update()` per plan C1.

## Files updated

- `src/mesh/cl_Mesh_PeriodicityFactory.cpp` (fixes above)
- `todo/periodic_bc_fix_plan.md` (N5 deferred subsection, logged 2026-06-02)
- `devlog/dl20260602_periodicity_factory_audit.md` (this file)
- `devlog/README.md` (index entry)
