# Step 5f — periodic seam edge sourcing fix (fix b)

**Date:** 2026-06-16
**Purpose:** Implement and (attempt to) validate the Step 5f fix for the periodic-edge
continuity bug; record the input-config blockers found along the way.
**Module:** src/mesh (PeriodicityFactory), with input.conf parse findings

---

## What was done

### Fix (b): hybrid edge source in `PeriodicityFactory::collect_edges()`

`collect_edges()` built the periodic seam edge list from each facet's **wrapper**
element (`tFacet->edge(e)`). After the CutFactory edge-wipe +
`create_edges_and_faces_on_mesh()` conductor-volume-only rebuild, the wrapper edge
containers go **stale**: `#DIAG 5f` showed 996/996 periodic facets edge-less, all
TRI3, all with `master_has_edges==1`. The old `has_edges()` skip would silently drop
every conductor seam edge.

**Implemented** (`src/mesh/cl_Mesh_PeriodicityFactory.cpp`): a `tCollectFacetEdges`
lambda used by all three flag-dedup passes —

- wrapper carries edges → use `tFacet->edge(e)` (unchanged, pass-1 behavior);
- else master carries edges → `master()->get_edges_of_facet(index_on_master())`;
- else → skip (genuine non-conductor: air / ferro / void / buffer).

The master branch is **gated on TRI3** (`number_of_nodes()==3`): admits TET4 (any
face) and PENTA triangular caps (3/4); a QUAD4 wrapper here would be a PENTA lateral
(0–2) that `PENTA6TS::get_edges_of_facet` cannot enumerate, so it raises a **loud
`BELFEM_ERROR`** pointing at the deferred thin-shell-wrapper edge-linking work rather
than throwing cryptically or dropping edges. The §5.1 invariant (shell normal never
in a periodic plane ⇒ 3D seam facets are always TRI3↔TET4) guarantees the gate is
never hit in valid meshes.

### Audit (both auditors)

Codex + Grok independently CONFIRMED Q1 (master `get_edges_of_facet` returns the same
canonical mesh `Edge*` the wrapper held), Q3 (flag dedup preserved), and REFUTED Q4
(no deliberate wrapper-only reason). Both flagged the same two boundary caveats:
- **Q2:** a material-asymmetric seam plane could desync source/target edge counts
  (debug-assert `:856/857`, silent in release). Not physical for periodic partners;
  absent in the uniform repro. Documented, not guarded.
- **Q4 (Codex):** must not become the general path for PENTA QUAD wrappers — handled
  by the TRI3 gate.

### Build

`make reset && mkdir -p bin lib && make hphirun -j 20` → success (192 MB binary).
`mpicxx -std=gnu++17 -fsyntax-only` clean in debug and release (`-Wall -Werror`).

---

## Validation: BLOCKED (config, not the fix)

The run could not validate fix (b) because the handed-over `input.conf` is a
**non-periodic** variant. Getting it to run at all required clearing a chain of
input/material parse failures, all caused by the active materials refactor (commit
`fd46a497`) being stricter than the config:

| # | Symptom | Site | Config fix applied |
|---|---------|------|--------------------|
| a | `Nodes undefined for gauge/bearing` | `cl_MaxwellBoundaryConditionFactory.cpp:51` | `bearing { node : 9 }` → `nodes : 9` |
| b | `materials must have custom section or builtin key` | `cl_MaterialFactory.cpp:203` | added `builtin : <name>` to copper/silver/ybco (ybco lost stray `v` label) |
| c | `Key hastelloy not found in map` | `cl_Map.hpp:230` (from `create_and_assign_materials`) | `hastelloy/magnesia : builtin ;` → brace + `builtin :` form |

After (a)–(c) the run reached the **solve** and ran cleanly to **timestep 33**;
fix (b)'s new `BELFEM_ERROR` never fired (good — fix doesn't break the non-periodic
path).

**But periodicity is entirely inactive in this config.** It requires
`topology { periodic { source : <3 node ids>; target : <3 node ids>; } }`
(`cl_MaxwellFactory.cpp:111-113`, `:381-399`) — **absent** from the handed-over file
(confirmed via backup). The run shows `#DIAG A/B periodic faces : 0`,
`#DIAG 4j asym/sym 0`, and **zero** `#DIAG 5f` calls. So fix (b)'s code path was
never exercised. The `out_grep5f.txt` run (996 periodic facets) used a config that
**had** the `periodic` subsection; the source/target plane node-IDs for `corc.msh`
are not recorded anywhere findable.

---

## State left for the user

- `src/mesh/cl_Mesh_PeriodicityFactory.cpp`: fix (b) + the `#DIAG 5f` probe (probe to
  be removed in Step 8 cleanup).
- `cmake-build-debug/input.conf`: edited to parse/run (a–c); original preserved at
  `cmake-build-debug/input.conf.orig_handoff` and `/tmp/input.conf.orig`.
- Docs updated: `src/mesh/doc/thin_shell_geometry_and_periodicity.md` §5.1 (normal /
  TRI3↔TET4 invariant).

## Next steps

1. **User:** supply the `topology{periodic{source;target}}` node-IDs that reproduce
   the 996-facet periodic case (or the `out_grep5f.txt` input.conf). Then re-run to
   validate fix (b) on the real periodic path (expect: clean periodicity update,
   past CutProcessor, into the solve — with `#DIAG 5f SUMMARY` now showing the
   edge-less facets handled, not skipped).
2. Decide canonical material/bearing input format vs the refactored code (the
   `name : builtin ;` colon-form not registering looks like an `fd46a497`
   regression — separate from the periodic work).
3. Step 6c (fragmentation extra dups) and Step 8 cleanup (remove probes) still pending.

---

## Addendum 2026-06-17 — fix (b) VALIDATED at runtime; new blocker exposed

The user added the periodic block (`topology { periodic { source : 1,2,3; target : 5,6,7 } }`).
Re-run with periodicity active (279/54/27/101 periodic faces):

- **Fix (b) is validated.** `#DIAG 5f` pass 2 shows all 996 thin-shell-pass wrappers
  master-sourced; the TRI3 gate never fired; the run advanced past `collect_edges`.
- **Decisive faithfulness test (`#DIAG 5h`):** `wrapper-vs-master-node-mismatch = 0/996`
  in both passes — `master->get_nodes_of_facet(index_on_master())` is **identical** to
  the wrapper element's nodes, so master-sourcing returns exactly the edges the wrapper
  would have. Master-source is faithful (**H1 refuted, H2 confirmed**); fix (b) is correct.

**New blocker (pre-existing, not caused by fix b).** The run now reaches
`match_edges` (`cl_Mesh_PeriodicityFactory.cpp:904`) in the **post-cut** periodicity
pass and asserts on unequal seam-edge counts:

| pass | nodes | edges |
|------|-------|-------|
| 1 (pre-cut) | 536/536 ✓ | 1534/1534 ✓ |
| 2 (post-cut) | 846/846 ✓ | 1782/1797 ✗ |

`#DIAG 5g` (paired-index edge-key set diff): **src-only 565, tgt-only 580** (net +15).
So ~1145 seam edges have no periodic partner — a *large* cut-induced restructuring of
the seam edge topology, not a few jump edges. Pass 1 is perfect, so it is entirely
cut-induced. Because fix (b) is faithful (5h), the original wrapper-based code would
hit the **same** wall — it just crashed earlier on edge-less wrappers.

This is the **edge-level** manifestation of the period-direction cut jump (handled at
the node/CutSet level in Step 6b), now in periodicity edge-matching. It is a Step 6/7
**design decision** (next section of the plan), not a unilateral fix:
- (i) `match_edges` skips period-wrapping/cut-crossing seam edges with no partner
  (relax the equal-count assert + null-safe lookup), or
- (ii) seam duplication (6b) also pairs/symmetrizes seam **edges**, keeping the seams
  isomorphic.

Recommended next: diagnose the 565/580 edges (all cut-crossing? endpoints on in-plane
cuts 0/3 vs the period cut?) analogous to the node-level diagnosis, then Codex/Grok
review before implementing. Probes 5g/5h are in place (removal: plan 8i).
