# BFM Load Profiling — Speedup Options

**Date:** 2026-07-17
**Purpose:** Analyze `cmake-build-debug/corc/profiler.callgrind` of `hphirun` to find why loading `corc.bfm` dominates runtime; propose speedup options. Read-only investigation, no code changes.
**Module:** mesh (BfmFile, ConnectivityCalculator, ProtoMesh)

## Profile Facts

Sampling profile, 1,782 hits total for the whole run (debug build: `-Og -g`, **no `-DNDEBUG`** — all `BELFEM_ASSERT` and sanity loops active).

Subtree costs (share of whole run):

| Subtree | Hits | Share |
|---|---|---|
| `BfmFile::load` | 655 | ~37% |
| ├─ `Mesh::finalize` | 469 | ~26% |
| │ ├─ `connect_thin_shells_to_thin_shells` | 198 | ~11% |
| │ │ └─ `populate_element_neighbors` 117 + `extract_thin_shell_mesh` 46 | | |
| │ ├─ `finalize_edges` (→ `connect_edges_to_edges` 97) | 195 | ~11% |
| │ └─ `connect_nodes_to_nodes` | 60 | ~3% |
| ├─ `load_edge_data` (→ `reconstruct_edge_connectivity` 93) | 144 | ~8% |
| └─ raw HDF5 I/O (Dataset ctor etc.) | ~25 | <2% |
| `Kernel::distribute_mesh` (→ `connect_elements_to_elements` 265) | 298 | ~17% |
| `DofManager::set_equation` | 295 | ~17% |

Cross-cutting flat costs: **~22% of all samples inside `std::unordered_map`** internals (`belfem::Map` wraps it; every insert mallocs a node; `_M_rehash_aux` visible), **~13% in malloc/free**, ~3% OpenMP idle-thread noise (`omp_get_num_procs`).

**Key finding (confidence high):** reading the 157 MB HDF5 file is *not* the bottleneck (<2%). The time goes into re-deriving connectivity after load — edge→element, element↔element, thin-shell neighbors — mostly through hash-map lookups with per-insert allocation.

## Speedup Options (ranked)

1. **Release build for production runs** (confidence high, zero code change). `-Og` + active asserts inflate everything, especially libstdc++ hashtable code and the debug-only sanity loops (e.g. `ProtoMesh::reconstruct_edge_connectivity` re-walks all element edges under `#if !defined(NDEBUG)`).

2. **Kill the hash map in `connect_elements_to_elements_sub`** (`cl_Mesh_ConnectivityCalculator.cpp:1096`, confidence high, local/low-risk). The function already builds `tKeys`, sorts and uniques it — then copies it into a `Map<key128_t,index_t>` only to look keys up again. Binary search (`std::lower_bound`) over the retained sorted array replaces the map: no node mallocs, no 128-bit hashing. This function is the payload of `Kernel::distribute_mesh` (~15% of run) and the same pattern appears in `connect_nodes_to_edges` (facet_key2 + Map, `:427-470`).

3. **Persist connectivity in the .bfm file** (confidence medium, structural, biggest algorithmic win). The file stores edge→node topology only; element→edge is rebuilt via `Map` keyed `min(A,B)*N+max(A,B)` (`cl_ProtoMesh.cpp:1651`, ~5% of run). Storing element→edge index arrays makes this a direct read. Same idea for element↔element neighbors and ShellToShell — which would also let `distribute_mesh` skip its recompute. Note: `ConnectivityCalculator::save_connectivity_data` (`:1282`) already exists but has **no caller and no load counterpart** — orphaned infrastructure that could seed this. Cost: format version bump + stale-file handling.

4. **Avoid the temporary mesh in `connect_thin_shells_to_thin_shells`** (`:1225`, confidence medium, ~11% of run). It calls `extract_thin_shell_mesh`, which constructs and *fully finalizes* a second Mesh (recursive `Mesh::finalize` visible in profile) just to obtain neighbor relations, then throws it away. Computing neighbors directly on the shell elements — or persisting ShellToShell per option 3 — removes most of this.

5. **Minor:** expose `reserve()` on `belfem::Map` and pre-size the maps that survive options 2/3 (rehash was ~1.3% + malloc churn); set `OMP_NUM_THREADS=1` when profiling to remove the ~3% idle-thread noise.

Options 2–4 together target roughly 25–30 percentage points of the run (load path + distribute_mesh recompute), on top of whatever the release build recovers.

## Status

**Option 2 implemented (2026-07-17, same session).** Christian added
`find_index_in_unique_cell()` to `cl_Cell.hpp` (`std::lower_bound` over the
sorted-unique array that `unique()` already produces — the key's index IS its
position, so the `Map` was a memoization of information the array already
contains) and converted `connect_elements_to_elements_sub`. Claude then:

- added a debug assert to `find_index_in_unique_cell` (`cl_Cell.hpp:487`):
  `lower_bound` result must be an exact match (`!(aMember < *it)`, needs only
  `operator<`); catches absent members and, in practice, unsorted cells at O(1)
- fixed a conversion leftover: `tNumKeys = tKeys.size()` after `unique()` in
  `connect_elements_to_elements_sub` (`cl_Mesh_ConnectivityCalculator.cpp:1137`)
  — the deleted map-build loop had been what reset `tNumKeys` from the raw facet
  count to the unique count; without it `tConnectivity` was ~2× oversized
- converted the twin pattern in `connect_facets_to_facets` (2D branch, `:431-467`)
  — NOT `connect_nodes_to_edges` as the analysis above mislabeled it (that
  function is plain counting loops, no map)

Design decisions:
- **No debug-only "is sorted" flag on Cell:** mutations flow through
  `vector_data()` and non-const `operator()` references, so a flag can't observe
  them; and a debug-only member changes `sizeof(Cell)` between build types (ABI
  trap). The point-of-use assert is the right granularity.
- **Option 3 (persist connectivity in .bfm) REJECTED by Christian:** the arrays
  cost too much memory/file size. Speeding up reconstruction (this change) is
  the accepted approach. `save_connectivity_data` stays orphaned.

Codex audit (thread `tmp/ai_exchange/find_index_unique_cell_audit.md`, 2026-07-17):
**clean, no blocking findings, all points confidence high.** Verified: lower_bound
position ≡ old map value at both sites; no reachable missing-key path; the
`tNumKeys` restoration correct with no stale pre-unique uses; helper
const-correctness. Two notes: (1) the helper's `BELFEM_ASSERT` compiles out in
release, weaker than `Map::operator()`'s always-active `BELFEM_ERROR` — fine
here since membership is proven by construction, but if the helper is ever used
on externally supplied keys, prefer an always-active check or sentinel return;
(2) labeling correction — the converted `connect_facets_to_facets` block is the
`else` (non-2D) branch; the `dim==2` branch is the unchanged bitset path above it.
Build/run validation by Christian pending.

## Round 2 (same day)

Re-profile after round 1 (1,782 → 1,691 samples): hashtable flat cost 22% → 13%,
`connect_elements_to_elements` −29%, `lower_bound` replacement costs 2.5%.
`BfmFile::load` unchanged as expected (its map users were untouched). Christian
approved converting the remaining functions in the same fashion:

- `ProtoMesh::reconstruct_edge_connectivity` (60 hits) — edge pointers sorted by
  node-pair key (`tEdgeKey` lambda) + parallel key array. CORRECTED same day:
  node-pair keys are NOT unique on cut-enriched meshes (see below).
- `Mesh::populate_element_neighbors` (52 hits) — graph vertices created in
  sorted-unique-key order in a `Cell`, indexed via `find_index_in_unique_cell`;
  map removed, cleanup iterates the Cell.

Deliberately NOT converted: ProtoMesh id→pointer member maps (queried from
BfmFile/Distributor/PeriodicityFactory), `Mesh::create_maps`/`create_edge_map`,
DofManager maps — long-lived id maps, different pattern, wider surface.

Both files pass `-fsyntax-only` with production flags.

Codex round-2 audit (same exchange thread): **clean, no blocking findings.**
Verified: identical `Edge*` selection vs the old map; edge-key uniqueness holds
for BELFEM-created edges incl. cut/TSF (EdgeFactory dedups keys, node indices
unique); no map-iteration-order dependency in `populate_element_neighbors`
(map order only ever fed the cleanup loop); lambda captures/lifetimes sound;
`Cell` copy is pointer-shallow as intended. One suggestion adopted: debug-only
adjacent-duplicate-key assert after filling `tKeys` in
`reconstruct_edge_connectivity`, so a corrupted/foreign BFM with duplicate edge
records fails loudly instead of resolving to an arbitrary edge (old map was
silent last-write-wins there too). Re-passed `-fsyntax-only`.
Christian build/re-profile pending.

## Round 2 correction: node-pair keys are NOT unique on enriched meshes

The new assert fired on a genuine corc.bfm reload:
`Duplicate edge node pair ( edges 3057010 and 3110026 )`. Not overflow
(`key_t` = u64, N = 267,368 → N² ≈ 7.2e10) and not corruption: **twin edges
along cut boundary curves share both end nodes** when neither node was
duplicated by the cut. "Unique by construction" only holds for freshly computed
edges (EdgeFactory dedup); a SAVED cut-enriched mesh legitimately contains
duplicate node-pair edge records. The old `Map` silently resolved these
last-write-wins in container order — and reload was verified under that
semantics (D13–D25).

Fix in `reconstruct_edge_connectivity`: `std::stable_sort` (equal keys keep
container order), duplicate-key assert removed, lookup via local `tFindEdge`
lambda using `upper_bound(key) − 1` = last edge of the equal-key range = exactly
the old map's winner. The other three converted functions are unaffected —
their key arrays are `unique()`-ed and keys map to positions, not entities.
Codex correction-audit in the exchange thread.

**Christian validated 2026-07-17: corc.bfm reload works with the fix.**

Codex correction-audit (retry after an OpenAI capacity error): **clean, high
confidence** — `upper_bound − 1` on the stable-sorted array provably returns the
old map's last-write-wins winner; `tIndex > 0` guard sound incl. empty range;
no entity-keyed assumption leaked into the other three conversions. One caveat
adopted: the lookup guard upgraded `BELFEM_ASSERT` → `BELFEM_ERROR` so a
corrupted/foreign BFM fails cleanly in release too (file-I/O validation
category), matching the old `Map::operator()` release behavior.
