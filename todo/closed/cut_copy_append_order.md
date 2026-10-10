# Cut copies appended in ascending copy id, and the cut pipeline swept for hash-ordered loops

**Date:** 2026-10-09
**Purpose:** The copies that `CutProcessor` makes of the nodes on a thin cut get their ids in a sorted walk but are appended to the mesh in the iteration order of each cut set's hash map, so their positions, and every edge id and dof number built from positions, differ between platforms. This commit appends them in ascending copy id (the abstract nodes first, as now), and sweeps the rest of the cut pipeline so that the numbering promise in the docs is as wide as the code supports and no wider. Third commit of the DR-198 numbering line, after cfaf7e8 (interface copies) and the revert 99f8a8f.
**Module:** `src/homology` (`cl_CutProcessor`), `tests/homology`, `src/homology/doc`, `CHANGELOG.md`
**AIs involved:** Claude (sweep, plan, implementation, gates), Codex `gpt-6-astra` + Grok `grok-4.7` (plan jury, code jury); the Server session (blind delta audit and the pipeline on the pushed SHA); the Yufan session (all decks compared on the final SHA)
**Status:** ✅ COMPLETE 2026-10-09. On devel: `collect_duplicates` sorts the kept copies by id after the abstract nodes (`std::sort` on the suffix, debug-only ascending assert), the `mNodes` sort comment names what it decides; `tests/homology/test_CutProcessor.cpp` (the first bare `CutProcessor` test; red on 99f8a8f on the ordering assertion alone, green with the fix); guide and CHANGELOG (luna-swept). Gates: debug tree zero warnings, `test_homology` 58/58, `make check` 17/17; release tree, 21/21 decks through setup on the baseline and the fixed binary, dof banners identical, node order changed on the 19 decks with cuts and kept on the two disk decks. Juries: plan astra/high + grok-4.7/high (no P0/P1); code grok-4.7/high, then astra/high after the usage-limit reset (approve, two P3 on the test, folded in). Finding: `MaxwellFactory` renumbers the enriched mesh with `symrcm` after `create_cuts`, with `std::sort` on equal degrees, so the numbering is reproducible per toolchain and not platform-independent until the graph sort gets a tie-break (reported, `src/math/graph`, not done). Residue: that tie-break; the dead members `CutData::mCutCaseMap` and `CutFactory::mIndicesOfOriginalTerminalNodes`.

> **Scope guards (from the Server's brief):**
> - The fix is `CutProcessor::collect_duplicates` only. Anything else the sweep finds is classified in §2 as order-insensitive with its reason, or fixed here if it is in `src/homology`; a finding outside `src/homology` is reported, not widened into.
> - The abstract nodes stay first, before every cut copy.
> - "Ascending copy id" within a cut set is the allocation order of `CutSet::create_duplicates`: the non-periodic originals in ascending original id, then the periodic pairs in the order that id-sorted walk first meets them (`cl_CutSet.cpp:66-85`, `:100-162`; the lower partner is met first when it is in the list, which the plan jury noted is not unconditional). Across sets the order is the set order, the sorted hex of the pattern (`mCutSets` is an `OrderedMap`, `cl_CutProcessor.hpp:70`).
> - The `.bfm` caveat stays: a reused file replays the numbering it was written with.
> - The docs say what O3 measures about the decks, not more (plan jury, Codex 2).
> - devel only (tip 99f8a8f); never main. No register row.
> - All line numbers below are those of the private worktree at origin/devel 99f8a8f (the shared checkout carries a peer's unpushed commits and its guide sits at other lines).

---

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| Cut-copy ids are deterministic | `CutSet::create_duplicates` walks `mNodeOriginals`, which `CutProcessor::collect_nodes` sorted by id (`cl_CutProcessor.cpp:700`), and allocates `new Node( ++aMaxNodeID, … )` in that walk; non-periodic originals first (`cl_CutSet.cpp:73-85`), periodic pairs after (`:100-162`); `duplicate_nodes` runs the sets in `mCutSets` order, the sorted hex (`cl_CutProcessor.cpp:1087-1093`) | source trace |
| Cut-copy positions follow hash order | `collect_duplicates` counts and then appends each set's copies by iterating `duplicate_map()`, a `belfem::Map` = `std::unordered_map` (`cl_CutSet.hpp:40`; `cl_CutProcessor.cpp:1196-1204`, `:1214-1228`), then `append( mMesh->nodes(), tNodes )` (`:1230`). libstdc++ iterates by bucket, so the position a copy gets depends on the hash table. The abstract nodes go first (`:1209-1212`) | source trace; O2 reproduces it (three inversions on the two-hole strip) |
| Consequence | the interface copies come after the cut copies (`cl_CutFactory.cpp:431-443`), so a model with cuts has platform-dependent positions for every copy; edge ids and the dof numbering follow positions | source trace |

**Bottom line:** sort the kept copies by id after the abstract nodes; the appended order is then the allocation order, a function of the node order that enters the cut processor.

## 2. The Sweep: every container walk in the cut pipeline that can reach mesh state

Mesh state means: the order of `mMesh->nodes()`, `edges()`, `facets()`, `sidesets()` and `curves()`, an id handed out, an element's node slot, a periodic pair, a backup entry, a sideset's facet order. Walks that only delete, count, look up or set a flag do not reach it.

| # | Site | Container | Reaches mesh state? | Verdict |
|---|---|---|---|---|
| 1 | `CutProcessor::check_node_bitsets` fills `mHashToBitsets` (`cl_CutProcessor.cpp:1029-1082`); `determine_cut_sets` walks it assigning `set_index( tCount++ )` and pushing `to_hex()` into `mPatterns` (`:778-789`) | `Map< size_t, Cell< DynamicBitset* > >` (`.hpp:88`) | The pattern index is consumed only by `pattern_index` (`:793-806`), whose result is used as `mPatterns( index )` to get the hex, which keys `mCutSets` (`:889`, `:930`, `:973`, `:1020`); `create_cut_sets` inserts `mPatterns` into the `OrderedMap` (`:824-833`), whose order is the sorted hex regardless of insertion order | insensitive |
| 2 | `CutProcessor::~CutProcessor` walks `mCutSets` and `mHashToBitsets` (`:105-125`) | both | deletes only | insensitive |
| 3 | `collect_nodes` (`:667-706`) | `mMesh->nodes()` with flags, then `sort` by id (`:700`) | `mNodes` and the node indices it sets are id-sorted | insensitive |
| 4 | `collect_elements` (`:233-278`) | `unique` + `sort` by id (`:269-271`) | `mElements` id-sorted; the element indices set from it | insensitive |
| 5 | `collect_edges` (`:143-187`), `compute_edge_bitsets` (`:194-226`) | `mMesh->edges()` with flags; `mCutData` in cut order | mesh order | insensitive |
| 6 | `collect_facets` (`:284-600`) | `mMesh->edges()` / `faces()` with flags, `mCutData` in cut order, sidesets by id list | the thin-cut edge/face lists are mesh-ordered subsets | insensitive |
| 7 | `determine_cut_sets`, `compute_node_bitsets`, `relink_elements` (`:713-770`, `:602-660`, `:1098-1160`) | `mElements` (id-sorted) | bitset flips and `insert_node` per element; `relink_element` looks the copy up by hex and original (`:1165-1178`) | insensitive |
| 8 | `duplicate_nodes` (`:1087-1093`) → `CutSet::create_duplicates` | `mCutSets` (`OrderedMap`), `mNodeOriginals` (id-sorted) | ids allocated in set order then walk order; `add_node_pair_to_backup` called in that order (`cl_CutSet.cpp:143`) | insensitive |
| 9 | **`collect_duplicates` (`:1184-1231`)** | **`duplicate_map()` per set, a `Map`** | **appends copies in hash order** | **FIXED (R1)** |
| 10 | `create_thin_cut_sidesets` (`:837-856`) → `CutData::add_thin_cut_sidesets_to_mesh` (`cl_CutData.cpp:55-180`) | `mCutData` in cut order; `mThinCutEdges` / `mThinCutFaces` (mesh-ordered subsets from #6) | sideset ids `++aMaxSideSetID` per cut, facet ids `++aMaxElementID` per thin-cut edge/face, in that order | insensitive |
| 11 | `CutData::collect_edges`, `collect_coefficients` (`cl_CutData.cpp:185-280`) | `OrderedMap< index_t, int >` of the cochain (`std::map`) | sorted by edge index | insensitive |
| 12 | `CutData::collect_elements` (`:283-375`) | `mMesh->elements()` with flags | mesh order | insensitive |
| 13 | `CutData::mCutCaseMap` (`cl_CutData.hpp:52`) | `Map< id_t, int >` | declared, never referenced | dead member; left alone (out of the brief) |
| 14 | `CutFactory::~CutFactory` walks `mIndicesOfOriginalTerminalNodes` (`cl_CutFactory.cpp:72`) | `Map< id_t, Vector< index_t > * >` (`.hpp:80`) | deletes only; never filled anywhere in `src/` | insensitive (dead member) |
| 15 | `CutFactory::link_node_duplicates_and_originals` (`:2387-2473`) | `mMesh->nodes()` with flags | the per-original duplicate containers are filled walking the nodes in mesh order, deterministic after R1 | insensitive once R1 lands |
| 16 | `CutFactory::collect_orphan_nodes` (`:2476-2507`) | `mMesh->nodes()` with flags | mesh order | insensitive |
| 17 | `CutFactory::duplicate_nodes_on_face_sidesets` (`:1181-1402`), the thin-shell path | `DynamicBitset::where` → ascending node index | thin-shell copies allocated and later appended (`append_move`, `:1167`) in ascending index of the original; periodic backup in that order (`:1352-1378`) | insensitive |
| 18 | `relink_slave_elements_with_duplicate_nodes`, `duplicate_and_relink_facets` (`:1406-1545`) | sideset facet lists, `DynamicBitset::where` | facet ids `++tFacetIdA/B` in facet order | insensitive |
| 19 | `close_terminal_loops` / `close_terminal_loop` (`:1575-1700`), `create_side_curves_for_thinshells_3d` (`:2125-2157`) | `mProtoshells`, curve segment lists | segment ids and curve labels in shell and segment order | insensitive |
| 20 | `CutFactory::save_edges` (`:524-591`) | `OrderedMap`, a local `Map< id_t, Node* >` used for lookup | debug VTK only, deleted | insensitive |
| 21 | `Topology::run`, `update_block_map`, destructor (`cl_Topology.cpp:31`, `:100-119`, `:546`) | `Map< DomainType, Vector< id_t >* >` | the broadcast packs groups in hash order, but each group's `Vector` keeps its build order and is read by type tag, never by position in the buffer; `mBlockTypes` is a lookup map | insensitive |
| 22 | `PeriodicityFactory::update_periodicity` (`cl_Mesh_PeriodicityFactory.cpp:367-430`) | the backup `Cell` of pairs | `restore_node_pairs` rebuilds `mMasterNodes/mSlaveNodes` in backup order (`cl_Mesh_Periodicity.cpp:277-301`), which #8 and #17 fill in sorted order | insensitive |
| 23 | `PeriodicityFactory::map_facets` (`:498-591`) | two local `Map< id_t, Facet* >` | lookup only; the facet pairing is sorted by kd-tree index (`:576`) | insensitive |
| 24 | `PeriodicityFactory::match_facets_and_faces` (`:1321-1381`) | walks `tSources`, a `Map< key128_t, Facet* >` (`:1341`) | sets `A->set_periodic( B )` per pair: each pair independent, the result the same in any order; the face lists are then built walking `aSourceFacets` (`:1356`), a `Cell` | insensitive |
| 25 | `InterfaceProcessor::add_duplicate_nodes_to_mesh` (`cl_InterfaceProcessor.cpp:791-867`), the interface append | the sets' hash maps, then `sort( tDuplicates, tByID )` (`:826-827`) | ascending copy id since cfaf7e8 | insensitive |
| 26 | `CutProcessorManual` (`cl_CutProcessorManual.cpp:65-74`, `:151-174`) appends in walk order | — | never constructed: `grep -rn CutProcessorManual src/ tests/` finds only its own files and an `#include` at `cl_MaxwellFactory.cpp:26` | unreachable; not touched |

**Upstream of the sweep, reported and not widened into (plan jury, Codex + Grok):** `CutFactory` renumbers `mMesh->nodes()` before the cohomology, with reverse Cuthill-McKee by default (`cl_CutFactory.cpp:139-146`, `:178-201`; `symrcm` sorts neighbors by degree, `fn_Graph_symrcm.cpp:126`, with `std::sort`, so equal degrees have no tie-break) or the Poisson ordering (`:2632-2717`, sorted on a solved scalar). That order is what "enters the cut processor"; the docs state the contract relative to it. `src/math/graph` is outside `src/homology`.

What this sweep does **not** cover, and the docs must not claim: the renumbering above, the `.bfm` replay, anything downstream of `CutFactory` (the dof manager's own ordering), and partitioned runs.

## 3. Mechanism

### R1 — `collect_duplicates` appends in ascending copy id

Keep the structure (`:1184-1231`): flag the φ-block nodes, count the kept copies, abstract nodes first, the unkept copies deleted and the maps cleared as before. Then `std::sort( tNodes.begin() + mAbstractNodes.size(), tNodes.end(), tByID )` with `tByID = []( const Node * aA, const Node * aB ) { return aA->id() < aB->id(); }` (the interface processor's comparator, `cl_InterfaceProcessor.cpp:826`); the abstract-node prefix is excluded by the slice, not by the id inequality. A debug-only loop (`#if !defined( NDEBUG ) || defined( DEBUG )`) asserts strictly ascending ids over the suffix. No scratch container. The comment at the `mNodes` sort (`:699`, "not neccessary") now says that the sort decides the allocation order. `belfem::sort( Cell, C&, n )` is not used: it sorts a prefix, and a temporary lambda does not bind to `C&` (plan jury, Grok).

### R2 — the test

`tests/homology/test_CutProcessor.cpp`: the punched strip of `test_CutData.cpp` with k = 2 (builder copied; both sit in anonymous namespaces), coreduced, `Homology` → `Cohomology::updatekGeneratorsFromHomology` → `clean()`, then a bare `CutProcessor( mesh, cohomology, phiBlocks = {1}, nonPhi = {}, phiInterfaces = {}, phiBoundaries = {}, phiPeriodic = {} )`. Oracle on the container order and `Node::id()` (never `index()`, the mesh is finalized before the processor's trailing `finalize()`): originals keep their slots; the abstract nodes come next in cut order with ids `max+1`, `max+2`; every later node is a hanging copy on exactly one original plus the abstract nodes of its pattern, in strictly ascending id; at least two patterns; at least one pattern with two or more copies (a one-entry map iterates in order whatever its buckets).

**Measured:** against the unfixed code the test fails on the ordering assertion alone (copies 42, 41 / 45, 44 / 46, 43 at slots 40-45); the fixture, pattern and multiplicity assertions pass there, so the fixture risk of the plan round did not materialize.

### R3 — docs

- Guide item 4 (`homology_usage_guide.md:1627-1630`): the cut copies are appended in ascending copy id by `collect_duplicates`; the contract is relative to the node order entering the cut processor (after the RCM renumbering); `.bfm`, dof manager and partitioned runs excluded.
- Guide "Cut Processing / Internal Workflow" step 4 (`:1319`) and the `CutSet` algorithm (`:1708-1711`): allocation order; the CutSet text does not say "append" (it never touches the mesh), it points at `collect_duplicates`.
- `CHANGELOG.md` Changed (`:57-59`): the sentence replaced; the deck statement written from O3's count.
- Luna sweep (`gpt-5.6-luna` / medium) on the guide and CHANGELOG text.

### R4 — gates

- O1: debug tree, build zero warnings, `test_homology` 58/58, `make check`.
- O2: the new test against the unfixed code: red on the ordering assertion (done, above).
- O3: release tree, the 21 decks through setup with no `.bfm`, one deck at a time, on the baseline binary (99f8a8f) and on the fixed binary; per deck the dof banner compared (`$SP/decks/dof_diff.sh`, a scratchpad script of this session, not in the repo) and the node order read from each run's `.bfm` (`/nodes/ids`) compared, to count the decks whose node order changed. Not before 19:00: the machine is shared with the Yufan session's reference runs.
- O4: the Yufan session compares the decks on the pushed SHA (their call).

## 4. Steps

- [x] plan jury (astra/high + grok-4.7/high): no P0/P1; the sort API, the periodic wording, the oracle strength, the RCM caveat, the comparator naming and the `:699` comment folded in
- [x] R1 `collect_duplicates` sort + comment
- [x] R2 `test_CutProcessor.cpp` + `tests/homology/CMakeLists.txt`
- [x] O2 red on the unfixed code, O1 build + `test_homology` green
- [x] O1 `make check` 17/17
- [x] code jury on the own-files diff: Grok (text, assert message, macro corrections) and, after the usage-limit reset, Codex astra (approve, two P3 on the test), all folded in
- [x] O3: 21/21 setup-ok on both binaries, dof banners identical, node order changed on 19 decks (all with cuts), identical on the two disk decks (no cuts), id sets identical
- [x] R3 docs + luna sweep
- [x] devlog, this file closed, push devel, SHA to the Server; O4 handed to Yufan

## 5. Open Questions

- O-a: whole-range suffix sort over per-set collection — RESOLVED 2026-10-09 → suffix sort (both auditors).
- O-b: `mCutCaseMap` and `mIndicesOfOriginalTerminalNodes` are dead members (a data member and, for the second, a destructor walk); removing them is outside the brief. Reported, not touched.
- O-c: the RCM tie order (Codex 1, Grok). O3 showed it is decisive: `MaxwellFactory` renumbers the enriched mesh with `symrcm` after `create_cuts` (`cl_MaxwellFactory.cpp:675`, `:725`), from container positions, with `std::sort` on equal degrees; helix's first changed slot is 715 of 37521 originals. The commit removes every hash-map order from the numbering; platform independence needs an id tie-break in `symrcm` (`src/math/graph`), reported to Christian, not done here. The docs say so.
