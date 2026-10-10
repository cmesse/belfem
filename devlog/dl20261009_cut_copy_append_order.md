# Cut copies appended in ascending copy id, and the cut pipeline swept for hash-ordered loops

**Date:** 2026-10-09
**Purpose:** Third commit of the DR-198 numbering line. `CutProcessor::collect_duplicates` appended
each cut set's copies in the iteration order of a hash map, so their positions in the mesh, and
the edge ids and dof numbers built from positions, depended on the platform. The copies are now
sorted by id after the abstract nodes, which is their allocation order. The rest of the cut
pipeline was read for other container walks that reach mesh state; every one is order-insensitive
for a stated reason. The deck gate then showed that the reverse Cuthill-McKee renumbering after the
enrichment, not any hash map, decides the final node order, so the docs promise no more than "no
hash-map iteration order reaches any node position".
**Module:** `src/homology` (`cl_CutProcessor`), `tests/homology`, `src/homology/doc`, `CHANGELOG.md`
**AIs:** Claude Fable 5.1 (DrClean session: sweep, plan, implementation, gates); Codex `gpt-6-astra`
and Grok `grok-4.7`, both at high, on the plan; Grok alone on the code round, since Codex hit its
usage limit mid-round and was rerun after the reset (approve, two P3 on the test, folded in); Codex `gpt-5.6-luna` on the four doc paragraphs; the Server session (Christian's relay); the Yufan session
(deck comparison on the pushed SHA)
**Status:** COMPLETE — on devel, rebased onto 7011f44 (DR-177 index keys) with the debug gate rerun there; gates: debug tree zero warnings, `test_homology` 58/58 (the new test red on the unfixed code, green with the fix), `make check` 17/17; release tree, 21 of 21 decks through setup on the baseline and on the fixed binary, dof banners identical, node order changed on the 19 decks with cuts. Residue for Christian: an id tie-break in `symrcm` (`src/math/graph`) if platform-independent numbering is wanted; two dead members in `CutData` and `CutFactory`

## Summary

- **The defect.** `CutSet::create_duplicates` hands out copy ids walking `mNodeOriginals`, which
  `CutProcessor::collect_nodes` sorts by id, set after set in the sorted order of the pattern hex
  strings; the ids were already deterministic. `collect_duplicates` then appended each set's copies
  by iterating `duplicate_map()`, a `belfem::Map` (`std::unordered_map`), in bucket order. On the
  two-hole punched strip of the new test that order is 42, 41, …, 45, 44, …, 46, 43 (three
  inversions, measured on the unfixed code).
- **The fix.** After the abstract nodes, `std::sort` on the suffix of the collected `Cell` by
  `Node::id()`, with the comparator the interface processor uses; a debug-only loop asserts the
  suffix is strictly ascending. No scratch container, no change to who is kept or deleted.
- **The sweep.** Twenty-six sites in `CutProcessor`, `CutSet`, `CutData`, `CutFactory`,
  `Topology`, `InterfaceProcessor` and the periodicity rebuild, each classified in the plan
  (`todo/closed/cut_copy_append_order.md` §2): the hash-to-bitset table reaches the mesh only
  through the hex key of an ordered map; every other walk is over an id-sorted `Cell`, a mesh
  container with flags, an ordered map, a `DynamicBitset::where` list, or does lookups and deletes
  only. Two dead members (`CutData::mCutCaseMap`, `CutFactory::mIndicesOfOriginalTerminalNodes`)
  and the never-constructed `CutProcessorManual` are reported, not touched.
- **The boundary.** `CutFactory` renumbers the mesh nodes with reverse Cuthill-McKee before the
  cohomology (`compute_rcm_problem`; `symrcm` sorts neighbors by degree with `std::sort`, no
  tie-break for equal degrees), or with the Poisson ordering. Both plan auditors raised it. The
  docs state the contract relative to the node order that enters the cut processor; whether to
  give `symrcm` an id tie-break is Christian's call, in `src/math/graph`.
- **The test.** `tests/homology/test_CutProcessor.cpp`, the first test to construct a bare
  `CutProcessor`: the punched strip with two holes, coreduced and cleaned as the k-hole tests do,
  then the processor with block 1 as the only φ block. Oracle on container order and ids:
  originals in place, abstract nodes next in cut order, every later node a hanging copy on one
  original plus the abstract nodes of its pattern, strictly ascending; at least two patterns; at
  least one pattern with two or more copies, since a one-entry hash map iterates in order whatever
  its buckets.
- **Gates.** Debug tree (SCLS debug flavor, gcc, Armadillo): the new test red on the unfixed
  99f8a8f code on the ordering assertion alone, green with the fix; build zero warnings,
  `test_homology` 58/58, `make check` 17/17. Release tree
  (SCLS `gcc`): the 21 example decks through setup on the baseline 99f8a8f binary and on the fixed
  binary, one deck at a time, no `.bfm`: 21 of 21 setup-ok on both, dof banners identical on all
  21; the node order read from each run's `.bfm` (`/nodes/ids`) changed on 19 decks and is
  identical on the two disk decks, which have no cuts; the id sets are identical everywhere.
- **Docs.** The guide's InterfaceProcessor item, the Cut Processing workflow step 4 and the CutSet
  algorithm, and the CHANGELOG Changed note describe the allocation and append order, name
  `collect_duplicates`, state the contract relative to the order entering the cut processor, and
  keep the `.bfm`, dof-manager and partitioned-run exclusions. Luna sweep applied.

## Key Findings

- The plan round (both seats) corrected the recipe: `belfem::sort( Cell, C&, n )` sorts a prefix
  and takes an lvalue comparator, so the suffix sort is `std::sort` on the container's iterators;
  the oracle needed the two-copies-in-one-set requirement; the periodic allocation order is "the
  order the id-sorted walk meets the pair", not unconditionally the lower original first; "every
  deck with cuts differs" is a measurement, not a consequence of the sort. Grok's doc line numbers
  were the shared checkout's (it carries a peer's unpushed commits); the plan cites the worktree at
  99f8a8f. Grok's claim that `CutProcessorManual` may be reachable is refuted: it is never
  constructed.
- The code round (Grok) found comment wording, the assert message (an inversion was reported as a
  shared id), the `BELFEM_ASSERTIONS_ACTIVE` macro instead of a restated `#if`, and two test
  comments; all folded in. Its fixture-abort risk did not materialize: the test runs to its
  assertions on both builds.
- The `.bfm` holds the node order as `/nodes/ids`, which made the per-deck node-order comparison a
  file diff between a baseline and a fixed setup run.
- **The numbering is not platform-independent after this commit, and the docs say so.** On
  helix the first changed position is slot 715 of 37521 original nodes: originals are permuted,
  not only copies. `MaxwellFactory` runs `mesh::symrcm( mMesh->nodes() )` after `create_cuts`
  (`cl_MaxwellFactory.cpp:675`, `:725`) on the enriched mesh, before the `.bfm` is saved;
  `fn_Mesh_symrcm_nodes.cpp` builds each node's adjacency from container positions, copies
  included, and `graph::symrcm` sorts neighbors by degree with `std::sort` and no tie-break
  (`fn_Graph_symrcm.cpp:126`, `op_Graph_Vertex_Degree.hpp:27`). So the final order of every
  deck with cuts is that RCM's output on an input this commit changed, and across toolchains it
  still depends on how the standard library's sort leaves equal degrees. What this commit
  removes is every hash-map iteration order from the numbering, which was the brief. An id
  tie-break in `symrcm` (`src/math/graph`, outside the brief) would close the gap; reported to
  Christian through the Server, not done here.

## Changes Made

- `src/homology/cl_CutProcessor.cpp` — `collect_duplicates` suffix sort and assert; the comment at
  the `mNodes` sort
- `tests/homology/test_CutProcessor.cpp` (new), `tests/homology/CMakeLists.txt`
- `src/homology/doc/homology_usage_guide.md`, `CHANGELOG.md`
- `todo/closed/cut_copy_append_order.md`, `todo/README.md`, this devlog, `devlog/README.md`
