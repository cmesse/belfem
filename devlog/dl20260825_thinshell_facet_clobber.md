# Thin-shell facet clobbering by update_facet_nodes (gantry crash)

**Date:** 2026-08-25
**Topic:** `hphirun` abort on the large 2D gantry deck (464 tapes, generalized-Pellikka cuts):
`DynamicBitset` assert "Index 53908 out of range (expect < 9744)" in
`ThinShellFactory::collect_nodes` (`cl_ThinShellFactory.cpp:815`).

## Diagnosis

Probe on the failing facet (tape sideset 215, facet 2641, at the tape tip):

```
#Probe f: 2641 k: 0 n: 53852 nidx: 53908 selforig: 1 oid: 53852 oidx: 53908 hang: 1 nsrc: 209
#Probe    src 0..208 : ids 51231-51439, all EntityType::NODE
```

Node 53852 is an unlinked hanging node from the thin-cut stage (self-original, 209 sources —
the slave duplicates of one whole tape; the MMF/jump carrier of the cut that terminates on
tape 215). Failure chain:

1. `CutFactory::create_thin_shell_cuts()` freezes the tape master-node list (9744 nodes) and
   stashes the pristine tape facets.
2. The thin-cut stage relinks the jump-side air element — master of tape facet 2641 —
   replacing real node 2427 with node 53852. Correct for the φ assembly.
3. `restore_thin_shell_sidesets()` restores the pristine facets, but `CutFactory::run()`
   ends with `unfinalize()/finalize()`, and `Mesh::finalize()` calls `update_facet_nodes()`
   (`cl_Mesh.cpp:790`), which re-derives every facet's nodes from its master element —
   silently copying the unlinked node into the restored tape facet.
4. Such nodes never get an `original()` link (`link_node_duplicates_and_originals`,
   `cl_CutFactory.cpp:2588-2592`), so `collect_nodes` sees a mesh-wide index and the
   debug-build assert fires. In a release build this was a raw out-of-bounds heap write.

Smaller decks survive because no cut terminus lands on a tape facet whose master is on the
jump side.

## Fix (audited: Codex + Grok plan round, Codex + Grok code round)

- `Mesh::update_facet_nodes()` (`cl_Mesh.cpp`): skip sidesets typed `DomainType::ThinShell`
  or `DomainType::GeometryOnly`. Tape facets and the ThinShellFactory's aggregate sideset
  are extrusion geometry; their node lists are restored deliberately by the cut factory and
  must not be re-derived from (relinked) masters. `GeometryOnly` is assigned to a mesh
  sideset in exactly one place (`cl_ThinShellFactory.cpp:79`), cannot come from the input
  parser, and covers the re-finalize inside `create_thinshells()`
  (`cl_MaxwellFactory.cpp:1172-1175`). The skip is also load-bearing for MPI: the
  distributor ships rank-0 facet node ids verbatim and worker finalize never runs
  `update_facet_nodes`, so rank 0's lists are what every rank gets.
  Ghost / Periodic / Cut / EnrichedInterface sidesets are deliberately NOT skipped.
- `ThinShellFactory::collect_nodes()`: `BELFEM_ERROR` (setup path, once per shell) naming
  facet id, node id, original index and master count, placed before the bitset write so a
  release build aborts loudly instead of corrupting the heap.

Audit trail: `tmp/ai_exchange/thinshell_abstract_node_facet.md` (plan, two audit rounds,
reconciliation). Codex found the post-create GeometryOnly gap in the original ThinShell-only
skip; Grok confirmed it independently and supplied the MPI argument.

## Residuals (follow-up material, not fixed here)

- `.bfm` facet save carries no node topology (`save_facet_data(false)`); reload reconstructs
  facet nodes from masters (`cl_ProtoMesh.cpp:1007-1020`). A save/reload of a realized
  thin-shell + cut mesh therefore re-derives GeometryOnly facet node identity from
  cut-relinked masters regardless of this fix. Needs its own task before the round-trip is
  relied on for such meshes.
- The exact pedigree of the injected node (origin-placed abstract node vs self-original cut
  duplicate created at the terminus) is unresolved; the fix does not depend on it.

Also repaired in this session: `.claude/scripts/ask_grok.sh` — `run_grok()` re-enabled
errexit before returning a nonzero grok exit, killing the wrapper at attempt 1 with no
retries and no error output; transient grok API failures looked like silent script deaths.

## Gates

- [x] User rebuild + gantry rerun past `create_thinshells` (the original reproducer) —
      passed same evening; the run reached BDF1 step 6 with clean convergence
- [x] `make check-fast` — **9/9 green** (2026-08-26, `cmake-build-claude`, the private
      build tree Christian approved so gates stop queueing behind him). Includes the
      homology suite, which is also the gate DR-107's cleanup owed.
- [x] **Fresh-build gantry at np=10 — PASSED** (2026-08-26). `gantry.bfm` deleted so the
      mesh is built from `gantry.msh` and the distributed path is genuinely exercised.
      **Zero "Could not find master" aborts** where previously all ten ranks died, so the
      connectivity skip is no longer unexecuted: it ran, on worker meshes, and worked.
      The run then built its thin shells, distributed, and integrated six BDF1 steps with
      clean convergence (−80.70 dB at step 5, −71.75 dB at step 6) before terminating in
      `Postprocessor::synch_node_indices` on `MPI_ERR_TRUNCATE` from ranks 3 and 7 at tags
      60/140 — that is DR-100, a pre-existing and unrelated defect, not a regression of
      this fix. Commit c6ef67ef's message declares both of these gates owed; they are
      cleared here rather than by amending the (unpushed but already-referenced) commit.

## Follow-up same evening: distribution regression and its fix

The gantry rerun surfaced a real regression of C1 on fresh-built parallel meshes:
`ConnectivityCalculator::connect_facets_to_elements()` re-derives every facet's master
by node-identity matching, ignoring masters already set — and post-C1 the restored tape
facet nodes legitimately diverge from cut-relinked master elements, so every worker rank
aborts with "Could not find master for facet N" during `Kernel::distribute_mesh`
(deterministic; masked whenever the mesh loads from `gantry.bfm`, whose loader rebuilds
facet nodes from masters). Both plan-audit rounds had missed this consumer.

Christian's ruling: the ownership-inheritance design (facet and slave inherit from the
master) is correct and the calculator should honor it. Fix (own plan + Codex/Grok rounds,
exchange `connect_facets_trust_masters.md`): facets with `has_master()` are trusted —
the skip wires the facet-wrapper element container (the calculator is the only place
that ever fills it) and deliberately does NOT call `set_master`, whose default node
relink would re-inject the cut duplicates on workers. Fresh gmsh reads still derive
(facets start masterless). Grok additionally established the skip is load-bearing in
release builds: the "no master" abort is an assert, and release would have called
`set_master( nullptr, ... )`.

Gate: fresh-build gantry at np=10 must distribute and run; serial run and check-fast
unchanged.
