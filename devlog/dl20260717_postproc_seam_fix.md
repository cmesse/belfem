# Maxwell Postprocessor: B-Field Seam at Partition Boundaries

**Date:** 2026-07-17
**Purpose:** Diagnose and fix the cosmetic seam in the recovered B-field along MPI partition boundaries
**Module:** fem/kernel (Postprocessor, DofMgr BlockData, Element), consumed by fem/maxwell

## Symptom

The nodal B-field recovered by the Maxwell postprocessors shows a visible
jump along the boundaries of the MPI domain decomposition
(`cmake-build-debug/corc/postproc_bug.png`). Purely cosmetic — the solve is
unaffected — but the aura (ghost layer) should, in theory, make the
recovery partition-independent.

## Root cause (confidence: high)

`Postprocessor::recover_fields()` performs a per-node weighted
least-squares patch fit: `compute_node_matrices()` builds a per-node
Vandermonde normal matrix over the node's element disc, and phase 1 of
`recover_fields()` accumulates the matching right-hand side element by
element. For a node on a partition boundary, the full disc is only
available if the aura elements participate. Three independent gaps
prevented that — the recovery was *consistently one-sided* (both the
Vandermonde and the rhs saw only the owned half-disc), which produces a
valid but side-biased fit, i.e. the seam:

1. **Selection ignored the aura.**
   `Postprocessor::select_elements_and_owned_nodes()` flagged only
   `mField->block(b)->elements()` — the owned FEM elements.
   `aura_elements()` were never flagged, so partition-boundary patches only
   contained the owned side. Primary cause.

2. **The FEM aura was facet-based, recovery needs node-based.**
   `BlockData::collect_element_indices()` expanded the aura via
   `element->element(k)`, which `connect_elements_to_elements()` builds
   from shared-facet keys. A recovery patch needs every element touching
   the node, including vertex-only neighbors. The mesh distributor already
   ghosts the full node-based disc (incl. duplicate-node discs,
   `Mesh_Distributor::select_entities()`), so the elements exist on the
   submesh — they just had no FEM wrapper.

3. **Aura elements had uninitialized edge orientations.**
   The dof-less aura `fem::Element` constructor never called
   `compute_edge_directions()` (only `link_dofs()` does). The Nédélec edge
   functions (`EF_TRI3::link` etc.) read `edge_directions()`, so aura
   elements entering a conductor-block recovery would have used garbage
   sign bits. Latent until gap 1 was fixed.

## Changes

- `src/fem/kernel/cl_FEM_Element.cpp` — aura ctor computes edge directions
  when `mElement->has_edges()`.
- `src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp` — aura selection walks the
  node discs (`node->original()->elements()` plus duplicates' elements) of
  every node of every owned element, still guarded by the candidate flag
  (non-owned elements of selected blocks / sideset masters+slaves only).
  Serial path unchanged.
- `src/fem/kernel/cl_FEM_Postprocessor.cpp` — both branches of
  `select_elements_and_owned_nodes()` iterate `elements()` **and**
  `aura_elements()` with the identical flag/claim body, so aura elements
  join the patches and the credited ⊆ claimed invariant of
  `recover_fields()` phase 1 is preserved.

Everything downstream was verified aura-safe before the change:
`_Volumes` is redistributed to non-owner copies
(`Kernel::compute_element_volumes`), `synch_source_field` covers all
entities derived from `mMyElementIndices`, `nedelec_data_*`/`node_data`
read mesh fields (no dof objects), and the block element maps include aura
wrappers. The only other consumers of `aura_elements()` are the element
maps, so the enlarged aura does not change any assembly path.

## Audit

Codex audit thread: `tmp/ai_exchange/postproc_seam_fix.md` (asks A1–A5:
claimed-node invariant, submesh entity completeness, `has_edges()` safety
across element types, other aura consumers, thermal-`q()` gap).
**Result: no blocking findings; A1–A4 confirmed with line citations
(confidence high).** Notable confirmations: the distributor packs edge/face
IDs for every ghosted element (`cl_Mesh_Distributor.cpp:1188-1258`), so no
aura element can reach the edge path without its edges present; every
factory-created mesh element is an `ElementTemplate` overriding
`has_edges()`, so the base-class BELFEM_ERROR is unreachable from the aura
ctor. A5 (thermal `q()` staleness) confirmed real; Codex recommends the
dof-free `node_data("T")` route before relying on thermal-coupled
superconducting postprocessing, and notes it additionally requires `"T"`
to be added to the postprocessor's synchronized source fields (currently
only `edge_h`/`face_h`).

## Regression found in first parallel run (hphirun) and fixed

First parallel run crashed with SIGSEGV in `EF_PENTA6TS::link`
(cl_EF_PENTA6TS.cpp:42, `aElement->facet()->node(0)` with null facet),
reached from `compute_node_matrices()` linking a **thin-shell aura
element**. Root cause: `BlockData::link_thin_shell_facets_parallel()` sent
the element→facet pairs to each element's *owner rank only*, so aura
thin-shell wrappers always had `mFacet == nullptr`. Latent before this
session because recovery never linked aura elements.

Fix (two parts):

- `link_thin_shell_facets_parallel()` now `share`s the full element→facet
  id table to all procs; every rank links every wrapper it holds (owned
  and aura), skipping pairs whose element/block/wrapper/facet is not on
  the local submesh. Owner-side completeness is still asserted.
- `Postprocessor::select_elements_and_owned_nodes()` skips thin-shell
  elements whose wrapper has no facet (possible when the facet of a
  far-aura shell element was not ghosted): they are excluded from the
  patch instead of crashing. Owned thin-shell elements always have their
  facet, so owned behavior is unchanged.

Design note (Christian's question, whether a calculator may run on a
non-owned element at all): yes for the recovery paths — the postproc
calculator uses only ghosted geometry, fem-wrapper edge orientations, and
mesh field values that `synch_source_fields()` distributes from rank 0's
gathered global fields to exactly the selected entity sets. That synch IS
the "request missing information" step, rank-0-mediated instead of
neighbor-to-neighbor. The known holes where data is *not* synched to aura
elements are (a) thermal T via `q()` (dof objects, see follow-up 1) and
(b) `phi_m`/`phi_s` in `compute_superconductor_ts`, which read phi on the
facet's master/slave *volume* elements whose off-facet nodes are outside
the synch set — pre-existing for owned shell elements with non-owned
masters/slaves too. Proper fix for both: extend the postproc source
fields / synch node sets, not a new neighbor exchange.

## Second regression: "calculator is not allocated" — fixed, and thin
## shells reverted to side-local recovery

Second parallel run hit `BELFEM_ERROR( mIsAllocated, ... )` in
`Calculator::link`. Cause: a rank can hold a block with **only aura
elements** (e.g. it owns air next to a conductor owned entirely by another
rank). `Calculator::allocate()` early-returned when
`number_of_elements() == 0` (owned count), so such a block's calculator
was never allocated — harmless before, fatal once recovery links it.

Fixes in `cl_FEM_Calculator.cpp`:

- `allocate()` proceeds when the group has aura elements; the dof-count
  probe `elements()(0)` falls back to 0 dofs for aura-only groups (the
  recovery paths never touch the dof-sized work arrays).
- Coupled-run guards: the thermal dof manager's aura is expanded from a
  smaller owned-element set than Maxwell's, so an aura-only Maxwell block
  (or a single aura element) may be absent on the thermal side.
  `MaxwellData` construction now requires `block_exists` on both peers;
  `link_element_maxwell_thermal` relinks the thermal peer only if the
  element exists in the thermal group; `compute_superconductor[_ts]`
  falls back to `gTbulk` when the thermal block is absent. All guards are
  no-ops on owned paths (owning elements implies both dof managers have
  the block).

**Thin shells reverted to side-local (owned-only) recovery** — the
facet-null skip is replaced by skipping the whole aura pass for
`DomainType::ThinShell` blocks. Reason: `compute_superconductor_ts` →
`get_normal_calculator()` reads `phi` directly from the facet's master
and slave *volume* elements (cl_FEM_Calculator.cpp:1436-1448), and `phi`
is not in the TS postprocessor's synch set — an aura TS element would
recover from stale phi. One-sided-but-correct beats full-disc-but-stale.
The TS aura path additionally depends on the sideset wrapper/calculator
of the facet existing locally. Making TS recovery partition-independent
needs the phi/T synch extension first (see follow-ups). The full
element→facet share in `link_thin_shell_facets_parallel()` is kept — it
removes the null-facet landmine on aura wrappers regardless.

## Run 3: crash persists — joint Codex+Grok audit reframes it as a
## latent .bfm-reload defect, not a seam-fix regression

Run 3 crashed identically (EF_PENTA6TS::link, null facet) despite the
thin-shell aura skip. Joint audit (thread in
`tmp/ai_exchange/postproc_seam_fix.md`, sections "Codex run-3 analysis"
and "Grok run-3 analysis") established:

- The debug runs are **serial** (`Debug.cmds` launches hphirun without
  mpirun → MPI singleton). In serial every seam-fix change is provably
  inert: aura containers empty, owned selection byte-identical,
  serial facet-linking path untouched (Grok, high confidence).
- File forensics: `corc.msh` regenerated 14:23, seam runs (working
  postproc) until 16:46, `corc.bfm` **minted 17:22 as a checksum cache
  by the first post-fix run**. Every run since silently takes the
  `.bfm` reload path (`MaxwellFactory` checksum match). The crash
  "appearing after the seam fix" is calendar correlation — the trigger
  is the fresh cache file, not the diff.
- The solve never touches block-wrapper facets (thin-shell physics
  assembles on sideset elements), so reload survives the solve and dies
  at the first postprocessor that links a PENTA6TS block element.
- Block 36 (magnesia layer) is `DomainType::Buffer` via
  `create_buffers` reclassification — the Air/phi postprocessor selects
  `Air || Buffer` blocks, which explains the crash appearing in the phi
  postprocessing.
- Static analysis of the reload path finds all facet-linking gates open
  (`/thinshells` present, tape sideset 32 has 34560 facets, blocks
  33-39 live); the exact null-facet mechanism on reload is unproven —
  needs runtime probes.
- Codex side-finding: `collect_thin_shell_facet_ids()` has a real
  count/populate asymmetry (counts selected shells, populates ALL
  shells) — buffer overrun on multi-shell partial selection; corc has
  one shell, so not this crash, but should be fixed.

Discriminating experiment for Christian: rename `corc.bfm` away and
re-run (forces the factory path — crash should vanish); optionally
`git stash` + re-run with the cache present (pre-session code should
crash the same way, proving latency).

## Run 4: cold start still crashes, on proc 2 — root cause found and
## fixed (Buffer-reclassified shell layer escaped the thin-shell guard)

The cold-start run (no .bfm cache) crashed identically **on proc 2** —
so the debug runs were parallel all along (per-rank lldb under mpirun;
`Debug.cmds` runs per rank). The serial verdict and the .bfm timeline
were red herrings for this crash.

With parallel confirmed, Grok's H3 audit finding is the mechanism:
`create_buffers()` (cl_ThinShellFactory.cpp:1732-1736) reclassifies
shell layers without `MaterialProperty::rho` — corc's magnesia layer,
block 36 — to `DomainType::Buffer`. The thin-shell aura skip keyed on
`mesh::Block::domain_type() == ThinShell` therefore did not fire for
block 36: its aura pass ran, and an aura PENTA6TS wrapper whose shell
facet was not ghosted onto proc 2 stayed facet-less (the share-based
linking correctly skips non-local facets), was flagged, entered a
recovery patch, and `EF_PENTA6TS::link` dereferenced null. This one
hole explains all four crash runs consistently (run 1: no TS handling
at all; run 2: facet-linked TS aura wrappers reached a then-unallocated
aura-only calculator; runs 3-4: six ThinShell-typed layers skipped, the
Buffer-typed one not).

**Fix:** the thin-shell test in `select_elements_and_owned_nodes()` now
keys on the **element type** (PENTA6TS / QUAD4TS / HEX8TS) in addition
to the domain type, in both branches — the element type is what
actually requires the facet link. Buffer-reclassified layers recover
side-locally like their ThinShell-typed siblings.

## Run 5: crash gone, residual artifacts at partition-boundary corners
## — aura widened from owned elements to owned nodes

With the crash fixed, the recovered field (seam.png,
seam_with_elements.png, owners.png) showed the large patch offsets
gone, but isolated single-element dimples remained exactly at the
reentrant corners of the ownership boundary. Christian's diagnosis
("aura not big enough") is correct in a precise sense: the FEM aura in
`BlockData::collect_element_indices()` expanded from the node discs of
**owned elements**, but a proc can own a node while owning none of the
elements around it (reentrant boundary corners). Such a node's disc
members only entered the aura if they happened to touch some other
owned element — at sharp notches they don't, leaving the node's
recovery patch incomplete.

Fix: a second expansion pass walks the discs (original + duplicates) of
every **owned node** and adds the flagged candidates. The distributor
already ghosts the complete disc of every owned node
(`Mesh_Distributor::select_entities()`), so all candidates exist on the
submesh. Serial remains a no-op (no flagged candidates). This makes the
FEM aura the exact closure needed by patch recovery: every element of
every owned node's disc has a local FEM wrapper.

## Slave-surface artifacts: T-matrix mirroring hypothesis tested and
## rejected (three-AI audit, 2026-07-19)

Christian observed the residual artifacts only on the periodic slave
surface and hypothesized that mirroring the master's T-matrices onto
the slave (instead of merging/averaging both sides) leaves the slave
"a bit off." Claude analyzed the chain and Codex + Grok audited
independently (thread: `tmp/ai_exchange/postproc_seam_fix.md`,
sections `t-matrix-mirroring`). All three converge, confidence high:

1. **The elimination is not variationally one-sided.** Element
   T-matrices project `TᵀKT` / `Tᵀr` in every assembly path
   (`cl_FEM_Tmatrix.cpp:101-152`, consumers in
   `cl_FEM_DofManager.cpp`), so the retained master dof's equation
   accumulates both discs' physics; slave elements "talk back" via Tᵀ.
2. **All periodic weights in the first-order corc chain are exact
   combinatorial constants.** Pure periodic node/edge ties carry
   weight +1.0 — `PeriodicityFactory::match_edges` swaps anti-aligned
   slave edge endpoints before tying
   (`cl_Mesh_PeriodicityFactory.cpp:1093-1116`) — and interface
   cascades multiply through ±1, 0.5, 1/3, 4/3. No geometric content;
   for a node-exact conforming pairing the constraint u_B = u_A is
   exact, and copying A's flattened chain onto B is the chain rule,
   not an approximation. Grok's four adversarial attacks on the
   Nedelec-sign (`mS`) × periodic-weight composition all failed.
3. **Averaging would violate the constraint.** If both sides'
   descriptions agree it is a no-op; if they disagree it enforces
   u_B ≠ u_A and would create a genuine physical seam. Mortar-style
   merged operators only make sense for non-conforming interfaces.

Residual non-exact classes flagged by the audit (neither supports
averaging): half-cut rim edges deliberately left untied
(`cl_Mesh_PeriodicityFactory.cpp:1131-1137`, closes through node
hanging); silent hanging→free demotion when no matching source dof
type exists (`cl_FEM_DofMgr_DofData.cpp:3538-3554`); higher-order
face/facet periodic weights are NaN-deferred (inert at first order,
`mFaceDofMultiplicity = 0`). The slave-only artifact asymmetry
therefore still points at parallel recovery-set completeness on the
slave face — the owned-node aura widening above targets exactly that
and awaits its first run.

## Root cause of the slave-surface seams: duplicate links never reach
## the workers (orphaned-original selection hole)

Christian probed `recover_fields` on the seam family 26517 (orphaned
original, periodic partner of 23624) / 110088 / 112743:

- serial: `110088 org:26517`, `112743 org:26517`, `23624 per:26517` —
  links intact
- 10 procs: `110088 org:110088`, `112743 org:112743`, no `per:` —
  original/duplicate and periodic links absent on the workers

The distributor does transmit duplicate links (packed keyed on the
**original's** id, `cl_Mesh_Distributor.cpp:1045-1071`; applied in
`ProtoMesh::create_nodes`, `cl_ProtoMesh.cpp:233-258`), but the packing
only emits records for nodes in the target's node table — and the
selection closure set the bits of an original's *duplicates* without
ever setting the original itself. An element-less original (orphaned
seam node, disc carried entirely by its duplicates) is selected by no
other rule, so its family record was never packed and the worker linked
nothing. With `original()` degenerate on the workers, the sibling
crediting in `recover_fields`, `compute_node_matrices`, and the
owned-node aura widening all silently no-op → each seam node recovers
alone from its half-disc. This also explains the slave-surface-only
asymmetry: the cut-duplicate families and orphaned originals live on
the slave surface; master-surface nodes are their own originals with
complete discs and need no links. Serial keeps the links → clean.

Fix: select the original into the node table inside the duplicate
closure (`cl_Mesh_Distributor.cpp:363-378`). The closure is
all-or-nothing per family, so the worker-side `mNodeMap` lookups stay
safe.

Related finding: `ProtoMesh::create_periodicitiy` is only invoked on
the .bfm reload path on rank 0 (`cl_Mesh_BfmFile.cpp:1769`) — the
distributor never rebuilds periodic node links on workers. Currently
benign (solve uses transmitted t-matrices, recovery is side-local),
but worth remembering for any future worker-side code that follows
`Node::periodic()`.

## Verification

- Syntax-only compile of all edited files with the real `flags.make`
  flags (`libbelfem_kernel`, maxwell target): clean.
- Run 1 exposed the thin-shell facet regression (parallel), run 2 the
  aura-only-block allocation gap (parallel); both fixed. Run 3's crash
  is attributed to the latent .bfm-reload defect above, pending the
  discriminating experiment.

## Open follow-ups (not fixed here)

1. **Thermal-coupled superconductor recovery:** `compute_superconductor*`
   fetches T via the thermal peer's `Calculator::q()`, which reads dof
   objects; on dof-less aura elements the loop no-ops and `mq` keeps the
   previous element's values → slightly stale T in `jc()` near partition
   boundaries (JJc diagnostic only, thermal-coupled runs only). Candidate
   fix: dof-free `node_data("T")`, after confirming the vector-map entry
   exists.
2. **Quadratic edge-field synch (pre-existing, unverified, confidence
   low):** `synch_source_field` transfers one value per edge index, while
   quadratic edge fields store two dofs per edge (`2*index`, `2*index+1`).
   Needs a layout check before calling it a bug.
