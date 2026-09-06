# Periodic BC Fix Plan

**Date:** 2026-05-13 (original) · **Consolidated resume plan added 2026-06-01**
**Purpose:** Step-by-step plan to fix periodic boundary condition handling across DOF conversion, MPI distribution, and input wiring.
**Severity:** High (periodic BCs are silently broken at the DOF level and on multi-rank runs; cohomology/input machinery is missing from `periodic_new` entirely)
**Affected Modules:** `src/fem/kernel/` (DOF manager), `src/mesh/` (input + distributor), `src/fem/maxwell/`, `src/homology/`, `src/io/`
**Related Audits:** `devlog/dl20260513_periodic_bc_audit.md`, `devlog/dl20260513_periodic_branch_gap_audit.md`, `devlog/dl20260601_periodicity_factory_plan.md`, `devlog/dl20260608_periodic_thick_to_thin_audit.md`, `devlog/dl20260608_periodic_cut_claim_audit.md`
**Sub-plan:** The cohomology thin-cut continuity defect (periodic cuts not realized correctly) is now broken out into its own focused plan: **`todo/periodic_thin_cut_continuity_fix.md`** (3-way confirmed 2026-06-08). It supersedes and absorbs the earlier `todo/periodic_thin_cut_strategy.md` (now removed; its useful diagnostics/acceptance criteria were folded into the sub-plan).

---

## Consolidated Resume Plan (2026-06-01) — AUTHORITATIVE

This section merges the original Steps 0–8 (below) with the Codex lifecycle trace
(`devlog/dl20260601_periodicity_factory_plan.md`) and a Claude re-verification of
the current `periodic_new` tree. **Start here; the Step 0–8 detail below is kept for
reference but the ordering is superseded by this section.**

> **⚠ Stale symbol names (Codex audit 2026-06-08).** Several function names used
> throughout this plan **no longer exist in the current source** and must be
> re-mapped before acting on the DOF-half steps:
> `create_dofwise_periodicities_master()`, `flag_periodic_entities_12()`,
> `collect_nodes_from_flags_12()`, and `Periodicity::match_faces()` are all gone.
> The live post-cut path is `mMesh->periodicity()->update(); set_entity_dependencies();`
> (`src/fem/maxwell/cl_MaxwellFactory.cpp:476`), and periodic DOFs are now consumed
> through the **ordinary hanging-DOF T-matrix path** (DofData has no periodic-named
> functions). The real matchers are `PeriodicityFactory::match_nodes`/`match_edges`/
> `match_facets_and_faces`. Treat every mention of the four dead names below as
> "verify against current source first."

### Verified current state (code-checked 2026-06-01)

| Layer | State | Evidence |
|---|---|---|
| Step 0 (cohomology/input port) | **Done** | commits `3940461`, `2ca7a04`, `cec4fb7`, `75d35a7` on `periodic_new` |
| DOF entanglement + sequencing (orig. Steps 1–2) | **Done** | `cl_FEM_DofMgr_DofData.cpp`: `entangle()` :3596, `unfold_row()` :3556, `mHangingDOFs` rebuilt after `create_dofwise_periodicities_master()` via the `tAfter` pass :2444-2481 |
| **`PeriodicityFactory::create_periodicity()`** | **Done (verified 2026-06-08)** | `exit(0)` removed; now builds `Periodicity`, sets planes, calls `update_periodicity()`, returns (`cl_Mesh_PeriodicityFactory.cpp:10-21`). Edge/face derivation lives in `update()` as intended. |
| **`Periodicity::match_faces()`** | **Done (verified 2026-06-08)** | hard-coded element IDs 26903/25950 and the `exit(0)` are gone; no `exit(` remains in either periodicity file. |
| DOF periodic path actually firing | **Can now fire** | `create_periodicity()` runs ⇒ `mMesh->has_periodicity()` true ⇒ the DofData periodic path is reachable. End-to-end DOF verification (orig. Steps 4–5) still pending. |
| Maxwell input wiring | Parses planes only | `MaxwellFactory::create_periodic()` :377 sets source/target planes; the factory is run later via `CutFactory::run()` :134, not from `create_periodic()` |

**Correction to the 2026-05 notes:** `create_periodicity()` *is* invoked — from
`CutFactory::run()` at `cl_homology/cl_CutFactory.cpp:134` (after full-mesh edges/faces
are created at :123-129, reset at :180-183). It is not dead-wired; it is unfinished.

### Intended two-phase lifecycle (the key architectural insight)

Periodic edge/face derivation must run **twice** against two different edge/face universes,
so that logic lives in `Periodicity::update()`, not in the factory:

1. **Pre-cohomology (full temporary mesh).** `CutFactory::run()` creates all edges/faces,
   then calls `PeriodicityFactory::create_periodicity()`. This must populate node pairs +
   enough edge/face pairs for `SimplicialComplex` / `Cohomology::clean()` / `CutData` to
   fold the periodic boundary. Cohomology then duplicates interface nodes; periodic
   pointers and role flags must propagate to the duplicates. Temporary edges/faces are
   reset at `CutFactory::run()` :180-183.
2. **Post-cohomology (final Maxwell universe).** `MaxwellFactory::create_magnetic_kernel()`
   recreates only conductor-needed edges/faces, then rebuilds periodic node lists via
   `collect_nodes_from_flags_12()` and re-derives edges/faces via `update()`
   (`cl_MaxwellFactory.cpp:447-450`).

### Consolidated step sequence

- **C1 — Finish `create_periodicity()` as a node/facet geometry matcher only.** ✅ **DONE
  (verified 2026-06-08).** `exit(0)` removed; `create_periodicity()` builds the `Periodicity`,
  sets planes, and delegates edge/face derivation to `update()` (`cl_Mesh_PeriodicityFactory.cpp:10-21`).

- **C2 — Clean the face matcher.** ✅ **DONE (verified 2026-06-08).** The old `match_faces()`
  has been folded into `PeriodicityFactory::match_facets_and_faces()` (`:820`); the hard-coded
  element IDs 26903/25950 and the `exit(0)` are gone, and no `exit(` remains in the periodicity files.

- **C3 — Fix role-flag timing.** `MaxwellFactory::create_cuts()` calls
  `flag_periodic_entities_12()` at :764-767 **before** `create_cuts_sub_master()` →
  `CutFactory::run()` → `create_periodicity()` has created the periodicity, so the guarded
  call is currently a no-op. Set flags 1/2 immediately after the first successful
  `create_periodicity()`/`update()` (inside the factory after `set_periodicity()`, or in
  `CutFactory::run()` right after :134). Keep/audit the later flag use for the
  post-cohomology rebuild.

- **C4 — Smoke-test mesh-level periodicity** (before any input-syntax work). Confirm
  `has_periodicity()`, node pairs, edge pairs, face pairs, no `exit()`, and that role flags
  survive cut-duplicate node creation. Use `corc/corc.msh` via `periodictest`.
  **Correction (Codex audit 2026-06-08):** the earlier note that `periodictest.cpp` is stale
  because it uses "removed" `set_master_plane`/`set_slave_plane` is **refuted** — those are the
  current API names (`cl_Mesh_PeriodicityFactory.hpp:51-75`). `periodictest.cpp` already calls
  `set_master_plane`/`set_slave_plane` + `create_periodicity()` and is usable as-is; verify it
  builds/runs rather than renaming.

- **C5 — Test the post-cohomology rebuild.** After `create_edges_and_faces_on_mesh()` +
  `create_thinshells()`, confirm `collect_nodes_from_flags_12()` includes original **plus**
  duplicated/thin-shell-layer nodes, and `update()` recreates valid final edge/face pairs.

- **C6 — Re-validate the DOF periodic path** (now that it can actually fire): serial
  regression that the entangled slaves enter `mHangingDOFs` and are condensed out.
  Then resolve the line-2354 TODO (orig. Step 3: 2nd-order face slave-facet T-matrices in
  the periodic case).

- **C7 — Maxwell input wiring** (orig. Steps 6–7): the explicit sideset-filter syntax in
  `todo/periodic_input_extension_plan.md`. Only after C1–C6 are green.

### Key open design question (must resolve during C1/C2)

`Periodicity::match_edges()` mode 0 collects edges only from `DomainType::Conductor`
blocks (`cl_Mesh_Periodicity.cpp:208-222`). That may be correct for the **final** Maxwell
edge universe, but cohomology cuts traverse **air**, so the **pre-cohomology** phase likely
needs periodic edge pairs outside conductor regions. If so, `update()`/`match_edges()` need
an explicit scope/mode argument (e.g. "all edges" for phase 1 vs. "conductor+thin-shell"
for phase 2). Verify against what `SimplicialComplex`/`CutData` actually consume before
committing C1.

### Deferred — thin-shell coincident-duplicate node matching (logged 2026-06-02)

`PeriodicityFactory::match_nodes()` (`cl_Mesh_PeriodicityFactory.cpp:318-355`) pairs nodes
facet-locally and takes the **first** unflagged target node `Q` whose in-plane coordinates
match the source node `P`. Where a thin-shell tape pierces a periodic plane there are
*coincident* duplicate nodes (top/bottom of the tape) at identical in-plane coordinates, so
`P_top` can bind to `Q_bottom`. The completeness assert at `:357` guarantees the pair *count*
is correct but **not** that each pairing is the topologically correct one; a wrong node pair
then propagates silently into the edge/face keys (correct keys, wrong geometry).

This is harmless on meshes without in-plane coincident duplicates and is **deferred** until
thin-shell handling is addressed generally. When picked up, the likely fix is a
side-aware / global-kdtree disambiguation that resolves which duplicate a given `P` belongs
to (e.g. by the tape-normal side), rather than first-match-wins. Confirm against the `corc`
test mesh whether it actually contains coincident in-plane duplicates before investing.

---

## Branch Context

Current branch: `periodic_new` (forked from `sideconnectors`, today). A sibling branch `periodic` (last touched 2026-03-11, `origin/periodic`) holds substantive periodic-cohomology fixes that were never merged forward through the sideconnectors line. Both branches diverged from common ancestor `e3f28ce` and are 6 commits / 24+ commits ahead respectively. Neither is merged into `main`.

The original Codex+Claude audits were performed against `periodic_new` only and incorrectly verdicted cohomology as "present and correct" because the periodic-aware code lives on `periodic`, not on the current branch. Section 1 of this plan's predecessor audit is corrected in `devlog/dl20260513_periodic_branch_gap_audit.md`.

## Canonical Findings (Claude + Codex agree)

| Area | State | Confidence |
|---|---|---|
| Cohomology **thick** cut (generator) | Present and periodicity-aware (quotient fold + `clean()`) | high |
| Cohomology **thin** cut realization | **Broken under periodicity** — cut faces trimmed at periodic seam + cut duplicates never tied across the period. See `todo/periodic_thin_cut_continuity_fix.md` (3-way confirmed 2026-06-08) | high |
| Thin-shell factory (base layers) | Present and correct; side connectors out of scope | high |
| Periodic mesh → DOF conversion | **Present but mis-sequenced**: `create_dofwise_periodicities_master()` runs after `mHangingDOFs` is frozen, so periodic slaves are never registered there | high |
| Parallel handling | Mesh-level distribution works; DOF-level fails because the hanging-DOF MPI exchange iterates `mHangingDOFs` and misses the periodic slaves | high |
| Input file wiring | **Missing**: only the procedural `PeriodicityFactory` demonstration in `src/mesh/main.cpp:45-50`; no input parsing, no domain-type string parsing, no Gmsh `$Periodic` reader | high |

**Key locations:**
- `cl_FEM_DofMgr_DofData.cpp:2416-2424` — `mHangingDOFs` populated (too early)
- `cl_FEM_DofMgr_DofData.cpp:2431` — `create_dofwise_periodicities_master()` called (after freeze)
- `cl_FEM_DofMgr_DofData.cpp:3436-3557` — the periodic conversion routine
- `cl_FEM_DofMgr_DofData.cpp:2353` — TODO marking that slave-facet T-matrix construction ignores the periodic case
- `cl_MaxwellFactory.cpp:410-414, 728-730` — cut + magnetic-kernel periodic touchpoints
- `cl_ThinShellFactory.cpp:1137-1179` — layer-node periodic propagation
- `cl_Mesh_Distributor.cpp:115-148, 321-349` — mesh-level periodic distribution

## Ordering Rationale

Codex originally proposed input wiring first; Claude proposed DOF sequencing first. **Canonical order: port the missing cohomology/input machinery first (Step 0), then DOF sequencing.**

Reason: the procedural `PeriodicityFactory` API would normally be sufficient to drive tests without input wiring, but on `periodic_new` the cohomology layer itself is missing periodic awareness — `Cohomology::clean()`, `SimplicialComplex` chain/cochain coupling, `CutData` periodic edge handling, and the `Input_Section` `unique()`-call bug all need fixing before any downstream test could be valid. DOF and MPI fixes (Steps 2-5) are real bugs but cannot be verified against a still-broken cohomology layer. Input parsing wiring (Steps 6-7) remains the last layer of work.

## Step-by-Step Plan

### Step 0 — Port periodic cohomology/input machinery from `periodic` to `periodic_new`

Manually re-apply the four substantive commits from the `periodic` branch:

1. `45bf086` — bulk periodic cohomology/input plumbing (`SimplicialComplex` aPeriodicity param, `CutData` periodic edges, `MaxwellFactory` PeriodicityFactory plumbing, en_DomainType parser entries)
2. `6c5c985` — `Cohomology::clean()` periodic refinement (slave-flag-aware coboundary redirection)
3. `b1371f8` — further `clean()` and `SimplicialComplex` refinement (mOriginalEdges tracking, clean-loop guards)
4. `be97d6b` — comment out the five `unique()` calls in `cl_Input_Section.cpp::get_ids()` that sort IDs and destroy periodic node associations

Skip the two `Merge branch 'postproc_fix' into periodic` commits — `periodic_new` already has its own postproc_fix history via sideconnectors.

**Conflict risk by file:**

| Risk | Files | Notes |
|---|---|---|
| Low (port cleanly) | `cl_Input_Section.cpp`, `cl_CutData.cpp`, `cl_CutProcessor.cpp`, `cl_SimplicialComplex.hpp`, `cl_Homology.cpp` | `periodic_new` either didn't touch these or only touched them minimally in unrelated regions |
| Medium | `cl_Cohomology.hpp`, `cl_CutFactory.hpp`, `cl_MaxwellFactory.hpp` | Both branches added small things in similar areas |
| **High (manual merge required)** | `cl_Cohomology.cpp`, `cl_CutFactory.cpp`, `cl_MaxwellFactory.cpp`, `en_DomainType.cpp` | Both branches added substantively; high overlap risk |

**Required preservations during port:**
- `cl_SimplicialComplex.cpp:281` element-selection predicate on `periodic_new` is `is_flagged() && dimension() == tDim` — a sideconnector-line fix. The `periodic` branch still has the older `is_flagged() and number_of_corner_nodes()==tDim+1`. Keep `periodic_new`'s version.

**Required exclusions (debug artifacts on `periodic` that must NOT be carried over):**
- `cl_MaxwellFactory.cpp:431` (on `periodic`): unconditional `mMesh->save("mesh.exo")`
- `cl_CutFactory.cpp:500` (on `periodic`): uncommented `write_debug_cohomology(tSuggestedHomology)` (line 511 has the commented-out form — only the active one is the issue)

**en_DomainType.cpp nuance:**
`periodic_new` already has `DomainType::Periodic = 12` plus `AirPeriodic`/`BufferPeriodic`/`ConductorPeriodic`/`FerroPeriodic` (15/18/21/24). The `periodic` branch's contribution is only the *string parser* mapping `"periodic"` → `DomainType::Periodic`. Port the parser, extending it to also recognize `"air periodic"`, `"buffer periodic"`, `"conductor periodic"`, `"ferro periodic"` and dispatch to the specific values. Do NOT touch the enum on `periodic_new` (it's already richer than `periodic`'s).

**Verification after port:**
- Build with `make reset && make hphirun -j 20` per the build-system note in CLAUDE.md (static libraries don't update on plain `make`).
- Re-run the original audit on the merged tree — confirm Section 1 (cohomology/cuts) now actually matches "present and correct."
- Run any existing periodic-using executable to catch obvious regressions before continuing.

**Effort estimate:** half-day to full day depending on `cl_Cohomology.cpp` and `cl_MaxwellFactory.cpp` conflict resolution.

### Step 1 — Slave-DOF source policy: smaller-ID rule (decided 2026-05-13)

**Decision:** For every periodic DOF pair, hang the DOF with the **larger** `Dof::id()` on the DOF with the smaller `id()`. The mesh-level master/slave designation in `mesh::Periodicity` is ignored for the purposes of DOF hanging direction.

**Rationale:**
- Deterministic and rank-stable: every rank arrives at the same answer without negotiating who's master.
- Decouples DOF-level hanging from mesh-level master/slave designation, which is set at mesh-construction time and may not reflect a useful global ordering.
- Avoids the "both sides claim to hang on the other" ambiguity if a future pipeline reverses designations.

**Implementation in Step 2:** Inside `create_dofwise_periodicities_master()` at `cl_FEM_DofMgr_DofData.cpp:3436-3557`, for each (`tM`, `tS`) DOF pair, compare `tM->id()` and `tS->id()`; the larger-ID DOF takes its source from the smaller-ID DOF (with weight 1.0). If the smaller-ID DOF is already hanging from mesh-basis logic, the larger-ID DOF inherits the same source chain.

**Edge cases:**
- Same ID: impossible — they're distinct DOFs.
- Both already hanging: assert it doesn't happen (case is undefined in current Maxwell pipelines). Promote to a documented merge policy only if a future use case requires it.
- Slave (larger-ID) DOF already has unrelated sources from elsewhere: assert, same reason.

**No MPI extension needed (decided 2026-05-13):** Investigation confirmed that `Mesh_Distributor::select_entities()` at `src/mesh/cl_Mesh_Distributor.cpp:531-542,646-650,691-694` explicitly ghosts periodic partners onto every rank that owns one half of a pair, for nodes, edges, and faces alike. Consequently, when `create_dofwise_periodicities_master()` runs on a rank holding only the slave, the partner master `Dof*` is already locally available via the ghost copy. The existing hanging-DOF MPI exchange handles the rest once Step 2 fixes the `mHangingDOFs` sequencing. The partitioner itself does **not** constrain periodic pairs to the same rank (`cl_Mesh_Partitioner.cpp` has no references to `periodicity()` or `has_periodicity()`), but the distributor's ghosting makes that constraint unnecessary.

### Step 2 — Fix the DOF sequencing

Refactor `cl_FEM_DofMgr_DofData.cpp:2400-2460` so the order becomes:

1. `create_dofwise_t_matrices_master()` — mesh-basis hanging DOF promotion (unchanged).
2. `create_dofwise_periodicities_master()` — periodic slaves promoted (unchanged).
3. **NEW**: Rebuild `mHangingDOFs` from a single pass over `mDOFs` filtering on `Dof::is_hanging()`.
4. Apply the existing `DynamicBitset` pruning (currently at 2436-2459) only if still needed; the rebuild may make it redundant.

Notes:
- The rebuild must run on **all ranks**, not just rank 0, because non-root ranks also need `mHangingDOFs` populated for their slice (mesh periodicity is already distributed at this point).
- Update `mNumberOfHangingDofs` from the rebuilt list.
- Add a `BELFEM_ASSERT` consistent with the Step 1 policy.

Deliverable: ~30-line refactor; same function entry points.

### Step 3 — Resolve TODO at `cl_FEM_DofMgr_DofData.cpp:2353`

Extend slave-facet T-matrix construction to consider the periodic case. Today the surrounding loop only handles hanging facets from mesh subdivision. Mirror the periodic-facet logic from Step 2 so face DOFs on periodic slave facets get the right T-matrix entries even when the master facet is itself non-hanging.

Deliverable: small extension to the existing slave-facet loop; covered by Step 4 test.

### Step 4 — Serial periodic regression test

Build a small Maxwell test case using the procedural `PeriodicityFactory` API (no input wiring yet). Required coverage:

1. Original periodic nodes participate.
2. Cohomology cut duplicates remain paired.
3. Thin-shell layer nodes propagate periodicity per `cl_ThinShellFactory.cpp:1137-1179`.
4. Derived periodic edges and faces appear in `mHangingDOFs` after Step 2.
5. Element consolidation eliminates the slave global DOFs.

Side connectors are out of scope here (per audit scope). Use an existing thin-shell test mesh; modify it programmatically.

Deliverable: new test under `src/fem/maxwell/` or wherever similar regressions live; runs in serial.

### Step 5 — MPI periodic regression test (verification only — no MPI extension expected)

Same case as Step 4, run on 2 ranks with the periodic pair deliberately spanning ranks. Asserts:

1. Non-root ranks receive periodic node/edge/face pairs (already works via `populate_periodicity_data`).
2. Non-root ranks have ghost copies of periodic partners via `Mesh_Distributor::select_entities()` (already implemented at `cl_Mesh_Distributor.cpp:531-542,646-650,691-694`).
3. Non-root ranks rebuild `mHangingDOFs` with the periodic slaves (works after Step 2).
4. Final solution satisfies the periodic constraint at the cross-rank pair.

Per the Step 1 MPI analysis, the existing distributor ghosts periodic partners onto every rank holding one half of a pair, so the cross-rank `Dof*` source lookup in `create_dofwise_periodicities_master()` works without modification. If this test fails despite Step 2, the regression is likely in the distributor's ghosting (e.g. for thin-shell duplicate nodes) and the fix lives there, not in a new partitioner constraint.

Deliverable: MPI test exercising the same scenario as Step 4.

### Step 6 — Input wiring: design

Decide the user-facing syntax. Options to compare with Christian:

- **A.** Explicit `Periodicity` section in the input file with master/slave sideset IDs or plane equations.
- **B.** Domain-type prefix on sidesets — extend `en_DomainType.{hpp,cpp}` parser so strings like `AirPeriodic` round-trip through `domain_type(string)` (`src/mesh/en_DomainType.cpp:103` currently rejects them).
- **C.** Gmsh `$Periodic` block parsing in the mesh reader.

Recommendation: A as the primary path (explicit, decoupled from mesh format); B as a complement so the `AirPeriodic`/`BufferPeriodic`/`ConductorPeriodic`/`FerroPeriodic` enum values at `en_DomainType.hpp:39-55` become reachable from input. C deferred unless needed for Gmsh-native workflows.

Deliverable: short design note appended to this file or a new `todo/periodic_input_design.md`.

Refinement after inspecting `cmake-build-debug/input.conf` on 2026-05-15:
use the existing plane-based `periodic { source : ... ; target : ... ; }`
syntax now, and implement explicit sideset filters as the first input extension.
Detailed plan: `todo/periodic_input_extension_plan.md`.

### Step 7 — Input wiring: implementation

Implement the chosen syntax. Must be parsed and `mesh::Periodicity` constructed **before** `MaxwellFactory::create_cuts()` runs (`cl_MaxwellFactory.cpp:728`). Update:

- Input parser (likely `cl_InputFile` and Maxwell-specific parameter classes).
- `domain_type(string)` if option B is in scope.
- Wherever `Periodicity` is currently nullable so that input-driven construction lands in the same lifecycle slot the procedural API uses today.

Deliverable: input parsing path with parity to the Step 4 procedural setup.

### Step 8 — End-to-end input regression

Identical scenario to Steps 4 and 5, but driven from an input file. Confirms the wiring lands `mesh::Periodicity` correctly and the rest of the pipeline (steps 2-3) handles it transparently.

Deliverable: input-file-driven test; CI green on serial and MPI.

## Open Questions

1. **Slave-DOF policy** (Step 1) — assert or merge?
2. **Partitioner weighting** (Step 5) — needed, or do current flags suffice?
3. **Input syntax** (Step 6) — option A, B, or both?
4. **Existing periodic Maxwell tests** — none found in the audit; is there hidden test surface that this work might break? Worth a quick scan before Step 2.

## Out of Scope

- Periodic BCs in non-Maxwell physics (heat, mechanical). The fix is in shared DOF-manager code, so it should benefit them automatically, but no targeted testing is planned.

## Effort Estimate

| Step | Effort |
|---|---|
| **0 — port periodic cohomology/input from `periodic` branch** | **half-day to full day** |
| 1 — slave-DOF policy | 30 min discussion |
| 2 — DOF sequencing fix | 1-2 h |
| 3 — slave-facet T-matrix TODO | 1-2 h |
| 4 — serial regression | half day |
| 5 — MPI regression | half day |
| 6 — input design | 1 h discussion |
| 7 — input implementation | 1-2 days |
| 8 — end-to-end test | half day |

Total: ~4-6 days of focused work, gated on the Step 1 and Step 6 decisions and on the Step 0 port landing cleanly. Strongly recommend re-auditing after Step 0 before continuing to Step 1, since the original audit baseline was wrong.
