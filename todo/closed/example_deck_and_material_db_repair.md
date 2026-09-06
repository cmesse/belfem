# Repair Plan: First-Run Failures in `examples/` and the Material rho Database

**Date:** 2026-08-10
**Purpose:** Make a fresh checkout of BELFEM runnable by someone who has never run it —
specifically on the two paths that fail today for anyone without pre-existing cached
artifacts: a parallel first run (aborts in METIS during material construction) and a stale
rho-database cache (aborts with an opaque HDF5 message).
**Module:** `src/physics/materials` (R1, R2), `examples/` (R3)
**AIs involved:** Claude (investigation + plan + implementation), Codex + Grok (audit of this
file before any source edit — Christian's instruction, 2026-08-10)
**Status:** **CLOSED 2026-08-11** (currentness sweep) — moved to `closed/`. What verified it:
every §6 definition-of-done item is ticked against a run, not a review — `mpirun -np 2` in a
clean helix dir builds the database in-run and marches timesteps with zero aborts (2026-08-10
01:22), the genuine pre-format-change `Copper_RRR30.hdf5` triggers the warning and a correct
regeneration, `examples/README.md`'s gmsh line was executed verbatim, and D4's replacement deck
was run from a clean directory by Christian (61,418 nodes / 341,712 elements, `create_periodic`,
both rho databases built in-run, 12 timesteps to t = 300 ms, warm restart from its own memdump).
All four register rows the campaign owned are struck: DR-54, DR-57, DR-58, and DR-59 — the last
by deletion of the retired parallel builders in `6ce93549`. **The one residual, O2, is not part
of this plan's scope and is carried forward as register row DR-65** so it is not lost with the
move. *(Superseded status, kept for the record:)* **ALL FOUR DEFECTS CLOSED (2026-08-10).** D1–D3 are done and committed in
`bc578b5e`; **D4 was closed by Christian the same day** by replacing the `examples/corc`
deck+mesh pair with his working 6-tape model, which also retired the mechanism debate this
plan carried (→ `debt_register.md` DR-54, closed). Remaining items are §8's residual
decisions, none of them defects. Implemented 2026-08-10; all three defects
acceptance-tested (see §7, three attempts recorded). Post-implementation jury round COMPLETE:
**three findings**, verified and fixed same session — worker-side null-`SpMatrix` reference
bind in the Projector made defined, work-mesh leak in both serial builders closed, probe
hardened to the raw `htri_t` tri-state; re-acceptance green on both tests. *(An earlier
version of this line read "2/2 findings" while enumerating three; the count is three, matching
DR-57.)* Remaining Christian decisions are in §8. Plan-stage jury audit completed 2026-08-10
(Codex + Grok, blind; thread `tmp/ai_exchange/review_example_deck_repair.md`). Verdict: no-go
for R1 as originally written; all findings were verified in source and folded in below.
Session record: `devlog/dl20260810_examples_and_rho_database_repair.md`.

> **The audit's load-bearing correction:** R1's original rationale was **false**. The call
> `partition( tCommSize, false, true, false )` passes `aSetProcOwners = false` in the *second*
> position, and `Partitioner` applies METIS's result to element/node owners only inside
> `if( aSetProcOwnerships )` (`cl_Mesh_Partitioner.cpp:45-57`). Entities default to owner 0
> (`cl_Mesh_Basis.cpp:22-25`), `Distributor` ships only `owner() == aTarget`
> (`cl_Mesh_Distributor.cpp:308`), and rho evaluation skips non-owned nodes
> (`cl_Material_Metal.cpp:864`). **So this path does all its work on rank 0 regardless of
> contiguity.** Dropping CONTIG stops an abort; it does not create parallelism, and the
> "each rank computes its own nodes" claim must not be repeated. Found independently by both
> auditors, verified by Claude at every cited line.

> **Scope guards:**
> - **Do NOT change the `Mesh::partition` overload defaults.** The `false` siblings are
>   deliberate (Christian, 2026-08-09). R1 changes one argument at two *call sites*, not a
>   default.
> - No redesign of the rho-database format, no new file version scheme beyond the minimum
>   needed to detect an old file. R2 is a guard, not a migration.
> - `tmp/examples/Tape_Quench/*` and the `Validation` set are OUT (see O2).
> - No change to what the database *contains* or how rho is computed. Purely the
>   build/load/partition plumbing.

---

## 1. Background — established by direct test, 2026-08-09/10

All six examples were surveyed and two were run at HEAD.

| example | ships .geo | ships .msh | ships rho DB |
|---|---|---|---|
| circuit, costheta, garber, helix, sidecoating | yes | no | **no** |
| corc | **no** | yes (gmsh **2.2**) | yes |

**The convention is: ship the geometry, generate the mesh.** `corc` is the sole anomaly — a
stale binary mesh with no geometry source. Its `input.conf:111-115` declares
`periodic { source : 1,2,3 ; target : 5,6,7 ; }` in **vertex IDs**, and the shipped mesh
carries only **one** vertex, so it aborts in `MaxwellFactory::create_periodic` →
`PeriodicityFactory::set_master_plane` → `Mesh::vertex(1)` before any assembly.
`examples/garber/corc.geo` is **not** a drop-in source: its own deck uses different periodic
vertex IDs (`2,3,4 / 116,118,123`).

> **MECHANISM CORRECTED 2026-08-10** (gdb backtrace + mesh dumps, after Christian pointed out
> that his corc model runs). This section originally blamed the gmsh **2.2** format — "no
> `$Entities` block and therefore no vertices". **That is wrong, and the format is
> irrelevant.** BELFEM builds vertices from **type-15 point elements**, taking the id from the
> geometry tag, in both format paths (`cl_Mesh_GmshReader.cpp:326-338` for 2.2, `:630-635` for
> 4.1). `cmake-build-debug/corc/corc.msh` is *also* 2.2 and holds **10** of them, so it runs.
> The example mesh holds exactly one — id 89005, the **last** element in `$Elements`, geometry
> tag 9, serving `bearing { nodes : 9 ; }` — appended by the Python tool that produced the
> mesh (gmsh writes points *first*, as ids 1-10 in the working file). The six periodic
> vertices were simply never emitted. The two directories are also **different models**:
> `examples/corc` is 3-tape / 78,878 elements, Christian's run dir is 6-tape / 390k, which is
> why "corc runs" and "the shipped corc aborts" are both true.

**Consequence of the survey (CORRECTED after audit — my first draft overstated this):**
of the five examples shipping no rho database, only those declaring a **PureMetal with an
`RRR` key and angle-dependent rho** reach D1 — `MaterialFactory` gates `set_RRR` on that
(`cl_MaterialFactory.cpp:90-97`) and `Metal::set_RRR` gates the database on
`depends( rho, angleBxJ ) && mComputeTables` (`cl_Material_Metal.cpp:135`). That is
**helix, circuit and sidecoating**; `garber` declares only `ybco` and `costheta` only builtin
iron, so neither is affected. D1 is still a first-contact defect for half the examples.

There is no `README` in `examples/`, no `Allrun` counterpart to the four `Allclean` scripts,
and `examples/` is referenced nowhere in `doc/` or the top-level `README.md`.

## 2. Defects

- [x] **D1 (P1) — a run on ≥2 ranks aborts if the rho database must be built.** *(fixed
  2026-08-10, committed in `bc578b5e`; see §7 — the fix went three layers deeper than R1)*
  `Metal::populate_rho_database()` (`cl_Material_Metal.cpp:685-703`) builds the cache only
  when `<label>_RRR<n>.hdf5` is absent; serial populates serially, `comm_size() >= 2` calls
  `populate_rho_database_parallel()`. That partitions the work mesh with
  `tMesh->partition( tCommSize, false, true, false )` — against the signature
  `partition( nParts, aSetProcOwners = true, aForceContinuousPartitions = true,
  aResetVertexContainers = true )` (`cl_Mesh.hpp:706-710`), the third argument is
  **`aForceContinuousPartitions`, explicitly `true`**. The mesh is
  `create_database_mesh()`'s structured 95×35×37 order-2 tensor grid over (T, log B, angle)
  (`fn_create_database_mesh.hpp:20-24`). METIS returns −4, *"A contiguous partition is
  requested for a non-contiguous input graph"* — note that BELFEM's own formatter reports
  status −4 as *"Unknown Error"* via `metis_status()`, so the quoted sentence is METIS's own
  stderr, not BELFEM's message (auditors, 2/2) — and the abort lands in
  `MaxwellFactory::create_materials()` **before the mesh is read**.
  **Reproduced and complement-tested 2026-08-09:** helix on 2 ranks aborts with no cache;
  the identical run reaches t = 5.8 s once a serially-built `Copper_RRR50.hdf5` is present.
  **The same call appears in `Material_Alloy`** (`cl_Material_Alloy.cpp:537`).
  *Note: an earlier `debt_register.md` DR-57 wording said this flag was "inherited from a
  default". That was wrong — it is passed explicitly. Corrected here and in the register.*

- [x] **D2 (P2) — a stale rho cache silently shadows regeneration, then fails opaquely.**
  *(fixed 2026-08-10, committed in `bc578b5e` — `fn_rho_database_is_current.hpp`)*
  The gate is `file_exists()` only; the file carries no version stamp. The format has changed
  — old files hold `rrr` (lowercase) plus a `lambda` dataset, current ones hold `RRR` and no
  `lambda` — so `load_rho_database` dies with *"Dataset RRR of type double does not exist."*
  Confirmed by dumping keys: `tmp/examples/Tape_Quench/BuiltinMat/Copper_RRR30.hdf5` =
  `[Bmax,Bmin,Tmax,Tmin,label,lambda,rho,rrr]`; `sidecoatings/`, `examples/corc/` and a freshly
  generated file = `[Bmax,Bmin,RRR,Tmax,Tmin,label,rho]`. The message names an HDF5 dataset
  rather than telling the user their cache is stale.

- [x] **D3 (P2, docs) — the mesh-generation step is undocumented.** Five examples require
  `gmsh -3 <name>.geo -o <name>.msh` and nothing says so. *(fixed 2026-08-10 —
  `examples/README.md`, committed in `bc578b5e`)*

- [x] ~~**D4 — `examples/corc` ships a mesh missing its six periodic vertices.**~~
  **CLOSED 2026-08-10 by Christian, and none of the proposed repairs was used.** He replaced
  the deck+mesh pair outright with his working **6-tape** model — `examples/corc/` is now
  byte-identical to `cmake-build-debug/corc/` (`cmp` clean on both files) and carries 10 point
  elements (geometry tags 1-10) against a deck asking for `2,3,4 / 7,8,9`. The 3-tape,
  one-vertex model is gone, and with it the whole restore-vs-retire-vs-repoint question.
  **Verified from a clean dir holding only `corc.msh` + `input.conf`:** reads (61,418 nodes /
  341,712 elements / 3 blocks / 9 sidesets), passes `create_periodic`, builds both rho
  databases in-run, runs 12 timesteps to t = 300 ms, and warm-restarts from its own memdump.
  The two shipped `*_RRR*.hdf5` caches were deleted in the same change — safe **because** R1's
  fix landed first, and now demonstrated on a second deck. History of this item: see §1's
  mechanism correction, which was wrong twice before it was right. This is the last defect of
  this plan and the one remaining public-facing item (`debt_register.md` DR-54); it also
  gates DR-02, whose §4.1 matrix names corc as a deck.

## 3. Ordered steps

- [x] **R1 — D1 fix (REFRAMED after audit).** *(done — kept, but necessary-not-sufficient;
  clearing the abort exposed two further layers, §7.)* Change the **third** argument from `true` to
  `false` at `cl_Material_Metal.cpp:755` and `cl_Material_Alloy.cpp:537` (note the corrected
  Metal line — `:731` in the first draft was the *serial* path), with a comment recording
  both facts. **Correct rationale:** *stop forcing CONTIG on a partition call whose ownership
  side-effects are deliberately disabled and which currently aborts METIS before any rho work
  happens.* Contiguity cannot buy locality here because ownership is never applied; the whole
  computation stays on rank 0 either way. This is the minimal change that clears the abort,
  and it leaves the overload defaults untouched as required.
  **Verification:** delete `*_RRR*.hdf5` from a helix run dir, `mpirun -np 2 hphirun`; expect
  the database to build and the run to proceed. Then `h5diff` against the serially-generated
  file. **Read that check correctly (auditors' caution):** a green h5diff proves "rank 0 still
  computes the whole field", **not** "partitioning did not perturb the physics" — the latter
  is untestable here because no distribution actually occurs.

- [x] **R2 — D2 guard (EXTENDED after audit).** *(done — MPI-symmetric probe-and-broadcast,
  warn-and-regenerate, verified against a genuine old-format file.)* In
  `Metal::populate_rho_database()` **and the
  `Alloy` twin — O3 is closed, Alloy caches identically** (`cl_Material_Alloy.cpp:454-482`,
  `:837-840`) — treat a cache lacking the `RRR` dataset as stale: warn at Default level naming
  the file and the reason, then regenerate. Probe with the existing `hdf5::dataset_exists`
  rather than adding a version stamp. **The stale decision must be MPI-symmetric** (auditors,
  2/2): every rank takes this branch, so probe on rank 0 and broadcast the boolean, or have
  all ranks probe consistently — never let ranks diverge into regenerate-vs-load.
  Regeneration overwrites safely: `FileMode::NEW` is `H5F_ACC_TRUNC`.
  **Verification:** drop the old-format `Copper_RRR30.hdf5` into a run dir, run serially,
  expect a warning plus a regenerated file with the current key set.

- [x] **R3 — D3 docs.** *(done — `examples/README.md` written; its gmsh line was executed
  verbatim on helix before it was documented.)* Add `examples/README.md`: the ship-geometry/generate-mesh convention,
  the exact `gmsh` line, the note that a first run generates material databases (and that a
  first *parallel* run needs R1), and the `corc` anomaly with its current status.
  No `Allrun` scripts — out of scope unless Christian wants them.

- [x] **R4 — report D4 to Christian.** No code. *(reported 2026-08-10; the ruling itself is
  still owed — D4 stays open.)*

## 4. Risks

- **R1 is a behaviour change in a parallel path I can only test on 2 ranks** (2 cores
  available). Non-contiguous partitions could in principle change *which* rank computes which
  node; the bit-identity check in R1 is designed to catch any resulting difference.
- **R1's premise is that this partition needs no locality.** If some downstream consumer of
  the database mesh assumes contiguous ownership, the change is wrong. I checked the
  immediate consumer (`Distributor` + gather-by-ID) and found none — but this is the single
  assumption most worth an auditor's attention.
- **R2 changes an abort into a warning + recompute.** If a cache is stale for a reason other
  than the format change, silently regenerating hides it. Mitigation: the warning names the
  file and the reason.
- Neither R1 nor R2 is exercised by any existing test; there is no test target for materials.

## 5. Open questions

- [x] ~~**O1** — should R2 warn-and-regenerate, or hard-error with a clear "delete this file"
  message?~~ **Resolved by implementation (2026-08-10): warn-and-regenerate**, on Christian's
  go-ahead for the plan as written. The warning names the file and the reason, so a
  hand-made cache is never discarded silently. Reopen only if a user reports losing one.
- [ ] **O2 (carried out of this plan on 2026-08-11 — now `debt_register.md` DR-65)** —
  `tmp/examples/Tape_Quench/BuiltinMat` has been dead since **2026-04-08**: its
  custom material is labelled `buffer`, reserved by `8eae5f00` ("add MgO"). Fixing it means
  renaming the label in the deck *and* in the `MatData` `.so` that registers it. Out of scope
  here; needs a ruling on whether that deck is still wanted (it is a §4.1 regression case).
- [x] ~~**O3** — does `Alloy` need the same D2 guard?~~ **CLOSED by the audit: yes.** Alloy
  mirrors Metal's cache/load/save including the top-level `RRR`
  (`cl_Material_Alloy.cpp:454-482`, `:837-840`). Folded into R2.

## 6. Definition of done

- [x] `mpirun -np 2` in a clean helix dir builds the database and runs. *(2026-08-10 01:22 —
  database built in-run, 0.7 s projection, timesteps marching, zero aborts.)*
- [x] ~~Parallel-built database is bit-identical to the serial one.~~ **Amended by evidence
  (§7):** unachievable — the MUMPS/MKL projection solve is multithread-nondeterministic, and
  two *serial* builds already differ at ~1e-7. The achievable statement, which holds:
  **single-threaded builds are bit-identical**; 2-rank-vs-serial (7.3e-8) is at the same
  scale as serial-vs-serial (8.7e-8).
- [x] Old-format cache triggers a warning and correct regeneration. *(verified against the
  genuine pre-format-change `Copper_RRR30.hdf5`.)*
- [x] `examples/README.md` exists and its gmsh line has been executed verbatim.
- [x] `debt_register.md` DR-57/DR-58 updated, including the DR-57 attribution correction.
- [x] Jury findings folded in with attribution; devlog written. *(post-implementation jury
  round: 3 verified findings fixed same session; `devlog/dl20260810_examples_and_rho_database_repair.md`.)*

**Not part of the DoD, recorded for traceability:** everything above is committed in
`bc578b5e` (working tree clean, verified 2026-08-10). Christian's build in his own tree is
the only remaining confirmation.

---

## 7. Implementation result (2026-08-10, on Christian's go-ahead; final)

**R2 and R3: done and verified** (stale-cache guard in Metal+Alloy via
`fn_rho_database_is_current.hpp`, MPI-symmetric probe-and-broadcast, warn-and-regenerate —
tested against the genuine old-format file; `examples/README.md` written).

**D1 took three attempts; the record matters more than the tidy ending:**

1. **R1 (contiguity `true`→`false`)** — necessary, kept, insufficient. Clearing the METIS
   abort exposed the Distributor crashing on the empty worker sets (rank 1,
   `ProtoMesh::create_t_matrices`, "Cell index out of bounds").
2. **Master-only build + all-ranks load** — WRONG, withdrawn. Deadlocked: the "serial"
   builder's `Database( tMesh, … )` ctor is collective at `comm_size > 1` (config broadcasts
   + the Projector's barrier/solve/share pair), so rank 0 and rank 1 entered mismatched
   communication sequences. Root cause of the design error: the pre-flight "no hidden
   collectives" check read the function bodies but not their callees.
3. **Lockstep build (final)** — reading the whole chain showed the architecture beneath was
   already symmetric: `Projector::project` is a designed collective pair (master assembles,
   workers join the MUMPS solve, share/receive of the result), the Projector ctor already
   guards its assembly on rank 0, and `populate_tensor_mesh` populates only the master — so a
   worker's `create_database_mesh()` yields an empty-but-valid mesh. The fix REMOVES the
   `comm_size` branch: every rank runs the former serial path in lockstep; the master owns all
   nodes and does all evaluation; workers no-op through the node loop and receive the
   projected table inside the ctor. One enabling guard: `TensorMeshFactory::create_bsplines`
   returns early on an element-less mesh (its loops iterate the config-sized grid, not the
   container). The entire `populate_rho_database_parallel` apparatus is retired as never
   having been necessary.

**Acceptance evidence (2026-08-10 01:22):**
- Fresh dir, 2 ranks, no cache: database built in-run (projection 0.7 s), file written,
  timesteps marching, zero aborts — the exact first-user scenario.
- Residual history digit-identical to the serial run.
- **DoD amendment (evidence-forced):** "bit-identical h5diff" was unachievable from the
  start — two builds differ at ~1e-7 relative **even serially, before any of these changes**,
  because the MUMPS/MKL projection solve is multithread-nondeterministic. Proven both ways:
  serial-rebuild-vs-reference max rel 8.7e-8; 2-rank-vs-serial 7.3e-8 (same scale); and two
  `OMP_NUM_THREADS=1` builds are **bit-identical**. The achievable determinism statement is
  the single-threaded one, and it holds.

## 8. Residual items (for the post-implementation jury + Christian)

- [x] ~~The retired `populate_rho_database_parallel` bodies (~100 lines × 2, with the
  stuck-counter and dense-packing landmines) are now unreachable. Delete or keep: Christian.~~
  **DELETED 2026-08-10** (Christian: "if we are sure it is never called, we can remove it"),
  after a three-AI jury round confirmed no caller, no virtual/base-pointer path and no
  orphaned helpers — 222 lines out of `Metal` and `Alloy`. Thread:
  `tmp/ai_exchange/review_dead_code_deletion.md`. **Note the split:** only the *gather* defect
  died with the code. The dense-pack-vs-`node->index()` mismatch lives in the **shared**
  two-arg `populate_rho_database`, which the serial path still uses — `debt_register.md`
  DR-59 stays open for that half.
- The `create_bsplines` empty-mesh guard changes behaviour for EVERY tensor-mesh constructed
  on a non-master rank, not just this path — flagged to the jury as the main blast-radius
  question.
- [x] ~~D4 (corc example) still needs Christian's ruling~~ — **closed 2026-08-10** by replacing
  the deck+mesh pair (see D4 in §2). O2 (Tape_Quench `buffer` label) unchanged.
