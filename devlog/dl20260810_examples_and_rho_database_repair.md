# Devlog 2026-08-10 — First-Run Repair: `examples/` and the Material rho Database

**Date:** 2026-08-10
**Topic:** Making a fresh checkout runnable for someone who has never run BELFEM — the
parallel first run (METIS abort during material construction) and the stale rho-database
cache (opaque HDF5 abort), plus the missing `examples/README.md`
**AIs involved:** Claude (investigation, plan, implementation), Codex + Grok (blind jury,
twice: plan-stage and post-implementation)
**Claude Confidence:** high on the fixes (each acceptance-tested), medium on the blast radius
of the `create_bsplines` empty-mesh guard (§5)
**Verification:** acceptance run — fresh directory, no cache, `mpirun -np 2 hphirun` on helix
builds the database in-run and marches, residuals digit-identical to serial; stale-cache
regeneration verified against a genuine pre-format-change file. Committed in `bc578b5e`.
**Plan:** `todo/example_deck_and_material_db_repair.md` (D1–D4 / R1–R4, checkboxes now current)
**Register rows:** `todo/debt_register.md` DR-57 (closed), DR-58 (closed), DR-59 (open,
latent), DR-54 (open — the ruling below)

## Summary

Two defects stood between a fresh checkout and a first successful run, and neither was
visible to anyone with a working run directory, because a working run directory already
carries the cached artifacts that hide them. Both are fixed and committed. A third — the
`examples/corc` deck — needed no code change at all: Christian closed it by replacing the
deck+mesh pair, after the diagnosis I had recorded for it turned out to be wrong (§8).

The header fact worth carrying forward: **the parallel rho-database build never worked, and
never needed to exist.** What looked like a partitioning bug was a path that computed
everything on rank 0 anyway.

The session's own lesson is in §8: a plausible mechanism, one failing run, and no check
against the working counterexample sitting in the tree produced a confident and wrong entry
in the debt register — corrected only because Christian said "but my corc runs".

## 1. D1 — parallel first run aborted in METIS

`Metal::populate_rho_database()` built the cache only when `<label>_RRR<n>.hdf5` was absent,
and split on rank count: serial populated serially, `comm_size() >= 2` called
`populate_rho_database_parallel()`. That partitioned the work mesh with
`tMesh->partition( tCommSize, false, true, false )`, whose **third** argument is
`aForceContinuousPartitions`, passed explicitly `true`. The work mesh is
`create_database_mesh()`'s structured 95×35×37 order-2 tensor grid over (T, log B, angle),
whose graph is non-contiguous; METIS returned −4 and the run aborted inside
`MaxwellFactory::create_materials()` — **before the mesh was even read**.

It took three attempts, and the record matters more than the ending:

1. **Un-force contiguity** (R1) — necessary, kept, insufficient. Clearing the abort exposed
   the Distributor crashing on empty worker sets (rank 1, `ProtoMesh::create_t_matrices`,
   "Cell index out of bounds").
2. **Master-only build + all-ranks load** — wrong, withdrawn. It deadlocked: the "serial"
   builder's `Database( tMesh, … )` constructor is collective at `comm_size > 1` (config
   broadcasts plus the Projector's barrier/solve/share pair), so the two ranks entered
   mismatched communication sequences. The design error has a clean root cause: the pre-flight
   "no hidden collectives" check read the function bodies but not their callees.
3. **Lockstep build (final)** — reading the whole chain showed the architecture underneath was
   already symmetric. `Projector::project` is a *designed* collective pair (master assembles,
   workers join the MUMPS solve, then share/receive of the result), the Projector constructor
   already guards its assembly on rank 0, and `populate_tensor_mesh` populates only the master
   — so a worker's `create_database_mesh()` yields an empty-but-valid mesh. The fix therefore
   **removes the `comm_size` branch entirely**: every rank runs the former serial path in
   lockstep, the master owns all nodes and does all evaluation, workers no-op through the node
   loop and receive the projected table inside the constructor.

One enabling guard was needed: `TensorMeshFactory::create_bsplines` now returns early on an
element-less mesh (`cl_TensorMeshFactory.cpp:208-222`), because its loops iterate the
config-sized grid rather than the actual container.

**The load-bearing correction, found independently by both plan-stage auditors and verified at
every cited line:** the original rationale for R1 was false. `Partitioner` applies METIS's
result to element and node owners only inside `if( aSetProcOwnerships )`
(`cl_Mesh_Partitioner.cpp:45-57`), and the *second* argument disables exactly that. Entities
keep owner 0 (`cl_Mesh_Basis.cpp:22-25`), `Distributor` ships only `owner() == aTarget`, and
rho evaluation skips non-owned nodes. So this path did all its work on rank 0 regardless of
contiguity — dropping CONTIG stops an abort, it does not create parallelism. The claim "each
rank computes its own nodes" must not be repeated.

## 2. D2 — a stale cache shadowed regeneration, then failed opaquely

The gate was `file_exists()` alone, and the file carries no version stamp. The format has since
changed: old files hold `rrr` (lowercase) plus a `lambda` dataset, current ones hold `RRR` and
no `lambda`. So an old cache was loaded rather than rebuilt, and `load_rho_database` died with
*"Dataset RRR of type double does not exist."* — a message naming an HDF5 dataset instead of
telling the user their cache predates the format.

Fixed with `fn_rho_database_is_current.hpp`: probe for the `RRR` marker, and on a miss warn at
Default level naming the file and the reason, then regenerate (`FileMode::NEW` is
`H5F_ACC_TRUNC`, so the overwrite is safe). Applied to **both** `Metal` and `Alloy` — the twin
was found by the audit (O3), Alloy mirrors Metal's cache/load/save including the top-level
`RRR`. The decision is **MPI-symmetric by construction**: every rank reaches this branch, so
the master probes and broadcasts the boolean, and ranks can never diverge into
regenerate-vs-load.

The probe returns the raw `htri_t` tri-state and treats only `> 0` as current, so a probe
*error* lands on the rebuild side rather than being swallowed by a bool conversion — that
hardening came from the post-implementation jury.

## 3. D3 — the mesh-generation step was undocumented

`examples/` had no README, no `Allrun` counterpart to its four `Allclean` scripts, and was
referenced nowhere in `doc/` or the top-level `README.md`. Five of six examples ship a `.geo`
and no `.msh`, and nothing said `gmsh -3 <name>.geo -o <name>.msh`. `examples/README.md` now
documents the ship-geometry/generate-mesh convention, the exact gmsh line (executed verbatim
before being written down), the note that a first run generates material databases, and the
`corc` anomaly. The format-4.x requirement is stated explicitly, because several decks address
geometry entities by vertex IDs that only a 4.x file carries.

## 4. Jury rounds

Two blind Codex + Grok rounds, per the frozen protocol.

**Plan-stage:** verdict *no-go for R1 as written*. Both auditors independently found the false
rationale of §1; both flagged the MPI-symmetry requirement on D2; Codex closed O3 (Alloy needs
the same guard). All findings verified in source before the plan was rewritten. Thread:
`tmp/ai_exchange/review_example_deck_repair.md`.

**Post-implementation:** three findings, all verified and fixed the same session — a worker-side
null-`SpMatrix` reference bind in the Projector (undefined behaviour, made defined), a work-mesh
leak in both serial builders (`delete tMesh` after the `Database` takes its own copy), and the
`htri_t` probe hazard above. Re-acceptance green on both tests.

## 5. What is deliberately left open

- ~~**DR-54 / D4 — `examples/corc` ships a mesh missing its six periodic vertices.**~~
  **CLOSED the same day by Christian** — see §8. He replaced the deck+mesh pair with his
  working 6-tape model rather than repairing the 3-tape one, which made the whole
  restore-vs-retire-vs-repoint question moot.
- **The retired `populate_rho_database_parallel` bodies** (~100 lines in each of `Metal` and
  `Alloy`) are now **unreferenced dead code** carrying two landmines: a receive loop whose read
  counter never advances (`cl_Material_Metal.cpp:828-833`, identical in Alloy), and a dense
  owned-node packing that the `Database` projector reads by `node->index()` — correct only while
  one rank owns everything. Delete or keep is Christian's call; if kept, both must be fixed
  before anyone adds a caller. Tracked as DR-59.
- **The `create_bsplines` empty-mesh guard changes behaviour for every tensor mesh constructed
  on a non-master rank**, not only this path. Flagged to the jury as the main blast-radius
  question; no other consumer was found, but this is the assumption most worth revisiting.
- **O2** — `tmp/examples/Tape_Quench/BuiltinMat` has been dead since 2026-04-08 (its custom
  material is labelled `buffer`, reserved by `8eae5f00`). Out of scope here; it is a §4.1
  regression case, so it needs a ruling on whether the deck is still wanted.

## 6. Determinism, stated correctly

The plan's definition of done asked for a **bit-identical** h5diff between the parallel-built
and serial-built database. That was unachievable from the start, and not because of anything
changed here: two builds differ at ~1e-7 relative **even serially**, because the MUMPS/MKL
projection solve is multithread-nondeterministic. Proven both ways —
serial-rebuild-vs-reference max rel 8.7e-8, 2-rank-vs-serial 7.3e-8 (same scale), and two
`OMP_NUM_THREADS=1` builds **bit-identical**. The achievable determinism statement is the
single-threaded one, and it holds. The DoD was amended rather than quietly dropped.

## 7. Commit note

The work is in `bc578b5e` (2026-08-10). That commit is **not** limited to this campaign — it
also carries the PID/ARPACK controller work, the DR-12 `BfmFile::load` mesh-checker trust flag,
and an unrelated ruling on `mRhoMin` (floor set to 0.0 so the power-law resistivity and
`drho_powerlaw_dJ` agree across the whole subcritical band, from the
`review_newton_tangent_derivatives` thread). Read the file list, not the commit message, when
attributing changes.

## 8. Correction: the corc diagnosis was wrong, and how it was caught

Christian read the DR-54 write-up and objected that his corc model runs. It does. The two
statements were reconciled by direct evidence, and the original diagnosis did not survive.

**What was claimed:** *"`corc.msh` is gmsh format 2.2, which has no `$Entities` block and
therefore no Point/vertex entities"*, aborting with *"Tried to access invalid vertex id: 1"*.

**What is true:**

- **The format is irrelevant.** BELFEM builds vertices from **type-15 point elements**, taking
  the id from the geometry tag, in *both* format paths (`cl_Mesh_GmshReader.cpp:326-338` for
  2.2, `:630-635` for 4.1). A 2.2 file carries them fine — `cmake-build-debug/corc/corc.msh`
  is **also 2.2** and holds ten (geometry tags 1-10, physical tag 0, the `-save_all`
  signature), which is why its deck's `2,3,4 / 7,8,9` resolves.
- **The shipped example mesh holds exactly one point element** — id **89005**, the *last*
  entry in `$Elements`, geometry tag 9, on node 1 — serving `bearing { nodes : 9 ; }`. gmsh
  writes points *first* (ids 1-10 in the working file), so this one was **appended after the
  fact**, matching Christian's account that a Python tool produced the mesh. The six periodic
  vertices were never emitted.
- **The two directories are different models.** `examples/corc` is 3-tape (`Volume_1..2`,
  `Tape_1..3`, 78,878 elements); the run dir is 6-tape (`Volume_1..3`, `Tape_1..6`, 390k).
  Both "corc runs" and "the shipped corc aborts" are true, of different meshes.
- **The quoted error message is not the one you get.** gdb on the shipped deck:
  `MaxwellFactory::create_periodic` (`cl_MaxwellFactory.cpp:123`) →
  `PeriodicityFactory::set_master_plane(A=1,B=2,C=3)` (`cl_Mesh_PeriodicityFactory.cpp:96`) →
  `Mesh::vertex(1)` → `Map<uint,Element*>::operator()` with `aKey = 1`. `BELFEM_ASSERT` is
  compiled out in this tree (NDEBUG), so the readable assert at `cl_Mesh.hpp:1433` never fires
  and the abort surfaces as `Map::operator()`'s always-on `BELFEM_ERROR` fallback:
  **"Key not found in map"** (`cl_Map.hpp:235`) — no id, no map name.

**Consequences.** The wrong mechanism had already propagated into `examples/README.md`
(user-facing), the plan §1, `debt_register.md` DR-54 and DR-02, and both index entries; all
corrected.

**Resolution, same day.** Christian did not repair the 3-tape mesh at all — he replaced the
pair. `examples/corc/` is now byte-identical to `cmake-build-debug/corc/` (`cmp` clean on
mesh and deck): the **6-tape** model, 10 point elements with geometry tags 1-10, deck asking
for `2,3,4 / 7,8,9`. He also deleted the two shipped `*_RRR*.hdf5` caches, which is safe
**because** the D1 fix above landed first.

Verified from a clean scratch dir holding only `corc.msh` and `input.conf`: reads (61,418
nodes / 341,712 elements / 3 blocks / 9 sidesets), passes `create_periodic`, **builds
`Copper_RRR50.hdf5` and `Silver_RRR20.hdf5` in-run from the empty directory**, runs 12
timesteps to t = 300 ms, writes exodus + `memdump.hdf5`, and a second invocation warm-restarts
from it at t = 300 ms and continues. That is the D1 first-run path demonstrated on a second,
larger deck — the two halves of this session's work meeting in one run.

**Two lessons worth keeping.** First, the original conclusion was reached from a plausible
format difference and a single failing run, without checking a *working* counterexample that
sat in the tree the whole time — one `head -3` on the other mesh would have refuted it.
Second, `Map::operator()`'s release-path message discards the diagnostic its debug-path
sibling carries; that asymmetry is what made a precise defect look like a generic map miss,
and it is a small always-on-error-message defect in its own right.
