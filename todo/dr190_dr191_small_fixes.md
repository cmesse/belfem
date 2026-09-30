# DR-190 and DR-191: the `DofManager` across a `create_field` refusal, and the unset transverse Kohler switch

**Date:** 2026-09-29
**Purpose:** two small fixes behind one plan jury and one code jury. DR-190: `Kernel::create_field`
leaks the whole `DofManager` on every refusal between `new` and the push; the register row's other
half (a rank-dependent gate before the collective thin-shell linker) is refuted by the source and is
struck with the evidence. DR-191: valgrind's uninitialised read in `Lead::kohler` and
`WhiteTin::kohler` is the transverse Kohler switch `mBSKohlerTransSwitch( 0 )`, sized in
`create_kohler` and never assigned; the fix is the one assignment the longitudinal branch already
makes, plus a continuity test across the switch.
**Module:** `src/fem/kernel` (DR-190), `src/physics/materials` (DR-191), `tests/physics`
**AIs involved:** Claude (plan + implementation), Codex `gpt-6-astra` + Grok (juries; the kernel puts the round on the astra row)
**Status:** IN PROGRESS — approved 2026-09-29 (Christian, quiz); R1-R5 done: plan jury and code jury both Codex astra/high + Grok 4.7/high, reconciled (no P0/P1 in the code round); Claude's gates on the private non-MKL tree green (10/10, 9/9, no `create_field` leak record, no `kohler` valgrind context; mutation red at `x/4`). Two residues ruled 2026-09-29 (Christian, quiz: file, commit as is): the dangling `IWG::mField` after a CAUGHT refusal (fix shape `IWG::detach_field()` + try/catch, no in-tree caller reuses an equation after a catch) and the test's `( x/4, 2x )` window (samples at `0.9x`/`1.1x` would pin the stored switch). R6 closed 2026-09-30 by the residue plan `dr113_dr190_dr191_residues.md` (both rows struck); `make check` in the shared tree, Christian, 2026-09-30: 17/17 (fast 8, mpi 4, 21.24 s real)

> **Scope guards:**
> - DR-190: ownership in `create_field` only. No change to what is refused, to `set_equation`'s
>   order, or to the production (MPI-abort) reaction — the fix serves the debug-throw path that
>   `make check` and valgrind exercise.
> - DR-191: the two assignments, the WhiteTin size, one test. No change to the Kohler polynomials
>   or to the alloy table format.
> - The MUMPS context valgrind reports on the same test (§2.3) is recorded and dispositioned, not
>   silently folded in.

---

## 1. DR-190

### 1.1 The gate half is refuted

Row text: `BlockData::create_blocks` "decides whether to enter the COLLECTIVE
`link_thin_shell_facets()` from `tAllSelectedThinShellFacets`, which is built from this rank's
`ThinShell` records ( `:134-137`, collector `:141-170` )".

Source at HEAD `31be8ee`: the gate reads `tAllSelectedThinShellFacets.length() > 0 || tHaveConnectorBlocks`
(`cl_FEM_DofMgr_BlockData.cpp:149-152`). The vector is filled by
`collect_thin_shell_facet_ids` (`:193-208`): rank 0 runs `selected_thin_shell_facet_ids` on its
mesh and `share`s the result, every other rank `receive`s it — so all ranks hold root's list, and
root holds the full mesh. the static helper `selected_thin_shell_facet_ids` (`:158-190`) is what the collector
calls on rank 0 only; the row read the helper without its caller, and its `:134-137` / `:141-170`
pair no longer points at either. The second half of the gate reads the IWG's block-type map, which
`IWG_Maxwell::initialize` replicates from rank 0 (`cl_IWG_Maxwell.cpp:144-203`); the worker
`ThinShell` subset of `cl_ProtoMesh.cpp:1416-1425` is not consulted by the gate at all. `git log -S'share( aAllSelectedThinShellFacets )'`
shows the share arrived in `f81dd84` (2026-09-05, the release commit) and was never removed, so it
predates the row (filed 2026-09-29). The second half of the gate (`tHaveConnectorBlocks` from the
IWG's broadcast block-type map, `:141-148`) was already rank-uniform. **No gate fix is owed**; the
row is rewritten to say so (confidence high: source trace plus history).

### 1.2 The leak

`Kernel::create_field( IWG * )`, `cl_FEM_Kernel.cpp:588-608`:

```
aEquation->initialize();                              // :590
DofManager * aField = new DofManager( this, n );      // :591
if ( mCommRank == 0 ) BELFEM_ERROR( check_mesh == 0 ) // :595-598   refusal A
comm_barrier();                                       // :601
aField->set_equation( aEquation );                    // :603       refusals B (inside)
mDofManagers.push( aField );                          // :605
```

A refusal at A, or anywhere inside `set_equation` (`cl_FEM_DofManager.cpp:85-130`: `create_dofs`,
`create_blocks` → the connector refusals in both linkers, `create_sidesets`, …), throws past the
push under the debug-throw reaction (policy `src/core/assert.hpp:89-101`, the throw `:190-192`),
and the object is never reachable again: `~Kernel` frees the managers it holds (`cl_FEM_Kernel.cpp:126-128`,
then IWGs, comm tables and the mesh `:131-147`), never an unpushed one. Found by the
DR-175 plan jury (Grok); the `ThinShellSelection.*Refused` cases are the reproducer.

### 1.3 Is `~DofManager` safe on an object whose `set_equation` half-ran?

The constructor allocates eight sub-objects and nothing else (`cl_FEM_DofManager.cpp:31-54`); the
destructor deletes the postprocessors and the eight (`:58-80`). Each sub-destructor is `reset()`
plus null-guarded deletes:

| sub-object | destructor | state after a half-run `set_equation` | safe? |
|---|---|---|---|
| `DofData` | `reset()`: deletes `mDOFs`, `mHangingDOFs`, clears maps (`cl_FEM_DofMgr_DofData.cpp:64-95`) | dofs may exist (`create_dofs` ran) | yes: the same path `set_equation` itself calls first (`DofManager::reset`, `:133-144`) |
| `BlockData` | `reset()`: guarded `mEmptyBlock`, deletes `mBlocks` (`cl_FEM_DofMgr_BlockData.cpp:40-64`) | blocks and `mEmptyBlock` built before the linker refuses (`:127-135`) | yes |
| `SideSetData`, `BearingData` | `reset()` with guarded empties (`…SideSetData.cpp:40-43`, `…BearingData.cpp:37-40`) | nothing built (refusals precede `create_sidesets`) | yes |
| `FieldData` | `reset()` clears counters and lists (`…FieldData.cpp:46-62`) | empty | yes |
| `SolverData` | `reset()`; `mSolver`, `mAndersonX/R`, `mSystemMatrix` null-initialized in the header (`…SolverData.hpp:102`, `:194`, `:218-219`) | none set (`set_solver` is a later call) | yes |
| `EigenValues` | four guarded deletes, members null-initialized (`…EigenValues.hpp:75`, `:204`, `:256`, `:265`) | none set | yes |
| `Parameters` | trivial | — | yes |

Two pointers the fix leaves behind, both already tolerated by the normal teardown: `IWG::mField`
(set by `aIWG->set_field( this )`, `cl_IWG.cpp:710-714`) dangles after the delete, and `~IWG` does
not read it (`:70-80`); `~Kernel` already deletes every field before every IWG (`:126` before
`:133`), so the IWG outliving its field is the existing contract. The `ThinShellSelection`
fixture deletes the kernel, then the parameters, then the mesh (`test_ThinShellSelection.cpp:351-358`).

Residues, established in the jury round: `Element::set_facet` writes a FEM element the `Block`
owns (`cl_FEM_Element.hpp:550-552`), so `BlockData::reset` takes it with the block; the element
and facet flags set in `collect_element_indices` (`:234`, `:252`, `:276`, `:289`) survive the
refusal and are cleared by the next `create_blocks` at entry (`:72-74`) — nothing in between reads
them in the refusal tests, which destroy the mesh with the fixture. A throw inside `create_dofs`
(its own `BELFEM_ERROR`s) was not audited for partial allocations; the gate exercises the
`create_blocks` refusals only.

### 1.4 Fix

`create_field` owns the object until the push, the pattern the September factories use
(`cl_MaterialFactory.cpp:149`, `cl_Mesh_PeriodicityFactory.cpp:158`, `cl_CurveFactory.cpp:613`):

```cpp
std::unique_ptr< DofManager > tOwner( new DofManager( this, mDofManagers.size() ) );
DofManager * tField = tOwner.get();
…
// the push may throw, so the owner keeps the field until the slot exists
mDofManagers.push( tField );
tOwner.release();
return tField ;
```

`Cell::push` is `std::vector::push_back` (`cl_Cell.hpp:265-267`), which can throw; releasing first
would reopen the leak on exactly that throw. The material factory's hand-off has the same order
(`cl_MaterialFactory.cpp:198-201`). The local is `tField`, not the argument-prefixed `aField` the
old line used.

`<memory>` added to `cl_FEM_Kernel.cpp`. A one-off setup allocation, so the smart pointer is inside
the coding philosophy's allowance.

## 2. DR-191

### 2.1 Localized

`valgrind --track-origins=yes test_physics --gtest_filter='AlloyRRR.ARhoDatabaseForAnotherRRRIsRefused'`
on Claude's private non-MKL debug tree (`/opt/scls/debug`, Armadillo, HEAD `31be8ee`, this session):
246051 errors from 3 contexts.

| context | at | origin |
|---|---|---|
| 1 | `Lead::kohler`, `cl_Material_Lead.hpp:114` — `if ( BS < mBSKohlerTransSwitch( 0 ) )` | the `Lead` object's heap block (`MaterialFactory::create_material`, `cl_MaterialFactory.cpp:638`, via `Alloy::create_component`) |
| 2 | `WhiteTin::kohler`, `cl_Material_WhiteTin.hpp:94` — same test | the `WhiteTin` object (`cl_MaterialFactory.cpp:637`) |
| 3 | `dmumps_` under `mumpstools_create_solver` ← `MUMPS::initialize` (`cl_SolverMUMPS.cpp:287`) ← `Projector::project` (`cl_DatabaseProjector.cpp:275`) ← `Database::Database` (`cl_Database.cpp:74`) | a `malloc` in `mumpstools_create_solver` |

Contexts 1 and 2 are the same defect: `create_kohler` sizes `mBSKohlerTransSwitch( 1 )`
(`cl_Material_Lead.cpp:207`, `cl_Material_WhiteTin.cpp:214`) and assigns nothing to it
(`grep -rn KohlerTransSwitch src/physics/materials`: the two `set_size` lines and the two reads are
the only hits). The longitudinal branch sets its switches at the Taylor matching points
(`mBSKohlerLongSwitch( 0 ) = std::exp( w1 )`, `:247` / `:263`); the transverse low-field cubic `q0`
is matched at `x = std::exp( w )` with `w = 3.0` (Lead, `:250-251`) and `w = 8.4` (WhiteTin,
`:266-267`) and the switch that should carry that `x` is never written. Copper's single-switch
variant does write it (`cl_Material_Copper.cpp:285`, `mBSKohlerSwitch( 0 ) = x`). So `Atrans`
picks its branch on garbage: for `beta ≠ 0` the alloy table's field dependence is wrong wherever
the garbage lands on the wrong side of the true switch — a wrong-answer defect in every deck that
builds a `Pb38Sn62` table (`pancake`, `tapestack3d`, `tapestack3d_quench`), invisible when the
garbage happens to be large (then `q0` is used everywhere, which is right below `x` and merely
inexact above it).

### 2.2 Fix

- `cl_Material_Lead.cpp`, after `q0( 3 ) = 0.0 ;` (`:267`): `mBSKohlerTransSwitch( 0 ) = x ;`
- `cl_Material_WhiteTin.cpp`, after `q0( 3 ) = 0.0 ;` (`:283`): `mBSKohlerTransSwitch( 0 ) = x ;`
- `cl_Material_WhiteTin.cpp:212`: `mBSKohlerLongSwitch.set_size( 3 )` → `( 2 )`; only entries 0 and
  1 are read (`cl_Material_WhiteTin.hpp:80-92`) or written (`:263-264`), so the third is dead.

Regression, `tests/physics/test_MetalFieldDependence.cpp`: `TEST( MetalFieldDependence,
KohlerTransverseBranchesMeetAtTheStoredSwitch )` — for `Lead` and `WhiteTin` built through
`MaterialFactory`, at `beta = π/2` (pure transverse, `c² = 0`) and `S = 1`:
- **far-point branch checks (the deterministic oracle):** `kohler( x/4 )` equals the low-field cubic
  and `kohler( 2x )` equals `exp( q1( ln 2x ) )`, each to `1e-10` relative, with `x = e^{3.0}` /
  `e^{8.4}` and the test's own oracle built from the hardcoded `q1` (`{ 1.458, −9.7936 }`,
  `{ 1.8073, −12.545 }`) by the same value/slope/curvature match. The two branches differ by several
  percent at `x/4` (Codex's probe: ~10.7 % Lead, ~7.4 % WhiteTin), so a zero, a huge, or a removed
  switch is red on one of the two points on either backend;
- **continuity at the switch:** `kohler( x( 1 ∓ δ ) )` agree to `1e-6` relative at `δ = 1e-4`
  (the C2 match makes the true gap `O(δ³)`); a control, not the oracle;
- at `beta = 0` the same far-point shape over the longitudinal switches as the control.

A continuity-only test is not a mutation gate: a switch outside the sweep selects one smooth branch
throughout, and both branches are increasing there (jury round, Codex's numeric probe).

**Stale caches.** `rho_database_is_current` accepts any file carrying an `RRR` link
(`fn_rho_database_is_current.hpp:47`) and an accepted file is loaded, not rebuilt
(`cl_Material_Alloy.cpp:512`, `:542`), so every `Pb*Sn*_RRR*.hdf5` written before the fix keeps the
wrong transverse branch. Operational disposition, no format change: delete those caches in every
run directory (the shipped examples write them beside the deck); said in the devlog and the
register row.

### 2.3 The MUMPS context — disposition open

The block is the solver pool `gSolvers( gMaxNumSolvers )`, a Fortran `allocate` of
`DMUMPS_STRUC` in `mumpstools_create_solver` (`src/sparse/mumpstools.f90:65-70`), which Fortran does
not zero. Before `DMUMPS` the wrapper sets exactly the four fields the MUMPS interface documents
for `JOB = −1` — `COMM`, `PAR`, `SYM`, `JOB` (`:104-107`) — and `dmumps_` then branches on
another field of the same block. `libdmumps.so` in `/opt/scls/debug` carries no line info, so which
field is not known from here; the candidate is MUMPS's own "already initialized" tag (`KEEP(40)`),
which the driver reads on every call, including the first, to tell a fresh instance from a reused one
(recollection, not verified against MUMPS 5.9.1's `dmumps_driver.F`; low confidence). If that is it,
the read is MUMPS inspecting memory only the caller can have zeroed, and a wrapper-side clear of the
slot's `KEEP` before `JOB = −1` would make it deterministic — an L-21 wrapper matter, not a BELFEM
logic defect. Not fixed in this round; recorded on the register. Jury: name the field if you can, from
the MUMPS 5.9 user guide or driver source.

## 3. Gate

- **DR-190:** on the private non-MKL tree, `valgrind --leak-check=full test_fem --gtest_filter='ThinShellSelection.*Refused'`
  — before: definite-loss records with `DofManager` frames (positive control, this session);
  after: none. Plus `test_fem --gtest_filter='ThinShellSelection.*'` green.
- **DR-191:** `valgrind --track-origins=yes test_physics --gtest_filter='AlloyRRR.ARhoDatabaseForAnotherRRRIsRefused'`
  — after: the two `kohler` contexts gone (the MUMPS one per §2.3); the new test green, and red
  with the two assignments reverted (the far-point check fails on whatever the unset entry holds).
- **Christian:** `make check-fast` on the shared trees.

## 4. Steps

- [x] **R1** plan jury (`review_dr190_dr191_small_fixes`, Codex `gpt-6-astra`/`high`, Grok `grok-4.7`/`high`).
- [x] **R2** DR-190 ownership in `create_field`; register row rewritten (gate half struck with §1.1).
- [x] **R3** DR-191 assignments, WhiteTin size, the continuity test.
- [x] **R4** code jury on the own-files diff, same slug and depth.
- [x] **R5** gates of §3 on the private tree, mutation checks recorded; Christian's `make check-fast`.
- [x] **R6** close-out — done 2026-09-30 through `todo/dr113_dr190_dr191_residues.md`: both rows struck, devlog `dl20260930_dr113_dr190_dr191_residues`, `todo/README.md`.

Defects found on the way:
- **D1** (DR-191, this session) the unset switch, §2.1. Confidence high (valgrind origin + source).
- **D2** (DR-190, this session) the gate half is a misreading, §1.1. Confidence high.
- **D3** (open) the MUMPS context, §2.3.

## 5. Open questions

- **O1** — **RESOLVED 2026-09-29 → not this round, filed on the DR-191 residue line.** Should `create_kohler` in the nine metals that carry switches (`test_MetalFieldDependence.cpp:28-30`)
  get a `set_size( n, BELFEM_QUIET_NAN )` fill plus a `BELFEM_ASSERT` that every entry is finite at
  the end? The pair is a debug-build detector (Grok: the assert fires on an unwritten NaN; a NaN fill
  alone would silently select the exponential branch in release). Proposal: not in this round — nine
  files for a class the far-point test and valgrind now cover; filed as a follow-up line on the
  register row. Christian's call.
- **O2** — **RESOLVED 2026-09-29 → DR-190 stays open as `[W]`, narrowed to the dangling-field window; the gate half is struck as refuted and the leak as fixed.** Original: with the gate refuted and the leak fixed the row closes; proposal was to strike it.

## 6. Audit trail

`tmp/ai_exchange/review_dr190_dr191_small_fixes.md`; valgrind logs in the session scratchpad
(`nomkl/vg/alloy_rrr.vg`), the origin frames quoted in §2.1.
