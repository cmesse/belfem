# Error-Path Hygiene Cluster: DR-165, DR-167, DR-175 as One Plan Round

**Date:** 2026-09-29
**Purpose:** Three debt rows of one class, a refusal that throws in throw mode leaves what the
refused path had already built: the HDF5 loaders leak their dataset handles (DR-165), the
periodicity factory leaks its centroid nodes and the object under construction (DR-167), and
the material and boundary-condition factories leak the material, condition or source function
built before the refusal (DR-175). Each row keeps its own commit; the cluster shares the plan jury
and the code jury because the diffs are small, disjoint, and use the same two mechanisms: check
before allocating where the check needs nothing built, a scope owner where it does.
**Module:** `src/io` (DR-165), `src/mesh` (DR-167), `src/physics/materials`, `src/fem/maxwell`,
`src/fem/thermal` (DR-175); `tests/io`, `tests/physics`, `tests/fem`
**AIs involved:** Claude Fable 5.1 (plan, orchestration, commits), Codex `gpt-6-astra` + Grok
`grok-4.7` (plan and code juries; the astra row because DR-167 is `src/mesh`), Opus (implementation)
**Status:** ✅ COMPLETE (2026-09-29). Three commits on `devel` (`cd71661` DR-165, `84df884` DR-167, the
DR-175 commit carrying `dl20260929_error_path_hygiene_cluster`) behind one plan jury (astra/high +
grok-4.7/high, 18 findings, none refuted), one code jury (astra/xhigh + grok-4.7/xhigh, one P1 settled by a
probe, four P2, two P3) and a terra/medium quick check on the fix delta F1-F6. Verified on a private
non-MKL tree against `/opt/scls/debug`: valgrind clean of every BELFEM frame on the three refusal filters
after a positive control at `HEAD` showed the leaks; the DR-165 assertion red at `HEAD` and green on the
fix; io 67/67, mesh 130/130, fem 299/299, physics 158 + 1 pre-existing skip. Owed: `make check` in the
shared tree (Christian). Residue: DR-190 widened with the `DofManager` ownership gap; DR-191 filed for the
uninitialised read in `Lead::kohler` / `WhiteTin::kohler` the new table-building case exposed

> **Scope guards:**
> - Three rows, three commits, one plan jury, one code jury on the own-files diff of all three.
> - No behavior change on any accepted path: every accepted read, mesh and deck produces what it
>   produces today. Only what a REFUSED path leaves behind changes. No `input.conf` key moves.
> - The leaks matter only in throw mode. `gThrowOnError` defaults to `BELFEM_ASSERTIONS_ACTIVE`
>   (`src/core/assert.cpp:310`), so every DEBUG build throws by default and every serial test
>   `main` sets it besides; a release build aborts inside `error()` before any unwind, so nothing
>   here changes a production run.
> - The cohomology core is not touched. `src/fem/kernel` is not touched: the `DofManager` that a
>   thin-shell selection refusal leaves half-built (`Kernel::create_field`, see O4) is kernel work
>   and stays out.
> - Of the save paths only `save_strings_to_file` is converted: it refuses at `:891` and `:907`
>   with its dataspace, type and dataset open (Codex 6, Grok 1). The scalar, bool, string, array,
>   vector and matrix savers check before `H5Dcreate` and assert after their closes, so they are
>   left alone.

---

## 1. Current Behavior and How It Fails

| Row | Site | What is built before the refusal | Refusal | Evidence |
|---|---|---|---|---|
| DR-165 | `hdf5_tools.hpp` `load_scalar_from_file` `:306-312` | `H5Dopen1` handle, `H5Dget_type` handle | `check_read_datatype` (`:206-216`, always-active) throws with both open; the closes sit at `:324-326` | source trace |
| DR-165 | `load_array_from_file` `:462-476` | dataset, type, then dataspace | class check `:468`; the length check `:476` (`tDims[0] == aLength`) throws with THREE handles open | source trace |
| DR-165 | `load_vector_from_file` `:597-603`, `load_matrix_from_file` `:776-782`, `load_bool_from_file` (`hdf5_tools.cpp:87-93`) | dataset, type | class check | source trace |
| DR-165 | `load_strings_from_file` `:927-969` | dataset, dataspace (`:927-929`), then a copied string type, then `new char *[ n ]` (`:952`) | the extent refusal `:932` leaks dataset and dataspace; the `H5Dread` refusal `:954` leaks three handles, the pointer array, and whatever variable-length strings HDF5 allocated before it failed (Codex probed the installed 1.14.6: a read that fails on the third element leaves two allocated) | source trace; the partial-allocation claim is Codex's probe, not re-run here |
| DR-165 | `save_strings_to_file` `:877-911` | dataspace, type, dataset | `BELFEM_ERROR` at `:891` (dataset creation) and `:907` (write) before the closes at `:909-911` | source trace (Codex 6, Grok 1) |
| DR-165 | `load_string_from_file` (`hdf5_tools.cpp:180-221`) | dataset, dataspace, type | NO refusal between open and close; converted for uniformity only, not as a leak fix (Grok 6) | source trace |
| DR-165 | `get_array_size` (`hdf5_tools.cpp:232-249`) | dataset, dataspace | rank check `:245` | source trace |
| DR-165 | consequence | HDF5's default close degree is WEAK: `H5Fclose` (`cl_HDF5.cpp:128`) returns but the file stays open inside the library until the dangling handles close, so a process that refuses a read and then opens the same path with `FileMode::NEW` (`H5F_ACC_TRUNC`) hits a still-open file | HDF5 documented semantics; the test suite never re-opens a refused path, so no case shows it today (medium) |
| DR-167 | `PeriodicityFactory::map_facets` `:465-566` | one centroid `Node` per facet on both planes (`create_temporary_nodes` `:1517`), kd-tree child containers inside them (`kdsort` `:1626-1631`, freed by `Node::~Node` → `delete_containers`) | count mismatch `:513`, no partner within tolerance `:532`; the deletion loops are at `:558-566` | source trace; `Periodicity.FacetPerturbationIsRefused` exercises `:532` |
| DR-167 | `create_periodicity` `:131` | `new Periodicity( mMesh )` | anything `update_periodicity` refuses: NaN planes `:338-339`, the two above, `match_nodes` `:746`, `match_edges` `:1262`, `:1285` | source trace; `Periodicity.NodeMismatchIsRefusedInEveryBuild` exercises `:746` |
| DR-167 | `from_proto` `:241-332` | `new Periodicity( mMesh )` | it never enters `update_periodicity` (Grok 4); its refusals are `create_hesse` (`:1731`) and every `mProtoMesh->node/edge/face/facet( id )` lookup, which goes through `Map::operator()` and refuses a missing id (`cl_Map.hpp:224-232`) | source trace |
| DR-175 | `Alloy::set_components( names, X, RRR )` `:125-199` | `delete_components` then one `Metal` per component through a nested `MaterialFactory` (`create_component` `:236-256`) | RRR gate `:184`, AFTER the components exist; the alloy constructor is still running, so `~Alloy` (`:75-83`) never frees them | source trace; `AlloyRRR.AlloyWithoutValidRRRIsRefused` |
| DR-175 | same function, after the components: `create_tables()` `:195` and `populate_rho_database()` `:198` | the components, plus `n+8` raw `Tensor` objects in a local `Cell` inside `create_tables` (`:286-289`, deleted at `:414-418`) | `self_consistent_KG` refuses on non-convergence (`:720`) with the tensors alive; `load_rho_database` refuses a stale database in the cwd whose label or RRR differs (`:758`, `:763`), a deck-reachable path (rerun with another `RRR` in the same directory) | source trace (Codex 2, Grok 3); no test reaches either |
| DR-175 | `Alloy::create_component` `:243` | the `Material` just created | "not a pure metal" throws holding it; and the components built for earlier names stay in `mComponents` | source trace; no test reaches it |
| DR-175 | `MaterialFactory::MaterialFactory` `:116-521` | the material of the current section (`new material::Iron` at `:154`, `create_material` at `:231` and `:320`; inserts at `:175` then `continue` at `:176`, and `:521`) | every `BELFEM_ERROR` between creation and the map insert: missing jc/n `:246`, defect keys `:266-267`, `:327-328`, resistivity type `:290`, `:351`, riva without n `:361`, the three `critical temperature` refusals `:406-418`, heating keys `:435-436`, the allow-list refusals in `check_unused_input`, the duplicate label `:516` | source trace; `AllowList.UnusedMaterialKeysAreRefused` (three decks), `CurveShape.KnownNonIronBuiltinBesideACurveIsRefused` fires BEFORE creation and is clean |
| DR-175 | same constructor, second layer | the materials of the sections BEFORE the refused one, already in `mMaterialsMap` | a constructor that throws runs no destructor, and `Map< string, Material * >` does not delete pointees; the factory has `~MaterialFactory() = default` (`.hpp:95`) because the consumer takes the map by copy (`cl_MaxwellFactory.cpp:2781`) | source trace; no shipped test has two materials with the second refused |
| DR-175 | `MaxwellBoundaryConditionFactory` `:375-560`, `ThermalBoundaryConditionFactory` `:165-359` | `new PhysicalBoundaryCondition` pushed into the member Cell (Maxwell `:72`, `:182`, `:197`; thermal `:73`, `:82`); `new SourceFunction()` at `:375` (`:165` thermal) | the `expk` refusal `:397` (`:183`), every key refusal after it, the two `default` refusals; `set_function` is at `:560` (`:359`) | source trace; `AdmissionRefusals.ExpkOn{Maxwell,Thermal}BcIsRefused`, whose header comment states the leak (`test_AdmissionRefusals.cpp:16-20`) |
| DR-175 | `JcFunctionDatabase` window guard | its `Database` | **ALREADY CLOSED**: `open_checked` holds the `Database` in a `std::unique_ptr` across the dimension, window, slope and floor refusals and releases it into the member (`cl_JcFunction_Database.hpp:170-211`), landed in `bc1ca8c` on 2026-09-23, the day before the row was widened with the Tier 2 finding. No code; the clause is struck in the register with that SHA | source trace of the current file and of `git show bc1ca8c` |

**Bottom line:** every site is a refusal placed after an allocation it does not own yet. Where
the refusal needs nothing that was built, it moves before the allocation (the alloy's RRR gate and
its component names). Where it needs the built thing (a dataset's type, a matched node set, a
material's properties), the built thing gets a scope owner whose destructor runs on the throw.

## 2. Mechanisms

Two, and only two, so the reader of any site recognizes the shape:

1. **Check before allocate.** The DR-171 K-1 pattern (`cl_MaterialFactory.cpp:213-228`: the RRR
   contradiction is checked on the name, before `create_material`). Used where the predicate reads
   only inputs.
2. **A scope owner at setup scope.** `std::unique_ptr` is sanctioned outside the assembly, solve
   and communication loops by the rule of thumb at `doc/coding_philosophy.md:524` (precedents
   `cl_CurveFactory.cpp:613`, `:729`, `cl_MaxwellFactory.cpp:3151`, `cl_JcFunction_Database.hpp:173`).
   For a resource that is not a single heap object (three HDF5 identifiers, a Cell of nodes, a
   Cell of tensors) a small named owner with a destructor, the same idea spelled out. Every such
   owner declares or deletes all five special members (`coding_philosophy.md:509`, Codex 5, Grok):
   copy and move deleted, since a copy would double-close or double-delete. The named owner also
   keeps the census sentence at `:518` true: `src/io` still holds no smart pointer.

For the **second layer**, the objects a factory constructor has already handed to a member
container when a later section is refused, the destructor of the half-built factory never runs.
Two house precedents exist: `ElectricalCircuitFactory` holds its products in `unique_ptr`
MEMBERS for exactly this reason (`cl_ElectricalCircuitFactory.hpp:39-45`, `.cpp:54-56`), and
`UserDefinedMaterial` wraps the rest of its constructor in `try { … } catch ( ... ) { dlclose;
throw ; }` (`cl_Material_UserDefined.cpp:48-62`). The member form would change the type the
consumers read (`materials()` returns the raw map, `boundary_conditions()` a `Cell< PBC * > &`
that `ThermalFactory` iterates three times after the handover, `:211-256`), so this plan takes the
`try`/`catch`-rethrow form for the three factories, with the comment that the block is live in
throw mode only. **O2 puts that choice to Christian.**

Rejected: a Cell that deletes its pointees (an owning `Cell< T * >` contradicts the container
contract `cl_Cell.hpp:34-36` and every consumer); an `mOwn` flag on the factory (the factory
principle in CLAUDE.md forbids it).

## 3. Consumer Inventory

| Row | Consumers of the changed path | Impact |
|---|---|---|
| DR-165 | every `HDF5::load_data` overload (`cl_HDF5.hpp:162-386`) forwards to a loader; callers are the mesh reader, `Database`, the restart loader, `load_bh_curve`. The three closes move from explicit lines to the owner's destructor at scope exit, the same point in the accepted path | none on an accepted read (high). The close ORDER changes to dataspace, datatype, dataset; HDF5 imposes no order between these three identifiers (`H5Dclose` does not require its derived dataspace or datatype to be closed first, and the strings loader already closes in that order at `:967-969`) |
| DR-165 | `OPEN_RDONLY_PARALLEL` | same loaders; unchanged (high) |
| DR-167 | `CutFactory` `:131` and `ProtoMesh::from_proto` `:1464` receive the pointer through `Mesh::set_periodicity`; `Periodicity::update()` (`cl_Mesh_Periodicity.cpp:38-46`) runs `update_periodicity` on an existing object and allocates nothing | none on an accepted mesh (high) |
| DR-175 | `MaxwellFactory::create_materials` copies the map (`:2781`); `Alloy::create_component` builds a nested `MaterialFactory` with the default constructor (`:240`); the test helper `refusal()` deletes the map on acceptance (`test_MaterialFactoryDeck.cpp`) | none on an accepted deck (high); `materials()` keeps its signature |
| DR-175 | `ThermalFactory` and `MaxwellFactory` hand `boundary_conditions()` to `Kernel::add_boundary_condition` and keep iterating the factory's Cell (`cl_ThermalFactory.cpp:211-256`) | none; the Cell and the accessor are unchanged (high) |
| DR-175 | the three shipped alloy sections (`examples/pancake`, `tapestack3d`, `tapestack3d_quench`, `Pb38Sn62` with `RRR : 1.5`) pass the hoisted RRR gate and the name check exactly as they pass today | none (high) |

## 4. Ordered Steps

- [ ] **R1 (DR-165) — the handle owner.** In `hdf5_tools.hpp`, inside `#ifdef BELFEM_HDF5`, a
  non-copyable `struct DatasetHandles { hid_t mDataSet = H5I_INVALID_HID ; hid_t mDataType = … ;
  hid_t mDataSpace = … ; ~DatasetHandles() }` whose destructor closes each member that is
  `>= 0` (`H5Sclose`, `H5Tclose`, `H5Dclose`). `H5I_INVALID_HID` is `(-1)` in the linked
  1.14.6 (`/opt/scls/mkl/include/H5Ipublic.h:75`). Every loader assigns its three
  identifiers into one local owner and loses its explicit close lines: `load_scalar_from_file`,
  `load_array_from_file`, `load_vector_from_file`, `load_matrix_from_file`,
  `load_strings_from_file` (hpp) and `load_bool_from_file`, `load_string_from_file`,
  `get_array_size` (cpp), and `save_strings_to_file` (hpp, the one saver that refuses with handles
  open). Every identifier is assigned into the owner at the call that opens it, so the extent
  refusal at `:932` is covered as well as the read refusal. The strings loader's `new char *[ n ]`
  becomes a zero-initialized `Cell< char * > tBuffer( tNumStrings, nullptr )` passed as
  `tBuffer.data()`, the read status is kept in a local, and `H5Dvlen_reclaim` runs on BOTH paths
  before the refusal is raised: a failed read can leave elements allocated (Codex 1), and over a
  null-initialized buffer the reclaim frees exactly those. The owner closes the type and the
  dataspace only after the reclaim, which needs both. The duplicate `BELFEM_ASSERT` class check in
  `load_matrix_from_file` `:784-788`, which restates the always-on `check_read_datatype` at
  `:782`, is deleted in the same hunk.
- [ ] **R2 (DR-165) — the wrapper accessor and the gate.** `HDF5::number_of_open_objects()` in
  `cl_HDF5.hpp/.cpp`: `H5Fget_obj_count( mFile, H5F_OBJ_DATASET | H5F_OBJ_GROUP | H5F_OBJ_DATATYPE
  | H5F_OBJ_ATTR )`, documented as counting the objects open IN this file other than the file
  itself (dataspaces are not file objects and are not counted; a selected group counts as one).
  **That count sees datasets only** (Codex 3, Grok 5): a transient datatype from `H5Dget_type`
  is counted by `H5Fget_obj_count` only with `H5F_OBJ_ALL` as the file argument
  (`/opt/scls/mkl/include/H5Fpublic.h:786-796`), and a dataspace is not a file object at all. So
  the wrapper gains two more static accessors over `H5Inmembers` (`H5Ipublic.h:640`):
  `HDF5::number_of_open_datatypes()` and `HDF5::number_of_open_dataspaces()`, library-wide
  counts whose absolute value includes the predefined types, so the test reads them as DELTAS
  across a refused read. The wrapper is extended rather than the test calling the C API because
  L-21 is a hard rule for library calls; the wrapper-policy checker scans `src/` only
  (`scripts/check_wrapper_policy.py:207`) and would not have caught a direct call in the test,
  which is not a reason (Grok 10). In `tests/io/test_HDF5.cpp`,
  `ClassMismatchIsRefusedInEveryBuild` asserts `number_of_open_objects() == 0` and zero deltas of
  the other two after each of its five refused reads, and a new `RefusedReadsLeaveNoHandleBehind`
  covers the refusals the class check does not: the raw-array length mismatch through
  `load_data( label, real *, n )` (`cl_HDF5.hpp:386-396`) and the rank refusal of
  `get_array_size` through `load_data( label, Cell< T > & )` on a 2-D dataset
  (`cl_HDF5.hpp:283-296`, a live path, Grok 10). **Discriminating direction:** on the pre-fix
  binary each refused read leaves one dataset open on the still-open file, so the count is at
  least 1 and grows with every refusal; the assertion is red there. R7(b) records the observed
  values instead of asserting a predicted number (Grok 5).
- [ ] **R3 (DR-167) — the centroid-node owner and the object owner.** In
  `cl_Mesh_PeriodicityFactory.cpp`, a file-local `struct TemporaryNodes { Cell< Node * > mNodes ;
  ~TemporaryNodes() }` deleting its pointees (null slots included: `create_temporary_nodes` sizes
  the Cell with `nullptr` before filling it, `:1497-1519`; `delete nullptr` is a no-op), copy and
  move deleted; `map_facets` declares `tSourceNodes` and `tTargetNodes` as two of them, passes
  `.mNodes` to `create_temporary_nodes`, `create_kdtree` and the loops, and loses the deletion
  loops at `:558-566`. The in-place sort at `:540` is harmless to the owner. `create_periodicity`
  and `from_proto` hold their `Periodicity` in a `std::unique_ptr` and `return tOwner.release()`;
  in `from_proto` the owner covers the `create_hesse` and `Map` refusals named in §1, not
  `update_periodicity`, which that path never calls. The destructor of a half-built `Periodicity`
  is safe: each `reset_*` is gated on the entity flag its successful stage sets
  (`cl_Mesh_Periodicity.cpp:173-251`), the backup clear is a `Cell::clear()`.
  `Periodicity::update()` (`cl_Mesh_Periodicity.cpp:38-46`) runs on an object the mesh already
  owns and deletes (`cl_Mesh.cpp:196-200`), so only the centroid owner applies there. Comment above
  the owners: they exist for the refusal path; the accepted path releases or deletes at the same
  point as before.
- [ ] **R4 (DR-175) — the alloy checks before it allocates.** In `Alloy::set_components( names,
  X, RRR )` (`:125-197`): the RRR resolution and the `> 1` gate (`:176-186`) move to the top of
  the function, before `delete_components()` and `mComponents.set_size`; both read only the
  argument and the constant already set by the pair overload. Stated consequence on a refused
  programmatic RE-call on a live alloy: the previous component set now survives the refusal
  instead of being deleted first (Grok), the stronger guarantee. Then every component NAME is
  checked before the loop through a new public static
  `MaterialFactory::is_pure_metal_label( const string & )`, the string form of the file-local
  `is_pure_metal( builtin_name( … ) )` that `builtin_kind` uses (`cl_MaterialFactory.cpp:689`). The
  name differs from the file-local predicate on purpose: a member named `is_pure_metal` would hide
  the anonymous-namespace overload inside every member function and break the `BuiltinName` call
  at `:689` (Codex 4, confirmed by a syntax probe). The message names the alloy and the offending
  component. `create_component`'s own type check (`:243`) stays as the contract of a function
  other callers may reach, and its created `Material` is held in a `unique_ptr` until the check
  passes. Result: the RRR and name refusals allocate nothing. The refusals AFTER the components
  (`create_tables` and the rho-database reload, §1) are **O6**.
- [ ] **R4b (DR-175, pending O6) — the alloy's post-component refusals.** A file-local owner for
  the tensor scratch of `create_tables` (`Cell< Tensor< real > * >` `:286-289`, the `TemporaryNodes`
  shape; the deletion loop at `:414-418` goes), and `set_components( names, X, RRR )` wraps
  `create_tables()` and `populate_rho_database()` in `try { … } catch ( ... ) {
  this->delete_components() ; throw ; }` so a refusal there leaves the alloy without components,
  the state `~Alloy` would produce. No test reaches either refusal; the rho-database mismatch is
  the deck-reachable one and gets a case: build an alloy with tables into a scratch directory,
  then build the same label with another RRR in that directory and expect the `:763` refusal
  (under valgrind, it is the leak witness).
- [ ] **R5 (DR-175) — the material factory's two layers.** (i) In the constructor loop the
  current material is a `std::unique_ptr< Material > tMat`, assigned at the three creation sites
  (`:154`, `:231`, `:320`) and released into the map at the two inserts (`:175`, before the
  `continue` at `:176`; `:521`); `tMat->` and `tMat.get()` replace the raw uses. (ii) The loop
  body moves to a private `create_materials( aSection, aBuildTables )`; the constructor is
  `try { this->create_materials( … ) ; } catch ( ... ) { delete every pointee of mMaterialsMap ;
  throw ; }` with the comment that a throwing constructor runs no destructor and the block is
  reached in throw mode only (release aborts inside `BELFEM_ERROR`). The two layers do not
  overlap (Grok): on a throw the local owner dies first and frees the current material, which is
  not in the map yet; the `catch` deletes only what the map holds. Precedent
  `cl_Material_UserDefined.cpp:41-57`. **Pending O2.** A new case
  `AllowList.ARefusedSecondMaterialFreesTheFirst` writes a deck with an accepted copper followed by
  a refused ybco; it can only turn red under valgrind, as the DR-169 test does, and says so.
- [ ] **R6 (DR-175) — the boundary-condition factories.** In both constructors `tFunction`
  becomes `std::unique_ptr< SourceFunction >`, and `release()` IS the argument of `set_function`
  (`cl_FEM_PhysicalBoundaryCondition.cpp:43-46` is a bare store, so the store must be the
  handover, Grok); the loop body is wrapped as in R5(ii) with `catch ( ... ) { delete every
  pointee of mPhysicalBoundaryConditions ; throw ; }` (`~PhysicalBoundaryCondition` deletes a
  null `mFunction` safely, `.hpp:52`, `.cpp:27-30`; a condition whose function was already stored
  deletes it). Two new cases in `test_AdmissionRefusals.cpp`, one per factory: an accepted first
  section followed by a refused second one, so the `catch` runs over a condition that already
  owns its function (Codex 7). The stale sentence in `test_AdmissionRefusals.cpp:16-20` ("On a
  refused deck the condition and its source function created before the throw leak") is rewritten
  to state the new contract. **Pending O2.**
- [ ] **R7 — gates (Christian runs builds and binaries).** (a) `make check` in
  `cmake-build-debug` (release and debug configurations as the previous cluster did; the
  Tier 2 MPI suites are in `check`). (b) The R2 assertions once against a binary WITH the R2
  accessors and WITHOUT the R1 owner (the accessors are what the test calls, so `src/io` is not
  stashed whole, Codex 3), the observed counts recorded in the devlog, then green with R1. (c) A non-MKL tree
  against `/opt/scls/debug` (the 2026-09-23 precedent: a private out-of-tree build in the
  scratchpad, 22 min on one core), `valgrind --leak-check=full` over `test_mesh
  --gtest_filter='Periodicity.*Refused*'`, `test_physics --gtest_filter='AlloyRRR.*:AllowList.*:
  CurveShape.*:ConstantLaws.*'` and `test_fem --gtest_filter='AdmissionRefusals.*'`, FIRST on
  `HEAD` as the positive control (the pre-existing records with `PeriodicityFactory`, `Alloy`,
  `MaterialFactory` and `*BoundaryConditionFactory` frames must appear), then on the diff, where
  no leak record may carry a frame in those units; the R4b and R6 two-section cases are in the
  filters. The gate is "no BELFEM-owned allocation in any
  leak record", never `ERROR SUMMARY: 0` (Open MPI's `MPI_Init` records are permanent).
  `ThinShellSelection.*` is deliberately outside the filter (O4).

### 4.0 Implementation Progress (updated 2026-09-29)

**Implemented (three Opus agents, every diff read by Fable before the code jury; syntax-gated
with the tree's own `flags.make`, `-Wall -Werror -pedantic-errors`, rc 0 on nine TUs; nothing
built or run at this point):**

- [x] R1 — `DatasetHandles` in `hdf5_tools.hpp` (five special members, closes ids `>= 0`);
  nine sites converted including `save_strings_to_file`; the strings loader reads into a
  null-initialized `Cell< char * >` and reclaims on both paths, guarded on a non-empty dataset
  because an empty Cell may hand out a null `data()` that the reclaim refuses (agent's reading of
  HDF5 1.14.0 `H5T.c:3048`; the tree links 1.14.6). The duplicate matrix assert is gone.
- [x] R2 — `HDF5::number_of_open_objects()` plus the static `number_of_open_datatypes()` and
  `number_of_open_dataspaces()`; `ClassMismatchIsRefusedInEveryBuild` asserts all three after
  each refused read; `RefusedReadsLeaveNoHandleBehind` covers the length and rank refusals.
- [x] R3 — `TemporaryNodes` owner in `map_facets`, `unique_ptr` in `create_periodicity` and
  `from_proto` with `aPeriodicity` kept as a raw alias.
- [x] R4 — RRR gate and the component-name check (`MaterialFactory::is_pure_metal_label`,
  case-insensitive) hoisted above `delete_components()`; `create_component` holds its material in
  a `unique_ptr` until the type check passes.
- [x] R4b — `TensorScratch` owner in `create_tables`; `try`/`catch` around `create_tables` and
  `populate_rho_database` deleting the components; `AlloyRRR.ARhoDatabaseForAnotherRRRIsRefused`
  builds RRR 5.0 then 5.5 in a scratch working directory (both map to `Pb40Sn60_RRR5.hdf5`, since
  the cache name carries the RRR as an integer; RRR 7 would have written another file).
- [x] R5 — `unique_ptr< Material > tMat` released at both inserts; the constructor body renamed
  in place to the private `create_materials`, the new constructor wraps it;
  `AllowList.ARefusedSecondMaterialFreesTheFirst`.
- [x] R6 — both BC factories: `unique_ptr< SourceFunction >` released into `set_function`; the
  bodies renamed in place to the private `create_conditions`; two two-section cases; the stale
  header sentence rewritten.
- [x] R7 — gates. **(c) positive control ran 2026-09-29** on a detached worktree at `HEAD`
  (`2fc3ba3`) built against `/opt/scls/debug` (`USE_MKL=OFF`, `USE_DEBUG=ON`, zero `mkl` symbols
  in the four test binaries), `valgrind --leak-check=full` over the three filters, all cases
  passing:

  | suite | filter | definite-lost records | BELFEM frames in the leak records |
  |---|---|---|---|
  | `test_mesh` | `Periodicity.*Refused*` | 33 | `PeriodicityFactory::create_periodicity` (`:138`), `update_periodicity` (`:353`), both refusal tests |
  | `test_physics` | `AlloyRRR.*:AllowList.*:CurveShape.*:ConstantLaws.*` | 43 | `Alloy::create_component` (`:241`), `Alloy::set_components` (`:159`, `:113`), `Alloy::Alloy` (`:72`), `MaterialFactory::MaterialFactory` (`:231`) |
  | `test_fem` | `AdmissionRefusals.*` | 41 | `SourceFunction::SourceFunction` from both factories (`Maxwell:375`, `Thermal:165`), `PhysicalBoundaryCondition::set_label` from the pushes (`Maxwell:197`, `Thermal:82`) |

  The Open MPI `Communicator::init_thread` records are the permanent baseline. **(c) clean run on the
  diff 2026-09-29:** 29 / 29 / 29 records, every frame `Communicator::init_thread`; the four suites
  natively io 67/67, mesh 130/130, fem 299/299, physics 158 + 1 pre-existing skip. **(b) ran
  2026-09-29:** the `HEAD` worktree with only the accessors and the test rebuilt, `test_io` red with open
  datasets 1, 2, 3, 4 and datatypes 1, 2, 3 across the refused reads; green on the diff. **(a) `make
  check` in the shared tree: owed (Christian).**

**Defects and residues found during implementation:**

- **D1 (LOW, residue, not fixed).** `Alloy::populate_rho_database` deletes `mRhoData` without
  nulling it (`cl_Material_Alloy.cpp:456`); a throw in `save_rho_database` after
  `populate_rho_database_serial` stored the new `Database` would leak it from the constructor,
  and a refused programmatic re-call on a live alloy would leave a dangling pointer for `~Alloy`.
  No test reaches either. Recorded in DR-175's status cell.
- **D2 (LOW, residue, by design).** The component loop in `set_components` is outside the
  `try`; after the name check it can fail only on `bad_alloc` or inside a pure-metal constructor
  built with tables off.
- **D3 (LOW, behavior on a refused deck).** An unknown alloy component is now refused by the
  name check with the alloy-named message instead of `create_material`'s "Material zn is not
  known"; when both the RRR and a component are bad, the RRR refusal fires first. No test
  pinned the old message.
- **Code jury 2026-09-29 (Codex `gpt-6-astra`/xhigh, Grok `grok-4.7`/xhigh, blind):** one P1,
  four P2, two P3, none refuted; the fix set F1-F6 ruled by Christian and applied the same day.
- **D5 (HIGH, fixed 2026-09-29).** `H5Inmembers` refuses the library-defined identifier types
  (Codex 1; confirmed by a probe against the linked 1.14.6: −1 for both), so the census test
  would have thrown before its guarded read. The datatype census is now
  `H5Fget_obj_count( H5F_OBJ_ALL, H5F_OBJ_DATATYPE )` (probe: 0, 1 after `H5Dget_type`, 0 after
  the close); dataspaces have no census in the library and the test says so.
- **D6 (LOW, fixed 2026-09-29, was D1).** `populate_rho_database` nulls `mRhoData` after the
  delete, and the alloy `catch` frees the database (Codex 2, Grok 1).
- **D7 (LOW, fixed 2026-09-29).** The cache test skips without `BELFEM_HDF5` (Codex 3).
- **D8 (LOW, fixed 2026-09-29).** Both map inserts bind the slot reference first and release
  into it, since C++17 sequences `release()` before a throwing `operator[]` (Codex 4, Grok 5).
- **D9 (LOW, fixed 2026-09-29).** `number_of_open_objects` refuses a closed file (`mFile` is
  unassigned on a non-root rank in serial mode) and its comment states the root and stack
  counts correctly (Grok 3, 4).
- **D10 (HIGH, pre-existing, fixed 2026-09-29 as the F6 widening).** The array, vector and
  strings loaders read the extents into a one-element buffer with no rank check, so a dataset of
  another rank overwrote the stack before any owner ran; each now refuses `rank != 1` the way
  `get_array_size` does (Grok 2).
- **D4 (closed 2026-09-29).** `ARhoDatabaseForAnotherRRRIsRefused` costs 6.7 s in a non-MKL debug
  build; `test_physics` carries no `fast` label, so `check-fast` is untouched. Under valgrind it
  exposed DR-191.

## 5. Open Design Questions (not silently decided)

- **O1 — RESOLVED 2026-09-29 by the dispatch.** One plan jury at `gpt-6-astra`/high +
  `grok-4.7`/high; DR-167 is `src/mesh` and DR-175 touches a `src/fem/maxwell` factory, both on
  the astra row (`doc/ai_collaboration_protocol.md:466`), so the premise that DR-175 alone would
  be terra was wrong (Grok 7). The code jury is **O7**.
- **O2 — RESOLVED 2026-09-29 (Christian): the `try`/`catch`-rethrow form.** The second-layer mechanism in the three factories. `try`/`catch`-rethrow in the
  constructor (R5(ii), R6; precedent `UserDefinedMaterial`), or `unique_ptr` members with the
  consumers' signatures kept by a raw view (precedent `ElectricalCircuitFactory`), or leave the
  second layer open (only the current object gets its owner; the k−1 earlier ones still leak when
  the k-th section is refused, which no shipped test exercises). Recommendation: the
  `try`/`catch` form; it is six lines per factory, it is explicit, and the valgrind gate then has
  no known residue to explain.
- **O3 — partial entity links after a match refusal.** `match_nodes` sets `set_periodic`,
  flags and indexes the pairs it found (`:731-736`) before it refuses at `:746`, and `match_edges`
  does the same (`:1135-1136`, `:1240-1241`) before `:1262` and `:1285`; the `Periodicity`
  destructor resets an entity type only under its flag, which the caller sets after a SUCCESSFUL
  match (`:390-396`, `:423-431`), so a refused mesh keeps those links (Grok). The
  `create_periodicity` callers stop on the refusal and the fixtures destroy their mesh; a caught
  `Periodicity::update()` leaves them on a live mesh. Recommendation: record in the row as a
  known residue, no code in this round.
- **O4 — RESOLVED 2026-09-29 (Christian): out of this round; recorded in DR-175's status cell.** The whole `DofManager` on a thin-shell selection refusal. Mechanism corrected
  (Grok 2): the `DofManager` constructor only allocates its sub-objects
  (`cl_FEM_DofManager.cpp:31-54`); `create_blocks` runs from `set_equation` (`:99`), which
  `Kernel::create_field` calls at `cl_FEM_Kernel.cpp:603` AFTER `new` returned (`:591`) and after
  the `check_mesh` refusal (`:597`). Every refusal on that path therefore leaks the entire object,
  not yet pushed (`:605`). A `unique_ptr` in `create_field` released at the push is three lines,
  but whether `~DofManager` is safe on an object whose `set_equation` half-ran is not
  established, and the file is `src/fem/kernel`. Recommendation: out of this round; recorded in
  DR-175's status cell as the residue that keeps `ThinShellSelection.*` out of the valgrind
  filter, or as its own `[W]` row, Christian's choice.
- **O6 — RESOLVED 2026-09-29 (Christian): R4b is in.** The alloy's post-component refusals (R4b). Include R4b in this round (a tensor
  owner plus a `try`/`catch` in `set_components`, one new case), or record the two refusals as a
  residue in the row. Recommendation: include; the rho-database reload refusal is deck-reachable
  and the mechanisms are the round's own two.
- **O7 — RESOLVED 2026-09-29 (Christian): astra/xhigh + grok-4.7/xhigh.** The code jury's depth. The code round is round 2 on an ownership subject, which the
  protocol puts at `xhigh` for both seats (`:465`, astra for the mesh and maxwell files, `:466`);
  the previous cluster ran its code jury at astra/high. Recommendation: astra/xhigh + grok
  4.7/xhigh, cost noted, or high if Christian prefers the previous cluster's depth.
- **O5 — RESOLVED 2026-09-29 (Christian): Claude configures and builds the scratchpad tree against `/opt/scls/debug`; the positive control is built from a git worktree at `HEAD` so the working-tree edits cannot reach it.** The valgrind tree. No non-MKL tree exists in `belfem/` (both `cmake-build-debug` and
  `build` are MKL+debug); the 2026-09-23 run used a scratchpad build against `/opt/scls/debug`
  that this session would configure and build on Christian's word. Recommendation: the same.

## 6. Definition-of-Done Checklist

- [x] Every row in §1 is mapped to a step or an open question.
- [x] Each claimed leak is backed by a citation, not assumption.
- [x] R7(b): the DR-165 assertion has been seen red on the pre-fix binary once.
- [x] R7(c): the valgrind positive control ran on `HEAD`, then clean on the diff.
- [x] Register: DR-165, DR-167, DR-175 struck; the `JcFunctionDatabase` clause closed with `bc1ca8c`;
  the O4 residue on DR-190 (Christian's ruling); DR-191 filed.

## 7. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_error_path_hygiene.md` (plan round, code round, in
  one thread).
