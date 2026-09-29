# Error-path hygiene cluster: DR-165, DR-167, DR-175 landed as three commits behind one plan jury, one code jury and a valgrind gate with a positive control

**Date:** 2026-09-29
**Purpose:** close the three register rows of one class, a refusal that throws leaves what the refused
path had already built: the HDF5 loaders' dataset handles (DR-165), the periodicity factory's centroid
nodes and the object under construction (DR-167), the material and boundary-condition factories'
products (DR-175). Two mechanisms only: check before allocate where the predicate reads inputs, a scope
owner at setup scope where it does not.
**Module:** `src/io` (`hdf5_tools`, `cl_HDF5`), `src/mesh` (`cl_Mesh_PeriodicityFactory`),
`src/physics/materials` (`cl_MaterialFactory`, `cl_Material_Alloy`), `src/fem/maxwell` and
`src/fem/thermal` (the two boundary-condition factories), `tests/io`, `tests/mesh`, `tests/physics`, `tests/fem`
**AIs involved:** Claude Fable 5.1 (plan, orchestration, the fix step, gates, commits); Codex `gpt-6-astra`/high
and Grok `grok-4.7`/high (plan jury, blind); Codex `gpt-6-astra`/xhigh and Grok `grok-4.7`/xhigh (code jury,
blind, the ownership row of the depth table); Codex `gpt-5.6-terra`/medium (quick check on the fix delta);
Claude Opus (three implementation agents, every diff read by Fable before the code jury)
**Status:** COMPLETE — three commits on `devel` (`cd71661` DR-165, `84df884` DR-167, the DR-175 commit carries
this entry); the three rows struck; DR-190 widened with the `DofManager` residue; DR-191 filed
**Verification:** on a private non-MKL tree against `/opt/scls/debug` (`USE_DEBUG=ON`, `USE_MKL=OFF`, zero
`mkl` symbols): `test_io` 67/67, `test_mesh` 130/130, `test_fem` 299/299, `test_physics` 158 passed + 1
pre-existing skip; `valgrind --leak-check=full` over the refusal filters clean of every BELFEM frame, after
a positive control on a worktree at `HEAD` (`2fc3ba3`) that showed the leaks at the filed sites; the DR-165
object-count assertion red on that `HEAD` worktree (1, 2, 3, 4 open datasets across the refused reads) and
green on the diff. **`make check` in the shared tree, Christian, 2026-09-29: 17/17** (fast 8 tests, 12.63 s·proc; mpi 4 tests, 3.13 s·proc; 21.39 s real). Not run: the other shared tree
**Exchange:** `tmp/ai_exchange/review_error_path_hygiene.md` (plan round, code round, fix step, quick check in
one thread), `review_error_path_hygiene{,_fix,_v2}.diff`

## Summary

Christian named DR-165, DR-167 and DR-175 as the next register cluster. The plan
(`todo/closed/error_path_hygiene_cluster.md`) found one clause already closed: the `JcFunctionDatabase`
window guard had held its `Database` in a `unique_ptr` since `bc1ca8c` (2026-09-23), the day before the
row was widened with it. The register's fix sketch for DR-165, "check the class before `H5Dget_type`",
is impossible, since the file type IS that handle; the fix is an owner. Five rulings were taken by quiz:
the second layer of the three factories is a `try`/`catch`-rethrow in the constructor (O2, precedent
`UserDefinedMaterial`); the alloy's post-component refusals are in scope (O6); the half-built
`DofManager` of a thin-shell selection refusal is out (O4); Claude builds the non-MKL scratchpad tree
(O5); the code jury runs at xhigh (O7, the protocol's ownership row).

## What landed

**DR-165 (`cd71661`).** `hdf5::DatasetHandles`, a non-copyable owner of the dataset, datatype and
dataspace identifiers, assigned at the call that opens each and closed at scope exit, in the eight
loaders and in `save_strings_to_file`, the one saver that refuses with handles open (Codex 6, Grok 1 of
the plan round). The strings loader reads into a null-initialized `Cell< char * >` and reclaims the
variable-length payload on both paths before the refusal, since a failed `H5Dread` can leave elements
allocated (Codex probed the installed library). The array, vector and strings loaders refuse a dataset of
another rank before reading its extents into their one-element buffer, the `get_array_size` shape (Grok
2 of the code round; pre-existing stack overwrite, widened in by ruling). Gate: `HDF5::number_of_open_objects()`
and the static `number_of_open_datatypes()`. The first draft counted datatypes and dataspaces with
`H5Inmembers`; Codex's code-round P1 said the library refuses its own identifier types, and a 30-line
probe against the linked 1.14.6 confirmed it (−1 for both). The supported census is
`H5Fget_obj_count( H5F_OBJ_ALL, H5F_OBJ_DATATYPE )` (probe: 0, 1 after `H5Dget_type`, 0 after the
close); dataspaces have no census and the test says so.

**DR-167 (`84df884`).** A file-local `TemporaryNodes` owner for the centroid nodes of `map_facets`,
`unique_ptr` for the `Periodicity` in `create_periodicity` and `from_proto`. Grok corrected the plan:
`from_proto` never enters `update_periodicity`; its refusals are `create_hesse` and the `Map` misses of
the proto lookups. The `Periodicity` destructor is safe on a half-built object because each reset is
gated on the flag its successful stage sets.

**DR-175 (this commit).** `Alloy::set_components` checks the RRR and every component name
(`MaterialFactory::is_pure_metal_label`, case-insensitive; named so because a member `is_pure_metal`
would hide the file-local overload inside `builtin_kind`, Codex 4 of the plan round) before it allocates;
`create_component` holds its material until the type check passes; a `TensorScratch` owner in
`create_tables`; a `try`/`catch` around the table build that frees the components and the rho database
(the jury made the `mRhoData` residue a fix: `populate_rho_database` deleted it without nulling).
`MaterialFactory` and both boundary-condition factories: the product under construction in a
`unique_ptr`, released into the map by a bound slot reference (C++17 sequences `release()` before a
throwing `operator[]`, Codex 4 / Grok 5), or into `set_function`; the constructor bodies renamed in place
to a private `create_materials` / `create_conditions` and wrapped by a new constructor whose `catch`
deletes what the member container holds. New cases: a stale rho database for another RRR in a scratch
working directory (RRR 5.0 then 5.5, both `Pb40Sn60_RRR5.hdf5`; RRR 7 would have written another file),
a refused second material after an accepted first, and the same shape for both boundary-condition
factories.

## Juries

Plan round: 18 findings, none refuted. Beyond those above: the object count is blind to transient
datatypes (Codex 3, Grok 5); `save_strings_to_file` leaks (both); the alloy refuses after its components
exist (both); Rule of Five on both owners (Codex 5); the whole `DofManager` leaks on every
`create_field` refusal because `create_blocks` runs from `set_equation` after `new` returned (Grok 2,
the plan had the mechanism wrong). Code round: one P1 (the census), four P2, two P3, none refuted; the
fix set F1-F6 applied by Fable and quick-checked (one P3, a comment about what the census counts, fixed).
Every line citation of both rounds was re-opened; my own pre-registration was wrong twice (the pre-fix
count is not 2 but one dataset per refused read, accumulating; the strings loader's close order is not
evidence for the no-order claim).

## Gates

| gate | tree | result |
|---|---|---|
| positive control, `valgrind --leak-check=full` on the refusal filters | worktree at `HEAD`, non-MKL | 33 / 43 / 41 definite records in mesh / physics / fem; frames at `create_periodicity`, `update_periodicity`, `Alloy::create_component`, `Alloy::set_components`, `MaterialFactory::MaterialFactory`, both factories' `SourceFunction` and condition pushes |
| same filters on the diff | working tree, non-MKL | 29 / 29 / 29 records, every frame `Communicator::init_thread` (the Open MPI baseline) |
| the DR-165 assertion at `HEAD` with only the accessors and the test | worktree at `HEAD` | red: open datasets 1, 2, 3, 4 and datatypes 1, 2, 3 across the refused reads |
| the four suites natively on the diff | working tree, non-MKL | io 67/67, mesh 130/130, fem 299/299, physics 158 + 1 pre-existing skip |
| `make check` | shared, Christian | 17/17, fast 8 + mpi 4, 21.39 s real |

The new rho-database case costs 6.7 s in that debug build; `test_physics` carries no `fast` label, so
`check-fast` is untouched.

## Open

- **DR-191** (new, `[P]`): the diff's physics run under valgrind reports a conditional jump on an
  uninitialised value in `Lead::kohler` and `WhiteTin::kohler`, reached from `Alloy::populate_rho_database`
  (246079 hits), plus one context inside `dmumps_` under the rho-database projection. The `HEAD` run had
  none because nothing built tables there; every deck with `Pb38Sn62` builds them.
- **DR-190** widened: every refusal on the `Kernel::create_field` path leaks the whole `DofManager`
  (`ThinShellSelection.*` is outside the valgrind filter for that reason).
- Residues recorded in the plan: the component loop of `set_components` is outside the `try` (by
  design, only `bad_alloc` or a pure-metal constructor can throw there); partial `set_periodic` links on a
  refused mesh (O3); the unknown-component message changed.

## Files

`src/io/hdf5_tools.hpp/.cpp`, `cl_HDF5.hpp/.cpp`; `src/mesh/cl_Mesh_PeriodicityFactory.cpp`;
`src/physics/materials/cl_MaterialFactory.cpp/.hpp`, `cl_Material_Alloy.cpp`;
`src/fem/maxwell/cl_MaxwellBoundaryConditionFactory.cpp/.hpp`, `src/fem/thermal/cl_ThermalBoundaryConditionFactory.cpp/.hpp`;
`tests/io/test_HDF5.cpp`, `tests/physics/test_MaterialFactoryDeck.cpp`, `tests/fem/test_AdmissionRefusals.cpp`;
`todo/debt_register.md` (live 10, archive 176), `todo/debt_register_closed.md`,
`todo/closed/error_path_hygiene_cluster.md`, `todo/README.md`, `devlog/README.md`.
