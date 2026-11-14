# Usage-guide jury sweep: all 22 module usage guides audited against the source

**Date:** 2026-10-22
**Purpose:** Overnight staleness sweep requested by Christian: run the blind two-auditor jury over
every `src/**/doc/*usage_guide*.md`, one guide at a time, and verify every finding against the
source before reporting. No guide edited; fixes await approval.
**Module:** all modules with a usage guide (circuit, comm, containers, core, fem/interpolation,
fem/iwg, fem/kernel (dof manager), fem/maxwell, homology, io, linalg (lapack, linalg),
math/graph, math/quaternion, math/tensor, mesh, numerics/spline, physics/database,
physics/gasmodels, physics/gastables, physics/materials, sparse)
**AIs involved:** Claude Opus 5.5 (driver, mechanical pre-registration, verification and
reconciliation via seven verification subagents, spot checks); Codex `gpt-5.6-terra`/high and
Grok `grok-4.7`/high (blind jury, one round per guide)
**Status:** COMPLETE — 22 of 22 guides audited, verified, and fixed (one commit per guide,
3ce9e6c … 31a5a64); source defects handed to a parallel session
**Verification:** source trace for every finding (each citation re-opened at HEAD e9670bc /
76918d0); one compile-only probe (`clang++ -fsyntax-only` against the real linalg headers, both
`-DBELFEM_ARMADILLO` and `-DBELFEM_BLAZE`) for linalg F1/F3/F4; installed `cholmod.h` for sparse
F21. No build, no test, no reproducer. Every finding below is **reviewed**, not verified
**Exchange:** `tmp/ai_exchange/review_<stem>_usage_guide.md` (22 files: pre-registration, CODEX,
GROK, verification, reconciliation); driver and subjects in `tmp/ai_exchange/usage_guide_sweep/`

## Summary

Every usage guide has drifted from the code, most of them badly. The common shape: the guides were
written against an older API, so their examples call methods that no longer exist or have
changed arity, and several examples have ownership errors (double free, use-after-free, leaks)
or MPI recipes that deadlock. About 131 confirmed P1 findings (following the guide gives wrong
behaviour, an abort, or code that does not compile) across 21 guides; `database` has none.
Roughly 25 auditor claims were refuted outright, among them five auditor P1s. The sweep also
surfaced about a dozen defects in the **source**, not the docs; those are routed to Christian.
None of the fixes touches the cohomology core.

## Process

- Driver `usage_guide_sweep/run_sweep.sh`, sequential. Per guide: a mechanical `# CLAUDE`
  pre-registration (`mech_check.py`: every cited path and `file:line` exists), a staleness
  brief (`subject_<stem>.md`: APIs and signatures, behaviour, citations, CMake and input keys,
  omissions; P1/P2 with confidence), then `cross_review.sh --jury`.
- Quota: the first night stopped at `io` with both seats out (Codex usage limit, Grok 402). The
  driver was fixed to skip a guide only when a non-failed CODEX entry exists and to stop when
  both legs fail; the 13 failed files were moved to `failed_quota/`. The rerun from 01:46 audited
  the remaining 13; Codex ran out again on `sparse`, and its leg was rerun blind at 06:48 under a
  temp slug (`rerun_sparse_codex.sh`), so it did not see Grok's entry.
- Verification: seven subagents, each finding labelled CONFIRMED/REFUTED/UNVERIFIABLE with quoted
  evidence, severity re-graded where warranted, reconciliation table appended per file. Claude
  spot-checked at least two P1s per batch against the source; none failed.

## Counts (from the reconciliation tables)

| guide | findings | confirmed | refuted (whole or part) | confirmed P1 |
|---|---|---|---|---|
| circuit | 13 | 13 | 0 | 4 |
| comm | 20 | 16 | 5 | 4 |
| container | 16 | 14 | 2 | 2 |
| core | 32 | 30 | 1 | 9 |
| database | 12 | 10 | 1 | 0 |
| dof_manager | 27 | 26 | 1 | 13 |
| gasmodels | 21 | 20 | 1 | 7 |
| gastables | 23 | 19 | 4 | 6 |
| graph | 28 | 24 | 2 | 8 |
| homology | 26 | 25 | 1 | 18 |
| interpolation | 34 | 32 | 1 | 5 |
| io | 20 | 20 | 0 | 5 |
| iwg | 27 | 26 | 0 | 12 |
| lapack | 13 | 13 | 0 | 3 |
| linalg | 17 | 16 | 1 | 3 |
| materials | 18 | 17 | 1 | 3 |
| maxwell | 18 | 18 | 1 | 10 |
| mesh | 24 | 24 | 0 | 10 |
| quaternion | 9 | 9 | 1 | 1 |
| sparse | 21 | 20 | 0 (4 partial) | 5 |
| spline | 13 | 13 | 0 | 3 |
| tensor | 13 | 12 | 2 | 1 |

Counts are greps over the tables, so they are approximate where a row carries a mixed verdict.

## Confirmed P1s, per guide (detail and evidence in each exchange file)

- **circuit:** hybrid decks take node names, not indices; `CircuitResults.txt` IZ lags one step;
  the pitfall row has the build order wrong (factory before kernel); `userdefined` writes to
  `BC(tCount-i)`.
- **comm:** `send("literal",1)` binds the scalar overload and aborts; `-DUSE_MPI=OFF` fails to
  configure; `omp barrier` inside `omp master`; `gComm` is not defined by the header.
- **container:** GA workflow null dereference and leak; `Cell<int>(100)` only reserves.
- **core:** serial `rand()` can return 1.0; Timer pitfall recommends `stop()`; `filename()` and
  `format_with_leading_zeros` examples wrong (the latter returns a format string,
  `stringtools.cpp:1147-1149`); Example 4 `determinant()`; Example 7 assert fires; Example 1
  `gComm`; `message()` undeclared; `x() +=` on a by-value accessor.
- **dof_manager:** initialize/set_solver pitfall; E/I matrices null by default; D stored negated;
  `SpMatrix::set_value` missing; `owned_elements`, `number_of_integration_points`,
  `compute_jacobian(Element*)`, `owned_nodes`/`aura_nodes`, `InfoLevel::Info` do not exist; the
  `heat_flux_x` layout is interleaved; Pattern 1 use-after-free, Pattern 3 double free; the
  quick start is not a transient.
- **gasmodels:** equilibrium example passes a fresh vector as the element inventory
  (`cl_Gas.cpp:1896`); the formation table is normalized, not atom counts; `remix(x,false,false)`
  leaves state stale; HELMHOLTZ aborts above `mPmax`; the Mach snippet redeclares `M`;
  `dα_dδ`/`dα_dτ` do not exist; a cubic EoS never selects the liquid root.
- **gastables:** reserve-only `Cell` indexed (twice) and `Cell::length()`; the missing-species
  fallback is unreachable; Pattern 5 reaction enthalpy wrong under the per-species H shift;
  `heat_spline->deval` is molar Cp; POLY mode zeroes synthesized transport.
- **graph:** two-pass pattern dereferences null; one-argument `bfs` and the pseudo-peripheral
  search assume indices that are never set; `InfoLevel::Warning` and
  `DofManager::apply_permutation` do not exist; tagged `send`/`receive` arity; ordering wrappers
  turn a leftover flag into a self-loop; METIS debug loop indexes unset owners.
- **homology (18):** `CutFactory` arity and `{}` temporary; rank-0-only workflow deadlocks;
  `mPeriodicFactory` uninitialized; Topology/DomainType, Protoshell, CutProcessor, CutData APIs;
  `thin_shell_facets` cleared by `run()`; `save_debug_meshes`, `add_abstract_dof` absent;
  terminal ids vs indices; `run()` overwrites index and phi; double frees; Smith s/t labels and the
  Betti example. All doc-side or outside the cohomology core.
- **interpolation:** `Map::contains/count`; SideSet/Facet methods; Hermite example at ξ=0;
  `compute_facet_orientation` absent; `cout << Matrix`.
- **io:** `BELFEM_ASSERT` described as aborting in release; `temperature = 300 C` without `;`
  creates no key; `max iterations` vs `max_iterations`; one HDF5 file per thread on a
  non-thread-safe HDF5 (`H5pubconf.h:347`); Pattern 2 solver API does not compile.
- **iwg (12):** `set_omega` before `create_field`; `compute_mkf` on StaticHeat; per-element
  `reset_flags()`; leak and double free; nonexistent IWG types; TransientHeatConduction listed as
  safe with 1.0/111.0 placeholders; Example 4 without a material; `Mesh(comm)`; BDF ramp off by
  one; `rho()` used as density; conductivity tangent; braced temporary in `collect_node_data`.
- **lapack:** `gees` hides `sort='N'` (`fn_gees.hpp:421`), a callback alone does not sort;
  `gesvd` and `geev` take two job characters, so `(...,'A',false)` binds `false` to the second.
- **linalg:** copy-init from an expression fails on Armadillo (probe); Blaze leaves `Vector(3)`
  uninitialized; `Matrix / scalar` does not exist on either backend (probe).
- **materials:** label `"YBCO"` loads the table as an n table (`cl_JcFunction_Database.hpp:148`
  `mIsN = aLabel != "jc"`; guide `:299`); `const` evaluation writes scratch, so not thread-safe;
  Hastelloy χ(T) not routed (μ = μ0).
- **maxwell (10):** `IWG_MaxwellPostproc` arity and `kernel->solve()`; `set_tolerance`;
  `create_facets`/`create_elements`; `Block::set_type`; block selection order; `set_solver` after
  `initialize`; manual cut API; `PartitionMethod::METIS`; the §4.2 driver vs `belfem.cpp`;
  `rho(T,B,beta)` is Kohler, not E-J.
- **mesh (10):** MPI recipes deadlock (rank-0-only constructor); `partition` is a no-op off the
  master (`cl_Mesh.cpp:1185`); order-3 3-D tensor example aborts; LINE5/TRI15 abort `finalize()`;
  connectivity not computed on access; `flag_curved_elements` aborts on pyramids, TET20, HEX64;
  fields lost in `.hdf5`; duplicate nodes leak; edge/face fields zero length; `save_fields` omits
  globals and time.
- **quaternion:** `constant::pi` without `constants.hpp`.
- **sparse:** `set_mumps_blr` before initialize is clobbered (BLR at ε=0); CSR/CSC abort under
  PETSc at >1 rank and CSC under distributed STRUMPACK; STRUMPACK rejects CG, CGS, TFQMR, IBCGS
  although STRUMPACK is the default solver; Pattern 1 and the OpenMP example call
  `build_dof_graph`, `compute_element_stiffness`, `get_element_dofs`, which exist nowhere;
  `rcond()` result read as rcond(A) (see F21 below).
- **spline:** the help matrix must also match the BCs; no `ExtraMode` constructor;
  `update_data` keeps a stale extra mode.
- **tensor:** orthotropic Poisson convention (provisional on Christian).

## Code defects (the source, not the guide) — routed to Christian

Confidence high for each unless stated; all are source traces.

1. `cl_ElectricalCircuitFactory.cpp:65,301-302` — `userdefined` writes to
   `mPhysicalBoundaryConditions( tCount-i )`; `read_circuit` pushes no BC, so i==0 is out of
   range and other i hit the wrong BC.
2. `src/homology/cl_CutFactory.hpp:45` — `mPeriodicFactory` never initialized; latent (only
   non-root ranks skip `set_periodicity`, and the read is master-guarded).
3. `core/random.hpp:92` — `std::rand()/RAND_MAX` can return 1.0 (fix code or doc).
4. `IWG_TransientHeatConduction` — hard-coded 1.0/111.0 placeholders.
5. `src/sparse/fn_rcond.cpp:36` — `stype = 0` makes CHOLMOD factor A·A' (`cholmod.h:845`), so
   `rcond` returns ≈ rcond(A)²; the "LU" comment is wrong too.
6. `src/io/hdf5_tools.hpp:615,797` — HDF5 write status checked by `BELFEM_ASSERT` (silent in
   release); should likely be `BELFEM_ERROR`.
7. `src/io/cl_Ascii.cpp:88,94` — write guard rejects only `OPEN_RDONLY`; a parallel `save()`
   truncates one shared file from every rank.
8. `meshtools.cpp:322` — `interpolation_order_numeric` has no LINE5/TRI15 arm though the factory
   and Gmsh map build both; `interpolation_order()` lacks PYRA*, LINE5, TRI15, TET35, so a linear
   PYRA5 starts flagged curved.
9. `cl_Spline.cpp:599` / `cl_Spline.hpp:606` — `update_data` builds 4 rows without resetting
   `mExtraMode`; `integrate()` then reads row 4 (no caller hits it today); comment at
   `cl_Spline.hpp:223` false. `compute_tmatrix_quad16/hex64` call `bspline_quad`, which fills 3 of
   4 entries (`cl_TensorMeshFactory.cpp:1006-1017`).
10. `cl_TensorMeshConfig.hpp:212-213` — `element_ijk` casts a possibly negative value to unsigned
    before clamping.
11. gastables POLY mode — synthesized viscosity/conductivity evaluate to 0
    (`cl_GT_TransportPolyEmpty.cpp:57-59`) while `has_viscosity()` stays true (medium: the
    binding was spot-checked, the full chain was not).
12. materials — `buffer` reservation enforced only for `usermat` (`cl_MaterialFactory.cpp:336`);
    `UserDefinedMaterial::jc_custom`/`n_custom` call null slots; `have(lambda)` true with a null
    routing pointer; `mt_thermal_phi.cpp:28` ignores `density_correction`.
13. graph — `metis_ndp(g,1)` reaches `METIS_NodeNDP` because the redirect at
    `fn_Graph_METIS.cpp:111-115` is commented out (failure unproven); `Vertex` breaks the Rule of
    Five.

## Physics and design questions — Christian, not settled by vote

- circuit: is the one-step IZ lag intended (`cl_FEM_Controller.cpp:3061-3063`, save before
  compute)? is solving the circuit once per attempt (not per iteration) intended?
- core: `rand` range — code or doc?
- iwg: `impose_alpha` T_inf; BDF ramp off-by-one; conductivity tangent; `rho` as density.
- maxwell: quasi-Newton literature claim; μ0-only conductor equation; XFEM naming.
- homology: `CutSet::test`.
- tensor: Poisson-ratio naming convention in `fn_compliance_matrix.hpp:37-42` (ν12/ν13 mean the
  textbook ν21/ν31, ν23 means ν23).
- gasmodels: should a cubic `Gas` set its own liquid flag?; two-phase wording.
- gastables: correct reaction-enthalpy recipe under the hybrid enthalpy convention.
- materials: route Hastelloy χ?; iron ceiling 860 K; Debye-curve provenance.
- sparse: keep, cite or drop the solver size bands, "50–90 %" BLR memory, "10–100×" fill numbers.
- io: should mutable `line()` mark the buffer changed; should HDF5 in a no-HDF5 build refuse like
  `XML`; `Ascii::save()` always returns false.
- quaternion: `inv()` refuses norms below ~4.7e-8 that `normalize()` accepts.
- linalg: should the Armadillo wrapper gain the expression constructors Blaze has; should
  `Matrix / scalar` exist.
- lapack: only `gesvd` floors the `lwork` query and refuses to shrink Work.

## Side findings

- `literature/books/index.md` does not exist, but `CLAUDE.md` routes to it in three places
  (papers/* and lectures indexes do exist). Fix belongs in the literature repo or `CLAUDE.md`.
- `scripts/cross_review.sh` stops cleanly on quota; the sweep driver now distinguishes a failed leg
  from an audited guide.

## Fix round (same day)

Christian approved fixing the guides; source defects went to a parallel session (belfem-46) under
plan+jury, and the CLAUDE.md drift found on the way went to a third (belfem-f1: a610b15, fde2ba5).
One fixer subagent per guide applied every CONFIRMED finding, re-checked each edit against the
source, stated Christian-routed items as current behaviour without a ruling, ran a
`clang++ -fsyntax-only` probe of the examples against the real headers where code changed (both
backends for most guides), and ran a Codex `gpt-5.6-luna`/medium language sweep over the changed
sections. Claude read every diff, spot-checked the load-bearing claims, and committed one guide per
commit: 3ce9e6c quaternion, a63bf82 lapack, eeaec33 spline, 676f129 linalg, 3b441b3 tensor,
c0056b8 database, e77300a materials, 32bc757 io, aff9a49 gasmodels, 6df029d gastables,
a7ade12 graph, a4b84fd container, d99ed12 comm, d4411ca circuit, 382506e core, b3c49cc sparse,
183492a interpolation, 126d239 + 06da4bb mesh, 4df4d3a iwg, fc236d3 maxwell, 9d00f9d homology,
31a5a64 dof_manager; plus 1b3d7de (mesh README and contracts) and 46e2e20/b433677 (homology
README), which repeated the same stale recipes.

Verdicts the fix round overturned (each re-traced in source):

- iwg W15: no BDF off-by-one under `FEM::Controller`, which calls `shift_fields()` at the start of
  a step; only the guide's own loops shifted late. Removed from Christian's list.
- linalg F10: `dot(Matrix, Vector)` compiles only under Armadillo (probe).
- core C29: the MPI RNG is seeded per rank in `gComm.init()` (`random_device + rank`), not
  identical until seeded.
- graph F17: `parmetis_nd`/`ptscotch_nd` unflag first; only the serial wrappers turn a flag into a
  self-loop.
- container F12 (refuted as "UB"): a locked-bitset write runs in release and leaves the hash stale.
- A fixer's new claim that Helmholtz `dcpdT()` returns entropy was refuted (central difference of
  `cp`, at HEAD and in the tree).

New source defects found while fixing, forwarded to belfem-46: `clean_string("")` calls
`front()`/`back()` on an empty string (`stringtools.cpp:98`, confirmed); `CutFactory::cuts()` is
always empty because `mCuts` is never filled, so the `set_domain_type(Cut)` loop at
`cl_MaxwellFactory.cpp:1195` never runs (confirmed); `SolverData::assemble_jacobian` binds null
E/I pointers to references and the `cl_FEM_DofManager.hpp` comment has the D sign wrong (both
unverified). For Christian: `belfem`'s segregated thermal loop saves to a fixed
`hphi_results.e-s` (`belfem.cpp:303`) while the normal path uses `label()`.

The parallel mesh path, established while fixing mesh and reused in maxwell, homology and the mesh
README: every rank constructs the `Mesh` (the master reads), the `Kernel` constructor partitions on
the master (`partition_mesh` → `graph::metis_partition`) and runs the `Distributor` on every rank;
nothing in `src/` calls `Mesh::partition()`. `staticheat.cpp` is not built and does not compile.

Follow-ups after belfem-46's first plan round (78a3455, f50317c): `metis_ndp( g, 1 )` works
(standalone probe against SCLS libmetis: NodeNDP with npes = 1 returns the NodeND permutation), so
the graph guide no longer only traces the path; a cubic-EoS `Gas` takes the liquid root only after
`set_liquid_flag( true )`, and the flag is per object (per-point phase selection filed as a design
change, `todo/cubic_eos_phase_selection.md`). **Christian ruled that the `cut_<n>` sidesets stay
`DomainType::Default`** so the Maxwell solver does not select them: the only reader of
`DomainType::Cut` is `cl_Maxwell_FieldList.cpp:383` (thin-cut dof table), while thick cuts carry
the jump through the abstract nodes. The fix round had added `set_domain_type( Cut )` to three
homology examples; f50317c removes it and documents the reasoning. The `MaxwellFactory:1195` loop
over the always-empty `cuts()` is dead code; deleting it is an optional cleanup for belfem-46.

Left as is: closing "Revision History" lines in linalg and container that credit reviews (history,
Christian's call); a pseudo-code block in the homology guide labelled as such.

## Guides brought in line with the defect fixes

belfem-46 landed 28 one-defect commits (e55aee3..1109752). Ten fixer agents (brief: commit diff
plus current source, Codex luna/medium sweep over changed passages) rewrote every passage those
commits made false; each diff was reviewed against source before commit:
1eebaf6 tensor (cyclic nu23/nu31/nu12, trap, example_udfiber numbers), 904879a sparse (rcond = SPD
Cholesky of A, BELFEM_ERROR cases, indexing restored, direct thresholds), dd73c0d circuit (terminal
columns same-step, userdefined units via set_user_scale), c4d137c homology (cuts() gone, find
cut_<n> by label, set_periodicity optional), f75a069 core+container (rand [0,1) with LWG 2524
caveat, clean_string no-op cases), 72c9728 io (HDF5 status checks always-on, Ascii save refuses
both read-only modes), c1819d8 spline+database (bspline_cub used, test coverage stated narrowly,
element_ijk clamps before converting, update_data resets extra mode), 82455a1 mesh (pyramid and
quartic orders, OrderConverter refusal), 71e162d materials (have(lambda), buffer name, jc_custom),
371547b maxwell (segregated loop outputs, Cut row), 5384f20 dof_manager (D stored as -K_fd, late
use_jedi_force), 2192f84 graph (Vertex non-copyable), 8199d24 gastables (POLY reads synthesized
spline). Reviewed, not verified; no probe needed (no example code changed beyond identifiers).

Residuals the fixers reported (not acted on):
- Terminology, for Christian: the homology guide calls the `cut_<n>` sidesets thick cuts, but they
  are made by `CutData::add_thin_cut_sidesets_to_mesh` and `cl_CutProcessor.cpp` reports thick cuts
  as not implemented. Medium confidence it is an inconsistency.
- clean_string: a one-character string that needs cleaning (" ", "#") comes back unchanged
  (`tLength > 1` branch, stringtools.cpp).
- Genome::inherit: mutation index can still reach N*B on an LWG-2524 library in MPI builds.
- element_ijk with a NaN coordinate: clamp passes NaN, the cast is still UB.
- OrderConverter::check_input_mesh: the all-linear check is now unreachable.
- hdf5_tools.hpp comments near check_read_datatype carry history wording.
- doc/README.md:47 still says Nickel's E carries the ΔE dip (retired in materials guide rev 2.1).
- Header Version lines lag their revision tables in the materials (2.0 vs 2.2) and dof_manager
  (1.0 vs 1.3) guides.

## Open

- The code defects (with belfem-46) and the design questions above. When a fix changes behaviour
  a guide documents (rand range, rcond, Ascii parallel save, LINE5/TRI15, circuit lag, cuts()),
  the guide is updated to match.

## Files Updated

- All 22 `src/**/doc/*usage_guide*.md`, `src/mesh/doc/README.md`,
  `src/mesh/doc/mesh_contracts_and_invariants.md`, `src/homology/doc/README.md` (commits above)
- `devlog/dl20261022_usage_guide_jury_sweep.md` (this file), `devlog/README.md` index line
