# Changelog

All notable changes to BELFEM are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **New example deck `examples/rutherford`.** REBCO tapes wound around a flat
  Rutherford-type copper core with rounded edges, the flat flexible cable of Romancew (2026):
  a Python generator builds the tapes as thin shells on a stadium-shaped mandrel at a given
  lay angle, closes a fraction of a pitch with a translated periodic condition, and writes
  the matching `topology` section. Tapes per layer, criss-crossed layers and the pitch count are
  parameters; the copper core itself is meshed as air. Twenty-one decks now ship.

- **The dependency stack is now documented.** `doc/scls_toolchain.md` names SCLS — the
  Scientific Core Library Stack BELFEM is built, tested and run against — and covers what it
  supplies, how to install and activate a flavor, and the four separate mechanisms that read
  `$SCLS` during the configure: the build defaults derived from the flavor name, the
  link/include/rpath directories, the package search order, and the compiler pin. It also states
  plainly what the CMake does not state: SCLS is not enforced, the search falls back to `/usr`,
  and no other configuration is tested. `README.md`, `doc/getting_started.md` and `CLAUDE.md`
  now name the stack where a reader is about to run `cmake`; previously SCLS appeared in the
  documentation only as a mechanism for presetting `USE_MKL`, `USE_PARDISO` and `USE_DEBUG`,
  and was never introduced. The page separates SCLS's delivery paths — the RPM and DEB
  repositories, a source build into an isolated prefix on any other Linux, and the native
  macOS build — so the binary repositories are not read as the limit of where the stack runs,
  and it names Apple Silicon as a supported platform.

### Changed

These results intentionally differ from 0.9.2. Each is a correction, so a deck that exercises
the affected path now produces a different, correct answer.

- **Surface integrals on the quadrilateral side faces of a prism are twice their 0.9.2 value.**
  `normal_penta` scaled the side-face increment by 0.25 for a quadrilateral rule whose weights
  sum to four, so each side face integrated to half its area. The three side faces now integrate
  to their face areas; the triangular caps were already correct. Any sideset term on a PENTA
  side face, including Neumann and background-field terms, changes accordingly.
- **The ghost-penalty consistency term between stacked thin-shell layers is the weighted average
  used by the SWIP method.** The flux operator used the jump in the two layers' slopes with
  twice the harmonic weight, which reproduced the interface flux only when one resistivity was
  three times the other. It is now ρ_w (∂h_m + ∂h_s), with ρ_w = ρ_m ρ_s / (ρ_m + ρ_s) (Ern,
  Stephansen & Zunino 2009, Eqs. 7, 8, 15). Results change for any deck that enables the ghost
  term; the shipped decks leave it off.
- **Thin-shell layers retain their through-thickness jump in nodal output on curved shells.**
  Nodal recovery now fits in coordinates centered and scaled at the node, and in curvilinear
  shell coordinates for linear thin-shell layers. On the bends of the `rutherford` example,
  0.9.2 placed the YBCO current sheet in the silver and buffer layers of the Exodus output and
  in the printed J/Jc_max. The solve itself never read these fields.
- **REBCO jc and n below the tables' lowest field, 10 mT, now continue smoothly** instead of
  being held at the table edge: n is a parabola with zero slope at B = 0, and jc is a cubic
  Hermite interpolation to the measured self-field level. Both are C1 and non-increasing in |B|,
  with exact derivatives. BSCCO jc retains the clamp. The HDF5 tables are unchanged.
- **The specific heat of the pure metals is splined on a 1 K grid** instead of a 4 K grid, so
  the T³ start below 4 K is no longer undershot; lead's c_p was negative there. c_p changes
  between the old knots for every pure metal, most noticeably below 20 K.
- **Magnetoresistance of silver, lead, and tin at low field.** Silver's low-field Kohler cubics
  were placed in the wrong coefficient slots. The transverse switch point for lead and tin was
  never written, so branch selection read uninitialized memory. All three now use the matched
  cubics below the switch point.
- **YBCO's thermal conductivity is finite below 1 K.** The Debye and Callaway integrands are now
  evaluated in an overflow-free form. Values that were finite in 0.9.2 are unchanged.
- **Newton tangents.** The modified Kim jc law now carries its exact field derivative; it was
  previously zero. The piecewise power law's tangent below j1 is now the derivative of the value
  it returns. Converged answers are unchanged, but iteration counts and timestep histories may
  differ.
- **Element facets.** Facet 1 of the cubic QUAD16 now lists node 7 instead of node 8, so
  boundary integrals and Dirichlet nodes on that facet change.
- **Floor retries.** When the magnetic stagnation guard cuts an attempt already at the timestep
  floor, its Picard verdict is suspended until the next accepted step. The retry is bounded by
  the raised iteration budget and the watchdog instead of replaying the attempt. A step rejected
  by 0.9.2 at the floor can now be accepted.
- **The iteration line** now prints the residual where the head is measured, so the column names
  a state and, on the success path, ends with the accepted iterate. Scripts that parse the log
  read different numbers.
- **`examples/corc_periodic_bc` ships `tolerance : 1e-10` with Anderson mixing enabled.**
- **Cached resistivity tables carry a content revision.** The `<label>_RRR<n>.hdf5` files that
  metals and alloys write into the run directory now store a `revision`. A file without it, or
  with another value, is rebuilt once with a warning. Every 0.9.2 cache is rebuilt on the first
  run, so the corrected silver, lead, and tin laws reach alloys such as PbSn solder. A 0.9.2
  binary reading a 0.9.3 cache does not check the revision.

### Refusals

Inputs that 0.9.2 accepted and handled silently or incorrectly now stop during setup with a
message naming the problem. None of the shipped decks is affected.

- A sideset with a domain type on a periodic plane; it used to lose its kernel.
- A Dirichlet condition on a non-duplicate dof condensed onto more than one source;
  reconstruction overwrote the value.
- An `rrr<N>` token in an alloy name that disagrees with the `RRR` key, or a `builtin` value
  naming a non-iron material beside a b-h `curve`.
- A thin-shell side connector selected without its layer block.
- An h-conductor beside a thin shell when the conductor is not linear, or when its constant
  permeability is not μ0 within a relative 1e-6.
- An incomplete node match on a periodic pair, or a non-unit cochain coefficient on a cut.
- A periodic node pair whose interface copies lie on tie-weight interface sets of different
  master domain type (a conductor|air set against a buffer|air, ferro|air or air-master set) or
  on sets of different treatment (a coil touching the interface on one side of the plane only),
  when no compatible set holds the partner's copy. 0.9.2 tied such copies regardless. The
  message names both nodes and both domain types or treatments; the remedy is to move the
  periodic plane into a section of one domain type.
- An HDF5 array, vector, or string dataset whose rank is not 1.

A constant `jc` without an `ec` now gets the default 1e-4 V/m, as a jc function always did,
instead of producing a NaN power law.

### Fixed

- **A 3D thin-shell sideset can no longer leave setup with a mixed master side.**
  `MaxwellFactory::fix_facet_masters` orients facets that have the same domain type on both
  sides by a walk over shared edges. That walk also passes through faces that meet a thin shell
  at a T-junction, such as a block face split at the tape seam or the coil's top face. It could
  therefore reach one shell sideset from two directions and give part of it the opposite
  winding. No deck sign can repair that, and building the shell from it corrupted the heap
  (crash in `ThinShellFactory::update_node_facet_tables`). The new pass
  `MaxwellFactory::align_thin_shell_facets` runs at the end of the 3D walk and makes every
  thin-shell sideset consistent in three steps:
  - each connected patch of the sideset takes one winding by a vote: the cross-type facets vote
    when the patch has any (they carry the deterministic domain-type orientation), otherwise all
    facets do, and a tie keeps the winding of the patch's first facet;
  - patches of one sideset that have the same two blocks on their sides are then given the same
    master block, taken from the first patch whose vote had cross-type facets, or else from the
    patch with the smallest facet id; patches whose block pair no other patch shares, or whose
    facets do not share one pair, keep their own winding, and the log counts each kind;
  - a sideset with a facet that lacks an element on one side, or with one block on both of its
    sides, is refused with a `BELFEM_ERROR` that names the shell, the sideset and the block, since
    the cut opens a shell block by block. This also refuses a uniformly wound sheet that has the
    same block on both of its sides, which earlier versions accepted. A thin shell must separate
    two different volume blocks: split the volume block so that each facet has one block on each
    side, and split the sideset as well only where it genuinely lies between two block pairs.

  A one-patch sideset that was already uniform and has different blocks on its two sides is
  unchanged, so existing decks keep their signs. A sideset that cannot be oriented consistently
  stops setup with a `BELFEM_ERROR` that names it. 2D is unchanged. A processed mesh (`.bfm`)
  written by an earlier build is reloaded without these passes and keeps the windings it was
  saved with; delete it to have them applied.
- **A 3D thin shell whose outer edge is a closed loop is accepted when that loop is isolated.** When no part of the shell reaches the domain boundary (for example a tape
  fed by conductor leads instead of tape leads), its side curve is one closed loop.
  `CurveFactory::identify_subchains` gave an edge its chain index only when it led to an
  unvisited node, so the edge that closes the loop kept no index; every edge reached from a chain
  now gets that chain's index, and `CurveFactory::thin_shell_side_curves` builds the loop as one
  closed curve instead of refusing it; a side curve that branches is refused there by name. The
  nodes of a closed rim stay single-valued, while the nodes of a terminal curve are duplicated
  and a rim that touches another shell's side curve is duplicated as an open one, so a closed
  rim that shares a node, even at a single corner, with a terminal curve of any shell or with a
  side curve of another shell is refused in `CutFactory` before any node is duplicated, with a
  message that names the shell, the rim, the shared node and the other curve. Edge coating is
  not supported on a closed rim and stops setup with a message that says so.
- **A refusal during field setup no longer frees hanging dofs twice** in debug builds and tests.
  A release build aborts before unwinding and was not affected.
- Every error path in the material, boundary-condition, periodicity, and HDF5 loaders now frees
  what it built or opened before refusing the input.

## [0.9.2] — 2026-09-18

A maintenance release: memory ownership on the setup and teardown paths, a configure-time
MPI link check that catches the rpath defect a 0.9.1 build hit on another machine, and the
temperature-dependent elastic constants of the pure metals. Input files and output formats are
unchanged.

### Added

- **New example deck `examples/tapestack3d_quench`.** The `tapestack3d` geometry — eight
  soldered REBCO tapes, periodic, coupled h-φ/T — designed to be driven into a quench by a user
  plugin instead of a prescribed source: a linear 200 A/s transport ramp and a local Ic defect on
  the topmost tape, both reading one time program (`src/ramps.hpp`) so they cannot drift
  apart. It is the worked example for combining a `userdefined` current boundary condition
  with a material `defect { }` plugin in a single shared object. Nineteen decks now ship.

- **MATLAB derivation scripts for the interpolation tables**, recovered into
  `src/fem/interpolation/doc/matlab/`: the TET10 Nédélec generator with its zero-tests and
  DefElement cross-check, the triangle derivation notes, the PENTA18 second-derivative
  zero-test and the PENTA6 facet notebooks, each with a purpose comment and license header.
  `compare_tables.py` checks the transcriptions against the live C++ tables (300 entries,
  0 differences); the two asserting drivers run under Octave 9. Excluded from Doxygen.
- **`doc/commenting_guidelines.md`**, the rule set the source comments now follow, with the
  sweep tooling `scripts/check_comment_only.sh` and `scripts/strip_comments.py` that gate a
  comment-only commit by diffing the stripped sources.
- **Regression tests** for the fixes below: `tests/fem/test_SideSetRebuild.cpp` (edge
  functions after an integration-order change, debug builds), `tests/fem/test_BlockRebuild.cpp`
  (a block enriched twice and destroyed; the leak is what Valgrind checks on that fixture) and
  an identity check on the calculator's workspace vectors across a second
  `set_integration_order()` call.

### Changed

- **CMake 4.0 is now the minimum**, raised from 3.11 and applied everywhere the project states
  a floor: the top-level `CMakeLists.txt`, both plugin templates
  (`src/physics/materials/User{Material,Library}Template.cmake`, previously 3.13) and the five
  example-deck plugin builds. The project is developed against a toolchain that ships CMake 4,
  and `scripts/scls_env.sh` already required `cmake` to resolve there; the declared floors were
  the last place still claiming otherwise. Every policy up to 4.0 now defaults to NEW —
  a full configure was re-run and produces no warning and no policy or deprecation message. A
  distribution CMake older than 4.0 can no longer configure BELFEM.

- **Third-party headers are system headers for diagnostics.** The include directories of the
  linked libraries (SCLS, MKL, STRUMPACK, PETSc, HDF5, …) are passed to the C++ compiler with
  `-isystem` instead of `-I`. A warning raised inside one of these headers no longer appears in
  a BELFEM build and cannot become an error under `-Wall -Werror`. This is preventative: the
  0.9.0 build failure on Ubuntu came from two BELFEM headers and was fixed in 0.9.1, but the
  same mechanism (Ubuntu's GCC injecting `-Wformat-security`) would turn the next vendor-header
  warning into a build failure. BELFEM's own headers keep ordinary `-I` and are diagnosed as
  before; a directory the compiler already searches by default also stays `-I`. The per-test
  `BEFORE` include of `$SCLS/include` was removed as redundant, so the tests now search
  BELFEM's headers before SCLS, as the library always did. Plugins built from the
  `User*Template.cmake` templates are not affected; they still receive ordinary include paths.
- **`MaxwellFactory` and `ThermalFactory` are neither copyable nor movable**; all four
  operations are deleted, because each factory owns what it creates until it hands it over.
- **`Metal::create_mech` (protected, for derived metals) takes `(E2, nu2, T2, deltaK, deltaG)`**
  instead of `(E0, b, T1, T2, nu2)`. The five argument types are unchanged, so a derived class
  outside the tree still compiles and must be updated by hand to the quasi-harmonic closure.
- **`spline::create_helpmatrix` takes its size as `const index_t` and its spacing as
  `const real`**, both by value (previously `const real &` for both). Every caller in the tree
  already passed an integer count. A separately built consumer must be recompiled.

- **Validity ceilings of copper, silver, chromium and white tin lowered** (`T_max` 1358 → 1000 K,
  1235 → 900 K, 2180 → 570 K, 505 → 400 K): above these temperatures the new elastic curves depart
  by more than 5 % in shape from Blanke's dynamic-modulus compilation, the only high-temperature
  reference available. The declared validity range and the generated property tables of these
  metals, and of the formula alloys that contain them, now end there; iron (860 K) and aluminum
  (933 K) are unchanged.

- **The `material` report prints the Grüneisen parameter at 298.15 K**, the temperature of the
  printed density, as a diagnostic computed from the served bulk modulus, α, ρ and c_p. It is no
  longer a model input.

- **Configure now links an MPI probe and keeps the system library directories off the rpath.**
  A 0.9.1 build on another machine failed at every executable link with 38 undefined `PMIx_*`
  references out of `libmpi.so`: BELFEM recorded the system compiler's `/usr/lib64` first on
  every executable rpath, and GNU ld resolves a shared library's own `DT_NEEDED` entries
  through `-rpath` in order, so a distribution PMIx beat the SCLS one. `detect_gcc.cmake` no
  longer appends the compiler libdir for a `/usr` prefix, and `/usr/lib64`, `/usr/lib`,
  `/lib64`, `/lib` and the multiarch directory are pruned before the build and install rpaths
  are set. `find_mpi.cmake` then links a two-call MPI program with the wrapper right after the
  flavor probe; a failed link stops the configure with the linker output, a PMIx-shaped failure
  names the package to install, and under SCLS on Linux the summary prints the PMIx the loader
  resolves. Behaviour change: a compiler wrapper that cannot link `MPI_Init` is a configure
  error instead of a failure at the first executable.
- **The test launcher `mpirun` is taken from the MPI installation the wrapper belongs to.**
  `Add_Test.cmake` searches `MPI_HOME/bin` alone, then the wrapper's own directory together with
  its symlink target's, then CMake's default search, and prints the choice once.
  Before, CMake's default search visited its own install prefix first and could name a foreign
  launcher. A launcher from `MPI_HOME` that is not beside the wrapper, or one found by the
  fallback, is reported as uncertain.
- **Homology, 2-D terminal-curve orientation (`CutFactory`):** the preconditions of the
  master-side sense test (a first segment of non-zero length, a facet of the sideset that
  carries both of its nodes, that facet's master element carrying the same nodes, a non-zero
  in-plane facet length) are now always-active `BELFEM_ERROR` checks where before a debug-only
  assertion guarded a dereference. The sense is read off the master's
  edge nodes and the contract is stated once at the declaration. The orientation decision itself
  is unchanged; `reorient_generators` is untouched. `src/fem/maxwell/doc/thin_shell_facet_orientation.md` now
  says that BFS propagation is 3-D only and that the uniform 2-D master per tape is an input
  assumption (lower element id, contiguous ids per block).
- **Source comments brought under the commenting guideline; license headers everywhere.**
  Two comment-only sweeps over `src/` (gated by `check_comment_only.sh`, `make check` on both
  backends): wrong comments corrected first, dated and credited lines and about 270 lines of
  commented-out code removed, TODOs owned with exit conditions or deleted, narration cut, and
  the ownership, MPI-collectiveness, buffer-layout and unit contracts that the code cannot state
  added on the public entry points (`commtools.hpp`, the containers, the linalg dispatch headers,
  every `Mesh` accessor and adopting setter, the dof manager's collective calls, the solver
  wrappers). 121 IDE bylines removed; the license block is on every source file. The sweeps did
  not touch the closed cohomology core (`cl_Homology`, `cl_Cohomology`, `cl_SimplicialComplex`,
  `fn_Smith`).
- **Homology: the chain coefficient handed to `addChainToChain` is an `int` literal** (`±1`)
  instead of a `real` converted at the call. This preserves behaviour exactly and clears the
  last `-Wfloat-conversion` warning in the closed core; reviewed by two independent auditors.

### Fixed

- **`pardisotools.f90` calls MKL PARDISO through MKL's own Fortran interface** (`include 'mkl_pardiso.f90'`,
  the pattern `mumpstools.f90` uses for MUMPS). This removes the two gfortran argument-mismatch warnings and
  exposed two defects the implicit interface had hidden: the solve phase passed a seventeenth argument that
  exists only in the Panua/Schenk PARDISO API, and `PARDISO::get_determinant()` read a value MKL never
  writes. The extra argument is gone and `get_determinant()` now raises an error for PARDISO (MKL computes
  no determinant; use MUMPS). At the Fortran boundary `pardisotools_solve` now takes the right-hand
  side as mutable storage and `pardisotools_get_determinant` is removed. `USE_MKL_64BIT_API` together with `USE_PARDISO` is refused at compile time,
  because MKL's Fortran interface declares default `INTEGER` and the ILP64 library would need
  `-fdefault-integer-8`; before, that combination compiled with mismatched integer widths.

- **The remaining `-Wfloat-conversion` warnings in the library and the test headers are
  cleared.** The flag was added in 0.9.1; the implicit floating-point-to-integer conversions it
  reported — `std::ceil`/`std::round` results assigned to integers in the integration-point
  tables, `Section::get_int`, the `Alloy`, `Metal` and `SplineLookupTable` spline setups,
  `Genome::encode`, and two string helpers — are now explicit casts of the same expressions,
  or integer arithmetic where the string helpers only needed integers. Results are unchanged
  for the inputs the code is given in the tree. `<cmath>` is included
  directly in the three files that previously used `ceil`/`round` unqualified.
- **Input-file integers are range-checked.** `Section::get_int` rounded the stored real and
  converted it to `int` unchecked, which is undefined behaviour for a non-finite or huge value;
  it now stops with an error naming the key and section. A unit exponent (`m^2`, `s^-1`) was
  parsed with `std::stod` and truncated to `int` unchecked, so `m^1.5` silently became `m^1` and
  `m^abc` escaped as an uncaught exception; a non-integer, non-finite or out-of-range exponent
  (magnitude above 32) now stops with an error naming the unit. Rounding in `get_int` is
  unchanged. The exponent check runs while the unit is parsed, the range check on every
  `get_int` call; both are active in every build.

- **Poisson's ratio of the pure metals fell with temperature; it now rises, as measured.** The
  elastic closure of `Metal::create_mech` held the Grüneisen parameter constant and inverted the
  Grüneisen relation for the bulk modulus, which copied the unknown temperature dependence of γ
  into K: copper's bulk modulus softened 12 % between 0 and 293 K where measurement gives 3 %, and
  ν read 0.363 → 0.344 where Ledbetter 1981 (10.1002/pssa.2210660209) gives 0.340 → 0.347. Seven of
  the nine metals were inverted. The bulk and shear moduli are now exponentials of the logarithmic
  thermal strain the expansion curve already integrates (quasi-harmonic, Garai & Laugier 2007,
  10.1063/1.2424535), fitted per metal to measured single-crystal or polycrystal constants between
  4 K and room temperature (sources with DOIs in `material_property_sources.md`), with E and ν
  derived. The served moduli are isothermal, converted from the ultrasonic (adiabatic) data with the
  material's own α, ρ and c_p. Room-temperature levels move to the literature polycrystal values,
  including Hill averages of single-crystal measurements:
  copper E 133 → 128 GPa, silver 84 → 81, aluminum 73 → 70, tin 55 → 49, nickel 184 → 223 (the old
  value was the demagnetized state; BELFEM's nickel is magnetically saturated, so the ΔE dip below
  the Curie point is deliberately not represented, see the class header), lead 16.5 → 16.25 at the
  static level (Blanke's E with the crystal's bulk modulus; the dynamic Hill average would be 24 GPa).
  Chromium serves a constant ν = 0.237 because its spin-density-wave anomalies collapse the bulk
  modulus by 20 % toward the Néel point and no smooth closure follows that. Iron's isothermal ν is
  flat within 1e-4. Construction refuses constants that let ν fall by more than 1e-3 between
  sampled temperatures.
  New regression test `tests/physics/test_MetalElastic.cpp`.

- **`integrate_scalar_over_sidesets` broadcast the master rank instead of the value**, so
  every non-master rank returned zero. It now broadcasts the value; the contract states that the
  integral is computed on the master from its mesh and then broadcast, not distributed. No
  in-tree caller.
- **`IWG::collect_node_coords` wrote `nDim + 1` columns into the caller's matrix.** It now
  writes `nDim` columns into a self-sizing, asserted matrix, the dimension being the one the
  derived IWG's constructor sets.
- **Memory ownership on the setup and teardown paths**, found by two independent deep scans and
  Valgrind gates on `dipole` and on `tapestack3d` at two ranks, fixed in six audited tiers:
  - a kernel chained onto a parent allocated a placeholder mesh on every worker and never owned
    it; the mesh `MaxwellFactory` creates is handed to the kernel through the new
    `Kernel::claim_mesh_ownership()` and deleted by the factory only if no kernel was built;
  - the thermal equation was never handed to its kernel and so never freed; `ThermalFactory`
    now registers it with `add_equation()` before `create_field()` links it, and the kernel
    deletes it. The thermal boundary-condition factory, an empty sideset's calculator, the
    element's local-dof array and the distributor's control-point staging data were never freed
    either;
  - on the master of a distributed run the dof manager moved the hanging dofs out of the list
    the disconnect walked, and the next kernel on the same mesh allocated over the live basis
    dof container; the disconnect now walks the hanging list too and the container is freed
    before it is replaced;
  - the line elements of a 3-D Gmsh mesh (`mBoundaryEdges`) were never deleted; `memory()`
    counts them now;
  - the mesh-backed `SideSet` constructor built its integration lookup tables and then let its
    calculator's `set_integration_order()` build them again over the live pointers; one
    `delete_lookup_tables()` now runs before every rebuild and from the destructor. The
    calculator also precomputed a sideset's edge functions before the tables were rebuilt at the
    new order, latent because every caller passed the order the tables already held; the rebuild
    now runs first;
  - `~Block` never freed the last generation of its enrichment tables;
  - `Calculator::allocate_memory()` cleared its non-owning vector list on every
    `set_integration_order()` after the first, orphaning a generation of workspace vectors;
    vectors are now found by label and reused, as the matrices already were;
  - `CutSet::~CutSet` deleted `mBitset` but not `mNodeBitset`.

  The tables are rebuilt with the same arguments and the reused vectors keep their contents,
  so nothing the solver reads is meant to change; the `dipole` result file is byte-identical
  to the tier-1 run after every tier. Dead declarations removed
  alongside: `IWG::N`, the never-assigned `Wrapper` coordinates and their null-dereferencing
  accessors, `SideSet::mSideSetIntegrationData`.
- **Homology, 2-D terminal-curve orientation:** five 3-vectors were constructed without
  initialization, leaving components unwritten under the Blaze backend; they are now
  initialized to zero. A Blaze build could read those components before, which is undefined
  behaviour, so a Blaze result on the 2-D thin-shell path may differ from 0.9.1; the Armadillo
  path is traced unchanged.

## [0.9.1] — 2026-09-13

A portability release. BELFEM 0.9.0 did not build on Ubuntu 24.04 (GCC 13.3), and once it
built, some floating-point convergence and threshold tests used the wrong `abs` overload.
Both problems are fixed. Input files and output formats are unchanged.

### Fixed

- **Floating-point `abs` resolved to the integer C `abs`.** Eleven unqualified `abs()` calls on
  `double` values could bind to `int abs(int)` from `<stdlib.h>`. This happens with libstdc++
  when no header in the include chain declares `::abs(double)`. The argument was then truncated
  to an integer, so any residual below 1 counted as zero. All eleven calls now use `std::abs`.
  Affected code:
  - `Metal::set_RRR` (`cl_Material_Metal.cpp`): the secant loop for the residual resistivity
    `rho_0` never ran, so metals using this path kept the initial guess. In the verified YBCO case
    the thermal conductivity was off by −0.94 % at 20 K. This is what failed
    `YBCOThermalConductivity.RepresentativeTemperaturesAreFinitePositive`.
  - `Metal` Debye-temperature inversion, the YBCO thermal-expansion plateau and c<sub>p</sub>
    switch searches, and the root searches in `Alloy` and `HastelloyC276`.
  - `Tmatrix` (`cl_FEM_Tmatrix.cpp`): matrix entries with magnitude below 1 were dropped from
    the sparsity pattern.
  - `fn_Mesh_ratio.cpp`: the secant loop for the mesh grading ratio.
  - `belfem` (`belfem.cpp`): the thermal/electromagnetic time-synchronisation loop could exit
    while the two clocks still differed by almost 1 s.

  Builds that already found the `double` overload produce the same results as before. On
  affected builds, material properties and coupled thermal/electromagnetic transients may
  change, because the truncation defect is removed.
- **Build failure with `-Werror=format-security`.** Ubuntu's GCC driver injects
  `-Wformat-security` by default, which turned two warnings into errors under BELFEM's
  `-Wall -Werror`:
  - `Logger` (`cl_Logger.hpp`) passed an already formatted message to `fprintf` as the format
    string. A `%` in a message, such as a file path or `50%`, was read as a conversion
    specifier, which is undefined behaviour on every platform. The message is now printed
    with `"%s"`.
  - `sprint()` (`fn_sprint.hpp`) forwards a caller-supplied format by design, and its GCC
    pragma block now also suppresses `-Wformat-security`.
- **`Matrix<T>` fill constructor took the fill value as `double`.** For integer matrices the
  value made a `double` round trip. That was exact for every in-tree caller on default builds,
  but with `USE_MKL_64BIT_API` (`index_t` = `uint64_t`) filling with the sentinel `gNoIndex`
  (2⁶⁴−1) went through a `double` that cannot hold it, which is undefined behaviour. The
  constructor now takes `T` and passes the value through unchanged.

### Changed

- GCC and Clang builds now warn on implicit floating-point to integer narrowing
  (`-Wfloat-conversion`). It is a warning only, not an error, so it cannot fail a build; it
  is there to catch the `abs` class of defect above at compile time.
- Homology: the Smith normal form no longer routes integer divisions through `floor()` (a
  no-op on `int` operands), and its integer `abs` calls are spelled `std::abs`. Reviewed by
  two independent auditors as behaviour-preserving; the coefficient arithmetic is unchanged.

## [0.9.0] — 2026-09-06

First public release, presented at ASC 2026.

[0.9.2]: https://github.com/cmesse/belfem/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/cmesse/belfem/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/cmesse/belfem/releases/tag/v0.9.0
