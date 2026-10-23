# Guide-Sweep Code Defects — Module Part

**Date:** 2026-10-22
**Purpose:** Fix the code defects outside the critical framework that the 2026-10-22 usage-guide
jury sweep found (`devlog/dl20261022_usage_guide_jury_sweep.md`), plus the four items Christian
ruled on the same day. One change, and one commit, per defect.
**Module:** `src/circuit`, `src/numerics`, `src/homology`, `src/sparse`, `src/io`, `src/physics`,
`src/core`, `src/fem/iwg`, `src/fem/thermal`, `src/math/tensor`
**AIs involved:** Claude (trace + plan), Codex `gpt-5.6-terra` (audit), Grok `grok-4.7` (audit)
**Status:** ✅ COMPLETE (2026-10-22), B11 deferred — B1 4a9bd57, B2 c908dc1, B3 7b0af9d, B4 dca88b4, B5 473b910, B6 d996b87, B7a 0373ed5, B7b 344fb76, B7c 3225f01, B7d c0bc009, B8 fedd37a, B9 ef6d45a, B10 dc72f59, B12 635ba53; rulings C1 14c8f60 (segregated Exodus name), C2 1109752 (cuts() removed). Reviewed, not verified: plan jury, code jury and three follow-up rounds (Codex gpt-5.6-terra/high|medium, Grok grok-4.7/high), `-fsyntax-only`. **Owed gate (Christian):** `make check` (io and physics are not in `check-fast`), with `-DUSE_GASMODELS=ON` for B8, `example_udfiber` unchanged output, `test_MaterialFactoryDeck` (BufferLabel.*).

> **Scope guards:**
> - Only the defects listed here; framework items are in `todo/closed/guide_sweep_fixes_framework.md`.
> - No change to any usage guide (`src/**/doc/*usage_guide*.md`), which another session is editing.
>   Behaviour changes are reported to that session with the commit sha.
> - `homology/cl_CutFactory` is outside the cohomology core (§7.1 lists `cl_Cohomology`,
>   `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`); B2 touches only
>   `cl_CutFactory.hpp`.
> - Evidence so far: source trace only. Nothing is built.
> - **Physics questions are not settled by vote** (B12, B11): the jury checks them, and Christian
>   decides.

---

## Steps

- [x] **B1 — circuit `userdefined` source writes into a boundary condition it does not own.**
  `cl_ElectricalCircuitFactory.cpp:301-302` was copied from the Maxwell BC factory
  (`cl_MaxwellBoundaryConditionFactory.cpp:573-574`). It calls
  `mPhysicalBoundaryConditions( tCount-i )->set_units / scale()*=`, but `read_circuit` pushes no BC,
  so `tCount - i` indexes past the end at `i == 0` and an unrelated BC otherwise. Circuit sources have
  no BC to carry a unit factor (`cl_VoltageSource.cpp:56`, `cl_CurrentSource.cpp:52` call
  `mFunction->compute( t )` directly).
  **Fix:** give `SourceFunction` a rank-local `real mUserScale = 1.0` with `set_user_scale( real )`,
  applied only in `function_userdefined` (`cl_SourceFunction.hpp:450-453`:
  `return mUserScale * mUserFunction( aTime );`). In the circuit branch, delete the two BC lines and
  call `tFunction->set_user_scale( tValue.first )` after `read_user_defined`. The Maxwell and thermal
  BC paths never call it, so their numbers are unchanged. It is rank-local because the circuit lives
  on rank 0 (`cl_FEM_Controller.cpp:3059`). Alternative rejected: reusing the amplitude slot, which
  defaults to NaN (`cl_SourceFunction.cpp:87`) and is broadcast with `mValues`.
  **Gate:** compile; a circuit deck with a `userdefined` source in mA reproduces the same source in A.
- [x] **B2 — `CutFactory::mPeriodicFactory` is uninitialized** (`cl_CutFactory.hpp:45`; the
  constructor at `cl_CutFactory.cpp:42-59` never sets it; it is read at `:127`). Fix:
  `PeriodicityFactory * mPeriodicFactory = nullptr ;`. **Gate:** compile.
- [x] **B3 — `rcond` factors A·Aᵀ, not A.** `fn_rcond.cpp:36` sets `stype = 0`; CHOLMOD then works on
  A·Aᵀ (`/opt/scls/include/suitesparse/cholmod.h:845`, `:2664-2665`), so the result is about
  rcond(A)². The "LU-factorization" comment is also wrong (CHOLMOD is Cholesky). `rcond` has no
  caller in `src/`, `tests/` or `nonfree/`.
  **Fix:** make the contract SPD-only and say so in a header Doxygen comment. Set `stype = 1` (use the
  upper triangle of a symmetric matrix stored in full; CSR read as CSC is the transpose, which is the
  same matrix when A is symmetric). After `cholmod_factorize`, `BELFEM_ERROR( tCommon.status ==
  CHOLMOD_OK, … )`, so an indefinite matrix (`CHOLMOD_NOT_POSDEF`) is refused rather than returning a
  number for a partial factor. Restore the caller's indexing base before returning (sibling finding
  F16, same function). The comment says Cholesky, and that `cholmod_rcond` is the cheap
  (min/max diag L)² estimate, not a 1-norm estimate.
  **Behaviour change:** the meaning of the returned value; report to belfem-85.
  **Gate:** compile under `USE_SUITESPARSE`; a 1-D Laplacian whose rcond is known in closed form.
- [x] **B4 — HDF5 status checks vanish in release builds.** Every
  `BELFEM_ASSERT( aStatus == 0, … )` after `H5Dwrite`/`H5Dread` in `hdf5_tools.hpp` (`:321`, `:368`,
  `:466`, `:543`, `:615`, `:693`, `:797`, `:902`) becomes
  `BELFEM_ERROR( aStatus >= 0, … )`, matching `:197` and `:951`. This is I/O, a setup/output tier
  call rate. **Trap found in planning:** for an empty vector or matrix the write is skipped
  (`if( tLength > 0 )`, `:590`; the matrix path likewise), so `aStatus` keeps whatever the caller
  passed. With an always-on check, a garbage `aStatus` would abort an empty write in release. Set
  `aStatus = 0` before each skippable branch. **Gate:** compile + `make check-fast` (the HDF5
  round-trip tests).
- [x] **B5 — `Ascii::save` in `OPEN_RDONLY_PARALLEL` mode truncates the shared file from every rank.**
  The guard rejects only `OPEN_RDONLY` (`cl_Ascii.cpp:88`). Fix: reject both read-only modes.
  **Behaviour change:** a save in parallel read-only mode is now refused; report to belfem-85.
  **Gate:** compile.
- [x] **B6 — `Spline` keeps a stale extra mode.** Both `create_coeffs` overloads
  (`cl_Spline.cpp:≈560`, `:594-626`) shrink `mData` to 4 rows but leave `mExtraMode`, so after an
  `update_data` on a spline that had an `Integral` table, `integrate()` reads row 4 of a 4-row matrix
  (`cl_Spline.hpp:597-606`). No current caller does this (`cl_Gas.cpp:877-878`,
  `cl_GM_EoS_Cubic.cpp:631`). Fix: `mExtraMode = spline::ExtraMode::None;` in both `create_coeffs`.
  `create_entropy` and `create_integral` set it again afterwards. This also makes the header comment
  at `cl_Spline.hpp:223` true. **Gate:** compile + `make check-fast`.
- [x] **B7 — materials, four defects (one commit each).**
  - **B7a** `buffer` is reserved only for `usermat` (`cl_MaterialFactory.cpp:336`); the `curve`
    branch (`:159-205`) and an explicit `builtin: copper` (`:212-215` → `:260`) register a
    non-Magnesia material under `buffer`. Fix: lift the refusal so that a `buffer` label must
    resolve to `BuiltinName::Magnesia` (through the same name map as `:96-98`): refuse it in the curve
    branch, and in the builtin branch unless the resolved name is Magnesia. **Behaviour change:**
    `buffer { builtin: copper ; }` is refused. A deck that relied on that has to rename the section.
  - **B7b** `UserDefinedMaterial::jc_custom` / `n_custom` (`cl_Material_UserDefined.hpp:378-389`) call
    slots that are always null (`cl_Material_UserDefined.cpp:30`), because `set_user_defined_function`
    now refuses jc/n (`:65-69`). Fix: delete both overrides, so the base `BELFEM_ERROR`
    (`cl_Material.cpp:581-591`) fires instead of a segfault.
  - **B7c** `have( lambda )` is true with a null routing pointer for a pure metal without RRR:
    `set_lambda_coefficients` raises the flag (`cl_Material_Metal.cpp:643`), and the spline that
    routes λ is built only in `set_RRR` (`:275-280`), which Copper skips for a NaN RRR. Fix: the flag
    stops serving as the "has coefficients" marker. `set_lambda_coefficients` stores the
    coefficients only; `set_RRR` tests `mLambdaCoefficients.length() > 0`, raises the flag, and then
    builds the spline (`create_spline` asserts `have()` first, `cl_Material.cpp:1047`).
    **What could break:** any code between `set_lambda_coefficients` and `set_RRR` that reads
    `have( lambda )`; to be grepped at the code jury. A no-RRR metal now reports `have( lambda ) ==
    false`, so a thermal deck using one fails loudly instead of calling a null pointer.
  - **B7d** `T_phi` (`mt_thermal_phi.cpp:28`) uses `density( gTroom )`, which skips both the
    `density_correction` factor and the `ref_density` fallback that `MaxwellData` applies
    (`cl_FEM_Calculator.cpp:142-157`). Fix: the same three-step rule (density, else ref_density, then
    × density_correction), computed once before the integration loop. B11 needs the same rule; if
    the jury prefers it, factor it into one helper beside the calculator's.
  **Gate:** compile + `make check-fast` (`test_MaterialFactoryDeck`, `test_MaterialDerivatives`).
- [x] **B8 — gastables POLY mode returns zero for synthesized transport.** Chain (source trace,
  high): with no transport data but critical data, `create_tables` synthesizes μ and λ into the
  splines and raises `mHaveViscosity`/`mHaveConductivity` (`cl_GT_RefGas.cpp:1242-1251`,
  `:1261-1276`). The polynomial lists hold only the `TransportPolyEmpty` placeholder (`:297`, `:313`),
  whose `eval` returns 0 (`cl_GT_TransportPolyEmpty.cpp:57-59`). `set_mode( POLY )` binds
  `poly_Mu`/`poly_Lambda` unconditionally (`:841-846`), so a user who switches to POLY gets μ = λ = 0
  while `has_viscosity()` is true. Fix: two members, `mViscositySynthesized` and
  `mConductivitySynthesized`, set in the synthesis branches. `set_mode( POLY )` binds the spline
  functions for a synthesized property (there is no polynomial to be exact about).
  **Gate:** requires `-DUSE_GASMODELS=ON`; a species without transport data reports the same μ(T)
  in POLY and SPLINE modes.
- [x] **B9 — `clean_string("")` calls `front()`/`back()` on an empty string** (`stringtools.cpp:98`).
  The `front()/back() != ' '` test never decides anything: a space is ≤ 32, so the scan loop has
  already set `needs_cleaning`. Fix: `if ( !needs_cleaning ) return aString;`. The empty string
  returns early, and non-empty strings behave as before. **Gate:** compile + `make check-fast`.
- [x] **B10 — serial `rand()` can return 1.0** (`random.hpp:92`: `std::rand() / RAND_MAX`). The MPI
  branch already uses `std::uniform_real_distribution( 0, 1 )`, which ISO C++ defines as [0, 1).
  `cl_Genome.hpp:187,210` then flips bit `tNumBits`, one past the end, on a 1-in-2³¹ draw. Fix:
  `( real ) std::rand() / ( ( real ) RAND_MAX + 1.0 )` and document [0, 1) in the Doxygen comment.
  **Behaviour change:** report to belfem-85. **Gate:** compile.
- [ ] ~~**B11 — `IWG_TransientHeatConduction` placeholders** (`cl_IWG_TransientHeatConduction.cpp:101`
  mass coefficient `1.0`; `:105` conductivity `111.0`). Christian ruled 2026-10-22: use the material
  data, and assert that the material is not null as the Maxwell code does. Fix: take
  `const Material * tMaterial = aCalc->material()` (which carries the
  `BELFEM_ASSERT( mGroup->material() != nullptr )`, `cl_FEM_Calculator.cpp:2990`). Interpolate T at
  each integration point from `node_data( "T" )`, as `IWG_StaticHeatConduction` does (`:96-106`).
  Then M += N'·(ρ·c_p(T))·N and K += B'·λ(T)·B, with ρ from the B7d rule.
  **Out of scope, stated:** this IWG never calls `mTimeStepMatrices->initialize` and adds no load to
  `f` (review_iwg W6), so it is still not usable end to end after this change. That is a separate
  item. **Gate:** compile.~~ Deferred by Christian 2026-10-22 ("leave it for now"); see `todo/transient_heat_conduction_iwg.md`.
- [x] **B12 — compliance matrix Poisson convention (physics check by the jury, decision Christian's).**
  Today `compliance_matrix` (`fn_compliance_matrix.hpp:37-42`) writes `S(1,0) = -aNu12/aE2`, which is
  the textbook −ν21/E2, and `S(2,0) = -aNu13/aE3` (= −ν31/E3), but `S(2,1) = -aNu23/aE2` (= −ν23/E2).
  Convention here: ν_ij = −ε_j/ε_i under uniaxial stress along i, with ν_ij/E_i = ν_ji/E_j.
  Christian ruled 2026-10-22: use the **cyclic** set {ν23, ν31, ν12}, which matches the G23/G31/G12
  order and his Hesse-matrix derivation, and document the trap.
  **Fix:** parameters `aNu23, aNu31, aNu12` (and `aG23, aG31, aG12`, unchanged).
  S(1,0) = −ν12/E1, S(2,0) = −ν31/E3, S(2,1) = −ν23/E2, then symmetric. Only S(1,0) changes value
  for the same input. `Tensor::fill_orthotropic_elasticity` renames `aPoisson13 → aPoisson31` and
  `aShear13 → aShear31` to match (`cl_Tensor.hpp:275-299`); it is the only caller. The Doxygen block
  states the convention and the trap: for a fibre along axis 1, ν12 is the major ratio and ν31 is the
  minor one (= ν21 under transverse isotropy).
  **Same numbers:** the only user, `example_udfiber.cpp:46-71`, passes `nu12 = 0.013` into both the
  12 and 13 slots. That value acts as ν21 today (it matches T300 carbon fibre's ν21; medium
  confidence on the provenance). To keep every printed number identical, the example names
  `nu21 = 0.013`, derives `nu12 = nu21 * E1 / E2`, and passes `nu31 = nu21` into the 31 slot. Its
  output line `nu12 = -S(1,0)*E1` (`:115`) already computes the major ratio, so input and output then
  agree. **Gate:** `example_udfiber` prints identical numbers before and after (Christian runs it;
  it is built only under the coding-examples option).

## Comment-only items (no jury, landed on approval)

- [x] `cl_FEM_DofManager.hpp:388-392`: the block comment said `J x = f - D y` with `[ J D ]`; D is
  assembled with the sign flipped (`D -= aJacobian`, `cl_FEM_DofMgr_SolverData.cpp:1517`) and added
  to the RHS (`:1979`), so it now reads `[ J -D ]` and `J x = f + D y`.

## Amendment after plan jury round 1 (2026-10-22)

Record: `tmp/ai_exchange/review_guide_sweep_fixes_modules.md`. Codex terra/high; Grok rerun alone
after a max-turns failure (it saw the Codex entry, so that seat is less independent).

- **B1:** the factory runs on every rank (it is not rank-guarded), so every rank sets its own
  `mUserScale`; `synch()` does not carry it. Setter is `set_user_scale( const real aScale )`.
- **B3 hardened:** refuse non-square input before touching the indexing base; refuse a build with
  64-bit `int_t` (`CHOLMOD_INT` is int32). Check `cholmod_analyze`'s factor for null and the status
  after analyze and after factorize. Refuse `CHOLMOD_NOT_POSDEF`; accept `CHOLMOD_DSMALL`, which is
  a warning. Free the factor, finish CHOLMOD and restore the caller's indexing base before any
  error. The contract: SPD, both triangles stored.
- **B4 extended** to `hdf5_tools.cpp:65, 111, 169, 214`.
- **B7c reworked:** `set_lambda_coefficients` stops raising the flag. In `set_RRR`, the spline
  branch (`:275`) tests the coefficient count, raises the flag there, and the unconditional
  `set_have( lambda )` at `:293` goes. For release builds, `mFunctionLambda` defaults to a
  `Material::lambda_undefined` that calls `BELFEM_ERROR` naming the material, so an unrouted λ
  fails loudly in release too, instead of calling null. Gate: `Copper( NaN )` reports
  `have( lambda ) == false`, and calling `lambda()` on it throws under `throw_on_error`.
- **B7d:** if a helper is factored it keeps the NaN guard; B11 is deferred, so `T_phi` uses the
  rule inline.
- **B8:** also rebind the derivative slots (`dMudT`, `d2MudT2`, `dLambdadT`, `d2LambdadT2`). The
  internal `set_mode( POLY )` at `:1217` runs before the synthesis flags are set, so it is unaffected.
- **B11 DEFERRED — needs Christian.** Both auditors: the IWG sets no matrix flags, so `M` and `K` are
  never sized in production. `tests/fem/test_Pipette.cpp:408-439` calls `compute_mkf` with no
  material and asserts `Σ M = area`, which holds only for the unit coefficient. Doing it properly
  means: flags + `initialize`, a material and T field in the pipette fixture, and a new mass oracle.
  Expand B11 to that, or leave it?
- **B12:** also update `tests/math/test_TensorKernels.cpp:568` to the cyclic argument names (values
  unchanged) and add one assertion with E1 ≠ E2: `S(1,0) = −ν12/E1`, `S(2,0) = −ν31/E3`,
  `S(2,1) = −ν23/E2`. Guide lines to hand belfem-85: `tensor_usage_guide.md:157-175`.
