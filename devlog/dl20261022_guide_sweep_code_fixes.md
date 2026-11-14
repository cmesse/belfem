# Guide-sweep code defects: 26 fixes, two rulings, one new tensor-mesh defect

**Date:** 2026-10-22
**Purpose:** Fix the code defects that the usage-guide jury sweep (`dl20261022_usage_guide_jury_sweep.md`,
session belfem-85) found, one change per defect through plan jury and code jury, plus the items
Christian ruled on during the session.
**Module:** `src/mesh`, `src/fem/kernel`, `src/math/graph`, `src/circuit`, `src/numerics`, `src/homology`
(outside the cohomology core), `src/sparse`, `src/io`, `src/physics`, `src/core`, `src/fem/thermal`,
`src/math/tensor`, `src/executables`
**AIs involved:** Claude Opus 5.5 (primary). Plan juries: Codex `gpt-6-astra`/high + Grok `grok-4.7`/high
(framework), Codex `gpt-5.6-terra`/high + Grok `grok-4.7`/high (modules; Grok rerun alone at 200 turns
after a max-turns failure). Code juries: Grok `grok-4.7`/high on both diffs, Codex seat rerun with
`gpt-5.6-terra`/high on both after a usage-limit failure (astra waived by Christian). Follow-up rounds 3–5:
Codex `gpt-5.6-terra`/high, /medium, /medium + Grok `grok-4.7`/high. Codex `gpt-5.6-luna`/medium for the
language sweeps.
**Status:** COMPLETE — 28 commits e55aee3..1109752; B11 deferred; owed gates in §5
**Verification:** static only. `-fsyntax-only` with the debug tree's flags on every touched TU and on all
392 TUs of `cmake-build-debug/compile_commands.json` (Vertex copy ban); `fn_rcond.cpp` also with
`-DBELFEM_SUITESPARSE`; `cl_GT_RefGas.cpp` with the full include set. A standalone METIS probe
(scratchpad, against `/opt/scls/lib/libmetis`). A 1-D Python model of the cubic reproduction
oracle (max error 1.8e-15). `check_citations.py`: 526 resolved, 0 unresolved. Nothing built, linked or
run in the tree: **reviewed, not verified.**
**Literature References:** Zienkiewicz et al. 2013, §17.2.1, Eq. 17.4 (Cox–de Boor recursion)
**Exchange:** `tmp/ai_exchange/review_guide_sweep_fixes_framework.md`, `review_guide_sweep_fixes_modules.md`,
`review_guide_sweep_followups.md`, `review_guide_sweep_round3.md`, `review_guide_sweep_round4.md`,
`review_guide_sweep_round5.md`, `prose_guide_sweep_fixes.md`

## Summary

belfem-85 relayed the defect list; Christian approved source edits "as soon as the jury agrees" and ruled
on the five questions the list carried. Two plans (`todo/closed/guide_sweep_fixes_framework.md`,
`todo/closed/guide_sweep_fixes_modules.md`) went through plan juries; the code went through code juries
and three narrower follow-up rounds, each verified against source before anything was applied. Every
defect is its own commit. The juries changed the work materially: a second tensor-mesh defect (QUAD16
links), a third one in every order (`tPairs`), the empty-path HDF5 trap, the release-build λ null call,
the `rcond` hardening, and a wrong premise behind the `buffer` refusal.

## What landed (commit, defect, confidence of the trace)

| Commit | Change | Conf. |
|---|---|---|
| e55aee3 | `interpolation_order` gets PYRA5/13/14, LINE5, TRI15, TET35; `interpolation_order_numeric` LINE5/TRI15; `geometry_type` LINE5. `auto_integration_order` aborted on every pyramid and TET35 (`fn_intpoints_auto_integration_order.cpp:49-53`); the disabled IntegrationAuto regression is enabled and extended | high |
| d0af680 | `OrderConverter::check_input_mesh` refuses non-upgradable types by name before facet work | high |
| 19f4411 | `bspline_cub` was right-to-left (entry 0 = (1+ξ)³/48) and uncalled; reordered, cited, used by quad16/hex64 (the quadratic call left entry 3 uninitialized: `Vector(n)` does not fill, `cl_AR_Vector.hpp:55-56`) | high |
| 1e23a2b | `link_bsplines_quad16` slots 4/5 held (p,q+1),(p,q+2) instead of (p+1,q),(p+2,q): bottom-edge control points never linked. Found by me while checking the ordering, independently by both auditors; script check: quad9, hex27, hex64 match their idx tables | high |
| 5053447 | **D10**: `tPairs` declared outside the node loop in `connect_nodes_to_bsplines` and never cleared, so every node's sources included all earlier nodes' pairs, every order (Codex, follow-up round). Projector unaffected (element control points + `mT`) | high |
| 9abf660 | `tests/mesh/test_TensorMeshBSpline.cpp`: QUAD16 (isotropic, translated anisotropic), HEX64 (translated anisotropic): control-point count, partition of unity, coverage, exact x, x², x³ with polar-form control values at the Greville abscissae; stored coordinates deliberately not read | high (algebra), test not run |
| 89fe9e9 | `element_ijk` floor+clamp in real before the unsigned cast ([conv.fpint]) | high |
| 4287445 | `SolverData::assemble_jacobian` binds E/I only inside the JEDI branch | high |
| 17eab8f | Controller: `compute_circuit_current` before `save_timestep`; only terminal-current columns change (Christian: reorder; both auditors traced that node voltages, lumped columns and the memdump are unaffected) | high |
| abe679b | `graph::Vertex` deletes copy/move (closes `Dof`, `Component`, `ElectricNode`; mesh entities already inherit `Basis`'s deletes) | high |
| 975f94f | dead NodeNDP redirect removed: probe shows `METIS_NodeNDP(npes=1)` returns `METIS_OK` with NodeND's permutation, no write past `sizes[]` | high |
| a90d2cf | DofManager block comment: D stored as −K_fd, `J x = f + D y` (`cl_FEM_DofMgr_SolverData.cpp:1517`, `:1979`) | high |
| 4a9bd57 | circuit `userdefined` wrote into `mPhysicalBoundaryConditions( tCount − i )` (one past the end at i = 0); `SourceFunction::set_user_scale`, rank-local, every rank builds the factory | high |
| c908dc1 | `CutFactory::mPeriodicFactory = nullptr` | high |
| 7b0af9d | `rcond`: `stype = 0` factored A·Aᵀ; now SPD contract (CSC upper / CSR lower), non-square / non-CSC-CSR / 64-bit `int_t` refused, status after analyze and factorize (NOT_POSDEF refused, DSMALL accepted), base restored and CHOLMOD released before the error | high (header semantics medium) |
| dca88b4 | HDF5 status checks `BELFEM_ERROR( aStatus >= 0 )` (12 sites); empty save paths set `aStatus = 0` and check `H5Dcreate` after closing handles; load paths keep the extent-query status (a reset there would hide a failed `H5Sget_simple_extent_dims` — caught by me before the jury) | high |
| 473b910 | `Ascii::save` refuses `OPEN_RDONLY_PARALLEL` | high |
| d996b87 | `create_coeffs` clears `mExtraMode` | high |
| 0373ed5 | material named `buffer` must build Magnesia; contract docs + `BufferLabel` tests | high |
| 344fb76 | `UserDefinedMaterial::jc_custom/n_custom` removed (null slots) | high |
| 3225f01 | `have( lambda )` raised where the spline is built (`cl_Material_Metal.cpp:275`, the unconditional `:293` gone); `mFunctionLambda`/`mFunctiondLambdadT` default to `lambda_undefined` (BELFEM_ERROR) | high |
| c0bc009 | `T_phi` uses MaxwellData's reference density rule | high |
| fedd37a | RefGas POLY mode binds spline value and derivative functions for synthesized transport | high (static; needs USE_GASMODELS) |
| ef6d45a | `clean_string("")` no longer calls `front()`/`back()` | high |
| dc72f59 | serial `rand()` returns [0, 1) | high |
| 635ba53 | compliance matrix cyclic {ν23, ν31, ν12}; only S(1,0) changes; `example_udfiber` remapped (ν12 = ν21·E1/E2, ν31 = ν21) so its numbers are unchanged; new unequal-modulus test. Algebra checked by both auditors | high |
| 14c8f60 | segregated loop writes `<label>.e-s` (Christian: bug) | high |
| 1109752 | `CutFactory::cuts()`/`mCuts` and MaxwellFactory's empty retyping loop removed (Christian: cut_<n> stay Default, remove dead code) | high |

## Rulings by Christian (2026-10-22)

- `rand()`: ISO `uniform_real_distribution` is [a, b); the serial branch could return 1.0 and overrun
  `Genome::inherit`; fixed rather than documented.
- IWG_TransientHeatConduction: use material data with a null assert — then, after both auditors showed
  the IWG never sizes M/K and `test_Pipette.cpp:439` relies on the unit mass: **leave it for now**
  (`todo/transient_heat_conduction_iwg.md`).
- Circuit lag: reorder. Compliance: cyclic set, "a well known trap", document it.
- Cubic EoS liquid root: phase depends on the input pair and on history; a design change, not a fix
  (`todo/cubic_eos_phase_selection.md`).
- Follow-ups F1–F5 from the code juries: apply. `buffer` refusal: keep — then, after Grok showed the
  tagging is by missing `rho` (`cl_ThinShellFactory.cpp:1966`) and not by name: keep, with the real
  reason. Spline B6 guard: keep `BELFEM_ASSERT` (Codex had called it P1). Cuts: remove dead code.
  hphi name: bug.
- New standing rule: decisions go through the Claude quiz (memory updated).
- Tensor control-point coordinates: Christian challenged Codex's claim and gave his convention
  (ξ ∈ [-1, 1]; control points at -2, 0, 2 and -3, -1, 1, 3). That convention is the Greville placement
  and agrees with the basis; the factories space by `step()` (node step), which matches only for order 1.
  Cosmetic for the projector. Filed as `todo/tensor_mesh_control_point_coordinates.md`; Christian
  expects nothing to fix and deferred the discussion.

## What the juries refuted or redirected

- Codex (framework plan): the geometric linear-reproduction gate was wrong (control coordinates never
  enter the projection); replaced by the polar-form weight test.
- My pre-registration's CHOLMOD header citations (`:845`, `:2664`) were weaker than Grok's
  (`:1191-1196`, `:2762`).
- Grok on the code: the Eq. 17.4 citation over-attributed the closed forms (reworded twice); the
  `element_ijk` comment overstated the UB window.
- Refuted by me: `aCond` naming (house rule: returned results take the `a` prefix); `tTol` at file scope
  (matches `test_Mesh.cpp`'s `tEps`); "46 hex64 slots unchecked" (my script compared all 64).

## Process notes

- `scripts/cross_review.sh` uses `mktemp …XXXXXX.md`; BSD mktemp does not randomize X's before a suffix,
  so two juries started together collide. Worked around with a per-run `TMPDIR`; not patched (tooling
  needs its own round).
- The Codex usage limit hit the astra seat; Christian waived astra for the rerun. The rerun Codex could
  read the Grok entries, so those code rounds were less independent than blind ones (recorded in the
  exchange files).

## Open / owed

- **Gates (Christian runs builds):** full build incl. tests and nonfree (Vertex copy ban); `make check`
  with `-DUSE_GASMODELS=ON` (io and physics are not in `check-fast`); `test_mesh` TensorMeshBSpline.*,
  `test_fem` IntegrationAuto.*, `test_physics` BufferLabel.*, `test_math` ComplianceMatrixCyclicPoissonConvention;
  `examples/rlc_circuit` (`L1` column vs FEM current); `example_udfiber` (unchanged output).
- Guide updates handed to belfem-85 with shas (rand, rcond, Ascii, element orders, circuit columns,
  userdefined units, buffer, compliance, spline T-matrix, cuts(), hphi name in the Maxwell guide).
- Residuals recorded, not fixed: `element_ijk` with NaN or a zero count; `load_matrix_from_file` rank
  check; `save_scalar_to_file` create check; `create_integral` early return before the mode assignment;
  JEDI late-enable is debug-checked only; acceptance `BufferLabel` test does not check the built type.

## Files Updated

Source and tests: see the 28 commits. Docs: `doc/getting_started.md`, `doc/input_file_reference.md`,
`doc/input_schema.yaml`, `src/homology/doc/README.md`, `src/physics/materials/doc/materials_contracts_and_invariants.md`
(language-swept). Todo: `todo/closed/guide_sweep_fixes_{framework,modules}.md`,
`todo/cubic_eos_phase_selection.md`, `todo/tensor_mesh_control_point_coordinates.md`,
`todo/transient_heat_conduction_iwg.md`, `todo/README.md`.

## Addendum 2026-10-23: residuals from the guide fixers

belfem-85 reported five residual items; Christian chose all five by quiz. One jury round (Codex
`gpt-5.6-terra`/medium, Grok `grok-4.7`/high; exchange `review_guide_sweep_residuals.md`).

- 96e5d51 `clean_string`: the `tLength > 1` guard returned a one-character string unchanged even when
  it needed cleaning; " ", "\t", "#" now give "" (test rewritten). Behaviour change.
- 8b60e1f `Genome::inherit`: `tMutate` clamped below N·B. The jury refuted my LWG 2524 comment (its
  resolution keeps [0, 1)); the reachable case is `rand() * tNumBits` rounding up for a draw just below 1.
- 5040629 `element_ijk`: both auditors showed the NaN assert vanishes in release with the cast still
  undefined; `! ( tS > 0.0 )` now sends NaN and negatives to element 0 at no release cost.
- c872dbc `OrderConverter`: the all-linear check was unreachable after the allow-list; removed.
- 4be9084 `hdf5_tools.hpp`: history wording removed from the `check_read_datatype` comment (comment-only).
- Cut naming: Christian first answered "thick is right"; shown `thick_thin_cuts_and_conjugate_edges.md`
  and the CutProcessor error texts, he ruled the `cut_<n>` sidesets are thin cuts derived from the
  thick cohomology cut. No code change; belfem-85 fixed the homology guide (721e9b7).

Reviewed, not verified (`-fsyntax-only` only). Owed with the rest: `test_core` (StringTools.*).
