# Task Planning

Active, deferred, and archived implementation plans for BELFEM development.

**Layout:**
- **Active** (this directory) — live plans and open bugs, forward work.
- **[`deferred/`](deferred/)** — sound ideas parked behind the current production path.
- **[`closed/`](closed/)** — completed / resolved / done plans, kept for the historical record.
- **[`closed/tests/`](closed/tests/)** — the test-suite planning campaign (modules 00–13).
- **[`cancelled/`](cancelled/)** — plans whose design was rejected, replaced by a different
  mechanism, or abandoned by ruling. Distinct from `closed/`: the work described here was
  **not** done and is not going to be.

**A plan belongs in `todo/` only while code or tests are genuinely missing from the tree.**
Once the substance has landed and the sole remnant is a run gate, a build, or a sign-off, the
plan moves to `closed/` with that gate named in its Status line — an owed gate is tracked in the
debt register, not by keeping a finished plan on the active list. Unchecked boxes are **not** the
test: they systematically overstate remaining work, because a plan is written when work is
planned and rarely updated when it lands.

**Retirement at session close (standing rule, 2026-10-05):** the session close (devlog + commit)
that finishes a plan's substance moves it to `closed/` (or `cancelled/`) in the same commit and
moves its line below; a plan that cannot close gets its finished steps ticked, its obsolete steps
struck, and its Status line rewritten (`doc/ai_collaboration_protocol.md` §6).

**Convention:** new substantial plans follow **[plan_template.md](plan_template.md)** (structure
extracted from `closed/meshfile_refactor_plan.md`): header block with living Status line, scope
guards, gap table with (a)/(b)/(c) classes, ordered `Rn` steps, `Dn` defect tracker with live
checkboxes, `On` open questions, definition-of-done, audit trail.

**Registers:** [debt_register.md](debt_register.md) (live rows) and
[debt_register_closed.md](debt_register_closed.md) (struck rows).

This index was rebuilt on 2026-10-05 as one line per plan. The earlier chronological index — dated
"Added" entries, the currentness sweeps of 2026-08-09, 08-11, 08-30 and 09-03 with their per-file
rationale — is in git: `git show 18e57fe:todo/README.md`.

---

## Active (6)

- [mkl_num_threads_allrun](mkl_num_threads_allrun.md) — Allrun pins `MKL_NUM_THREADS` to `NTHREADS` beside `OMP_NUM_THREADS`, and `MKL_DYNAMIC` gets decided
- [pancake_inner_leadout](pancake_inner_leadout.md) — Replace the inner lead of the single-tape pancake generator (`examples/pancake_tape/python/pancake/`), currently an L-shaped hard-way corner in the bore wall, by a lead that first separates from the inner turn on a smaller …
- [periodic_postprocessor_fields](periodic_postprocessor_fields.md) — make the recovered nodal fields periodic
- [test_gap_backlog](test_gap_backlog.md) — the ranked list of tests the suite is missing, distilled from the incident record
- [test_gap_followups](test_gap_followups.md) — the rulings of the Tier 1 test-gap run turned into twelve small rows
- [test_gap_tier2_execution](test_gap_tier2_execution.md) — Run Tier 2 of `tmp/devlog_sweep/triage.md` (§2, 48 rows) plus N-49 … N-53 reopened by the exclusion re-read (`todo/test_gap_backlog.md` §3b; N-54 stays Tier 3) through the standing loop, batched by suite: one plan per …

## Deferred (20) — [`deferred/`](deferred/)

- [2d_thinshell_todo](deferred/2d_thinshell_todo.md) — Dependency-ordered work items to bring the 2D thin-shell path up to the 3D reference. Derived from `todo/2d_thinshell_gap_analysis.md` (file:line evidence and confidence labels live there; v2 line references …
- [buffer_cut_implementation_plan_v2](deferred/buffer_cut_implementation_plan_v2.md) — Replace `todo/buffer_cut_implementation_plan.md` with a version that matches the current BELFEM thin-shell, cut, kernel, and thermal pipelines.
- [coreduce_complexPellikkaGeneralized_performance_findings](deferred/coreduce_complexPellikkaGeneralized_performance_findings.md) — Merged performance analysis, evaluated against actual code
- [cut_representative_options](deferred/cut_representative_options.md) — non-tight cohomology generators: tighten, go cut-free, or multi-sheet
- [edge_function_quadratic_shells](deferred/edge_function_quadratic_shells.md) — Quadratic Thin-Shell Nédélec Edge Functions (EF_QUAD9TS, EF_PENTA18TS)
- [gauging_context_from_the_quench_run](deferred/gauging_context_from_the_quench_run.md) — measured κ, the tolerance result, and the one test that decides whether a gauge term helps here.
- [handoff_scls_strumpack_openmp](deferred/handoff_scls_strumpack_openmp.md) — Handoff: SCLS ships STRUMPACK with OpenMP tasking disabled
- [input_conf_configurator_plan](deferred/input_conf_configurator_plan.md) — Author and validate `input.conf` decks from a GUI instead of by hand. The mechanism: a machine-readable schema of the input contract drives a headless validator, a mesh-introspection layer that *derives* the `topology` …
- [jacobian_init_remaining_bottlenecks](deferred/jacobian_init_remaining_bottlenecks.md) — Track the findings from the 2026-08-18 Jacobian-init investigation that are still open after the `DynamicBitset` and `SpMatrix` work landed
- [maxwell_postprocessor_gpu_acceleration](deferred/maxwell_postprocessor_gpu_acceleration.md) — Discussion of potential GPU acceleration for extreme-scale problems
- [mumps_symmetric_triangle_extraction](deferred/mumps_symmetric_triangle_extraction.md) — Design — ALGORITHM ONLY — for supplying MUMPS one triangle instead of the full matrix when `SYM != 0`, so that `SymmetryMode::GeneralSymmetric` / `PositiveDefiniteSymmetric` become supported rather than refused. Halves …
- [powerlaw_jc_n_field_derivatives](deferred/powerlaw_jc_n_field_derivatives.md) — Extend the Newton tangent of the collapsed h-kernels when critical current density `jc` and power-law exponent `n` come from lookup tables or other field-dependent `JcFunction`s: `jc = jc(|B|, β[, T])` and `n = n(|B| …
- [practice_drift_mutation_suite](deferred/practice_drift_mutation_suite.md) — follow-ups for the practice-drift checker's mutation suite
- [side_edge_fusing_cut_aware_plan](deferred/side_edge_fusing_cut_aware_plan.md) — Make the thin-shell rim fuse ("fuse-to-g") cut-consistent: today, different interface levels of the SAME rim station are constrained to different branches of the multivalued air trace g (with/without the free-cut λ) …
- [slave_side_enrichment](deferred/slave_side_enrichment.md) — Analysis of what it takes to enrich both master and slave sides of thin shell sidesets
- [slepc_eigensolver_integration](deferred/slepc_eigensolver_integration.md) — Interface SLEPc as a second eigenvalue backend alongside the (now repaired) ARPACK-ng driver, so BELFEM can resolve the *small* end of the spectrum — which mode-1 Arnoldi structurally cannot — via a spectral transform …
- [thinshell_automatic_sublayering](deferred/thinshell_automatic_sublayering.md) — Future enhancement — automatic geometric subdivision of thick layers
- [thinshell_hphi_formulation](deferred/thinshell_hphi_formulation.md) — Capture the motivation, theory, and implementation sketch for adding a per-layer H–φ formulation switch to BELFEM's thin-shell stack, so that insulator and high-contact-resistance interlayers can be modelled without the …
- [thinshell_selective_nitsche_coupling](deferred/thinshell_selective_nitsche_coupling.md) — Long-term architecture improvement for high-contrast thin-shell layers
- [user_api_header_install](deferred/user_api_header_install.md) — replace the 628-header install with a flat 14-header user API

## Cancelled (11) — [`cancelled/`](cancelled/)

- [arpack_small_end_configuration](cancelled/arpack_small_end_configuration.md) — Make the `compute conditioning` eigen fallback actually produce a number for the thermal system, and stop it burning wall clock when it cannot. Mechanism: size `ncv` / `maxit` from `( n, nev, which end )` instead of one …
- [conditioning_diagnostic_backends](cancelled/conditioning_diagnostic_backends.md) — Give the `compute conditioning` diagnostic a solver-native backend beyond MUMPS, using the Arioli/Demmel/Duff error analysis those libraries already implement, so it stops depending on an eigenvalue estimate that cannot …
- [conditioning_petsc_native](cancelled/conditioning_petsc_native.md) — Let a PETSc-backed field report a conditioning number from the Krylov solve it is already doing, instead of paying for a shift-invert eigen round. Mechanism: `KSPSetComputeSingularValues()` before the solve …
- [controller_picard_tapestack_regression](cancelled/controller_picard_tapestack_regression.md) — Gregory's 2D tapestack case (magnetic-only, `algorithm : Picard`) ran fine on `matfix` and "struggles a lot" on `sideconnectors`, while quench cases improved. Mechanism in one sentence: quench-tuned timestep-cut guards …
- [eigen_shift_invert_solver_reuse](cancelled/eigen_shift_invert_solver_reuse.md) — cut the eigen diagnostic's second resident MUMPS factorization
- [iterate_refactor_plan](cancelled/iterate_refactor_plan.md) — De-duplicate the three `iterate*()` functions and fix the drift between them, without changing convergence behavior (except where explicitly flagged).
- [rho_lambda_b_first_reorder](cancelled/rho_lambda_b_first_reorder.md) — reorder the field-dependent `rho`/`lambda` API to B-first
- [sch04_history_purge](cancelled/sch04_history_purge.md) — purge the Fujikura `sch04` tables from git history
- [thin_shell_overhang_bug_analysis](cancelled/thin_shell_overhang_bug_analysis.md) — Analyze why cohomology cuts fail when thin shells extend beyond a solid conductor block
- [timestep_collapse_mitigation_design](cancelled/timestep_collapse_mitigation_design.md) — Concrete, code-level mitigation design for the Δt-collapse failure analysed in `todo/timestep_collapse_residual_floor_plan.md` (the `sidecoatings` t = 13.8762 ms abort and the `greg5` regression). Mechanism being fixed …
- [timestep_collapse_residual_floor_plan](cancelled/timestep_collapse_residual_floor_plan.md) — Both `greg5` (2D tapestack, magnetic-only Picard) and `sidecoatings` (3D coated conductor, coupled Newton) ratchet Δt down to the deck minimum and abort. The mechanism, in one sentence: **the controller's promote/demote …

## Closed (114) — [`closed/`](closed/)

- [yufan_ni_coil_fixes](closed/yufan_ni_coil_fixes.md) — DONE 2026-10-08: Yufan's 3-D thin-shell alignment landed with a block-consistent patch relation and a block-on-both-sides refusal; closed side curves accepted when isolated (DR-188 policy narrowed on Alves 2022a §II-B); `pancake_tape` ships with its solder/air sideset; make check 17/17 ×3; copper-lead gate owed on DR-196
- [2d_thinshell_gap_analysis](closed/2d_thinshell_gap_analysis.md) — Procedure-by-procedure comparison of the 2D h-φ thin-shell pipeline against the 3D reference path; inventory of what is MISSING, STUBBED, or WRONG in 2D.
- [ac_loss_postprocessing](closed/ac_loss_postprocessing.md) — energy, decomposition and honesty on top of the power global
- [ai_wrapper_tier_selection](closed/ai_wrapper_tier_selection.md) — explicit model + effort selection for `ask_codex.sh` / `ask_grok.sh`
- [alloy_database_al_cr_sae310](closed/alloy_database_al_cr_sae310.md) — Extend the elemental material database with Aluminum and Chromium so that the `Alloy` homogenization tool can synthesize HTS-relevant alloys (first target: SAE 310 stainless), with lambda and rho data completed wherever …
- [anderson_picard_acceleration_plan](closed/anderson_picard_acceleration_plan.md) — Stop ts17-class nonlinear stalls at their mechanism: first repair the controller's relaxation/switching state machine (the E/D/F consensus bundle from the ts17 campaign), then add opt-in Anderson(m) mixing to the …
- [bdf_jacobian_scaling_bug](closed/bdf_jacobian_scaling_bug.md) — Document critical bug in Newton Jacobian assembly for BDF2-BDF5 time-stepping schemes
- [bdf_nonlinear_mass_verification](closed/bdf_nonlinear_mass_verification.md) — The two forward-looking items left open after the BDF Jacobian scaling bug was implemented and closed (`closed/bdf_jacobian_scaling_bug.md`). The α-scaling, coefficient wiring, startup ramp, and the reset/savepoint …
- [beta_angle_normal_sign](closed/beta_angle_normal_sign.md) — Document potential bug in β angle computation for HTS thin shells
- [bfm_stale_cache_detection](closed/bfm_stale_cache_detection.md) — A `.bfm` file caches an *enriched* mesh — cuts, thin-shell layers, edge-coating walls, periodicity, hanging entities. The reuse test compares only a checksum of the **raw** geometry, so changing an input option that …
- [circuit_demo_harvest_plan](closed/circuit_demo_harvest_plan.md) — harvest the `electricalCircuit` demo circuits into transient tests
- [citation_index_completion](closed/citation_index_completion.md) — closing the gap between the works the tree cites and the works it can identify
- [code_findings_from_comment_sweeps](closed/code_findings_from_comment_sweeps.md) — defects the comment auditors found while reading for contracts, parked for a code session
- [comment_cleanup_sweep_1](closed/comment_cleanup_sweep_1.md) — bring the existing `src/` comments under `doc/commenting_guidelines.md`
- [comment_cleanup_sweep_2](closed/comment_cleanup_sweep_2.md) — add the missing contracts on the public entry points, correct the comments the code contradicts, finish the narration tail
- [conditioning_shift_invert_fallback](closed/conditioning_shift_invert_fallback.md) — Implement Christian's ruling of 2026-08-28: the conditioning diagnostic reports the metric the active solver can supply natively, and falls back to a **shift-invert** eigenvalue computation when it cannot. Memory cost …
- [corner_edge_transition_problem_2026-03-13](closed/corner_edge_transition_problem_2026-03-13.md) — Handoff note for continuing investigation tomorrow
- [coulomb_gauge_penalty_report](closed/coulomb_gauge_penalty_report.md) — Closing report of the overnight task "Coulomb Gauge Penalty for the h-phi Formulation (theory + G-operator on the Nédélec edge functions)". Records what was derived, what was audited and by whom, the disagreements and …
- [cross_review_tooling](closed/cross_review_tooling.md) — Turn the existing manual Claude + Codex + Grok review protocol into a `/cross-review` slash command, an opt-in commit-triggered light review, and a public `METHODOLOGY.md`. Freeze rather than redesign: the scripts …
- [current_sign_2d_fix](closed/current_sign_2d_fix.md) — 2-D current-BC sign inversion (declared +I runs as -I)
- [curvefactory_dr188](closed/curvefactory_dr188.md) — CurveFactory: pairs from facets on every path, closed curves oriented, closed terminals refused
- [cut_pocket_removal_rules](closed/cut_pocket_removal_rules.md) — Agree on the rules for a pocket-removal / representative-canonicalization post-pass on rectified cohomology generators, before any code is written. The pass removes closed, support-only cut fragments by firing exact …
- [doc_currentness_fixes](closed/doc_currentness_fixes.md) — fix list from the 96-file documentation sweep
- [doxygen_site_structure](closed/doxygen_site_structure.md) — Make the generated Doxygen site navigable for three readers — a new user who wants to run a simulation, a new developer who wants to extend a module, and a returning developer looking something up. Wave 1 (the three …
- [dr105_picard_stall_plan](closed/dr105_picard_stall_plan.md) — DR-105: the pure-Picard T5 quench deck aborts on the floor-retry limit because the stagnation guard ends every attempt at iterate 6 and the floor retry raises budgets the guard never consumes while resetting ω (a replay …
- [dr106_dr144_joint_session](closed/dr106_dr144_joint_session.md) — close DR-106's `-9` retry together with DR-144's deferred R12/R13
- [dr111_quad4ts_winding](closed/dr111_quad4ts_winding.md) — DR-111: the 2-D thin-shell layer elements are wound clockwise
- [dr113_dr190_dr191_residues](closed/dr113_dr190_dr191_residues.md) — the residues of the small-fixes cluster: `IWG::detach_field()` across a caught `create_field` refusal (DR-190), the MUMPS pool's `KEEP(40)` zeroed at allocation because the 5.9.1 driver reads it before `JOB = -1` …
- [dr128_spap_low_field_axis](closed/dr128_spap_low_field_axis.md) — extend the sp-ap B axis below 10 mT
- [dr129_kernel_controller_init](closed/dr129_kernel_controller_init.md) — DR-129: `Kernel::mController` is never initialized
- [dr144_mumps_lifecycle_plan](closed/dr144_mumps_lifecycle_plan.md) — close DR-144's seven MUMPS instance-lifecycle gaps
- [dr147_plugin_path_resolution](closed/dr147_plugin_path_resolution.md) — DR-147: one search order for every plugin `.so`
- [dr18x_cluster](closed/dr18x_cluster.md) — DR-184 / DR-182 / DR-186 / DR-187 / DR-166 as five small rounds
- [dr190_dr191_small_fixes](closed/dr190_dr191_small_fixes.md) — DR-190: the `DofManager` owned across a `create_field` refusal (the rank-dependent-gate half is refuted by the source); DR-191: `mBSKohlerTransSwitch( 0 )` never assigned in Lead and WhiteTin, localized by valgrind
- [dr92_restart_step_ramp](closed/dr92_restart_step_ramp.md) — Theory of the DR-92 restart anomaly and the proposed re-entry cap. [ Title and framing corrected 2026-08-19: this is NOT a growth limit. The working hypothesis, not proof, is which matrix STRUMPACK's solver-entry state …
- [elastic_moduli_quasiharmonic](closed/elastic_moduli_quasiharmonic.md) — ν(T) and E(T) from the thermal strain instead of the constant-γ inversion
- [error_path_hygiene_cluster](closed/error_path_hygiene_cluster.md) — DR-165, DR-167, DR-175: the HDF5 loaders' dataset handles, the periodicity factory's centroid nodes and object under construction, the material and boundary-condition factories' products, each owned across its refusal
- [example_deck_and_material_db_repair](closed/example_deck_and_material_db_repair.md) — Make a fresh checkout of BELFEM runnable by someone who has never run it — specifically on the two paths that fail today for anyone without pre-existing cached artifacts: a parallel first run (aborts in METIS during …
- [falsification_tooling](closed/falsification_tooling.md) — Institutionalize the falsification side of the workflow: a permanent interface/orientation regression battery (`make check-fast`, both backends), an evidence hierarchy wired into the collaboration protocol and …
- [ferro_cut_flux_island](closed/ferro_cut_flux_island.md) — Ferro-cut flux island (gantry dead yoke arm) — investigation brief
- [ferro_undulator_mu_and_interface_bugs](closed/ferro_undulator_mu_and_interface_bugs.md) — Two defects reported by Gregory from a ferro undulator benchmark, with his two workarounds assessed. The first is a real debug-only crash whose proposed fix introduces undefined behaviour that silently turns the iron …
- [fix_piecewise_degenerate_window](closed/fix_piecewise_degenerate_window.md) — Eliminate the NaN that killed `tapestack3d_coarse` at timestep 122: for table-driven jc/n materials, the flux-flow Bézier machinery in the `rho_piecewise` family produces NaN when the table's n(T,B,θ) approaches 1 …
- [gas_correctness_fixes](closed/gas_correctness_fixes.md) — Checklist of the defects found by the two-round literature verification of `src/physics/{gastables,gasmodels}` against the primary sources (devlog `dl20260806_gas_literature_verification.md`). Every constant and formula …
- [gasmodels_open_source_migration](closed/gasmodels_open_source_migration.md) — Move `nonfree/physics/{gasmodels,gastables,atmosphere}` into the open-source `src/physics/` tree without redistributing third-party data.
- [gauge_newton_tangent_plan](closed/gauge_newton_tangent_plan.md) — consistent Newton tangent for the Coulomb gauge penalty
- [ghost_penalty_parameter_provenance](closed/ghost_penalty_parameter_provenance.md) — Trace where the two hardcoded ghost-stabilization constants in `maxwell::h_ghost()` come from, whether they are literature-grounded, and how to choose better-motivated values.
- [ghost_swip_average](closed/ghost_swip_average.md) — DR-168: the SWIP average flux weight replaces the slope jump in `h_ghost`, the coefficient tests follow, the manufactured-flux fixture (R7) lands, the doc and the kernel comment stop disagreeing on the α floor
- [handoff_20260831_open_defects](closed/handoff_20260831_open_defects.md) — the non-documentation findings of that campaign
- [handoff_double_corc_d1_session](closed/handoff_double_corc_d1_session.md) — Hand off to a fresh session. The periodic-seam campaign is complete; the double-layer corc now reaches its first Jacobian assembly and dies at the KNOWN D1 defect of the Maxwell kernel-collapse plan. This file contains …
- [handoff_thermal_iwg_collapse_session](closed/handoff_thermal_iwg_collapse_session.md) — Self-contained brief for a FRESH Fable session to tidy up `cl_IWG_MaxwellThermal.{hpp,cpp}` and `mt_thermal_h.{hpp,cpp}` the same way the Maxwell side was collapsed on 2026-07-21: dispatch tree -> collapsed-kernel pick …
- [hex8tb_phase2_fem_wiring](closed/hex8tb_phase2_fem_wiring.md) — updated 2026-08-09 — not swept
- [instruction_doc_currency](closed/instruction_doc_currency.md) — The 2026-08-11 `CLAUDE.md` sweep found nine claims that were false against the tree, four of which would have produced wrong code if followed. All nine are fixed. This plan records what was found, closes the residual …
- [interface_node_duplication_coil_ferro](closed/interface_node_duplication_coil_ferro.md) — Extend the `InterfaceProcessor` node-duplication scheme (currently exercised for the φ/h conductor interfaces) to ferro-air and air-coil interfaces, so ParaView output shows the physical field jump at the iron surface …
- [l21_wrapper_card_repair](closed/l21_wrapper_card_repair.md) — the wrapper-policy card cites a contract that moved, and cannot enforce itself
- [linalg_backend_verification](closed/linalg_backend_verification.md) — two halves of one problem, both settled on 2026-08-30.
- [make_doc_repair](closed/make_doc_repair.md) — repair `make doc`: Doxygen front-end breaks and the 1.18-vs-1.9.1 skew
- [materials_tangent_cluster](closed/materials_tangent_cluster.md) — debt cluster A: DR-163, DR-164, DR-170, DR-173, DR-174 in `src/physics/materials`
- [matfix_controller_port_plan](closed/matfix_controller_port_plan.md) — Harvest the Jul 27–28 controller/solver hardening that exists only on `matfix` (commits `8c52161f`, `7e5df34c`, `2371339c`) onto `sideconnectors`, reconciled with the ts17/Anderson work that exists only there. Root …
- [matlab_derivations_recovery](closed/matlab_derivations_recovery.md) — recover the MATLAB generator and verifier scripts from the pre-reset `tmp/` into the module `doc/` trees
- [matvec_csc_omp_stack_overflow](closed/matvec_csc_omp_stack_overflow.md) — `matvec_csc` overflows the OpenMP worker stack above 262,144 rows
- [maxwell_background_ds_dr181](closed/maxwell_background_ds_dr181.md) — `background_phi` read the point-0 surface increment at every point
- [maxwell_kernel_collapse_plan](closed/maxwell_kernel_collapse_plan.md) — Collapse the 38 near-identical matrix-assembly variants (25 `maxwell::h_*` in `mt_maxwell_h.{hpp,cpp}`, 13 `T_h_*` in `mt_thermal_h.{hpp,cpp}`) into a small set of generic kernels that delegate all per-integration-point …
- [maxwell_postprocessor_element_caching](closed/maxwell_postprocessor_element_caching.md) — CPU-only optimization to eliminate redundant physics computations
- [maxwell_postprocessor_review](closed/maxwell_postprocessor_review.md) — Document review findings for Maxwell postprocessor performance and architecture
- [mesh_oneliners_dr161_169_176_179_180](closed/mesh_oneliners_dr161_169_176_179_180.md) — five `src/mesh` one-liners from the debt register
- [meshfile_refactor_plan](closed/meshfile_refactor_plan.md) — Plan to repair and refactor the `.bfm` (HDF5) mesh save/load around the `ProtoMesh`/`Protoshell` reconstruction engine so a Maxwell run can stop, change parameters, and relaunch from the last stored timestep — loading a …
- [mpi_link_probe_plan](closed/mpi_link_probe_plan.md) — keep the system libdirs off the rpath so the SCLS PMIx wins; link a minimal MPI program at configure time so a broken MPI fails at cmake
- [mumps_error_analysis_opt_in](closed/mumps_error_analysis_opt_in.md) — split the MUMPS ADD diagnostic out of `compute conditioning`, make it opt-in
- [mumps_memory_budget_plan](closed/mumps_memory_budget_plan.md) — probe the machine, hand MUMPS ICNTL(23), stop guessing ICNTL(14)
- [mumps_workspace_cap_raise](closed/mumps_workspace_cap_raise.md) — The `-9` (out-of-workspace) retry ladder tops out at ICNTL(14) = 240 after three doublings. On the `tape_quench_usermat` deck the ladder is exhausted while the missing workspace is still only ~2.4 MB, the step is …
- [nedelec_edge_function_defects](closed/nedelec_edge_function_defects.md) — Track two edge-function defects surfaced as by-catch of the Nédélec LaTeX extraction audit round (`closed/nedelec_tex_extraction_plan.md`, §4.0/§4.0b). Both are in the `E` interpolation operators; the `C` curl operators …
- [nedelec_tex_extraction_plan](closed/nedelec_tex_extraction_plan.md) — The pre-BELFEM theory notes in `tmp/nedelec/` (LaTeX, Lagrange-multiplier era) contain a weak-form derivation and a Nédélec-element derivation that are still correct and worth keeping. Extract those two parts into the …
- [ngspice_parser_plan](closed/ngspice_parser_plan.md) — promoted from `deferred/`
- [nitrogen_eos_completion](closed/nitrogen_eos_completion.md) — Everything still needed to finish `EoS_Nitrogen` and wire it into `Gas`. Written as reconstruction insurance after an editor save failed to reach disk — the three outstanding derivative bodies are given in full below …
- [nonlinear_iteration_strategy_near_quench](closed/nonlinear_iteration_strategy_near_quench.md) — Theory analysis and recommendations for improving nonlinear solver performance near critical material transitions (quench)
- [orientation_claude](closed/orientation_claude.md) — Validate that facet integration points and orientation tables are correct
- [parmetis_ptscotch_wiring](closed/parmetis_ptscotch_wiring.md) — wire the parallel nested dissection (F2)
- [penalty_opt_in_plan](closed/penalty_opt_in_plan.md) — honour `eta : 0` (no ghost facets), make the gauge — and, as a second step, the ghost — opt-in
- [periodic_bc_fix_plan](closed/periodic_bc_fix_plan.md) — Step-by-step plan to fix periodic boundary condition handling across DOF conversion, MPI distribution, and input wiring.
- [periodic_cap_cut_emission](closed/periodic_cap_cut_emission.md) — Fix the solve-phase periodic rebuild abort — "number of edges on source and target periodic surfaces do not match" — root-caused (three-AI, 2026-07-16) to one-sided emission of thin-cut geometry on the periodic …
- [periodic_dofs_as_hanging_dofs](closed/periodic_dofs_as_hanging_dofs.md) — Capture the design decision that periodic DOFs are a special case of hanging DOFs, and what that implies for parallel (MPI) handling. Draft for review; to be folded into a `doc/` section once agreed.
- [periodic_input_extension_plan](closed/periodic_input_extension_plan.md) — Add explicit input-file syntax for declaring periodic sideset pairs while preserving the existing plane-based periodic definition.
- [periodic_parallel_source_closure_fix](closed/periodic_parallel_source_closure_fix.md) — Fix the distributed-mesh inconsistency that causes the z-periodic CORC model to abort in `MUMPS::solve()` on a worker rank, while the same model solves correctly in serial.
- [periodic_thin_cut_continuity_fix](closed/periodic_thin_cut_continuity_fix.md) — Implementation plan for the `phi`-continuity error that appears only when cohomology thin cuts are used with periodic boundary conditions.
- [periodic_volume_edge_seam_fix](closed/periodic_volume_edge_seam_fix.md) — Fix the silently missing periodic constraints on bulk-conductor edge dofs at periodic planes (helix example: J pattern rotated 90° across the seam, r = −0.96, while net current and air field are correct).
- [picard_line_head_residual](closed/picard_line_head_residual.md) — DR-192 item 2: the per-iterate `Magnetic Picard` line prints the head residual ( what the accept reads ) beside the lagged one ( what the guards act on ), so a −25 dB line before a converged `timestep succeeded` stops …
- [pid_timestep_controller_plan](closed/pid_timestep_controller_plan.md) — authored 2026-08-09 — not swept (it moved DRAFT → IMPLEMENTED while this sweep was running; entry below refreshed to match)
- [piecewise_bezier_root](closed/piecewise_bezier_root.md) — one clamped Bézier root for the piecewise power law's value and both tangents
- [practice_drift_checker](closed/practice_drift_checker.md) — a checker for convention documents that have drifted from practice
- [release093_precheck_fixes](closed/release093_precheck_fixes.md) — the three code findings of the 0.9.3 release pre-check: a content revision on the cached rho tables so the corrected Kohler laws are not shadowed by a 0.9.2 cache, `DofData::reset` deleting each dof once after a refusal …
- [release093_ultra_followups](closed/release093_ultra_followups.md) — the two findings of the `/code-review ultra` pass over the 0.9.3 branch that Christian ruled to fix: `DofData::reset` releases the mesh dof containers before deleting the dofs of a refused field (the normal path …
- [restart_circuit_state](closed/restart_circuit_state.md) — Persist the coupled electrical circuit's dynamic state in the `memdump.hdf5` restart file so a CORC restart resumes with a **warm circuit** instead of re-creating it cold from `input.conf`. This is the one restart item …
- [restart_circuit_verification](closed/restart_circuit_verification.md) — Close out the three verification/documentation items left after the coupled-circuit dynamic-state restart was **implemented** (`{time, delta_time, x, prev_x}` in the `/circuit` group of `memdump.hdf5`, scatter on load …
- [restart_timestep_consistency](closed/restart_timestep_consistency.md) — OBSOLETE. Both
- [rho_lambda_argument_convention](closed/rho_lambda_argument_convention.md) — Decide how to handle the deliberate but footgun-prone opposite argument conventions of the magneto-dependent `rho` and `lambda` families.
- [ruled_refusals_cluster](closed/ruled_refusals_cluster.md) — DR-162, DR-171, DR-172: the side-connector guard and refusal in both linkers, the material-deck contradiction refusals with both input artifacts, the narrow multi-source Dirichlet refusal
- [run_gate_batches](closed/run_gate_batches.md) — OBSOLETE. Derived from a classification of
- [shared_library_and_install_plan](closed/shared_library_and_install_plan.md) — Turn BELFEM from a build-tree-only, static-archive project into something that can be installed and run from `$CMAKE_INSTALL_PREFIX`. Three coupled changes: build the project as a shared library, add `install()` rules …
- [side_connector_effective_resistivity](closed/side_connector_effective_resistivity.md) — Self-contained problem statement and modeling recommendation for copper side conductivity of HTS coated-conductor tapes (current escape around the insulating buffer layer), replacing three abandoned geometric/enrichment …
- [sideconnector_bridge_plan](closed/sideconnector_bridge_plan.md) — Implement the collapsed edge-wall impedance bridge of `side_connector_effective_resistivity.md` (rev. 4) with O1 resolved to option (c) — the minimal degenerate wall element (N=1): corner side edges stay anchored to the …
- [smith_inverse_grouped_terminals](closed/smith_inverse_grouped_terminals.md) — branch `Claude` experiment: unit-normalized Smith diagonal + weighted cut jump for bracket groups whose curves cross one conductor k times (Rutherford 5×); evidence for Gregory, never merged. ✅ COMPLETE 2026-09-19 …
- [tapestack3d_J_asymmetry](closed/tapestack3d_J_asymmetry.md) — Localize the source of the jagged, z-periodic |J| streaks seen in one outer copper layer of the tapestack3d example, while the other outer copper layer shows smooth z-invariant bands.
- [tapestack3d_jjc_noise_2125ms](closed/tapestack3d_jjc_noise_2125ms.md) — j/jc noise onset at t = 2.125 s: verdict + mitigation
- [test_gap_tier1_execution](closed/test_gap_tier1_execution.md) — the Tier 1 test-gap overnight run
- [test_normals_and_pipette_coverage](closed/test_normals_and_pipette_coverage.md) — Rebuild the facet-normal and volume-integral (pipette) tests that were lost when the `geometry::` namespace was retired, using the surviving reference fixture.
- [thermal_matrices_cleanup_and_newton_plan](closed/thermal_matrices_cleanup_and_newton_plan.md) — Restructure `src/fem/thermal/matrices/mt_thermal_h.cpp` (13 copy-paste producers, 1314 lines) into single-sourced formulas, and add the missing temperature derivatives (∂M/∂T, ∂K/∂T, ∂f/∂T) so the thermal Newton solver …
- [thermal_picard_freeze_dr08](closed/thermal_picard_freeze_dr08.md) — Record the first clean reproducer of the DR-08 coupled freeze, and propose a fix for both the freeze and the mechanism that exposed it
- [thin_cut_nonunit_rectification_implementation](closed/thin_cut_nonunit_rectification_implementation.md) — Replace the greedy `Cohomology::clean()` with a global feasibility solve that rectifies a non-unit thick-cut cochain (`|c(e)| ≥ 2`) to unit coefficients when a unit representative exists, or returns a constructive …
- [thinshell_conductor_normal_field](closed/thinshell_conductor_normal_field.md) — thin-shell h_n recovery next to an h-conductor volume (corc_solder)
- [thinshell_nodal_recovery](closed/thinshell_nodal_recovery.md) — correct nodal H, B, J in curved thin-shell layers
- [thinshell_recovery_unit_test](closed/thinshell_recovery_unit_test.md) — unit tests that run the node-local recovery on code-built meshes
- [timestep_controller_cleanup](closed/timestep_controller_cleanup.md) — Record the resolution of the disabled error-controller `#HACK` in `Controller::adjust_timestep()` and park the larger controller/integrator ideas that are out of scope for now.
- [user_material_backend_boundary_refactor](closed/user_material_backend_boundary_refactor.md) — User Material Backend Boundary Refactor
- [weekly_report_20260828](closed/weekly_report_20260828.md) — week in review, 22 to 28 August 2026

### Test-suite planning campaign (15) — [`closed/tests/`](closed/tests/)

- [existing_tests_triage](closed/tests/existing_tests_triage.md) — Assess existing tests in `./test/`, determine what to keep, salvage, or retire before implementing the new test plans (tests_00 through tests_13).
- [tests_00_strategy](closed/tests/tests_00_strategy.md) — Master test strategy for the BELFEM framework using Google Test
- [tests_01_containers](closed/tests/tests_01_containers.md) — Method-level test matrix for all containers in `src/containers/`
- [tests_02_linalg](closed/tests/tests_02_linalg.md) — Method-level test matrix for the `src/linalg/` module
- [tests_03_quaternion](closed/tests/tests_03_quaternion.md) — Method-level test matrix for the quaternion module
- [tests_04_tensor](closed/tests/tests_04_tensor.md) — Method-level test matrix for the tensor module (`src/math/` or `src/tensor/`)
- [tests_05_mathtools](closed/tests/tests_05_mathtools.md) — Method-level test matrix for the math tools module
- [tests_06_comm](closed/tests/tests_06_comm.md) — Method-level test matrix for the comm module (`src/comm/`)
- [tests_07_graph](closed/tests/tests_07_graph.md) — Method-level test matrix for the graph module (`src/math/graph/`)
- [tests_08_sparse](closed/tests/tests_08_sparse.md) — Method-level test matrix for the sparse module (`src/sparse/`)
- [tests_09_spline](closed/tests/tests_09_spline.md) — Method-level test matrix for the spline module (`src/numerics/spline/`)
- [tests_10_mesh](closed/tests/tests_10_mesh.md) — Method-level test matrix for the mesh module (`src/mesh/`)
- [tests_11_interpolation](closed/tests/tests_11_interpolation.md) — Method-level test matrix for the interpolation module (`src/fem/interpolation/`)
- [tests_12_integration](closed/tests/tests_12_integration.md) — Method-level test matrix for the integration module (quadrature rules)
- [tests_13_orientation](closed/tests/tests_13_orientation.md) — Test plan for facet orientation permutations at element interfaces
