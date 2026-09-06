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

**Convention:** new substantial plans follow **[plan_template.md](plan_template.md)** (structure
extracted from `closed/meshfile_refactor_plan.md`): header block with living Status line, scope
guards, gap table with (a)/(b)/(c) classes, ordered `Rn` steps, `Dn` defect tracker with live
checkboxes, `On` open questions, definition-of-done, audit trail.

---

## Currentness sweep 2026-09-03 — round 3: purpose-served and relevance

**Active plans: 4** (6 after round 3; `user_api_header_install` split and deferred 2026-09-04,
`linalg_backend_verification` closed 2026-09-04 — see their entries). Rounds 1 and 2 (below) took the directory from 57 to 29 by asking whether a
plan's work was in the tree. Christian's review of the 29 set a sharper standard, and round 3
applied it: **a plan is closed when its purpose is served in the tree, even if its own refined
scope is not** (`ac_loss_postprocessing` is the model case — `dotQ` is collected every step and
written to the exodus, so the energy/decomposition refinements do not earn a live plan), and **a
plan that is genuinely undone but that nobody is blocked on is cancelled, not carried.**

Round 3 moved 23. To `closed/` (12): `ac_loss_postprocessing`, `conditioning_shift_invert_fallback`,
`fix_piecewise_degenerate_window`, `gasmodels_open_source_migration`,
`interface_node_duplication_coil_ferro`, `mumps_memory_budget_plan`, `penalty_opt_in_plan`,
`restart_circuit_verification`, `shared_library_and_install_plan`,
`thin_cut_nonunit_rectification_implementation`, `thinshell_conductor_normal_field`, and
`test_normals_and_pipette_coverage` (folded into the test-hardening plan, since removed from the
tree). To
`cancelled/` (9): the two conditioning-backend plans, `controller_picard_tapestack_regression`
(overtaken by the rebuilt controller), `eigen_shift_invert_solver_reuse`, `iterate_refactor_plan`,
`thin_shell_overhang_bug_analysis` (February; no user case has pulled it), both
`timestep_collapse_*` plans (superseded by the κ ∝ Δt finding), and `rho_lambda_b_first_reorder`
(Christian's ruling: keep T-first). To `deferred/` (2): `input_conf_configurator_plan`,
`powerlaw_jc_n_field_derivatives`. Every moved file carries a dated banner under its title saying
which bucket it went to and why, so its stale Status line cannot mislead.

**Those that remained:** `circuit_demo_harvest_plan` (DR-39, open row → `closed/` on 2026-09-04:
executed end to end, tests RAN green, executable deleted, DR-39 STRUCK), `linalg_backend_verification`,
`parmetis_ptscotch_wiring` (→ `closed/` on 2026-09-04: wired, green at np 2/4, benchmark → DR-156),
`user_api_header_install` (→ `deferred/` on 2026-09-04 after its cheap items landed). The Doxygen
sweep residue, the test-hardening campaign plan and the cohomology handoff for Gregory Giard were
removed from the tree on 2026-09-05 as internal working material. `linalg_backend_verification` → `closed/` on
2026-09-04 after a scope jury: the flavor presets are the deliverable, now with `USE_PARDISO` in
parity with `USE_MKL` (a flavor name containing `mkl` presets both ON, anything else OFF); the
link probe (R2–R6) and any flavor-versus-knob refusal are declined, since every non-MKL vendor is
assumed to behave like reference BLAS/LAPACK and an explicit `-D` always wins.

By-catch fixed the same day: `config/globals.cmake` was still installing the retired `hphirun` and
`hphiTrun`; removed, `check_doc_claims.py` 37/37.

## Round 2 (superseded by round 3, kept for the record)

**Round 1 was wrong and was redone.** It triaged by reading each plan's Status line and counting
unchecked boxes, and moved 9 files. Christian rejected the result: far more plans were obsolete.
The bias was directional — a plan's text is written when work is *planned* and rarely updated
when it *lands*, so unchecked boxes and OPEN status lines both systematically overstate remaining
work. Round 2 asked the opposite question of every plan — *is the work already in the tree?* —
and required a grep result, not a quotation, for every verdict.

**Round 2 moved 19 more. Active fell 57 → 29.** To `closed/` (substance landed; remnant is a run
gate or a sign-off, named in each Status line): `matvec_csc_omp_stack_overflow`,
`nitrogen_eos_completion`, `current_sign_2d_fix`, `dr111_quad4ts_winding`,
`nedelec_edge_function_defects`, `pid_timestep_controller_plan`, `mumps_error_analysis_opt_in`,
`thermal_matrices_cleanup_and_newton_plan`, `ngspice_parser_plan`, `maxwell_kernel_collapse_plan`,
`hex8tb_phase2_fem_wiring`, `falsification_tooling`, `make_doc_repair`, `gauge_newton_tangent_plan`,
`ferro_undulator_mu_and_interface_bugs`, `dr128_spap_low_field_axis`,
`rho_lambda_argument_convention`. To the new `cancelled/`: `arpack_small_end_configuration` (the
spectral fold it designs was **replaced** by shift-invert — `cl_FEM_DofMgr_EigenValues.hpp:129-136`
says so in the source) and `sch04_history_purge` (Christian's ruling: the history rewrite is
abandoned; the working-tree removal is the whole of what shipped).

**Plans whose text was flatly contradicted by the tree** — the yield that box-counting could never
reach. `nedelec_edge_function_defects` said "No source modified" while both its fixes and a 22/22
circulation battery were in the tree. `dr128_spap_low_field_axis` was framed as unimplemented
while the low-field transition sits at `cl_JcFunction_Database.hpp:440-480`. `nitrogen_eos_completion`
claimed three unwritten bodies and a file absent from `SOURCES`; all nine exist, it is in the build
at `CMakeLists.txt:14`, and `cl_Gas.cpp:845` constructs it. `thermal_matrices_cleanup_and_newton_plan`
awaited `T_h_newton`, which is 77 lines at `mt_thermal_h.cpp:57`. `maxwell_kernel_collapse_plan`
asserted no documentation mentions its kernels; false. Four more Status lines were corrected in
place during round 1 (`mumps_error_analysis_opt_in` contradicted itself in one file;
`powerlaw_jc_n_field_derivatives` gated on struck DR-69; `gauge_newton_tangent_plan` recorded a
`chi` default reversed on 2026-09-01).

**The 29 that remain are LIVE by tree evidence** — each has a named thing that was searched for and
not found (`belfem_probe_link`, `test_Normals.cpp`, `parmetis_nd` callers, `drho_powerlaw_dbeta`,
`belfemConfig.cmake`, Phase-4 min-cost L1, the `stall slope` key, a PETSc-native metric).

Method and full evidence: `devlog/dl20260903_todo_currentness_sweep.md`; auditor rows in
`tmp/ai_exchange/todo_currentness_sweep.md` (Codex `gpt-5.6-terra`/high and Grok `grok-4.6`/high,
disjoint batches of 24).

## Round 1 (superseded, kept for the record)

All 57 active plans triaged (Claude, with Codex `gpt-5.6-terra`/medium on 29 and Grok
`grok-4.6`/high on 28, disjoint batches; record in `tmp/ai_exchange/todo_currentness_sweep.md`).
**Nine relocated, 48 remain active.** To `closed/`: `dr106_dr144_joint_session`,
`dr144_mumps_lifecycle_plan`, `doc_currentness_fixes`, `handoff_20260831_open_defects`,
`weekly_report_20260828`, and `mumps_workspace_cap_raise` (superseded — moved on the
instruction in its own successor, `mumps_memory_budget_plan.md`). To `deferred/`:
`cut_representative_options` (Gregory + Frederic: thick cut is the eventual target, not a
current priority), `side_edge_fusing_cut_aware_plan` (debt closed, remainder is research that
gates nothing in the ship configuration), `mumps_symmetric_triangle_extraction` (Christian:
"not today", no source touched).

**No `cancelled/` directory was created.** Exactly one file was a supersession and its own
successor plan directs it to `closed/`; a third bucket holding one document costs more in
lookup than it saves. Raise it again if a genuinely abandoned plan appears.

**Four Status lines were factually wrong and are corrected in place** — the sweep's real yield,
since a wrong Status is what makes a plan unreadable without opening it:
`nitrogen_eos_completion` (claimed unfinished, unbuilt and unwired; is in fact written, in
`SOURCES` and constructed at `cl_Gas.cpp:845`), `mumps_error_analysis_opt_in` (header said
"only R8b remains, blocked on O3" while both are ticked in the same file),
`powerlaw_jc_n_field_derivatives` (β gated on DR-69, which is struck) and
`gauge_newton_tangent_plan` (records a `chi = 1e-4` default reversed to opt-in on 2026-09-01).
Nine stale `todo/` cross-references were repaired across the live plans, the README and
`debt_register*.md`; every README link now resolves.

**Seven plans are finished except for one owed gate and need Christian's DR-42/DR-49 ruling**
rather than a move: `current_sign_2d_fix`, `dr111_quad4ts_winding`,
`nedelec_edge_function_defects`, `pid_timestep_controller_plan`, `falsification_tooling`,
`gauge_newton_tangent_plan`, `make_doc_repair`. And `rho_lambda_argument_convention` (T-first
doc pass) directly contradicts `rho_lambda_b_first_reorder` (Christian's B-first reversal) —
one of the two must be retired before either is worked.

---

## Added 2026-09-02

### [cut_representative_options.md](deferred/cut_representative_options.md) — non-tight cohomology generators: tighten, go cut-free, or multi-sheet

**PLAN — drafted 2026-09-02 with Codex (code inventory) and Grok (formulation); not scheduled.** The corc_solder
generator is a valid unit cocycle whose band admits no side labelling, so the thin-cut factory eats it
(guard added the same day). The plan explains the three ways out and works out the cut-free one in
BELFEM terms: h = −∇φ + Σ I_k ψ_k with ψ_k the Whitney interpolant of the cocycle — curl-free by δz = 0,
circulation = pairing regardless of tightness — as an extra dof in crossed air elements, an augmented φ
kernel, and a third hanging source (weight z_e) on interface edges. Option 1 (min-cost tension, TU
incidence) recommended alongside; option 3 rejected (reduces to option 2 after condensation).
**Reviewed 2026-09-03 by Gregory Giard:** option 2 also dissolves the non-unit-coefficient restriction
(no `clean_spfa` refine-your-mesh abort) and the tet-only restriction of the cut logic (the cohomology
engine is shape-agnostic, `CutProcessor` is not). His ruling, with Frederic: the thick cut is the
eventual target, not a current priority — the limits are known and guarded.

## Added 2026-09-01

### [thinshell_conductor_normal_field.md](closed/thinshell_conductor_normal_field.md) — thin-shell h_n recovery next to an h-conductor volume (corc_solder)

**IN PROGRESS — plan audited and approved 2026-09-01; R3+R4+R5 implemented the same day (Fable), syntax-checked, build + `make check-fast` + corc_solder rerun owed; jury round 1 on the diff: one P1 (linear-order guard) pending adjudication.** `compute_hn` takes −∇φ from the solder
master although the solder has no φ dofs → ~1500 A/m fictitious normal field on the tape sheets at 0.15 A, in the exodus
and in jc(B). Per-side dispatch in `compute_hn` and its postprocessor twin; plus the frozen transport current
(I₁ = I_ramp(0.02 s) at t = 0.1 s) and the unpinned bore φ.

### [penalty_opt_in_plan.md](closed/penalty_opt_in_plan.md) — honour `eta : 0` (no ghost facets), make the gauge — and, as a second step, the ghost — opt-in

**PLAN — drafted 2026-09-01, jury pending on Christian's O1 (flip the ghost default now or after a
deck sweep).** The old `mUseNitsche` never left: it is `ThinShellFactory::mCreateGhostFacets`,
hardcoded `true` with the `false` path intact. Plan exposes it from the deck (`eta : 0` ⇔ no
duplicate dofs, no ghost sideset), reverts the 2026-08-27 gauge default to off, and sequences the
ghost-default flip behind a sweep of every thin-shell deck. Evidence: on `tape_quench` κ ≈
2.1e15 × Δt[ms] regardless of either penalty; `chi = 1e-4` is below the estimator's noise;
`eta = 4e-6` converges like `eta = 4`.

### [mumps_memory_budget_plan.md](closed/mumps_memory_budget_plan.md) — probe the machine, hand MUMPS ICNTL(23), stop guessing ICNTL(14)

**IN PROGRESS — code landed and jury-audited 2026-09-01; R12 gate owed (build + cold `tape_quench` under MUMPS).** The `-9` retry ladder
(`ICNTL(14)`: 30 → 60 → 120 → 240) is a compile-time constant that never sees the machine's
memory, and MUMPS's own guide says the percentage bounds nothing once dynamic allocation is in
play. Plan: a platform probe for available memory (`MemAvailable`, not `freeram`), a rank-uniform
per-process budget (node share × safety, global MIN), handed to MUMPS as `ICNTL(23)` on the first
`-9`; ladder kept as the fallback when the budget is unknown; INFO/INFOG copy widened to 80 so the
BLR lower bound `INFOG(36)` is reachable; the soft path prints the shortfall again. Trigger: the
intrinsic-material `tape_quench` deck could not cold-start under MUMPS at 240 %.

### [matvec_csc_omp_stack_overflow.md](closed/matvec_csc_omp_stack_overflow.md) — `matvec_csc` overflows the OpenMP worker stack above 262,144 rows

**RESOLVED 2026-09-01 by removing the parallelism, not repairing it — build gate owed.**
`splinalg.f90`'s `!$omp reduction(+:y)` put each worker thread's private copy of the whole result
vector on that thread's **stack**; a Darwin worker gets 2 MiB, so a multithreaded CSC product
crashed with SIGBUS above n = 262,144 = 2 MiB / 8 B. Found when `examples/tape_quench_usermat`
(295,315 dofs) died in its first BDF1 step. DR-155.

A Codex + Grok jury round declined to approve the planned heap-buffer fix — it would have produced
*silently wrong* matvecs (D1), and two of the plan's scope guards were false (D2: MKL builds do
reach the kernel). Christian then closed it a different way: both matvecs are **master-only** and
are noise beside assembly and factorization, and no threaded-vs-serial measurement had ever been
made, so BELFEM's own `!$omp` directives are now gated on a new **`BELFEM_OMP`** define
(`USE_BELFEM_OPENMP`, **OFF** by default) and the kernel no longer threads at all. Nine third-party
`#ifdef OMP` sites deliberately stayed put — moving them would have deleted the oversubscription
warning, serialised PARDISO and hidden the banner's thread budget.

The plan keeps D1-D14 and the struck R-steps as the record of the design that was rejected.
**`USE_BELFEM_OPENMP=ON` re-arms the defect**, and `OMP_STACKSIZE` still matters for third-party
workers — both documented in `doc/parallel_execution.md`. Owed: build + `make check` in both token
states, and the deck on a default environment.

---

## Added 2026-08-31

### [l21_wrapper_card_repair.md](closed/l21_wrapper_card_repair.md) — the wrapper-policy card cites a contract that moved, and cannot enforce itself — **DONE 2026-08-31**, 11 sweep findings routed to Christian

`L-21` is the hard-policy tripwire that third-party libraries are never called directly. Its rule
stands; its supporting apparatus does not. The card's central citation,
`src/core/assert.cpp:294-320`, points at a five-line stub whose own comment says the MPI lifecycle
contract "moved verbatim into `comm_abort`… on 2026-08-30" — the guards and the
`MPI_COMM_WORLD`-over-`gComm.world()` choice now live at `src/comm/cl_Communicator.hpp:234-257`.
That is `INC-543`'s exact shape (a genuine line number carrying a statement that is no longer
true) **inside the card written to prevent it**, and `doc/coding_philosophy.md:631` is stale the
same way — plus it repeats one sentence near-verbatim.

Two more that matter for a reader acting on the card: the Rule names `commtools.hpp` as *the* MPI
wrapper, but `comm_abort` — its own worked example — is declared in `cl_Communicator.hpp` **on
purpose**, because `commtools.hpp` pulls linalg and core cannot include it; and the proposed
enforcement grep excludes `.f90`, so it is blind to the 33 direct `MPI_` lines in
`src/sparse/mumpstools.f90` and `parpacktools.f90` while flagging `typedef int MPI_Comm`, struct
members and comments as violations. **O1 blocks the enforcement step and needs Christian:** those
Fortran shims talk to MUMPS and PARPACK across a language boundary the C++ wrapper layer does not
cross, and whether that is an intended exception is written down nowhere — which is the real
defect, independent of the answer.

Documentation only; no `src/` change is proposed, and the plan explicitly forbids "fixing" the
Fortran. Found by the 2026-08-31 round-2 jury (Codex `gpt-5.6-sol`/xhigh + Grok `grok-4.6`/xhigh),
every citation re-verified firsthand.

**Closed 2026-08-31** (`devlog/dl20260831_l21_card_repair.md`). O1 ruled by Christian in the
direction that costs nothing: the wrapper requirement is specific to C++, and the Fortran drivers
for MUMPS, PARDISO and ARPACK **are** the wrappers the C++ tree uses — now a named, printed
exclusion in `scripts/check_wrapper_policy.py`, which replaces the card's grep with a call-level
sweep and was falsifier-gated before it was believed. INC-565/566 filed; the plan's assumption that
those rows already existed was wrong. **No residue.** The sweep's first run was red on 11 sites, none of them MPI, and Christian ruled
all 11 by design the same day: PETSc's scope reaches into the communicator (`PetscInitialize` is
paired with `MPI_Init` and cannot live behind a `Solver*`), and the materials HDF5 database probes
are sanctioned as written. Both are now owner entries that **print their reason** on every run,
scoped tightly enough that an `MPI_` call in `physics/materials/` or an `H5*` call in
`physics/gasmodels/` still goes red — re-proved by falsifier after the widening. Sweep green.


### [doc_currentness_fixes.md](closed/doc_currentness_fixes.md) — fix list from the 96-file documentation sweep — **R1–R7 COMPLETE 2026-08-31**, 10 boxes open (all routed to Christian or the out-of-scope R8)

The overnight Claude/Codex/Grok sweep of every markdown file in `doc/` and `src/*/doc/` found 17
defect groups. The headline is D1: **`hphirun` and `hphiTrun` are no longer built**
(`src/executables/CMakeLists.txt:40-47` comments both out; the build tree confirms) and eleven
documents still name them, six of them in runnable instructions — including `README.md:44` and
`doc/getting_started.md:45-46`, the two files a newcomer reads first.

Four taught snippets are broken in ways that do not announce themselves: `rho( norm(b), angle, T )`
against a `rho( T, B, beta )` signature (all-`real`, so it compiles and returns a wrong
resistivity), a `comm.send(data, target, tag)` MPI API that does not exist, and — the sharpest —
`heat_spline->eval(T)` documented as `cp` where the spline holds **enthalpy**, so every line of
that block is off by one derivative. Two defaults are documented backwards (`USE_DEBUG`, the
Coulomb-gauge `chi`), two thermophysical correlations are printed in the wrong functional form
(NASA-7 under a "NASA-9" heading; log₁₀ where the kernel uses ln), one cohomology formula has
`k+1` where Pellikka has `k-1` — while the same tree carries the correct form in a glossary four
files away — and `src/homology/doc/homology_usage_guide.md` needs a rewrite rather than line
edits, with its owner.

Two items are **ROUTED to Christian and must not be edited by an AI session**: whether the
static `∇·(μ∇φ)=0` air form in the Maxwell guide is stale or deliberate (formulation), and how to
reword the static-condensation attribution, given that `messe2023.txt:509-517` says the paper
implemented Lagrange multipliers and credits condensation to Alves et al.

Closes with `R1`–`R4`: `check_doc_claims.py` passes 37/37 and was structurally blind to every
class above.

### [handoff_20260831_open_defects.md](closed/handoff_20260831_open_defects.md) — the non-documentation findings of that campaign

What the documentation sweep turned up that was **not** a documentation problem, so its
documentation-and-comments-only scope could not close it. **F1** an unassigned owner sentinel indexing
`tCount` out of bounds in a release build (`graphtools.hpp:142-146`) — low today only because **F2**
finds `parmetis_nd`/`ptscotch_nd`/`build_pargraph_adjacency` are **dead code with no caller anywhere**,
which needs a wire/archive/delete ruling. **F3** records that the `mOwner` fix made failure *quieter*
(`2^64-1` → `2^31-1`) and corrects two campaign claims downward — no mesh entity ever saw the default,
because `mesh::Basis` sets owner 0. **F4** is my own error: a ticked D6 box hiding
`homology_usage_guide.md:1329-1350`, which still teaches `manifold_filter_3d()` (zero hits tree-wide)
as the live pipeline — **F4 CLOSED 2026-08-31**, and the residue was four sites, not two. **F5** measures
the `grep`/`.gitignore` false-negative exposure in shipped docs at exactly one claim, already fixed.
**Q1** asks whether the absent AC-loss postprocessing is deliberate. Reviewed, not verified.

**Rulings, 2026-08-31, Christian:** F2 → **wire it**, planned in
[`parmetis_ptscotch_wiring.md`](closed/parmetis_ptscotch_wiring.md) (which carries F1, since fixing it is the
precondition for wiring). Q1 → **a real gap**, planned in
[`ac_loss_postprocessing.md`](closed/ac_loss_postprocessing.md). Both scheduled for 2026-09-01.

### [ac_loss_postprocessing.md](closed/ac_loss_postprocessing.md) — energy, decomposition and honesty on top of the power global

**Scope reduced 2026-09-01.** Drafted on the finding that BELFEM computes no AC loss — true at
`bc6fb67d`, false two commits ahead: `7323d528` "saving heatlosses, currents and voltages to exodus"
lands a mesh global — since renamed by Christian from `heatloss` to **`dotQ`**, because the old name
did not say whether the number was a power or an energy — accumulating `rho*|j|^2*dV` over integration
points, elements and MPI ranks (`collect_dotQ()`, `cl_FEM_Controller.cpp:2329-2344`). So summation and reduction
are **done**; that commit is on Christian's Mac and on no branch here. What it does **not** do is the
part that makes the number publishable. The quantity is **power [W]**, reset before every assembly —
AC loss is **energy per cycle**, and nothing integrates over time. It is one scalar for the whole mesh,
so tape/stabilizer/former cannot be separated. It is fed through a **clamped** resistivity — but the clamp
turns out to be **unreachable with stock materials**, and two successive framings of that caveat in
this plan were wrong and are recorded as retracted. All three HTS law families put the power-law
channel *in parallel* with the normal channel (`powerlaws.hpp:193,214,242,269,290,309,331,348` and
the branch-stable `rho_riva`), so `rho <= rhon` always, however far over-critical an iterate goes;
the author says as much at `:2466-2468` and at `cl_Communicator.cpp:62-64` ("the defaults make the
clamp a no-op"). Verified in this checkout. It can still fire two ways, neither deck-reachable: a
plugin material with `rhon > 1e10`, or an executable narrowing the window. So the clamp flag drops to
low priority — cheap insurance, not a correctness fix. And it is undocumented: zero mentions in
`doc/`, zero in `examples/`, nothing telling a user it is watts — **since addressed upstream** by a
new §9.3 in the maxwell usage guide, unverified here. Remaining: time integration, per-block
decomposition, the clamp flag, a strip-solution gate, and the two-rank run the commit's own devlog
records as owed. Note `dotQ` is a deliberate name reuse — `cl_IWG_StaticHeatConduction.cpp:63,181`
already has a **nodal** `dotQ` on a different problem and a different mesh; separate namespaces,
examined not overlooked. Units independently verified here
(`mJ = C*q`, so A/m^2, giving watts).

### [parmetis_ptscotch_wiring.md](closed/parmetis_ptscotch_wiring.md) — wire the parallel nested dissection (F2)

**→ `closed/` 2026-09-04.** D1 fixed and verified, D2–D5 found and fixed on the way (D5: the matrix
graph's self-loops corrupted `METIS_NodeNDP` on the existing PETSc+metis path), the wiring landed as
`reordering scheme : parmetis | ptscotch`, and the new Tier-2 `sparsempi` suite is green at np 2/4
inside `make check` 19/19. Remnant is the benchmark (W-R5), filed as DR-156.

`parmetis_nd` and `ptscotch_nd` are written, declared, and called from nowhere; `build_pargraph_adjacency`
is called only by them. Christian's ruling is to wire them. **D1 must land first:** that adjacency builder
sizes its counter to the rank count and guards only the `comm_size()` marker, so a vertex still holding
`gNoOwner` — `INT_MAX`, 2147483647 (`typedefs.hpp:42,59`) — indexes it out of bounds; debug asserts,
**release performs an out-of-bounds write** (`graphtools.hpp:136,142-146`). It is harmless only while
dead, and wiring is precisely what ends that. Recorded with it: the `mOwner` fix made this site *quieter*
(`2^64-1` always faulted, `2^31-1` is a plausible offset) and must not be reverted, since the `std::min`
ownership sweeps need the sentinel to behave as +infinity. Build gates `BELFEM_PARMETIS`/`BELFEM_PTSCOTCH`
already exist. Note the test gap: nothing in `make check` enters this path, and the suite runs multi-rank
only for `tests/comm`.

### [sch04_history_purge.md](cancelled/sch04_history_purge.md) — purge the Fujikura `sch04` tables from git history

`share/material/fesc-sch04.hdf5` and `fysc-sch04.hdf5` were digitised from a copyrighted Fujikura
brochure and were deleted from the working tree on 2026-08-31 before the release
(`devlog/dl20260831_sch04_tables_removed.md`). The blobs remain reachable from `568a0fff` and
`b99ba890` on `main`, `devel` and `claude`, and all three branches on `origin` still carry them; the
GitHub `backup` mirror does not. R1–R10 cover the `git-filter-repo` rewrite, its four gates, the
force-push, and GitLab-side object expiry. **Nothing has been run** — the rewrite invalidates every
SHA from `568a0fff` forward and needs collaborator coordination. O1 decides whether it is needed at
all: a tarball or squashed-history release makes it moot.

### [dr147_plugin_path_resolution.md](closed/dr147_plugin_path_resolution.md) — DR-147: one search order for every plugin `.so` — **CLOSED 2026-08-31**

The source-function plugin was the only one of the three `.so` keys that never reached
`material::data_file()` — its deck string went straight to `dlopen`. The row had been parked as
"blocked on a layering choice" since 2026-08-29, because `src/numerics/sources` cannot include
`src/physics/materials`. The block dissolved on inspection: the resolver body needs only
`file_exists` and `gBelfemDataPath`, and every module already has `src/core` and `src/io` on its
include path, so the search order moved into `io/filetools` with **no `src/` CMake change** and
`material::data_file()` became a forwarder.

**Both audit rounds changed the work substantially, which is the point of running them.** At plan
stage both auditors refuted the proposed gate — a unit test on the new function cannot run red on a
tree where that function does not exist — and both found the scope guard ("no CMakeLists changes")
contradicted the test step. That refutation is what produced the suite's **first MODULE target and
first `dlopen`-based test**. Both also found a handle leak in `Material::read_defect` that the plan
had missed entirely, and Codex found a pre-existing false claim in two documents (an unset
`$BELFEM_DATA` does *not* leave only step 1 on an installed tree). Grok supplied the argument that
settled the ownership question: the Maxwell factory builds one `SourceFunction` per group member
from one `( file, label )`, so the leak was N mappings rather than one.

**✅ CLOSED 2026-08-31 — DR-147 STRUCK and archived, and unlike most strikes this one is VERIFIED
BY EXECUTION rather than reviewed.** `make check` green, and the 17 new cases were confirmed present
in the built binaries by name and run (`SourcePluginPath` 5/5, `SearchDataFile` 12/12) rather than
inferred from the suite result — the plugin test is `#ifdef`-wrapped and would otherwise have passed
vacuously, which was checked for explicitly. **Both halves of the red-then-green ran in one binary:**
`EnvironmentIsNotContaminated` — the hermeticity control both code auditors demanded — turns out to
exercise the exact pre-fix path (empty root → resolver returns the bare name → `dlopen` fails), so
the discrimination is executed, not asserted.

One deliberate deviation from the audited plan is recorded in §3.2: the subdirectory constant was
made bare against Grok's round-1 advice, to keep one source of truth. It was recorded BEFORE the
code round and the auditors were pointed at it, so it was tested rather than discovered; Codex
confirmed `data_path()` byte-identical in output and `data_file()` matching in strings and branch
order against HEAD.

Residual follow-ups, deliberately not closed with the row: the Codex prose sweep over the rewritten
`input_file_reference.md` sections (R15); a multi-config-generator hazard on the fixture's
`LIBRARY_OUTPUT_DIRECTORY`, harmless while the generator stays pinned to Makefiles; and
`material::data_path()`'s zero production callers, which predates this work.

## Sweep 2026-08-30 — register + todo currentness, active set 48 → 44

A currentness sweep — "has this week's recorded work made this row or file solved, obsolete, or
materially mis-stated?" — over the live debt register and every `todo/*.md`, run with Codex
(`gpt-5.6-terra`, high) and Grok (`grok-4.6`, high) against a pre-registered brief. Scope excluded
DR-106, DR-144 and the homology rows on Christian's instruction. **Four files moved, none deleted;
no row struck.**

**→ `closed/` (4)**

- [`ai_wrapper_tier_selection.md`](closed/ai_wrapper_tier_selection.md) — SOLVED. Thirteen steps
  landed, the R10 gate ran green, zero open checkboxes; §9.1 of the protocol carries the depth table.
- [`tapestack3d_jjc_noise_2125ms.md`](closed/tapestack3d_jjc_noise_2125ms.md) — SOLVED. Its own
  Status has read CLOSED since 2026-08-29, run-verified; the residue is DR-128's own file.
- [`run_gate_batches.md`](closed/run_gate_batches.md) — OBSOLETE. Derived from a classification of
  **44** live register rows; the table now holds **24**, and Batches A–D are all discharged or
  struck. Its §6 evidence-decay rule is *not* in `doc/lessons_learned.md` and is flagged in the
  file's closure banner for promotion — that is a `doc/` edit, so it was left for Christian.
- [`restart_timestep_consistency.md`](closed/restart_timestep_consistency.md) — OBSOLETE. Both
  halves are gone: (a) fixed, DR-112 struck 2026-08-30 after its hand-run restart gate; (b)
  retracted as a false positive 2026-08-29, and the archived DR-112 row says in terms that this
  file is superseded by that retraction. Grok's find; neither Claude nor Codex caught it.

**Register (17 rows in scope, no strike).** Seven rows carry the fix in the tree while their
description still states the defect in the present tense — DR-97, DR-105, DR-110, DR-119, DR-120,
DR-135, DR-139 — so each gained a dated **CURRENTNESS** clause naming what actually survives, and
DR-139's "sits inside the uncommitted circuit-cluster WIP" clause was struck (that WIP is committed
in `b99ba890`). Eight `file:line` citations were re-anchored after both auditors showed them now
naming *different functions*: DR-97 `:3434` (a BLR warning), DR-120 `:3368`/`:3542`, DR-135 `:288`,
DR-138 `:411`, DR-147 `cl_MaterialFactory.cpp:379`, DR-111's whole path (`ThinShellFactory` moved to
`src/fem/kernel/`), DR-110's caller list (no `add_source` caller survives in `cl_MaxwellFactory.cpp`
at all), and DR-39's four anchors. **DR-138 is the one strike candidate and it is Christian's, not
the sweep's:** its unit gate RAN green 2026-08-29 (15/15, 1e-3), but the row still names a coupled
reject-retry A/B, so it is the DR-42/DR-49 single-run exception — a ruling, not a sweep result.

**Not moved, flagged for a ruling.** `make_doc_repair.md` split the auditors: its R-steps are done
and its one open box is **O5**, a policy question ("should the gate be able to fail a build?")
routed to Christian — Grok called that SOLVED, Codex STALE-BUT-LIVE. It stays active until the
ruling. `current_sign_2d_fix.md` says DONE but still owes O5 (DR-133) and the O7 voltage-driven
regression run. `shared_library_and_install_plan.md` and `gauge_newton_tangent_plan.md` both carry
checkboxes their own prose says landed.

### [dr128_spap_low_field_axis.md](closed/dr128_spap_low_field_axis.md) — extend the sp-ap B axis below 10 mT

Registered 2026-08-30: the file was created 2026-08-29 and never indexed here. The sp-ap table's
10 mT B-axis floor clamps jc/n to B-independent values across the whole operating range of mT-scale
decks and puts a ∂/∂B tangent discontinuity exactly at the edge. Table data only, no C++ change.
Status: PLAN v2 — §3 (axis extension) is superseded by §7 on Christian's ruling, follow the Kohler
pattern the metals already use and build a transition function rather than clamping; §3 survives as
the fallback, and round 1's jury verdict must be read against §3 rather than §7.

---

## Added 2026-08-30

### [eigen_shift_invert_solver_reuse.md](cancelled/eigen_shift_invert_solver_reuse.md) — cut the eigen diagnostic's second resident MUMPS factorization

**OPEN — round-1 audited 2026-08-30 (Codex terra/xhigh + Grok grok-4.6/xhigh); BOTH said revise
before code, and the recommendation is now REVERSED from the draft.** Christian's proposal was to
borrow the production MUMPS solver when the field already uses MUMPS, instead of the diagnostic
building a second instance for the same matrix. **The draft's central premise was wrong: the
matrices are not the same under Newton.** `SolverData::jacobian()` returns **`mSystemMatrix`**
(`SolverData.hpp:679`) while a Newton body solves **`mJacobianMatrix`** (`SolverData.cpp:2428`) —
separate allocations (`:396`, `:405`) with `dJdx` added only to the Jacobian (`:1683`); found by both
auditors, verified independently. My error was trusting an accessor named `jacobian()` to return the
Jacobian and stopping one hop short. Three more draft claims fell: the "compute_conditioning()
reassembles" reasoning is **false on the live path** (`mMatrixFlag` is set only in
`compute_lambda_max` and cleared by `reset()` on every assembly), the DR-106 "created and freed per
call" quote is **stale** (the instance is retained for the `EigenValues` lifetime), and the
motivation was **misattributed** — on the STRUMPACK/PETSc run that prompted this, the borrow never
fires and the saving is exactly zero. **Grok proposed the better move and it is now Option A:
`free()` the dedicated diagnostic instance per call.** It keeps the isolation the current design
states in so many words, removes the *persistent* double residency (both instances hold factors
until process teardown — the real cost, not a transient peak), and needs none of the hazard
machinery. The decisive fact: the comment defending the retained instance justifies it as *"a
per-call free/re-init burns a fresh slot in the mumpstools registry each time ( observed as the
tapestack3d step-8 abort )"* — and **DR-143 fixed exactly that on 2026-08-29**, gate-proven at 43
conditioning calls with zero pool-exhaustion aborts where the old allocator died on the 8th. The
rationale is one day stale; the residual cost of Option A is a real JOB -1 + JOB 6 per timestep.
Borrowing survives as Option B behind a measurement, with **H4** (`mMatrix`/JOB-5-vs-6 selection)
and **H5** (Fortran workspace unreleased by a C++ `ICNTL(14)` restore) added — both missed by the
draft — H2 shrunk (the controller already sets production soft-fail and re-clears the latch), H1
corrected to restore the *saved pre-diagnostic* relaxation rather than a hard 30 (writing 30 would
undo production's own DR-106 ladder), and the borrow decision required to be **broadcast** because
`SolverParameters::synchronize()` does not synchronize solver type. **G2 is not runnable as drafted**
— production had already laddered to 240 before the diagnostic ran, so a borrowed diagnostic starts
at 240 and can never prove restoration; it needs fault injection. **New O4 outlives this plan:**
under Newton the eigen ratio and `MUMPS ADD COND1` are computed from *different matrices*, which
answers the open O3 of the `mumps error analysis` plan and deserves its own row — the conditioning
diagnostic may be reporting a matrix the solver never factors.

### [mumps_error_analysis_opt_in.md](closed/mumps_error_analysis_opt_in.md) — split the MUMPS ADD diagnostic out of `compute conditioning`, make it opt-in

**PLAN drafted and round-1 audited 2026-08-30 — BOTH auditors blocked; amended the same day; key
name pinned.** The audit (Codex terra/high + Grok grok-4.6/high, parallel, identical brief) returned
**the same P0 from both, found independently: D1** — `arm_conditioning_thermal` does *two* jobs
behind one early return, and the draft's R3 would have re-gated the whole function, silently turning
off `set_symmetric(true)` (`:3282`, the only write of `mSymmetric`, default `false`). That sends the
thermal field through the nonsymmetric ARPACK driver and changes its footer row from `κ₂ of Thermal
System` to `|λ|max/|λ|min` — a silent regression in the very diagnostic the split exists to leave
alone, in a change whose default state is what every existing deck gets. Grok proved it reachable
rather than theoretical: `examples/3D_tapestack/input.conf` ships magnetic MUMPS (`:19`) + thermal
PETSc (`:48`), both `compute conditioning : true`, and no gate in the draft would have seen it.
Four more: **D2** hoisting the `print_footer` tail is *not* sufficient — the ADD rows are nested
three deep, so the draft's R4 left `(eigen=off,mumps=on)` printing nothing and `(eigen=on,mumps=off)`
still printing ADD, i.e. both behaviours inverted; **D3** R5 could not have worked, since
`set_params` runs with `mKernel2 == nullptr` (`cl_MaxwellFactory.cpp:983`) — which also means the
*existing* thermal no-MUMPS warning is already dead on the coupled path; **D4** G1 was vacuous
(`check_doc_claims.py` reads neither input-contract artifact); **D5** R8 asserted the causal claim O3
forbids. Amended: R3 splits the two jobs, R4 spells out four predicate edits, R5→R5a/R5b (thermal
warning moves to `set_thermal_kernel`), R8→R8a/R8b (causal half blocked on O3), new R9 migration
notice / R10 example fix / R11 stale-comment sweep, five gap rows added, gates rebuilt as G1–G5 with
**G3b** (the shipped tapestack shape) as the one that catches D1. **O2 resolved → high, unreachable.**
**Christian ruled on D1 the same day: resolve by removal, not acceptance.** Two new steps —
**R3c** drops the κ₂ label branch and always prints `|λ|max/|λ|min` (correct for symmetric and
nonsymmetric alike, never over-claims, and does not ask the reader to know κ₂; it is also the label
**INC-213** already pointed at, `lessons_learned_evidence.md:305`, having recorded κ₂ as itself part
of a misreading), and **R3b** moves `set_symmetric( true )` out of `arm_conditioning_thermal`
entirely into `set_thermal_kernel`, beside the late-attach block that already re-runs the gates
`set_params` could not. With the footer no longer calling `is_symmetric()` and the arming functions
no longer carrying an eigen side effect, R3 collapses back to the straight predicate swap the draft
wanted — D1 is designed out rather than commented around. New **O5** flags the one thing R3b assumes
and has not proven: that the thermal `EigenValues` object is not rebuilt after attach, or the
one-shot write is lost where the per-timestep write survived. **Round 2 (both, `xhigh`) returned NO blocking defect** — Codex "block only for
three precise plan amendments", Grok "implementable for R1-R4... I would start R1". **O5 resolved
safe** by two independent traces: `mEigenValues` is assigned only in the `DofManager` ctor
(`cl_FEM_DofManager.cpp:60`) and nothing rebuilds it — `EigenValues::reset()` clears only
`mMatrixFlag`, `DofManager::reset()` never touches it, `link_matrix()` rebinds the matrix not the
object. The code already assumed this: `cl_FEM_DofMgr_EigenValues.cpp:770-772` says `mSymmetric` "is
set once at setup". Both auditors then found the same three gaps independently: the two-argument
`Controller( Kernel *, Kernel * = nullptr )` constructor (`hpp:466`) attaches a thermal kernel
*without* `set_thermal_kernel`, so R3b/R5b need a shared helper; R9's `!key_exists` provenance is
discarded by R2's bools and its thermal half hits D3 again (no thermal solver at parse time); and
**G3b could not prove what it claimed** — pinning the thermal value cannot detect a lost
`set_symmetric`, since both ARPACK drivers converge to the same ratio, which contradicted this
plan's own reasoning. New **R12** adds the real tripwire, a `BELFEM_ASSERT` on `is_symmetric()`, and
R3c now explicitly keeps that accessor for it. Grok additionally caught that **G3b fought R10** (R10
adds the key to the very file G3b needs it absent from → G3b runs on a pre-R10 fixture), that R5b's
deferred-message option has no call site at all, and that G1's "anywhere" wording would fail on this
plan file and on INC-213. Now R1-R12 / G1-G5, **pending Christian's approval to begin R1**. Original
entry follows. One flag pair
(`cl_FEM_Controller.hpp:134-135`) gates two unrelated diagnostics: the ARPACK eigen kappa
(`cl_FEM_DofMgr_EigenValues.cpp:929`, the expensive one) and the MUMPS Arioli-Demmel-Duff
`ICNTL(11)` numbers (`mumpstools.f90:404`, cheap — omega statistics plus a Hager estimator against
the existing factorization, no refactorization). Neither can be requested alone, so a user who
wants kappa is forced to take the `MUMPS ADD COND1 / COND2` rows, and a user who wants those rows
is forced to pay for an ARPACK run. Christian's call: the ADD numbers are noise for his purposes
while kappa is actionable, so the MUMPS half becomes a **new independent boolean, default false**.
Eight steps R1-R8 plus four gates G1-G4; the only structural edit is hoisting the duplicated tail
out of the `print_footer` `if/else` (`:3733` vs `:3834-3846`, verbatim copies). §3.1 flags the one
way to get this wrong: `arm_conditioning_*` and `capture_conditioning_*` must be re-gated **in the
same edit**, because the disarm lives only in the capture path (`:3315`, `:3333`) — an asymmetric
change leaves `ICNTL(11)` armed for the rest of the run. Three open questions logged not decided,
including **O3**, raised by an external suggestion and unresolved: with the Coulomb-gauge penalty on
by default at chi=1e-4, does the operator ARPACK iterates carry the same Dirichlet/penalty
contributions as the one MUMPS factors? It blocks nothing here but must be settled before R8
asserts why the two numbers differ. **Name pinned `mumps error analysis`** by a blind jury round
(`tmp/ai_exchange/mumps_cond_key_name.md`, Codex terra/high + Grok grok-4.6/high) that **split 1-1**:
both independently ranked the three offered names identically and both found `compute ` is a
**singleton, not a house prefix**, but Grok proposed a fourth name — bare `error analysis`, vendor
dropped — on the verified convention that the deck puts capability in the key and vendor in the
note (`matching` is STRUMPACK MC64 and is not `strumpack matching`; `metis nodendp` is the lone
exception). Rejected on a ground neither auditor weighed: the flag is a **no-op without MUMPS**
(`:3259-3262`, `:3272-3274`), so a bare `error analysis` would sit inert in any STRUMPACK deck with
nothing in its name to say why. Two auditor citations were corrected on re-verification — Codex put
`ICNTL(11)` at `mumpstools.f90:400` (it is `:404`) and called portability a red herring by citing
SuperLU `dgscon` while the same section says `dgsrfs` returns the ADD quantities proper
(`conditioning_diagnostic_backends.md:136-140`). Kept from the dissent: the R5 no-MUMPS parse
warning and the R7 constraint that the reference row name the footer strings verbatim. Also fixes
two documents that currently describe the two quantities as two routes to one kappa
(`doc/input_file_reference.md:219`, `src/fem/doc/timestepping_strategy.md:154-177`, whose footer
diagram is additionally stale against `:3736-3818`).

### [dr111_quad4ts_winding.md](closed/dr111_quad4ts_winding.md) — DR-111: the 2-D thin-shell layer elements are wound clockwise

**✅ CLOSED 2026-08-30 — DR-111 STRUCK ( Christian's ruling, DR-42/DR-49 pattern ); struck is not verified, one assertions-ON gate owed.** Fix landed in the CONSUMER, not the factory: `Calculator::dV_quad4ts` and `Pipette::measure_quad4ts`, `|·|` wrappers at one dispatch site each, shared bodies byte-identical. Two plan-audit rounds + a code round, unanimous functional CONFIRMED. `ThinShellFactory::process_nodes_line2` builds the LINE2 extrusion normal as `n = ( ty, -tx )` (`cl_ThinShellFactory.cpp:923-924`) — the tangent rotated **clockwise** — while `create_elements_on_blocks_line2` winds `( b0, b1, t1, t0 )` (`:1543-1546`) and `QUAD4TS` interpolates with the CCW-referenced `QUAD4` Lagrange functions. Every 2-D thin-shell layer element BELFEM has ever built therefore has `det J < 0`; numeric replication of the four functions reproduces the register's gdb value `-5e-12` for `+x`, `+y`, 37° **and the reversed facet**. The 3-D branch is correct by construction ( right-hand normal at `:1101`, `+2.0e-16` for CW/CCW/tilted alike ), which is why only the 2-D thermal path aborts. **This corrects DR-111's own root cause**, which says the winding is conditional on sideset orientation — it is not, `n` is derived from `t`. Nothing caught it because the magnetic path never computes a signed determinant: `EF_QUAD4TS` uses `mDetJ = 0.25*thickness*length` (`cl_EF_QUAD4TS.cpp:128`) and takes its tangent from the **facet**, not the element. Fix A = reverse the in-plane pair to `( b1, b0, t0, t1 )`, the only permutation that keeps nodes `{0,1}` on the bottom curve as `link_elements_with_edges` (`:1926-1958`) and `get_nodes_of_edge` (`cl_Element_QUAD4TS.hpp:104-134`) jointly require. **The round's substance:** A flips both `mS` entries ( the layer block element takes the node-order-based `compute_edge_directions()`, `cl_FEM_Element.cpp:181-183`, not the `physical_tag` thin-shell variant ) **and** swaps the hang sources ( `cl_MaxwellFactory.cpp:1776-1834` ). Codex called that a cross-wire and rejected the plan; Grok called it a compensating sign flip and endorsed it; **adjudicated in Grok's favour on a first-hand trace** — the positional label is a permutation applied to both sides, so for a 2-node edge the source *set* is unchanged and `q = φ_s0 − φ_s1` is exactly negated, cancelling the `mS` flip. Held at ~70 %, so **R4 ( magnetic-residual parity on `2D_Tapestack`, revert branch pre-registered ) is the tie-breaker, not a third review round.** Not assumed for LINE3/QUAD9TS — that edge has a third node ( O1 ). By-catch: layer thickness has **no positivity check** ( unit check only, `cl_MaxwellFactory.cpp:3083-3094`; `cl_FEM_Kernel.cpp:987` asserts only `abs(t) > 0` ), so a negative deck value silently inverts the stack; and `MeshChecker`'s swap table maps `QUAD4TS -> swap_quad4` ( `cl_MeshChecker.cpp:93-97` ), which would make the tangential pair the through-thickness pair and feed `EF_QUAD4TS` the tape **length** as its thickness — unreachable today, latent trap. The `Calculator::dV` tier demotion Christian ruled the same day is a **no-op in every build configuration** ( release has no check at all ) and is unordered relative to the winding fix.

### [dr129_kernel_controller_init.md](closed/dr129_kernel_controller_init.md) — DR-129: `Kernel::mController` is never initialized — **CLOSED 2026-08-30**

**✅ CLOSED 2026-08-30 — filed, fixed, audited twice and all three run gates green in a single day; the register row is STRUCK and archived, clearing the debt register's only P1.** Christian ruled O1(c) + all-three-members scope the same day; R1/R2/R3/R4/R5/R8 landed across 7 files and the code-audit round ( Codex terra/high + Grok grok-4.6/high ) returned **unanimous APPROVE, zero logic defects** ( two comment nits, applied ). Debug AND release gates BOTH RAN GREEN same day ( `GhostElementContract [ OK ]` in the rebuilt debug `test_fem`, then again in a `USE_DEBUG=OFF` / `-O2 -DNDEBUG` rebuild with `make check` 17/17 suites green — release-ness confirmed by an ASSERT-only string being absent from the linked binary, so the pass is not coming from a debug guard ); the coupled thermal smoke RAN too ( release build, 8 ranks, genuinely coupled h-ɸ/T; it carried the R8 construction site at startup then ran 483 converged coupled timesteps ( steps 203 to 686 ) with zero DR-129-class failures, closed early because that site executes once at startup ). **All three run gates green — ready to strike on Christian's ruling.** The plan-audit round preceding approval: Round 1 added D3 — the implicit copy of an indeterminate `Kernel` at `cl_ThermalFactory.cpp:138`, executing on every coupled run ( new step R8: construct in place ) — escalated D2 ( `mMyCommIndex` is copied there too, so it is live UB, not a dead member ), refuted the draft's 'silent wrong physics' severity story ( `have_thermal()` has zero consumers ) and its order-shuffle gate ( post-fix the member is deterministically null ), and produced a unanimous three-voice O1(c) recommendation. DR-129 is the only **P1** row
left on the debt register. `Controller * mController` at `cl_FEM_Kernel.hpp:70` carries no default
member initializer, is absent from the constructor's initializer list (`cl_FEM_Kernel.cpp:46-52`),
and is written only by `set_controller` (`:1127`) — so `Kernel::controller()` on a kernel that was
never handed a Controller returns garbage, and the `BELFEM_ASSERT` at `:1137` can neither see a
garbage-non-null value nor survive a release build. Proven live 2026-08-28: the DR-46 ghost fixture's
`link_to_group` passed solo and SEGV'd deterministically after 21 sibling tests, and the fixture now
carries a comment explaining why the call is skipped (`cl_TS_TestStack.hpp:892-899`) instead of a fix.
**Two findings widen the row rather than branching it** (extend-don't-branch, 2026-08-29): **D1** the
row names one consumer but there are **seven** reads of `Kernel::controller()` — `link_to_group`
(`cl_IWG_Maxwell.cpp:266`), the `MaxwellData` constructor's *member-initializer list*
(`cl_FEM_Calculator.cpp:78`), `cl_FEM_Calculator.cpp:1985`, and four in `cl_MaxwellPostprocessor.cpp`;
**D2** `Kernel::mMyCommIndex` (`cl_FEM_Kernel.hpp:76`) is uninitialized *and* dead — one grep hit in
`src/`, `tests/` and `nonfree/` combined, the declaration itself. Two claims from the register row
were checked and corrected: **`has_controller()` does not exist anywhere in `src/`**, so the fix
shape needs it added (R2), and the consumer lives in `src/fem/maxwell/`, not `src/fem/iwg/`. The
row's severity assessment *survives* re-verification: every `link_to_group` call comes from a
`DofManager` assembly routine, and `belfem.cpp` builds the kernel (`:132`) before the controller
(`:144`) whose constructor calls `set_controller` (`cl_FEM_Controller.cpp:87`), so no shipped deck
reaches a controller-less read. Seven steps R1-R7; O1 (is a controller-less `Kernel` legal, or does
`controller()` become a `BELFEM_ERROR`?) is Christian's ruling and is explicitly not decided in the
plan. Gate is the DR-46 fixture with `link_to_group` re-enabled on a debug build — must fail with a
**named** error or pass, never SEGV — plus the battery order-shuffle, since the original defect was
order-dependent and a plain suite pass cannot observe it.

### [ai_wrapper_tier_selection.md](closed/ai_wrapper_tier_selection.md) — explicit model + effort selection for `ask_codex.sh` / `ask_grok.sh` — **CLOSED 2026-08-30**

Every Codex audit this project has run inherited `~/.codex/config.toml` (terra/medium) because
`ask_codex.sh` passes neither `-m` nor `-c model_reasoning_effort`; every Grok audit inherited
`~/.grok/config.toml` (grok-4.6/xhigh) because `GROK_EFFORT` is never set. Neither wrapper records
which depth produced a finding, so a thin round-1 result is indistinguishable from a thin deep one.
The plan adds `CODEX_MODEL`/`CODEX_EFFORT` and `GROK_MODEL`, an allowlist validator, and — the
load-bearing part — a `(model=…, effort=…)` stamp on the exchange header line, plus a tier table in
the protocol keyed on subject scope and round number rather than on Claude's self-assessed
confidence. **Jury-audited 2026-08-30** (Codex "revise before implementation", Grok "approve with
corrections"); both endorsed the architecture. Grok's decisive catch was tested rather than agreed:
an unexported assignment in `cross_review.sh` never reaches `ask_*.sh`, so the plan's own pinning of
the unattended `--quick` hook would have silently kept running at the wrapper default. Both auditors
also overruled the draft twice — the second rung now raises effort rather than switching model, and
Grok's depth is recorded at its current xhigh rather than reduced in the same change that starts
recording it. ~~Status: PLAN — audited and revised, blocked on Christian's decisions O1-O6; no script
modified.~~ **Superseded: COMPLETE 2026-08-30, all thirteen steps landed and the R10 gate ran green.
Moved to `closed/`.**

---

## Added 2026-08-30

### [dr106_dr144_joint_session.md](closed/dr106_dr144_joint_session.md) — close DR-106's `-9` retry together with DR-144's deferred R12/R13

**COMPLETE — both rows STRUCK 2026-08-30** after the same-day rebuild + restart recovered two real `-9` events through the new ladder and `check-fast` ran green ( see the file's Status line and `devlog/dl20260830_mumps_minus9_retry.md` ).

DR-106's workspace retry and DR-144's R13 edit the same two functions in `cl_SolverMUMPS.cpp`, so
running them as separate rounds means two plan+audit cycles over the same forty lines. R12/R13 were
deferred "blocked on ownership, not on difficulty" because the block was another session's
uncommitted work; **that session has closed and the tree is clean, so the block is gone.** The
2026-08-30 jury round supplied the constraint DR-106's text was missing for six days: the shim
exports rank-local `INFO`, not `INFOG`, so a retry keyed on `mInfo( 0 ) == -9` fires on one rank
and strands the rest inside a collective call. Eight steps, four open questions, and a warning that
`make check` cannot gate the retry — only a deck rerun that reproduces `-9` can.

---

## Added 2026-08-29

### [dr144_mumps_lifecycle_plan.md](closed/dr144_mumps_lifecycle_plan.md) — close DR-144's seven MUMPS instance-lifecycle gaps

The C++ wrapper tracks a MUMPS instance's lifecycle in three members set at three different moments,
and only the Fortran occupancy table is authoritative. Seven filed gaps, all re-verified open against
the working tree after DR-143's first-free scan landed. The round's main finding is that gaps 3 and 4
collapse into ONE line: `MUMPS::mInitialized` is written at three solve sites and **read at exactly
one**, `cl_SolverMUMPS.cpp:284`, so it neither carries the frozen-factorization scope nor the
soft-fail retry the register believed it did. Highest value per line is gap 6, a missing
`EigenOutcome::SolverFailed` arm that reports every solver failure as a basis-budget refusal.
Scope approved by Christian 2026-08-29 ("DR-144 in full, one round"). **Jury-audited the same day:
Codex "revise before implementation", Grok "approve with the required corrections" — both confirmed
the gap inventory and both rejected the step shapes.** The seven gaps became **nine** (a double-init
double-leak, and an `MPI_BARRIER` status argument of the wrong integer kind under `BELFEM_INT64`),
and three fix shapes the register itself proposed were wrong, including one that would have enabled
a cross-wrapper destroy. **Christian approved and it is now IMPLEMENTED (9 of 11 steps, 2026-08-29)**
— all nine gaps landed with every audit correction, syntax-checked clean in C++ and Fortran (the
Fortran check confirmed real by a negative control), code-audit jury dispatched. O1 resolved to
*delete the dead setter* rather than guard it, because the shim carries no error idiom at all and
the API has zero callers. **`make check` RAN GREEN 2026-08-29**, discharging R8 and R9 and DR-140's
shared rerun — and confirmed to have actually executed the new test rather than inferred from the
suite result. DR-140 is struck; **DR-144 is not**, because G10 and G11 are live residue in other
sessions' uncommitted code. Five steps stay reviewed-not-verified regardless: no in-tree trigger
reaches them.

### [user_api_header_install.md](deferred/user_api_header_install.md) — replace the 628-header install with a flat 14-header user API

**→ `deferred/` 2026-09-04 (Christian's ruling: split).** The layout-independent items landed the
same day — `LICENSE` at the prefix root (R7/D4), `BELFEM_`-prefixed guards on the umbrella and
`cl_SourceFunction.hpp` (D1 half, D9), the prohibition note and `gTbulk` sentence on the umbrella
(R2 half), a correct compile line in the shipped `cl_Material_UserDefined.hpp` (D10); D8 was already
gone. The spine (flat `install( FILES )`, umbrella move, isolated staged gate, template rewrite,
end-to-end gate) is parked: the 2026-08-31 repairs made an installed prefix build plugins, and R1
would now rewrite the include blocks those repairs just fixed. O1 and O5 remain Christian's call;
the revision-2 re-audit is owed before any spine code. Original entry follows.


`make install` ships every `.hpp` under `src/` to serve a documented plugin API of 14. Christian's
2026-08-29 API change made all three entry headers backend-free, which removes both objections the
jury round raised against a curated set. Measured: all four shipped plugin examples compile against
one flat directory with one `-I` and no backend define. **Revision 2** — revision 1 failed its
Codex+Grok audit on three P0s (R1 landing before the template probes; R4 as a deletion list; R3 not
removing the inherited `-I`, which made the pinning gate vacuous), all in the step spec, none in
the architecture. Re-audit owed before code. Open: D1 (umbrella in `src/fem/kernel`), D2 (Undulator
CMake), D3 (templates not installed), D4 (`LICENSE`), D6 (library template demands a forbidden
backend), D8–D10; O5 (R8 depth) needs Christian.

---

## Added 2026-08-28

### [current_sign_2d_fix.md](closed/current_sign_2d_fix.md) — 2-D current-BC sign inversion (declared +I runs as -I)

Measured on gantry: all 464 thin-shell tapes carry exactly -340 A for declared +340 A; 3-D correct.
~~Plan+audit in flight; G-B1 bulk-2-D discriminator run pending, then a one-line fix in the suggested-homology sign.~~ **Superseded 2026-08-30: Fix B landed and gate-verified, input contract updated. Still open in the file: O5 ( DR-133 ) and O7 ( the G-V1 voltage-driven regression run ) — so it stays active.**

### [circuit_demo_harvest_plan.md](closed/circuit_demo_harvest_plan.md) — harvest the `electricalCircuit` demo circuits into transient tests

**Status: ✅ COMPLETE 2026-09-04.** R2/R4/R5/R6 landed in the existing `tests/circuit/test_ElectricalCircuit.cpp` and RAN green in Christian's `make check` (each case `[ OK ]` by name); the diode Newton converged at every step with ω ≡ 1 (worst case 11 iterations), so R1′ was not needed; the netlist twin held at 1e-12 over 1600 steps. The executable was then deleted (R7), fifteen sites retargeted (R8), DR-39 struck with (b) re-homed to the ngspice plan's Phase 6 (R9), devlog `dl20260904_circuit_demo_harvest.md` (R10). One gate still owed: `make check` on the post-deletion tree, proving the executable list links. Previous status: ADOPTED and RE-SCOPED 2026-08-30 (Christian's ruling on DR-39(a)).** The four demo circuits hardcoded in
`src/executables/electricalCircuit.cpp` are the tree's only exercise of the circuit module's
*transient* path, and three of the four are commented out against a `SourceType::SINE` API that no
longer exists (`:68,75,81`); the R-only one would be singular inside the shared
`ElectricalCircuit(4)`, and the diode bridge — a topologically correct bridge — has no load between
its DC terminals, so it rectifies into nothing. Plan: rebuild the circuits as fixed-Δt gtest cases
with analytic gates (exact for the resistive case; a two-regime switch gate 0.9980 A → 1.5056 A that
fails loudly if the switch mistimes; full-wave rectification with a 9.6426 V peak from the diode
fixed point), add one netlist-equivalence test tying `NgspiceCircuitFactory` to the switch gate, then
**delete the executable** — which leaves `Controller::solve_circuit()` as the single Newton loop,
discharging DR-39(a) without writing a `CircuitSolver` class. Christian approved the deletion, the
bridge completion, and the direct-API + one-equivalence-test layering on 2026-08-28, and **ruled the
`CircuitSolver` extraction refuted for good on 2026-08-30**.

**Re-scoped 2026-08-30 against the tree — the plan had been overtaken three times.** ~~a new
`tests/circuit/test_TransientCircuits.cpp`~~ and ~~taking all of DR-118 with it~~ are struck:
**R1 and R3 are already satisfied in-tree** (the fixed-Δt driver `solve_attempt`/`take_step` and the
series-RLC circuit both landed in `tests/circuit/test_ElectricalCircuit.cpp` after this plan was
written), so the harvest goes into that existing TU rather than a new one — a second TU would
duplicate the driver, which is a poor look on a de-duplication row; and **DR-118 was fixed in place
on 2026-08-29 and struck on its own gate**, so deletion no longer discharges it. **O3 resolved by
events** (DR-123 fixed and struck 2026-08-29; the rejected-step switch case is now a green test).
The deletion case is unweakened and rests on what always carried it: the binary takes no input file
and can solve only Frederic's circuit, yet ships installed, and three of its four circuits name a
`SourceType` enum that exists **nowhere** in the tree. **The hole the ordering protects:**
`create_diode` appears in the whole suite only in the structural adjacency test — the diode Newton
has no transient gate at all, and R5 is its first. Open: R1′/R2/R4/R5/R6, then R7-R10; O2 open.
Note R7/R8 change `BELFEM_INSTALL_EXECUTABLES` and four doc sites on freeze day.

### [make_doc_repair.md](closed/make_doc_repair.md) — repair `make doc`: Doxygen front-end breaks and the 1.18-vs-1.9.1 skew

**Status: ✅ COMPLETE 2026-08-28 (1.9.1 half).** `make doc` went from a hard abort to exit 0. One root
cause, four faces: every Doxygen input artifact in the tree was authored against 1.18.0 while the only
installed Doxygen is 1.9.1. **B1** `doxygen -w html` with no config argument still reads `./Doxyfile`
and validates the `HTML_HEADER` that same command creates — fatal on any fresh tree; **B2** (the
reported abort) the logo-patch needle carried 1.18's `$logosize`; **B3** 46 ignored tags, silently
dropping the site to the MathJax 2 CDN; **B4** the 1.18 layout dump produced the run's only `error:`.
The docs were never broken. Fix: the target now derives every Doxygen input from the *installed*
Doxygen — Doxyfile copied then `-u`'d (never in place), layout regenerated with `-l` and BELFEM's
delta re-applied read from `doc/DoxygenLayout.xml`, header generated from a config copy with
`HTML_HEADER` blanked, logo needle idempotent by construction. **Measured: console tag warnings 46 →
0; log 200/1 → 90/0, all 90 in one file and zero elsewhere in the tree.** Two audit rounds corrected
three of Claude's own errors (fabricated Bezier identifiers, a false `configure_file` premise, a
misleading baseline comment). Residue for Christian: **O5** the gate only `message(WARNING)`s so it
cannot fail a build; **O6** `make doc` writes tracked files; **O7** unverifiable on 1.18 here.
Devlog: `devlog/dl20260828_make_doc_repair.md`.

### [weekly_report_20260828.md](closed/weekly_report_20260828.md) — week in review, 22 to 28 August 2026

**Status: REPORT — not a plan, nothing owed from this file.** A one-page distillation of the
week's 68 devlog entries and 45 commits, written for the project owner: release engineering
(one `libbelfem`, `make install`, the unified `belfem` executable), the quench-front solver arc
(increment-form Picard, the reset seed hole, `chi` default ON, and the negative result that
Picard is structurally insufficient at the front while deck-level BDF1 outran BDF5+Newton), the
three silent production corruptions (`face_key_3d` at order 2, the `uint8_t` source counter,
the orphaned DR-100 MPI send), the 9-of-10 G-operator battery, materials and gas physics,
ngspice parser v1, and the method findings worth carrying. Closes with the gate ledger: much of
the week is reviewed, not verified. Register movement: 40 → 41 live, 48 → 88 retired, 41 new IDs
filed. Codex prose sweep applied; Grok fact/completeness audit refuted four first-draft claims
(live-row count, DR-87 drift direction, driver-collapse count, §1 summary) and caught the DR-100
omission, all re-verified against the tree before applying.

### [tapestack3d_jjc_noise_2125ms.md](closed/tapestack3d_jjc_noise_2125ms.md) — j/jc noise onset at t = 2.125 s: verdict + mitigation — **CLOSED 2026-08-30**

**Status: CLOSED 2026-08-29 — symptom resolved and run-verified as resolved (Addendum 3).** The
tapestack3d deck rerun on the fixed build (DR-126 gauge autopins + DR-127 certified exit) runs
137 steps to t = 3.3 s with the 2.125 s frame growing at ×1.047, the same factor as every other
frame in the window; largest I₁ frame ratio anywhere ×1.115; zero cuts or resets; conditioning
1.26e9. DR-127 was STRUCK the same day once its last branch — the `solve_thermal` reset return —
was executed under a forced-reset probe and retried cleanly at a halved sub-step. Residual work is
DR-128 alone (the 10 mT table floor, still clamping in these very frames). Verdict below as delivered
2026-08-28 (blind jury, 3 voices reconciled). The noise is a ratcheting current reorganization (max|B|
×3 in one 25 ms step, tape-1 |Jx| ×40, I₁ ×40), not display jitter. Material table and piecewise
law REFUTED as trigger (table clamped below its 10 mT B-origin; laws agree to ~1e-28 in the visited
regime; onset ran ρ-only Picard). Mechanism: ~9 unpinned φ gauge constants (8 buffer patches +
floating air level despite the bearing — Christian's buffer-unpinned hypothesis, confirmed by
probe) × resistively unpinned sub-critical stack × always-accept Picard semantics (pre-update
residual certifies x₁, unchecked x₂ is time-stepped). riva-vs-piecewise correlation confounded
(build/backend/algorithm/BDF/save-cadence). Ranked mitigation inside; formulation calls ⚑ Christian.

### [rho_lambda_b_first_reorder.md](cancelled/rho_lambda_b_first_reorder.md) — reorder the field-dependent `rho`/`lambda` API to B-first

**Status: PLAN — drafted 2026-08-28, Codex+Grok plan audit dispatched, pending approval; no source
modified.** Christian's decision 2026-08-27: uniform field-first, `rho( B, angle, T )` and
`lambda( B, beta, T )`, replacing today's `( T, B, beta )` — implementation deferred out of that
session. Inventory: 20 signatures, exactly four positional call sites (all `cl_FEM_Calculator.hpp`),
one function-pointer forward (`mFunctionRhoKohler`) that no compiler can catch, zero sites in
`nonfree/`. Method is a two-pass reorder whose intermediate state fails to compile at every
positional caller, so the change is compiler-enforced rather than trust-based. **Claude's
plan-stage recommendation is AGAINST** (Appendix B): the code is already uniform — the split is
between the code and two fem guides — and B-first would break the overload family's prefix
consistency at 17 one-argument call sites. Reverses the 2026-07-11 option (c) ruling, so
[rho_lambda_argument_convention.md](closed/rho_lambda_argument_convention.md) goes HISTORICAL if it lands.
`[F]`-class: it changes the user-material plugin API, so O4 sequences it against the design freeze.

## Added 2026-08-27

### [gauge_newton_tangent_plan.md](closed/gauge_newton_tangent_plan.md) — consistent Newton tangent for the Coulomb gauge penalty

~~**Status: DRAFT — plan audit (Codex+Grok jury) dispatched 2026-08-27 evening; no source touched.**~~ **Superseded 2026-08-30: IMPLEMENTED 2026-08-27 late evening and audit-corrected; four validation/open-question items ( D2, D3, R3 gates, O2 ) remain, so it stays active.**
The gauge penalty's Newton cross-term (`χ·(GᵀGq)⊗(dρ/dJ·Jᵀ C/|J|)`) was skipped while the penalty
was believed vacuous; at current sharing it is (n−1)× the kept term and the gauged
tapestack3d_coarse run walls on its signature at t = 4.17 s. Two HTS Newton kernels in
`mt_maxwell_h.cpp` get the block; side-connector kernels verified correct (metal, dρ/dJ = 0).
Control gate G1 (gauged Picard-only run) is Christian's this evening; blocks the pending
"gauging on by default" ruling (O3).

### [ngspice_parser_plan.md](closed/ngspice_parser_plan.md) — promoted from `deferred/`

**Status: OPEN — Phase-0 signed off AND audit round 1 complete (2026-08-27, Codex+Grok): three P0s found and folded into §12 as O5–O9; O5–O9 ALL signed off same day; **Phases 1–3 implemented, three two-vendor code-audit rounds run, and gate-passed in-tree the same day — 66 tests green** (`spice_number_to_si`, `NetlistParser`, `NgspiceCircuitFactory` with two solved analytic circuits; SPICE polarity conventions confirmed by execution). Remaining for v1: Phase 3b (hphirun wiring + showcase + Input Contract), gated on the DR-115 fix; by-catch DR-117/DR-118/DR-121 filed.** **Superseded 2026-08-30: v1 is COMPLETE ( shipped 2026-08-27, `examples/circuit` runs from `tapestack.cir` ). What remains is Phases 4-6, which is the surviving half of DR-39.**
The 2026-06-28 SPICE-netlist proposal is active again on Christian's request (registration entry
under "Circuit Module: ngspice Importer & Standalone Runner" below). Decisions: FEM linkage stays
in `input.conf` (6b — moot ngspice-loadability argument resolved: Option A comment directives keep
the `.cir` vendor-valid either way), Option A encoding for superconductor/timed switch, no
`.subckt` in v1, v1 scope = Phases 1-3 + `circuit { file : ... }` in `hphirun` + converting
`examples/circuit`. The plan header now folds in the DR-39 2026-08-13 verification (two Newton
copies, `tEpsilon0` top-vs-bottom drift) and gates the example conversion on DR-115.


### [fix_piecewise_degenerate_window.md](closed/fix_piecewise_degenerate_window.md)

**Status: PLAN — direction decided (new opt-in `resistivity type : riva`, Duron parallel form,
2026-08-27 Christian), audit round dispatched.**
The tapestack3d_coarse timestep-122 abort: `drho_piecewise_dT` returns NaN through its own
degeneracy guard when the sp-ap table's n(T,B,θ) approaches 1 near T_crit (`j2` overflow →
`tDisc = inf − inf`; the guard's fallback `tFrozen` is computed upstream of the guard, already
poisoned). Same overflow route exists in `rho_piecewise` and `drho_piecewise_dJ` behind passing
`n > 1` asserts. Plan: hoist/harden the fallbacks (R2-R4), probe sweep (R5), deck gate (R6).
Open: fallback value in the n→1 window (O1), Duron-parallel model as post-release default (O2),
PL-branch residual/tangent inconsistency (O3). Literature: Duron et al. 2004, Riva 2021 (EPFL
thesis 8754).

## Added 2026-08-26

### Sweep 2026-08-26 — active set 47 → 33

A currentness sweep over `todo/*.md`, classified from each file's own Status line and its
checkbox state, cross-referenced against the live and archived register rows. **Fourteen files moved, none deleted.**

**→ `closed/` (7)** — work finished, nothing forward left:
`dr92_restart_step_ramp.md` (COMPLETE, verified by execution, five-step round finished) ·
`thermal_picard_freeze_dr08.md` and `bdf_nonlinear_mass_verification.md` (both say CLOSED in
their own Status lines; DR-08, DR-33 and DR-64 are all struck) ·
`tapestack3d_J_asymmetry.md` (investigation complete, verdict delivered) ·
`coulomb_gauge_penalty_report.md` (a closing report by its own title) ·
`instruction_doc_currency.md` (R1–R4 all done, DR-60 struck; only the open *question* O3 remains) ·
`nonlinear_iteration_strategy_near_quench.md` (its own Status: "THEORY REFERENCE — keep it for its
literature argument, not as a work plan").

**→ `deferred/` (5)** — sound, parked, nothing scheduled:
`2d_thinshell_todo.md` (PAUSED, sprint calendar expired) ·
`ngspice_parser_plan.md` (pure proposal, no code ever written) ·
`jacobian_init_remaining_bottlenecks.md` ("Open, not scheduled … only worth attacking if that 23 s
starts to matter") · `handoff_scls_strumpack_openmp.md` (an SCLS build-system issue, not a BELFEM
one) · `gauging_context_from_the_quench_run.md` (an evidence pack for gauging work nobody has
scheduled).

All references inside `todo/` were rewritten to the new paths, and three **pre-existing** broken
links in this index — left by earlier moves that never updated it — were repaired
(`anderson_picard_acceleration_plan.md`, `2d_thinshell_gap_analysis.md`,
`coreduce_complexPellikkaGeneralized_performance_findings.md`). `todo/` now has zero broken
markdown links.

**→ `nonfree/todo/` (2), on Christian's ruling.** `fvm_implementation_plan.md` and
`fvm_module_next_steps.md` described `src/fvm`, which **no longer exists** — the module moved to
`nonfree/fvm` (DR-41, 2026-08-14). Open-source planning documents for a proprietary module are
what the nonfree convention exists to prevent, so both were moved into the nonfree repository and
removed from this tree. Each carries a relocation banner; their `src/fvm` paths were rewritten to
`nonfree/fvm`, and their pointer to `fvm_mpfa_o_vs_fem_recommendation.md` — which was **already
broken here**, since that file has always lived in `nonfree/todo/` — now resolves.

- **[run_gate_batches.md](closed/run_gate_batches.md)** — **CLOSED 2026-08-30 as obsolete** — the executable residue of the debt register,
  batched for an autopilot session. All 44 live rows now carry a waiting-on tag in the ID cell
  (`[RUN]` / `[RUN-BLOCKED]` / `[RULING]` / `[CODE]` / `[MIXED]`). **Only 2 rows are genuinely
  run-only**; 15 are waiting on a sentence from Christian, not on a machine. Three gates are
  `[RUN-BLOCKED]` because the dump or deck they name was deleted — DR-90's t=2300 ms memdump (only
  the t=7.1 s one survives), DR-105's `build/tape_quench` (`build/` is gone), and DR-79's pre-fix
  binary. Classification pass from status columns as written; no source re-verified.

- **[debt_register_closed.md](debt_register_closed.md)** — the register's archive; **83 rows by
  direct count late 2026-08-27**, after a sweep moved six struck rows out of the live register
  (DR-40, DR-79, DR-82, DR-83, DR-107, DR-114) and DR-103 joined the same evening (its stale
  second-order-weights todo ruled out, residue split into DR-122; the sweep's recorded "84" did
  not match a recount of 82 before DR-103 — one row of drift, unreconciled). ~~Live register: 38
  rows.~~ ~~**Recounted 2026-08-30: the archive holds 121 rows and the live register 24.**~~
  **Recounted 2026-08-31: the archive holds 136 rows and the live register 13** ( 9 `[P]`,
  4 `[W]`, the `[F]` list empty ), after the 2026-08-31 archival pass moved DR-109, DR-111,
  DR-151, DR-152 and DR-154 out, and DR-135, DR-136 and DR-137 were struck and archived the
  same night — the background-field / periodic-Dirichlet cluster, closed on evidence, a design
  ruling and an accepted-exposure ruling respectively. It was seeded when the 65 struck
  rows moved out of `debt_register.md` (108 rows in, 65 + 43 out, round-trip verified lossless),
  and the register's 233-line pass-by-pass preamble was replaced by a ~115-line **operating
  manual**: how to read a row, when a row may be struck, how a sweep gathers evidence, the
  closure-evidence rung ladder, five standing traps, and the run/ruling/code triage lens. Codex
  and Grok distilled the rules independently; Grok's check found two pieces of the old preamble
  that had become false — `[seeded — confirm]` is claimed of "every row" but appears in **zero**
  table rows, and the 08-23 triage census names rows that are since struck. Both dropped.
  Christian's rule for the round: "we are not interested in a changelog of the debt register,
  that's why we have git."

## Added 2026-08-25

- **[shared_library_and_install_plan.md](closed/shared_library_and_install_plan.md)** — build BELFEM as a
  shared library, add `make install`, and give both trees working RPATHs on Linux and Darwin.
  Three findings decide the shape: the 24 module archives carry **five mutual dependency cycles**
  (`comm<->core`, `core<->io`, `integration<->mesh`, `interpolation<->kernel`, `iwg<->kernel`,
  measured with `nm`), which ELF tolerates in a `.so` but Mach-O `ld64` does not — so the plan
  builds **one** aggregate `libbelfem` from per-module OBJECT libraries rather than per-module
  `.so`s; `-fPIC` currently reaches C++ only, so the three `.f90` files in `libbelfem_sparse` would
  fail to link into a shared object; and `BELFEM_RPATH`, collected by fifteen config files, is
  never emitted as `-Wl,-rpath` at all — it only feeds `LINK_DIRECTORIES`, so even build-tree
  binaries depend on `LD_LIBRARY_PATH`. `share/` and `examples/` install via `install(DIRECTORY)`,
  which globs at install time and so tracks content churn without a re-configure. Two items need
  Christian's sign-off before R1: **O2**, losing the `make libbelfem_<module>.a` per-module gate
  used across the devlogs and the probe workflow, and **R8**, the one `src/` change — a compiled-in
  `BELFEM_INSTALL_DATADIR` fallback, without which installed binaries find no gas tables or
  material databases. **Status: OPEN — drafted, not audited, no source modified.** Linux first,
  Darwin gate (R10) right after. **Amended 2026-08-26 with R8b, which is separable and worth
  landing on its own:** `add_test` sets no `WORKING_DIRECTORY`, so ctest runs each test from its
  own build directory, from which the three-level CWD climb in `fn_GT_data_path.cpp` resolves to
  `<build>/share/fluid` and can never exist — and `fn_material_data_path.cpp` has no CWD fallback
  at all. The data-dependent tests therefore self-skip unless the developer happens to have
  exported `BELFEM_DATA`, so **a green `make check` does not establish that they ran.** Fix is one
  ctest `ENVIRONMENT` property pointing at `${CMAKE_SOURCE_DIR}/share` — harness only, never a
  compiled-in define, or the configuring machine's source path ships inside `libbelfem`.
  **Status 2026-08-27: R1-R8b and R10 implemented and gated on Darwin, R4 and both halves of R9
  gated on Linux — `make install` produces a working tree on both platforms.** O2 and R8 got their
  sign-off; O7 (`gComm`/`gLog` in `main()`) and O8 (a moved install still needs `BELFEM_DATA`) are
  open follow-ups, and the Linux gate added **O9**: neither plugin template can compile against an
  install, because both put `${BELFEM_DIR}/include` alone on the include path while the headers go
  to `include/belfem/<module>/`. **Owed: R11 (flip the default, fix the templates) and R12.**

## Added 2026-08-22

- **[coulomb_gauge_penalty_report.md](closed/coulomb_gauge_penalty_report.md)** — closing report of the
  Coulomb-gauge / G-operator campaign. Verdict: element-local ∇h penalties are vacuous on
  simplex/thin-shell meshes (the curl null space is element-wise constant, cuts included);
  the consistent K-channel form is γ = χ·ρ*/μ²; measured κ on the quench deck is flat and
  anti-correlates with solver difficulty, so a κ-only remedy has no measured problem to
  solve. The theory itself now lives in
  `src/fem/maxwell/doc/coulomb_gauge_penalty_theory.md`. **Step A signed off
  2026-08-22; Step B (EdgeFunction::G interface + stubs) landed and build-verified the same
  night; Step C COMPLETE 2026-08-23 for all 9 implementable elements (TET4, TRI3, QUAD4TS,
  PENTA6TS, HEX8, TRI6, TET10, HEX8TB, HEX8TS), each through its own C.1→C.4 pipeline with
  blind two-vendor rounds; DR-98 (EF_HEX8::C negated curl) found and fixed along the way.
  Final gates: test_fem 145/145, make check 14/14.** Still open from the
  note: the ρ* averaging fork and the provenance of the reported condition numbers.
  Devlogs: dl20260822_coulomb_gauge_stepA.md, dl20260822_coulomb_gauge_stepBC.md,
  dl20260823_coulomb_gauge_stepC_elements.md, dl20260823_hex8_g_and_dr98.md,
  dl20260823_tri6_tet10_g.md, dl20260823_hex8tb_hex8ts_g.md.

## Added 2026-08-18

- **[closed/doxygen_site_structure.md](closed/doxygen_site_structure.md)** — ✅ **DONE
  2026-08-18.** Three-AI sweep of the generated Doxygen site, then R1–R17. Three defaults that
  changed after `Doxyfile.in` was written had severed the whole hand-written tier
  (`IMPLICIT_DIR_DOCS`), killed every usage-guide TOC (`MARKDOWN_ID_STYLE`) and left 1437 requested
  graphs unrendered. Result: **1471 errors + 232 warnings → 0/0**, 2737 diagrams, 30 module index
  pages, a 30-group API topic tree where there was none, class descriptions 26/318 → 146/318, and
  `src/linalg` documented for the first time. Two code defects fixed on the way (a backend-dependent
  `eigen` tolerance, a `Vector<T>`/`Vector<real>` mismatch in `fn_BZ_linspace.hpp`) plus one tested
  API change (`eigen`'s `aAbortOnComplex`, and a new `eigen_sym`). `make check` 12/12.
  Residual: the Armadillo paths of `eigen`/`eigen_sym`/`linspace` are still uncompiled here.

- **[jacobian_init_remaining_bottlenecks.md](deferred/jacobian_init_remaining_bottlenecks.md)** — what is
  left of the Jacobian-init investigation after the `DynamicBitset` fix took init from 290 s to
  23 s: the hash-lookup, serial-graph-union and `symrcm`-fallback findings, plus two phases that
  were never measured. Not scheduled.

---

## Currentness sweep — 2026-08-09

Every active file except the two authored that day was re-read and cross-checked against the
tree, the devlogs and `git log`. **Five plans were finished and moved to `closed/`; nine had
materially stale status lines.** The one-line triage:

| file | verdict |
|---|---|
| [pid_timestep_controller_plan.md](closed/pid_timestep_controller_plan.md) | authored 2026-08-09 — not swept (it moved DRAFT → IMPLEMENTED while this sweep was running; entry below refreshed to match) |
| [hex8tb_phase2_fem_wiring.md](closed/hex8tb_phase2_fem_wiring.md) | updated 2026-08-09 — not swept |
| [maxwell_kernel_collapse_plan.md](closed/maxwell_kernel_collapse_plan.md) | **stale → refreshed.** Code-complete; R6–R10 ticked as superseded, D13 closed. R12 `[◐]`: the ≥2-rank §4.1 gate **and** the never-written `src/fem/maxwell/doc/` architecture note |
| [thermal_matrices_cleanup_and_newton_plan.md](closed/thermal_matrices_cleanup_and_newton_plan.md) | **stale → refreshed.** "Newton half open" is wrong: `T_h_newton` fills all four blocks and is wired. What is left is proof, not code |
| [anderson_picard_acceleration_plan.md](closed/anderson_picard_acceleration_plan.md) | **stale → refreshed.** Committed in `6b2a2b98`, not uncommitted; R7 done; **new blocking defect D4 recorded** |
| [2d_thinshell_todo.md](deferred/2d_thinshell_todo.md) | **stale → refreshed.** "Execution has not started" is wrong — E5 and D2 landed 2026-08-04 |
| [iterate_refactor_plan.md](cancelled/iterate_refactor_plan.md) | **stale → refreshed.** Substantially overtaken by the matfix port; 3 of 5 bugs survive, all re-anchored |
| `fvm_implementation_plan.md` · `fvm_module_next_steps.md` (since moved to `nonfree/todo/`) | **stale → refreshed.** STALLED: no `src/fvm` code change since 2026-07-09 and the module is excluded from the build |
| [nonlinear_iteration_strategy_near_quench.md](closed/nonlinear_iteration_strategy_near_quench.md) | **stale → refreshed.** Δt control moved to the PID plan; one recommendation was refuted by events |
| [debt_register.md](debt_register.md) | **corrected.** 8 rows closed/downgraded, 3 reworded, 2 added (DR-52, DR-53) |
| [rho_lambda_argument_convention.md](closed/rho_lambda_argument_convention.md) · [coreduce_…](deferred/coreduce_complexPellikkaGeneralized_performance_findings.md) · [thin_shell_overhang_bug_analysis.md](cancelled/thin_shell_overhang_bug_analysis.md) · [2d_thinshell_gap_analysis.md](closed/2d_thinshell_gap_analysis.md) | findings intact, **line anchors re-baselined** (all had drifted) |
| [thin_cut_nonunit_rectification_implementation.md](closed/thin_cut_nonunit_rectification_implementation.md) · [gasmodels_open_source_migration.md](closed/gasmodels_open_source_migration.md) | **boxes reconciled** with what actually landed |
| [gas_correctness_fixes.md](closed/gas_correctness_fixes.md) · [nitrogen_eos_completion.md](closed/nitrogen_eos_completion.md) · [restart_circuit_verification.md](closed/restart_circuit_verification.md) · [test_normals_and_pipette_coverage.md](closed/test_normals_and_pipette_coverage.md) · [powerlaw_jc_n_field_derivatives.md](deferred/powerlaw_jc_n_field_derivatives.md) · [ngspice_parser_plan.md](closed/ngspice_parser_plan.md) · [interface_node_duplication_coil_ferro.md](closed/interface_node_duplication_coil_ferro.md) · [matfix_controller_port_plan.md](closed/matfix_controller_port_plan.md) · [falsification_tooling.md](closed/falsification_tooling.md) · [side_edge_fusing_cut_aware_plan.md](deferred/side_edge_fusing_cut_aware_plan.md) | **accurate** — stamped re-verified, with new cross-links where a later finding touches them |
| moved to `closed/` | `cross_review_tooling` · `cut_pocket_removal_rules` · `periodic_cap_cut_emission` · `handoff_double_corc_d1_session` · `handoff_thermal_iwg_collapse_session` |

**The finding worth reading even if you read nothing else — DR-52 (new row, P0, already
fixed).** The residual reported under Anderson mixing was the linear solver's own roundoff,
not a nonlinear convergence measure. It was fixed on 2026-08-07 (refresh removed, ε is the
pre-update force residual for both paths, defaults rolled back to BDF1 + Anderson opt-in),
and the fix **reverses D5 of the Anderson campaign** — but it is **uncommitted and not
compiled**, sharing a working tree with the Picard line-search retirement and the PID
timestep controller. Three plans now gate on that one bundle building and running: Anderson
R8, matfix R9, and the PID R5 A/B. That is the single highest-leverage next action in this
directory.

*(Sweep note: DR-52 was first written up here as an open P0 straight from the jury devlog,
then corrected after checking the tree. The lesson generalises — a devlog's headline finding
and its same-day addendum can tell opposite stories, and the addendum is the one that
matches the code.)*

**Disposition questions this sweep could not answer** (they are Christian's):
FVM in or out of 1.0 — the module has been uncompiled for a month;
`deferred/ngspice_parser_plan.md` — a feature proposal with no consumer pressure, a candidate for
`deferred/`; and the `[seeded — confirm]` correction pass on `debt_register.md` plus the six
`devlog/campaigns/` pages, at least one of which (`controller_anderson.md`) is now factually
wrong.

---

## Currentness sweep — 2026-08-11

Second sweep of this directory, two days after the first, run against the tree rather than
against the plans. **Two plans closed, two deferred, two register rows struck, one row found
worse than recorded, one added.** The headline is the one found worse, because it breaks the
pattern every prior sweep established.

**DR-45(a) is the first row a sweep has found to be *worse* than recorded.** Every earlier
correction in this directory went one way — the debt was already paid and the register had not
noticed (DR-43's four gasmodels defects, DR-19's already-ticked items, B6 sitting done-but-unticked
for four weeks, DR-52 written up as an open P0 and found already fixed). This one went the other
way. Commit `5a2ddf81` fixed the worker-warning gate in `MUMPS::free` and in the Vector-RHS
`MUMPS::solve` — but there are **three** `INFO` sites, not the two the row asserted for a week,
and the Matrix-RHS `solve` overload still gated its warning branch on `this->rank() == 0`. A
worker-raised warning was still discarded there — including the `INFO(1) = +1` that the open half
of DR-45 is entirely about. Grok raised it in the in-flight (b) jury round; this sweep confirmed
it by reading both overloads rather than trusting the vote, and the session that owns the row
closed it in `78ea534d` alongside clause (b) within the hour. *(Line numbers are deliberately
omitted: the anchors this sweep read were the **pre-`78ea534d`** ones and `78ea534d` rewrote that
region. In the current tree both overloads — `cl_SolverMUMPS.cpp:243` and `:401` — route positive
`INFO(1)` through the new `MUMPS::check_warnings()` at `:1225`, and no rank gate remains on a
warning path.)* **DR-45 is struck, on Christian's instruction and as a recorded exception** —
the same call he made for DR-42 and DR-49 earlier the same day, on the same grounds: the design
work is finished and the residual is one run. **Struck ≠ verified, with full force.** The fix
compiles clean under `-DNDEBUG` and `-DDEBUG` with `-Werror` and has never executed; its risk is
a *false abort on a healthy deck*, which only a run can rule out. The gate lives on in the row's
status column, which is why that column is never struck.

**A second, smaller break rode in with it — the same failure in miniature.** The DR-45 text, in the
row and in this sweep's first draft, attributed the "an always-active abort broke every parallel
run" precedent to **DR-24**, dated 2026-08-10. Both are wrong. The incident is the `set_currents`
guard reverted in **`61bbfa11`, dated 2026-08-11**, and it has **no register row at all**; DR-24 is
the periodic free-cut `mPrescribedCurrents` row, struck and closed on 2026-08-09. The association
was not arbitrary — DR-24 had asked for a one-time assert on `set_currents` and the reverted guard
was a failed attempt at exactly that — which is why it survived a jury round, a row edit and two
independent reads. A plausible ID copied forward without being followed is precisely the defect
class this register exists to catch. Cite `61bbfa11`.

The general lesson is the one worth keeping: the register's error is **not** one-directional. Four
prior corrections all found debt already paid, which quietly built an expectation that a stale row
is a pessimistic row. One counter-example retires that.

**And a third finding, which is really the first two seen together.** This sweep produced three
slips that look unrelated and are not: the stale MUMPS line anchors (plus ~15 drifted sites in the
rho/lambda doc inventory and `fix_facet_masters` at `:1114`), the DR-24 mislabel, and "both `INFO`
reads" naming two of three sites. **In every case the pointer was correct when written and was
never re-followed.** None is a reasoning error; none would have been caught by thinking harder;
all three survived attentive readers, because a plausible reference *reads* right at a glance and
confirming it costs a `grep` nobody spends. `CLAUDE.md` already carries the countermeasure —
**anchor by a searchable token, never a line number** — learned the expensive way from
`doc/input_file_reference.md`. Today suggests it was scoped too narrowly: a register ID, a commit
hash and a prose enumeration are the same object as a `file:line`, and the line number is only the
instance that rots fastest. **Whether to widen that convention is a call for Christian**, listed
below; it is recorded here, not acted on. (Observation contributed by the session that owns DR-45.)

The one-line triage:

| file | verdict |
|---|---|
| [bfm_stale_cache_detection.md](closed/bfm_stale_cache_detection.md) | **→ `closed/`.** R8's schema half — the only thing holding it open — landed in a parallel session; `doc/input_schema.yaml` carries the `mesh_config_tag` vocabulary and the `fn_mesh_config_tag.hpp` rule. Confirmed on disk (a static read, not a gate). The gate itself ran: `hphirun` round-trips the tag, a thickness edit rebuilds naming the changed value, an old stamp-less `.bfm` rebuilds once. DR-21 struck |
| [example_deck_and_material_db_repair.md](closed/example_deck_and_material_db_repair.md) | **→ `closed/`.** All four defects closed, every DoD item ticked *against a run*, and all four of its register rows (DR-54/57/58/59) struck. Its one residual, O2, is carried out as **DR-65** so the move does not lose it |
| [coreduce_…performance_findings.md](deferred/coreduce_complexPellikkaGeneralized_performance_findings.md) | **→ `deferred/`.** Its own status line has said "deferred" since 2026-07-03; the move records that. The remaining item is sequenced behind the rectification tests, which are only part done |
| [slepc_eigensolver_integration.md](deferred/slepc_eigensolver_integration.md) | **→ `deferred/`.** Deprioritized *by measurement* on 2026-08-10, before this sweep — bookkeeping on a decision already taken, with a stated revival condition (re-measure ARPACK at CORC scale) |
| [instruction_doc_currency.md](closed/instruction_doc_currency.md) | **stale → refreshed.** R4 was done and unticked: DR-63 was closed in `7432385a` by *amending the convention* — `doc/literature_references.md` is the single source of truth, module docs link to it, 1-of-23 is the expected state. Also corrected: the checker guards **32** claims, not 21 (re-run this sweep, 32/32 hold). Only O3 is left, and it is a ruling |
| [iterate_refactor_plan.md](cancelled/iterate_refactor_plan.md) | **anchors re-baselined; all three bugs re-confirmed live.** The ω-ordering asymmetry was read out of the source, not carried over: `iterate_magnetic` sets then clamps (`:1569`/`:1570`), both other paths clamp then set (`:875`/`:914`, `:1756`/`:1759`) |
| [controller_picard_tapestack_regression.md](cancelled/controller_picard_tapestack_regression.md) | **anchors re-baselined** (`try_escalate_to_newton` `:542` → `:600`, the retracted-D10 residual site `:2286` → `:2788`). The D10 retraction re-read at its new anchor and still holds |
| [rho_lambda_argument_convention.md](closed/rho_lambda_argument_convention.md) | **inventory re-baselined and it GREW.** Three of six files had drifted and the list was short by ~15 sites — including a full `Metal::rho(real B, real beta, real T)` signature in `materials_contracts_and_invariants.md:314`. One **source** site joins them: `cl_Material.hpp:981`'s doc-comment still says `rho(B,angle,T)` while the declaration 400 lines above is T-first. Comment text only; recorded, not fixed |
| [maxwell_kernel_collapse_plan.md](closed/maxwell_kernel_collapse_plan.md) | **refreshed.** Both R12 sub-items still open — the architecture note is still unwritten (re-grepped: nothing under `src/fem/maxwell/doc/` or `src/fem/kernel/doc/` documents the *collapsed kernel set* — `h_ghost` is named in three module docs, but only as the ghost-penalty stabilization, never as one of five kernels that replaced 38 variants). But DR-02's **Gate A instrumentation now exists**: the env-gated assembly dump is committed in `bc578b5e` at `cl_FEM_Controller.cpp:921`/`:1577`. It does not revive Gate A; it means a future equivalence check has its mechanism |
| [thin_shell_overhang_bug_analysis.md](cancelled/thin_shell_overhang_bug_analysis.md) | **anchor re-baselined.** `fix_facet_masters()` moved +129 lines in two days (`:1114` → `:1243`) — the standing argument for locating by symbol |
| [pid_timestep_controller_plan.md](closed/pid_timestep_controller_plan.md) | **index entry was stale, file was not.** The README still said "uncommitted, not compiled"; it has been committed in `1d6ef305` since 2026-08-09 and the plan file has said so since 2026-08-10. Corrected here |
| [debt_register.md](debt_register.md) | **DR-21 and DR-45 struck; DR-45(a) found to be a half fix first (fixed by its owning session in `78ea534d`; the strike is Christian's recorded exception and the run gate survives in the status column); DR-02 amended; DR-65 added; the DR-24/`61bbfa11` mislabel corrected.** 65 IDs now |
| [2d_thinshell_todo.md](deferred/2d_thinshell_todo.md) · [hex8tb_phase2_fem_wiring.md](closed/hex8tb_phase2_fem_wiring.md) · [thin_cut_nonunit_rectification_implementation.md](closed/thin_cut_nonunit_rectification_implementation.md) · [input_conf_configurator_plan.md](deferred/input_conf_configurator_plan.md) · [bdf_nonlinear_mass_verification.md](closed/bdf_nonlinear_mass_verification.md) | **authored or updated earlier the same day — not re-swept; three claims spot-checked against the tree rather than accepted.** A7 is correctly ticked against `ca703ae9`; hex8tb's B7 reminder is correctly struck by DR-21; T8 is correctly "part done" against `942faba9` |
| [conditioning_diagnostic_backends.md](cancelled/conditioning_diagnostic_backends.md) · [gas_correctness_fixes.md](closed/gas_correctness_fixes.md) · [nitrogen_eos_completion.md](closed/nitrogen_eos_completion.md) · [falsification_tooling.md](closed/falsification_tooling.md) · [matfix_controller_port_plan.md](closed/matfix_controller_port_plan.md) · [interface_node_duplication_coil_ferro.md](closed/interface_node_duplication_coil_ferro.md) · [thermal_matrices_cleanup_and_newton_plan.md](closed/thermal_matrices_cleanup_and_newton_plan.md) · [restart_circuit_verification.md](closed/restart_circuit_verification.md) · [test_normals_and_pipette_coverage.md](closed/test_normals_and_pipette_coverage.md) · [side_edge_fusing_cut_aware_plan.md](deferred/side_edge_fusing_cut_aware_plan.md) · [gasmodels_open_source_migration.md](closed/gasmodels_open_source_migration.md) · [timestep_collapse_residual_floor_plan.md](cancelled/timestep_collapse_residual_floor_plan.md) · [timestep_collapse_mitigation_design.md](cancelled/timestep_collapse_mitigation_design.md) · [powerlaw_jc_n_field_derivatives.md](deferred/powerlaw_jc_n_field_derivatives.md) · [nonlinear_iteration_strategy_near_quench.md](closed/nonlinear_iteration_strategy_near_quench.md) | **accurate — stay active.** Spot-checked in source where a check was possible: MUMPS `get_forward_error`/`get_backward_error` exist (`cl_SolverMUMPS.cpp:1358`/`:1371` after `78ea534d`); the three nitrogen derivative bodies genuinely do not compile (`* *`, an undeclared `kk`, a dangling `* -`); `mt_thermal_h.cpp` is 129 lines with exactly `T_h_picard`/`T_h_newton` |

**Three links in the 2026-08-09 table above are now dead, deliberately.** They point at
`anderson_picard_acceleration_plan.md`, `2d_thinshell_gap_analysis.md` and
`coreduce_…performance_findings.md` in this directory; all three have since moved to `closed/` or
`deferred/`. A past sweep's table is a record of what that sweep believed *and where the files were
when it believed it* — repairing the links would quietly rewrite that. The live index below is the
one that must match the directory, and it does. (Checked mechanically: 104 links, those three the
only misses.)

**Sample of what "check the tree" meant here.** Six register rows were re-read at their cited
sites and each was exactly as written: DR-07's `real dbdT = 0.0` placeholder
(`cl_FEM_Calculator.hpp:2879`), DR-13's `uint16_t` counters (`cl_Vertex.hpp:62-63`), DR-26's
`mPairVerdict` built at `cl_CutProcessor.cpp:785` and consumed nowhere, DR-46's
`DISABLED_GhostElementContract`, DR-53's two `false` fuse flags
(`cl_ThinShellFactory.hpp:347-348`), and DR-64's signed dot-product fix live at both producer
sites.

One claim did **not** survive its challenge intact. The FVM plans say "no `src/fvm` code change
since 2026-07-09", but `git log -- src/fvm` shows two later commits: `e4c02ac8` (2026-08-06) edits
`src/fvm/CMakeLists.txt` and `4f2c11cd` (2026-08-07) adds `src/fvm/doc/README.md` and four other
`.md` files. The plans' *conclusion* holds — no `.cpp`/`.hpp` under `src/fvm` has changed since
`1173454d` — but "code change" is the wrong words for it, and a reader running that `git log` would
conclude the plan was stale when it is only imprecise. Corrected in the disposition note below. The FVM "no code change since 2026-07-09" claim survived a challenge: `git log` shows two
later commits touching `src/fvm`, but both are documentation only.

### What this sweep could NOT settle

- **The `make check-fast` gate has still never run.** (`USE_TEST` defaults ON since 2026-08-14, but the shared tree's cache predates the flip and needs a cmake re-run.) At the time of this sweep `USE_TEST` was OFF in the shared tree; the
  spline suite, `Hex8TbUnitCirculation` and the new `tests/math/test_GraphSpfa.cpp` have never
  executed. Every plan whose residual is that gate stays active with the gate named. Builds are
  Christian's.
- **Every "pending run" row is beyond a static sweep by construction** — DR-02's Gate B thermal
  coverage, DR-52's greg3 A/B, DR-19's five R7 checks, DR-64's ferro A/B, DR-15/18/22/34/38/40.
  A bookkeeping pass can confirm the code is in and the anchors are right; it cannot produce
  evidence that only a run produces.
- **DR-45** stays open on a run, not on a decision. The jury round completed while this sweep was
  running and reached its stop condition, and clause (b) was implemented and committed
  (`78ea534d`) by the session that owns the row — but nothing has executed. The gate is helix
  serial and on ≥2 ranks completing *without* the new error firing.

### Dispositions handed back to Christian

Unchanged from the 2026-08-09 sweep, and deliberately not acted on:

- ~~**FVM in or out of 1.0.**~~ **Ruled OUT 2026-08-14: the module moved to `nonfree/fvm`** (paths below are pre-move). `fvm_implementation_plan.md` + `fvm_module_next_steps.md`. Re-checked,
  and stated more precisely than the plan states it: `add_subdirectory( fvm )` is still commented
  out (`src/CMakeLists.txt:12`), and the last change to any `.cpp`/`.hpp` under `src/fvm` is still
  2026-07-09 (`1173454d`). Two commits *do* touch the directory since — `e4c02ac8` edits
  `src/fvm/CMakeLists.txt` and `4f2c11cd` adds `src/fvm/doc/README.md` plus four other `.md` files
  — so the plan's blanket "no code change" is loose where "no implementation change" is exact.
  **Stalled is not deferred**, and a month of silence is a scheduling fact, not a verdict — the
  files stay where they are until you rule.
- **`deferred/ngspice_parser_plan.md`** — a feature proposal with no consumer pressure. Reads like a
  `deferred/` candidate; that is a product call, not a bookkeeping one.
- **The `[seeded — confirm]` correction pass on `debt_register.md`** — roughly half the unstruck
  rows still carry severities and blocking-1.0 flags that are *proposals* from the 2026-08-05
  seeding, never your judgement.
- **DR-65** (new) — is the `Tape_Quench/BuiltinMat` deck still wanted?
- **Whether to widen the anchoring convention** (new). `CLAUDE.md` says the input schema anchors
  every key by a searchable token and never a line number. Three slips in one day — stale
  `file:line`, a wrong register ID, a prose enumeration short by one — were all the same failure:
  a pointer correct when written and never re-followed. The rule may deserve to cover *any*
  reference a reader confirms by recognition, not just line numbers. Cheap to state, and it would
  have caught the DR-24 mislabel that a jury round did not. An hour's work if yes; if
  no, delete it and take DR-02's §4.1 deck matrix down by one case rather than leaving a
  permanently-failing entry.
- **DR-45's run gate** — a helix run, serial and on ≥2 ranks, confirming the new hard error does
  not fire on a healthy deck. Both clauses are in the tree; nobody has executed them.

---

## Active Tasks

> **Currentness note, 2026-08-30.** This section is a legacy digest, not a live index of
> `todo/*.md`. Many entries below link into `closed/` or `deferred/` and are kept for their
> summaries — a link that resolves into `closed/` means the work is finished, whatever the section
> heading says. Treat the per-file **Status** line as the authority and the dated "Added" and
> "Sweep" sections above as the current index.

- **[handoff_scls_strumpack_openmp.md](deferred/handoff_scls_strumpack_openmp.md)** — **(2026-08-17)** HANDOFF to the SCLS session: `/opt/scls/mkl` ships STRUMPACK with `OPENMP_TASKLOOP` and `OPENMP_TASK_DEPEND` compiled OUT although both probes succeed on this host — an RPM build-environment artifact. Factorization time is flat across a 5× thread-count change (21.2 / 20.5 / 23.5 s), consistent with either a bandwidth-bound kernel or the missing tasking; the two are not yet separated. A like-for-like local build plus the one-line A/B is ready at `/home/christian/codes/strumpack-local/`.

### tapestack3d |J| Asymmetry Investigation (2026-08-15, overnight)

- **[tapestack3d_J_asymmetry.md](closed/tapestack3d_J_asymmetry.md)** - **(Three-AI investigation, verdict delivered + corrected; A1 fix applied; verification plan V1-V6 open.)** The "two outer coppers should look alike" premise is refuted (rho-less magnesia buffer → `DomainType::Buffer` φ-region, ghost facets skip it → two conducting halves per tape); the jagged face is a standing solution mode in the YBCO-side layers, INSIDE the 1e-11-converged solution — the overnight "tolerance never reached" claim was retracted next morning (the log's dB is 10·log10, so −110 dB IS 1e-11). Element/orientation/T-matrix/postproc layer handling all cleared. **A1 applied 2026-08-15 on approval: the real `_Volumes` SPR-weight misindexing bug (`cl_FEM_Kernel.cpp`), inert on this mesh, fixed + asserted (rebuild gate owed).** A3 ghost-penalty calibration and the V5 discretization-knob sweep are the follow-ups; A2 is obsolete. Devlog: `dl20260815_tapestack3d_j_asymmetry.md`.

### Alloy Database: Aluminum, Chromium, and the Synthesized SAE 310 Gate

- **[alloy_database_al_cr_sae310.md](closed/alloy_database_al_cr_sae310.md)** - **(2026-08-15)** DRAFT, not yet audited. Add Al and Cr as `Metal` feedstock (Al: RRR-critical stabilizer, Woodcraft 2005 lambda; Cr: plain `Metal`, the 311 K Neel/SDW anomaly deliberately smoothed — pure Cr is alloy feedstock only), then gate the whole element pipeline by synthesizing SAE 310 through `Alloy::create_tables()` and diffing against NIST cryogenic fits 4-300 K (acceptance: cp ~5%, lambda ~10%, rho ~5%). The Alloy architecture is already right — disorder resistivity is an alloy-level `rho_0` parameter, not mixed — but one extension is required: a **lattice heat-conduction channel** (default-off `a·T^b`), since the purely electronic `1/(w0+w_i)` underestimates stainless lambda ~40% at RT and ~2x at 4 K. Open: lattice-term shape (O1), T_max cap for synthesized Fe-alloys vs elemental Curie-cp contamination above ~500 K (O2), mechanics feedstock — `create_tables` evaluates E/nu of every component unconditionally while Fe/Ni/Al/Cr mech is deferred, blocks the gate unless resolved (O3), 310 vs 310S target (O4). Kohler and mech for the ferromagnets are follow-on by scope guard. Groundwork devlog: `dl20260814_bezier_cp_copper_audit.md`.

### Nédélec Edge Function Defects (D1/D2/D3 from the LaTeX extraction audit + test battery)

- **[nedelec_edge_function_defects.md](closed/nedelec_edge_function_defects.md)** - **(2026-08-14, all three defects FIXED same day; R3 gate = probe-green 46/46, `make check` + 3D run open.)** **D3 (found by the new battery on its first run): rows 8-9 of the TET10 `mNeta`/`mNzeta` shape-derivative tables were crossed — η/ζ Jacobian columns identical, J singular at generic points; four entries fixed; the straight path had masked it because `link()` takes its Jacobian at evaluation point 0, which sat on an η=0 edge.** The battery (`tests/fem/test_EdgeFunctions.cpp` + `cl_EF_TestVolume` fixture: circulation identity, C-vs-FD-curl, EXODUS detJ guard, edge flips, curved-path equivalence over TRI3/TRI6/TET4/TET10) reproduces the D2 defect signature against the stale pre-fix library and runs 22/22 green against the fixed source. D1: `EF_TET4::E()` edge 2 used `∇ζ` where `∇η` belongs — found by Grok, confirmed by Claude, verified by an exact sympy probe cross-checked against DefElement's published degree-0 basis (coded circulation matrix = identity except +1/2 / −1/2 in column 2; corrected = exact identity); **Christian ruled and applied the three-line fix in-tree.** D2: the `EF_TET10::precompute()` tables carried η↔ζ exchanged table-wide (24 conformity violations, face dofs leaking ±4/3, ±8/3). Root cause = **one line**: the table generator `tmp/tet10/tet10_generate.m` used the naive node map (the EXODUS booby trap Christian named: node 1 must carry ζ or `det J < 0`). Fixed MATLAB-first on Christian's direction: generator corrected, pins re-pinned (Octave + Christian's parallel MATLAB, all residuals zero), 14 tables ported to C++, final gate parsed the edited file back into sympy — tables, derivatives, circulation identity and curl chain all exact; syntax check with build-tree flags clean. **Remaining R3: rebuild, circulation regression test (falsification-battery slice; two scratchpad probes are ready templates), 3D bulk-conductor run.** Devlogs: `dl20260814_nedelec_tex_extraction.md`, `dl20260814_tet10_table_fix.md`.

### Thermal Picard Freeze (DR-08)

- **[thermal_picard_freeze_dr08.md](closed/thermal_picard_freeze_dr08.md)** - **(2026-08-12, closed 2026-08-13)** CLOSED — DR-08 closed on Christian's ruling: R1 instrumentation verified live across the tapestack3d campaign, R2/R3 ran in the production binary, freeze never recurred; R1b/R4 retired unrun, R5/O1 scope notes stand in the file. Original entry: reproducer captured and fix applied the same evening — reproducer captured and fix applied the same evening, **reviewed not verified**. Two `tapestack3d` runs warm-started from the same `memdump.hdf5`, differing only in decomposition (8×2 vs 4×4), landed a ~2 % partition-order roundoff apart on opposite sides of `tolerance switch : 1e-4`: one promoted to Newton and converged, the other stayed Picard and froze bit-flat for four iterates — **the decomposition selected the algorithm**, and the promotion gate deadlocks (Newton is gated on progress only Newton can make). Fix, three parts: `KSPGetConvergedReason`/`KSPGetIterationNumber` after every `KSPSolve` (previously never read — a diverged or maxits solve was accepted silently; "0 iterations" printed at Detailed confirms the returns-its-input hypothesis for the freeze); `try_escalate_thermal_to_newton` — the thermal twin of the existing magnetic `try_escalate_to_newton`, fired from the stagnation detector at two flat iterates, once per attempt; and demotion hysteresis (Newton kept until the residual retreats a decade above the switch). Staggered `iterate_thermal` path deliberately untouched (O1). Gate: re-run the 4×4 reproducer — it must promote by the third flat iterate. Register: DR-08.

### PID Timestep Controller

- **[pid_timestep_controller_plan.md](closed/pid_timestep_controller_plan.md)** - **(2026-08-09)** IMPLEMENTED (R3+R4 landed) and **COMMITTED in `1d6ef305`** — the "uncommitted, not compiled" wording here was stale from 2026-08-09 and is corrected in the 2026-08-11 sweep; the plan file itself has been right since 2026-08-10. On Christian's execute order; O2/O3 resolved (abort at floor + C5 escape hatch; PID default). Remaining: round-3 diff jury (R6, in flight), R5 A/B validation (legacy vs PID, bdf1 vs bdf5 — Christian), O1 gain calibration. Replaces the memoryless iteration-target sqrt rule and the blind `Δt *= 0.5` cut with a log-space PID (the legacy rule ≡ pure-I, kI = 0.5, kept behind an internal member switch), enforces the `minimum timestep` floor with a loud abort, adds a post-cut growth hold, and persists Δt through memdump so a warm start does not re-enter the transient at the deck's initial timestep. Motivated by the sidecoatings Δt collapse. Devlog: `dl20260809_pid_timestep_controller.md`.

### Cut-Aware Side-Edge Fusing (D8)

- **[side_edge_fusing_cut_aware_plan.md](deferred/side_edge_fusing_cut_aware_plan.md)** - **(2026-08-06, R3 landed same day; D8 closed 2026-08-08)** IMPLEMENTED (committed), O1 REOPENED — the cut-aware single authority (`compute_side_authority` + rewired `connect_side_edges`/`connect_side_nodes`) is in the tree; the first fused Garber run confirms the D8 incoherence is fixed: the rim-strung, time-oscillating J blobs are gone. D8 itself was closed 2026-08-08 in the hex8tb tracker (verified end-to-end); the remaining R4 census and R5 physics gates are demoted to `mFuseEdges` re-enable preconditions — they gate neither the ship configuration nor the HEX8TB campaign. That evening, the run exposed a deeper issue: 90 of 348 rim stacks carry genuinely different λ branches above vs below the tape (Christian's insight, confirmed against the pre-fix dump), so single-branch fusing suppresses a real through-thickness MMF transition — rim edge-current depression where Norris peaking is expected. Decision (Christian): option (a), free rims for the release; `mFuseEdges` is switched OFF in source (`cl_ThinShellFactory.hpp:276` — the tree had carried `true` since b42044da, which is why Gregory's production runs were fused). Remaining: O1 disposition with Prof. Sirous (graded two-branch authority vs the HEX8TB wall element as the proper physics), the post-fix R4 probe dump (residual "mixing" at exactly the 90 λ-mismatch stacks is the O1 signal, not a bug), D9/D10 loose ends, and R5 physics gates before any re-enable. Original finding for the record: 61/151 rim stations constrained different interface levels to different cut branches of g (`Δφ` vs `Δφ ± λ`); row-level sign defects refuted (141,901/141,901 clean); Newton paralysis resolved separately (bearing gauge eigenmode doc).

### Matfix Controller Port

- **[matfix_controller_port_plan.md](closed/matfix_controller_port_plan.md)** - **(2026-08-06)** EXECUTED (Fable session, same day) — R1–R8 landed on the sideconnectors working tree, builds + port-relevant test suites green, theory doc updated; post-port amendment A1 (moved-baseline detector gated on an actually-moving thermal baseline, from the first R9 trace). Remaining: R9 trace gates (Christian), R10 thermal-tangent second wave (O5), D2-io resolved, O6 FVM pseudoinverse disposition. Original brief: harvest the Jul 27–28 controller/solver hardening that exists only on `matfix` (acceptance-hole fix, damped Newton entry, trust growth, stagnation-ω/latch hygiene, moved-baseline detector, solver soft-fail contract, absolute-tolerance semantics) onto `sideconnectors`, reconciled with the ts17/Anderson machinery. Root cause of Gregory's CORC convergence regression — two controller fix lines that never merged (devlog `dl20260806_greg_corc_convergence_regression.md`).

### Falsification Tooling

- **[falsification_tooling.md](closed/falsification_tooling.md)** - **(2026-08-05)** PLAN — awaiting Christian's approval: institutionalize the falsification side per the four-AI devlog review — D1 permanent interface/orientation regression battery (`tests/fem/test_InterfaceOrientation.cpp`: circulation identities for the four TS/TB edge functions, N=4 Ampère telescope test that would have caught greg2 in minutes, ghost 12-dof contract, 18-case orientation sweep; `make check-fast` on both backends from `env -i`; mutation acceptance = red on reintroduced `-mS[1]` flip), D2 evidence hierarchy wired into protocol + `/cross-review` + devlog template ("reviewed" ≠ "verified", evidence ladder, review stop condition), D3 campaign state pages + one-table debt register with the blocking-1.0 September lens. Physics expected values pend Christian sign-off; key open question O1 = EF fixture strategy around the `parent()->parent()->mesh()` thickness chain.

### Gas Models: Open-Source Migration

- **[nitrogen_eos_completion.md](closed/nitrogen_eos_completion.md)** - **(2026-08-06)** Finishing `EoS_Nitrogen` (Span et al. 2000) and wiring it into `Gas`. Six of nine functions are correct on disk (`phir_dt` rewritten and FD-verified — the old form was wrong by 100-650 %, invisible where the Gaussians underflow); `phir_dd`/`phir_t`/`phir_tt` do not compile. **The file carries the three outstanding bodies in full**, derived from the paper and sign-checked against the verified methane code, written as reconstruction insurance after an editor save failed to reach disk (only one copy and one worktree exist, so unsaved buffers are unrecoverable). Also: the FD gate to run before wiring, and the verified wiring checklist — the ctor takes no `Gas &` today, `init_tables()` does not exist (Tc = 126.192 K, pc = 3.3958 MPa, ρc = 11.1839 mol/dm³, Tt = 63.151 K, range to 1000 K / 2200 MPa), no vapour ancillary, plus the `HelmholtzModel` enum entry and the `cl_Gas.cpp:790-820` dispatch case.
- **[gas_correctness_fixes.md](closed/gas_correctness_fixes.md)** - **(2026-08-06)** Tick-off list from the two-round literature verification of `src/physics/{gastables,gasmodels}` (devlog `dl20260806_gas_literature_verification.md`). Live defects: methane Gaussian centers ψ↔γ swapped (EoS mechanically unstable at its own critical point, cv/cp/w ~2× off near-critical); methane λ_cr missing its χ^0.468 and crossover-F factors (helpers exist as dead code; λ +60-90 % gas-phase everywhere, no critical enhancement); methane pc 4.5922 → 4.5992 MPa (the paper's own front-matter typo, faithfully transcribed — adjudicated against the EoS's own critical-point evaluation, the NBP anchor, Friend 1989 and NIST); MC-SRK generalized c1 0.16054 → 1.60539 (Coquelet 2004); Chung β 1.368 → 1.3168 (Poling); two oxygen φ0-derivative slips. Latent: dead-code Soave/PR-78 digit errors, Lucas F_Q° omission, polar-β validity. Everything else verified digit-exact — hydrogen (all three variants), methane EoS + η/λ0/λex, oxygen (all 73 numbers), cubic machinery, Coquelet tables, `cubicalpha.inp`, NASA-9/CEA layer.
- **[gasmodels_open_source_migration.md](closed/gasmodels_open_source_migration.md)** - **(2026-08-03)** Phase 0 COMPLETE (under Codex+Grok audit), R1 next: move `nonfree/physics/{gasmodels,gastables,atmosphere}` into the open-source `src/physics/` tree, leaving finite-rate `combustion/` private. LBNL IPO approval obtained; chain of title cleared (code was Christian's own, `@dlr.de` stamps were cosmetic and are removed). The enabling finding: `gasdata.inp` has 471 records transcribed from a McGraw-Hill appendix, but they collapse to 263 addressable labels, only **32** carry complete critical data and only **18** reach full capability — so a clean rebuild from CoolProp/Cantera is small, and nothing needs to stay private. Phase 0 fixed six latent defects (D1 isobutane CAS `75-28-50`→`75-28-5`, never matched; D2 eight spreadsheet-mangled CAS numbers in the data **and** the factory that was written against them; D3 `uint` underflow in the fixed-column parser; D4 `std::string(getenv(nullptr))` UB on every `Gas` construction; D5 unreachable alpha-function fallback; D6 column overflow giving 2,5-lutidine M=257.65). Audited by Codex+Grok — Codex caught a missing `<cmath>` in the D5 fix and a miscount that surfaced **D9**: 73 duplicated labels, and placeholder CAS numbers (`800000-…`) winning last-in-wins for `H2`/`D2`, so hydrogen misses *both* its `cubicalpha` PM row and the factory's hardcoded CCR entry and silently falls back to a generic ω correlation. Carried: D7 label-nomenclature mismatch stranding records, D8 six-way `C7H9N` collision, D9, O4 (normal vs equilibrium hydrogen — physics call).

### Nonlinear Controller

- **[timestep_collapse_residual_floor_plan.md](cancelled/timestep_collapse_residual_floor_plan.md)** - **(2026-08-10)** OPEN, **no source modified, no run performed** — both `greg5` (2D tapestack) and `sidecoatings` (3D coated conductor) ratchet Δt to the deck floor and abort. The plan asserts the two decks have **different** dominant causes and must be tested separately (O3). `sidecoatings`: the t = 13.8762 ms log (timestep 694, 133 iterations) is reconstructed iteration by iteration with nothing left unexplained — ε plateaus *on top of* `mEpsilonSwitch = 1e-4`, so the Picard↔Newton handoff chatters (promotions at it 3/13/23/53/59/66, each costing 6–18 dB and 6–19 recovery iterations = 61 % of the budget, **D1**); the Newton ×2 "trust" growth doubles ω into the 1.0 clamp and regresses immediately, 2 of 2 times (**D3**); the escalation re-enters at ω₀ = 1.0 on a state observed unstable above ω ≈ 0.3 (**D4**). Once `mJustPicard` latched at it 82 the run converged monotonically −27.4 dB in 45 iterations to 1.2e-5 — and was then killed by the stall guard at a **measured MAD of 0.1993 dB against a 0.200 dB band** (**D2**), ~25 iterations short of the 1e-6 tolerance with 67 of 200 unspent. `greg5`: **D5, new** — `cmake-build-debug/greg5/input.conf:39` reads `scheme : bdf5`; `matfix` parsed only `method` so it silently ran BDF1, `sideconnectors` honours `scheme` so the same deck now runs **BDF5** with adaptive stepping (supersedes the companion file's §0.1 "no difference" row). **D6:** the `sidecoatings` current source linearly interpolates a 53 µs ADC table whose plateau dithers by one 0.0189 A count, handing the solver a ±357 A/s sign-alternating dI/dt where the physics is flat — which is why the failure appears after the ramp ends, not during it. **D7:** no Δt cut in any log says which of six guards fired. Ordered gates: R0 cross-review, R1 one-word `bdf5`→`bdf1` A/B, R2 the controller-off warm-start run (settles O1: was it convergence or a floor?), R3 the ε₁-vs-Δt slope the abort message asserts but nobody has measured, R4 dither falsification, R5 cut-reason instrumentation. Nine mitigations proposed (M1 trend-based stall test replacing the MAD scatter test; M3 handoff hysteresis + dwell; M9 last-resort fixed-ω mode before the floor), **none approved**.
- **[controller_picard_tapestack_regression.md](cancelled/controller_picard_tapestack_regression.md)** - **(2026-08-10)** OPEN, **no source modified** — Gregory's 2D tapestack (magnetic-only, `algorithm : Picard`) ran fine on `matfix` and struggles on `sideconnectors`, while the quench cases improved. Root cause candidate, found by the `--jury` round and agreed by all three voices: the new and retuned cut guards were tuned on hybrid Newton/thermal quench traces, and every one of them assumes a Newton rescue that a pure-Picard deck can never reach (`try_escalate_to_newton` bails at `cl_FEM_Controller.cpp:542`), so stall *detection* becomes a Δt *cut* instead of an algorithm demotion. Four live defects: D1 stall band 0.001 → 0.2 dB (cuts on any Picard phase slower than 3.8 %/iterate); D2 new progress watchdog, default ON at window 30 — shorter than the deck's own `target iterations : 50`; D3 `mOmegaNoiseBand = 0.05`, which suspends the strict damp-on-regression rule of Messe 2023 Eq. 14 and freezes the divergence counter, leaving a 0–5 % regression band with **no** damping mechanism (severity disputed: Claude high, Codex medium, Grok ~30 % — routed to Christian as O1); D4 PID + post-failure hold, ~3× slower Δt recovery per cut. D1/D2 are deck-recoverable (`stall tolerance`, `watchdog window`), D3/D4 have **no input key**. Plus side defects D6 (`mAndersonCommits` written 4×, never read, header describes behaviour the code lacks), D7/D8/D9 (validation + stale comment + over-claimed Eq. 14 citation). **D10 retracted:** the Picard line-search removal is *correct* — ε is the pre-update residual and exactly ω-independent (`cl_FEM_DofMgr_SolverData.cpp:2286-2290`), so the old reject loop re-solved a bit-identical system on all 8 backtracks. Everything is control-flow certainty, **nothing is measured**: R1 is the A/B gate (`stall tolerance : 0.001 ; watchdog window : 0 ;`). Exchange thread `tmp/ai_exchange/review_controller_tapestack_regression.md`.
- **[closed/anderson_picard_acceleration_plan.md](closed/anderson_picard_acceleration_plan.md)** - **(2026-07-30, CLOSED and moved 2026-08-11 on Christian's ruling)** Controller repair + opt-in Anderson(m) mixing on the Picard branch. R1–R7 landed, diff-audited by Codex+Grok, committed in `6b2a2b98`; the P0 residual defect this design introduced (ε was the linear solver's own roundoff) was found by a jury round and fixed in `4f2c11cd`, reversing the plan's D5 and returning Anderson to opt-**in** at depth 0/0. R8's end-to-end gate is closed on Christian's say-so rather than a recorded A/B — the run evidence question lives on as `debt_register.md` DR-52 (greg3 A/B).

### Side Connectors (HEX8TB)

- **[hex8tb_phase2_fem_wiring.md](closed/hex8tb_phase2_fem_wiring.md)** - **(2026-07-29, updated 2026-08-06)** OPEN — construction side is done; FEM plumbing is next. Landed: R1 edge function (3-AI consensus 2026-07-31: `E_k = s_k·F_k(η,ζ)·∇ξ`, exact-cuboid metric; campaign doc `src/fem/maxwell/doc/side_coating_wall_element.md`), R2/R3 nedelec dof count + `EF_HEX8TB` + factory case, the `h_side_connector` kernel scaffolding (Christian 2026-08-05, jury-reviewed with fixes applied), the `.bfm` persistence set B1–B6 (2026-08-04), and the ownership rule D6 (connector elements + recovery facets inherit the recovery-facet master's owner in `Kernel::partition_mesh`). NEXT: R5 plumbing checklist P1–P9: connector recovery-facet pass in `BlockData::link_thin_shell_facets_*` (P1), FieldList edge-h dof table + group activation in the MaxwellFactory type switch (P2), `link_to_group` dispatch to `h_side_connector` (P3), material assignment (P4), input keys for `mCreateSideConnectors`/`mConnectorWidth` (P5), then remove the WIP stop in `create_thin_shell` (P6); after that R4 interpolation-factory case, R6 Calculator/Postprocessor lists, R7 verification gate (numeric handedness = O1, wall current vs analytic r′, R0.3 probe, twisted-helix + cut-station regressions). Parallel deferreds (P9): ThinShellData connector record and block-width shipping in the Distributor. Open defects: D2 (LOW, side-facet edge table) and B7 (stale-`.bfm` trap: the checksum ignores connector settings — delete the cached mesh before the first connector smoke run); D8 closed 2026-08-08 (incoherence fixed and verified end-to-end; the residual census is a fusing-campaign re-enable precondition, and the connector path never goes through the fuse). Coated `h_in` variant (O2), closed-loop/threefold demo (O4), and selective layers stay behind scope guards.

### Maxwell Kernel Collapse — Follow-ups

- **[powerlaw_jc_n_field_derivatives.md](deferred/powerlaw_jc_n_field_derivatives.md)** - **(2026-07-21, ACTIVATED 2026-08-12)** PLAN, ~~future work behind a scope guard~~ — **the guard has been crossed by the case it named**: the `tapestack3d` deck attaches the measured SuperPower AP characteristic (`ybco { file : sp-ap.hdf5 ; }`), so jc and n are now `jc(T, log10|B|, θ)` splines and the `return_zero` B/β tangent is inconsistent rather than exact. Converged answers unaffected (residual is true, tangent approximate); the cost is Newton iterations where `d jc/d|B|` is steepest, i.e. the Ic crossing. Only the bottom layer is in place — `Database` differentiates its own spline analytically (`evaluate_derivx/_derivy/_derivz` over the table's three axes), so no finite differencing is needed there; the JcFunction hooks, the `drho_powerlaw_dB/_dbeta` family with its defect and piecewise overloads, the dispatch and the thermal jc(T)/n(T) leg are all still to be written. Analytic `drho_powerlaw_dB/_dbeta` for the Newton tangent once jc/n come from lookup tables. Full derivation (∂ρ_PL/∂jc = −n·ρ_PL/jc, ∂ρ_PL/∂n = ln(J/jc)·ρ_PL, shared parallel-combination chain factor), JcFunction `deval_dB/dbeta/dT` API (ModifiedKim analytic / Database table-backed / UserDefined fallback; the dT hook feeds the thermal solver's `drho_powerlaw_dT`, formula in O3), overload matrix incl. defect + piecewise branches, MaxwellData HTS-branch dispatch, R1-R5 steps + O1-O3 open questions (isotropized vs exact ∂β/∂q row; clamp-consistent derivatives; T-variant for the thermal collapse).

### Thin-Cut Generation for Non-Unit Coefficients (Coarse CCT Meshes)

- **[thin_cut_nonunit_rectification_implementation.md](closed/thin_cut_nonunit_rectification_implementation.md)** - **(2026-06-22)** Implementation plan: replace the greedy `Cohomology::clean()` body with `rectify_to_unit_or_certify` — a queue-based Bellman–Ford–Moore (SPFA) difference-constraint solve with negative-cycle extraction on the mesh 1-skeleton; periodic **quotient** via slave→master fold. Phased (0 diagnostic → 4 optional min-cost-flow). **Codex + Grok audited (2026-06-22):** adopted generator loop, `original_edges()` domain filter, master∪slave quotient adjacency, edge-ID certificate (+ folded self-loop = length-1 cert), Cochain side-state consistency, `BELFEM_ERROR` guards, `Cell`/ring-buffer state + int64 distances. **Analysis/theory home:** `src/homology/doc/thin_cut_nonunit_rectification.md` (geometry, existence obstruction `|⟨c,z⟩|>length(z)`, difference-constraint feasibility, the seven-cut-case enumeration, verified literature with DOIs); 3-AI analysis record `devlog/dl20260612_thin_cut_nonunit_algorithm.md`. (The earlier 2026-06-12 `thin_cut_nonunit_coefficient_algorithm.md` analysis was fully folded into those two and removed 2026-06-22.) **Status:** PLAN — audited, pending user approval; no source modified. **Update 2026-07-01:** second Codex+Grok audit round closed the QC (periodic quotient) and null-complex risks by static analysis; adopted shared `feasibility_solve()` core, `aPeriodicity == true` precondition, and corc-is-unit-coefficient reproducer correction. Decision (Christian): focus on the feasible-regime fix — greedy `clean()` replacement — first. **Ordered step tracker T1–T9 added** (solver core → Phase-0 diagnostic → decision gate → write-back → guards → messaging → diagnosability split → tests → validation); doc gained §1.1–1.3 on the non-unit pushed-object geometry. **Update 2026-08-25: IMPLEMENTED + PIVOTED 2026-07-16 (certify-then-greedy `clean_spfa()`, dense θ demoted to fallback, pocket census after); T1–T9 ALL CLOSED**: solver tests 2026-08-11 (16 SPFA cases), layer suite 2026-08-13, ctest gate 2026-08-23, forced-repair fixture + T9 corc discharge 2026-08-24, **DR-23 struck**. Only Phase 4 (min-cost L1 representative, unscheduled) remains. The 2026-08-25 jury on the double `clean_spfa()` run (constructor + post-`updatekGeneratorsFromHomology`) ruled both runs required; residues → DR-107.

### Maxwell / Thermal Matrix-Kernel Collapse onto MaxwellData

- **[maxwell_kernel_collapse_plan.md](closed/maxwell_kernel_collapse_plan.md)** - **(2026-07-13, swept 2026-08-09)** **CODE-COMPLETE — only the R12 verification gate remains.** `mt_maxwell_h.cpp` is 478 lines (`h_picard`, `h_newton_mu0`, `h_newton_mu`, `h_ghost` + the two side-connector kernels), `mt_thermal_h.cpp` is 129 (`T_h_picard`, `T_h_newton`), and both IWG dispatch trees are collapsed to algorithm-based picks. R6–R10 are ticked as *superseded by the 2026-07-21 leap*, not as individually verified — §4.1 was never run family-by-family, which is exactly what makes DR-02 load-bearing. **R12 stays `[◐]` on two counts:** the ≥2-rank §4.1 gate, and the `src/fem/maxwell/doc/` kernel architecture note, which was never written — no doc under `src/` mentions the collapsed kernel set, so the only prose about it is in devlogs. D13 closed (thin-shell layers are FEM Blocks; MaxwellData is built for them). Original brief: Collapse the 38 near-identical assembly variants (25 `maxwell::h_*` in `mt_maxwell_h`, 13 `T_h_*` in `mt_thermal_h`) into one generic Maxwell kernel + one thermal kernel that delegate all per-point material math to the new memoizing `calculator::MaxwellData` (`aCalc->maxwell()`), then shrink both IWG dispatch trees (Conductor/ThinShell cases → one assignment each). Full 38-row variant inventory with per-trap equivalence verdicts (§3.1 T1-T16: four β conventions, mu0 placement FP exception, 2D-j, HTS-only Newton term, clamp gaps, UserDefined overload parity, thin-shell alloy explanation, vector-aliasing scope). **Five new defects found in the current helper wiring (D1-D5):** D1 MaxwellData built for material-less air blocks → null-deref (`cl_MaxwellFactory.cpp:730-737` + `cl_FEM_Calculator.cpp:67,80`); D2 thermal-side helper binds a never-assigned (null) material instead of the Maxwell group's; D3 shared β cache slot aliases bj/bn conventions; D4 UserDefined reduced-dependency overload parity broken (8-arg defect powerlaw asserts full deps, `powerlaws.hpp:229-238`); D5 peer-element linking is nobody's job (legacy `get_thermal_calculator` contract). Steps R1-R12: helper fixes → shadow-compare harness → baselines → one dispatch family flipped per step (alloy → metal → HTS bulk → HTS thin-shell → UserDefined → thermal) → cutover. Regression matrix: helix (periodic), corc (thin-shell), Tape_Quench Builtin/Custom, 2D Validation set. Supersedes the R2-R3 collapse mechanism of `thermal_matrices_cleanup_and_newton_plan.md` (its Newton tangent work rebases onto the new `T_h`); hosts the decision point for `rho_lambda_argument_convention.md` (O10). **Status:** PLAN — drafted, pending Codex/Grok audit + Christian's approval; no source modified.

### Thermal Coupling: Producer Cleanup + Newton Tangent

- **[rho_lambda_argument_convention.md](closed/rho_lambda_argument_convention.md)** - **(2026-07-11, resolved in stages through 2026-07-14; swept 2026-08-09 — the doc pass is ALL that is left, and it is bigger than first scoped: 6 files, ~30 sites, re-baselined in the file)** Option (c) — **unify** — chosen and fully executed: rho/lambda/drhodT/dlambdadT are T-first `(T, normB, beta)` (2026-07-11, 3-AI audit closed ~7 stale sites); the private Kohler pointer family and user callbacks followed (2026-07-13); the last family, Jc/**powerlaw** (`rho_powerlaw`/`rho_piecewise`/`drho_*_dJ` T-bearing overloads), moved to `(normJ, T, normB, angleNxB [,x,y,z,t])` (2026-07-14, decided via `maxwell_kernel_collapse_plan.md` O10, executed by a parallel Opus session, Codex+Grok audit 0 defects; `JcFunction::eval(normB,angleNxB,T)` internal order deliberately kept). **Status:** RESOLVED — code done; remaining tail is the materials/dof/maxwell usage-guide doc pass (guides still show the old orders).
- **[thermal_matrices_cleanup_and_newton_plan.md](closed/thermal_matrices_cleanup_and_newton_plan.md)** - **(2026-07-09, swept 2026-08-09)** **BOTH HALVES IMPLEMENTED — verification open.** `T_h_newton` fills all four Newton blocks (∂M/∂T via dcp/dT against `collect_qhist`, the unsymmetric (Bᵀ∇T)⊗(dλ/dT·N) conductivity block, and the dρ/dT quench-feedback `dfdx`) and is live in dispatch; O1 was decided as option (c) — analytic `Material` derivatives through `MaxwellData::compute_drhodT`. What is left is **proof, not code**: the global FD tangent check (R5's verification clause), the `tape_hphiTrun` iteration comparison (R6) and the audit (R7), plus D7/O3 (volume-HTS β = π/2 physics call) and the D10 comment. The `T_phi` placeholder-properties task named in its scope guards is separately fixed. Original brief: Restructure `mt_thermal_h.cpp` (13 copy-paste producers / 1314 lines; only two real axes: volume-vs-thin-shell geometry × material law, which the generic `T_h` already dispatches at runtime) into single-sourced formulas, and fill the **never-implemented thermal Newton derivative matrices** — the ctor sets `dKdX_times_x`/`dMdX_times_x`/`dMdX_times_h` flags with a `//todo`, no producer fills them, so thermal "Newton" runs a frozen-coefficient `M + ΔtK` tangent that cannot see the quench feedback ∂ρ/∂T. Term recipes mapped onto the `assemble_dJdx` contract (incl. the missing `dFdX` source-derivative flag) per Bathe §7.2.2/Table 7.2, mirroring the magnetic `dmudH·H1/H0` and `Ctj⊗(jᵀC)` producer patterns. Defect ledger D1–D9 (incl. hidden per-int-point `Coords` allocation in all defect producers; volume-HTS β silently π/2). Steps R1–R7 (golden baseline → helpers → collapse to 2 producers → Material dT-derivative helper → Newton terms + global FD tangent check). Open: O1 FD-vs-analytic-vs-API derivatives (the existing `drho_*_dJ` family is the pattern for a dT sibling), O2 shared `compute_bn` home, O3 volume-HTS angle physics, O4 helpers-vs-template shape. **Status:** PLAN — Codex-audited 2026-07-09 (§3 math confirmed high-confidence, 8 corrections folded in), O-decisions with Christian; no source modified.

### High Priority: Solver / Time-Stepping Bugs

- **[bdf_nonlinear_mass_verification.md](closed/bdf_nonlinear_mass_verification.md)** - **(2026-07-03, closed 2026-08-13)** CLOSED — DR-33 closed on Christian's ruling: the tapestack3d campaign ran BDF1→BDF5 end to end (coupled, T-dependent ρcp, savepoint restores), retiring the parked V1 run half; the nine-case unit coverage stays in `check-fast`. Original entry: **ONE ITEM LEFT: V1.** B7 is done — `phi_ferro_newton` contracts `dMdx_times_h` against `IWG_Timestep::collect_qhist()` and the `#HACK` is gone. V1 (the BDF2-vs-BDF1 Newton-iteration regression on `thermalTest`) is what remains — still P2: the 2026-08-07 BDF5-by-default change that would have made it urgent was rolled back to BDF1 the same day (`debt_register.md` DR-33). Original brief: Verification/exactness tail after the BDF Jacobian scaling bug was **implemented and closed** (`closed/bdf_jacobian_scaling_bug.md`: coefficient wiring, startup ramp, α scaling, reset/savepoint — all in the tree and Codex-audited, initial run confirmed working). What remains: **V1** (was B5) BDF2-vs-BDF1 Newton-iteration regression test on `thermalTest` + optional HTS benchmark, and **B7** the `dMdX_times_h` producer contract — `phi_ferro()` (`mt_maxwell_phi.cpp:90`, `#HACK`) contracts the history derivative with `phi0` only, exact for BDF1 but an unweighted approximation of the β-weighted history sum for BDF2-5 (affects Newton *rate*, not converged correctness). **Status:** OPEN — core fix done; verification + one producer-side exactness change.
- **[nonlinear_iteration_strategy_near_quench.md](closed/nonlinear_iteration_strategy_near_quench.md)** - **(2026-02-03, swept 2026-08-09)** THEORY REFERENCE — partly implemented, partly superseded, one recommendation refuted by events. Δt adaptation exists (`adjust_timestep`) and its future belongs to **[pid_timestep_controller_plan.md](closed/pid_timestep_controller_plan.md)**; the **quench-detection trigger and the near-singularity iteration target are the only unclaimed ideas left**. Its "NOT line search / NOT Anderson" advice did not survive contact: both shipped and are load-bearing. Do not cite it as settled. Original brief: Theory-based approach to reducing iteration counts near critical material transitions. Literature analysis (Bathe §8.4.1, Zienkiewicz Vol 2) shows high iteration counts near quench are a **temporal-discretization-error** symptom, not a solver deficiency → adaptive Δt with quench detection (NOT line search / arc-length). **Status:** THEORY COMPLETE, ready for implementation. Expected 5-10× speedup.

### Eigenvalue Backends / Conditioning Diagnostic

- **[conditioning_diagnostic_backends.md](cancelled/conditioning_diagnostic_backends.md)** - **(2026-08-10)** Give `compute conditioning` a solver-native backend beyond MUMPS, so it stops leaning on an eigenvalue estimate that cannot resolve a badly conditioned mixed h-φ system. **The eigenvalue fallback stays live** (decided 2026-08-10 after a brief deactivation the same day): it is correct wherever κ is moderate and is the only option for solvers with no error analysis, so it now degrades to a warning plus NaN rather than aborting, and the deck parser warns that conditioning without MUMPS is not recommended. The reason is arithmetic, not tuning: ARPACK accepts a Ritz value once its error bound falls below `tol*|λ|`, and that bound cannot go below the matvec backward error `≈ ε·|λ_max|`, so the small end needs `tol > ε·κ`. A mixed h-φ formulation carries κ ~ 1e17 **by construction** (block scaling + constraint rows — characteristic of the formulation, NOT singularity), which asks for a relative tolerance above 1. Measured: exhausts the restart limit, converges nothing, at any tolerance/ncv. `EigenValues` is intact and still correct at the large end — it is the wrong instrument, not broken. **The replacement is already half-present:** BELFEM stores all 20 MUMPS `RINFOG` entries but exposes only COND1/COND2, leaving the componentwise backward errors and — more useful for "how hard is this to solve" — the **forward error bound RINFOG(9)** unread in memory. Both other libraries can produce the same family: **PETSc** forwards MUMPS's via `MatMumpsSetIcntl`/`MatMumpsGetRinfog` (`petscmat.h:2401,2409`, LU-with-MUMPS configurations only — its `KSPComputeExtremeSingularValues` is a different quantity and needs a Krylov method BELFEM does not use); **SuperLU** natively, and cheaply — `dgscon` (`slu_ddefs.h:210`) estimates rcond straight from the `mL`/`mU` factors the wrapper already retains from `dgstrf`, while `dgsrfs` (`:214`) returns the ADD `ferr`/`berr` pair. R1–R6; **O1 is the gate**: COND1/COND2, rcond, ferr/berr and the retired spectral ratio are four different objects and the footer prints one number under one name. Recommendation recorded: the forward error bound. Also records why we do not compute ADD ourselves — no transpose solve exists anywhere in the `Solver` API, and Hager-Higham needs `A^-T`. **Status:** OPEN — scoping done, nothing implemented.


### FVM Thermal Thin-Shell (MPFA-O — **moved to `nonfree/fvm` 2026-08-14**, out of the open 1.0 release)

**Both plans moved to `nonfree/todo/` on 2026-08-26** and are no longer in this repository.
`src/fvm` is gone; the module lives at `nonfree/fvm` (DR-41, Christian's 2026-08-14 ruling that
FVM is out of the open 1.0 release). Open-source planning documents for a proprietary module are
what the nonfree convention exists to prevent. Both carry a relocation banner recording the move
and flagging their stale `src/CMakeLists.txt` claims.

### Restart Follow-up: Coupled-Circuit Verification & Docs

- **[restart_circuit_verification.md](closed/restart_circuit_verification.md)** - **(2026-07-03)** Verification/documentation tail left after the coupled-circuit dynamic-state restart was **implemented and closed** (`closed/restart_circuit_state.md`). The `{time, delta_time, x, prev_x}` `/circuit` group, load-side scatter, R4 length guard, and Controller wiring are all in the tree and tri-AI audited; what remains is **R5** (confirm `IWG_Timestep` cold-starts cleanly on load), **R6** (end-to-end coupled restart test — none exists yet; must cover a source crossing the restart, a pre-dump fired switch, and a first-step reject), and publishing the `/circuit` format + the accepted standalone-L/C v1 limitation in `src/circuit/doc/circuit_usage_guide.md` (currently unmentioned). **Status:** OPEN — code in and audited; verification + doc only.

### Circuit Module: ngspice Importer & Standalone Runner

- **[ngspice_parser_plan.md](closed/ngspice_parser_plan.md)** - **(2026-06-28)** Proposal to read an ngspice/SPICE netlist (`.cir`) and build a `belfem::electronics::ElectricalCircuit` through the existing factory `create_*` calls; minimise the `circuit{}` block in `input.conf` (hybrid `file:` reference); and add a FEM-free `circuitrun` executable that mirrors `Controller::solve_circuit()`. Covers the SPICE→BELFEM element mapping (R/L/C/V/I/D map cleanly; `terminal pair` + `superconductor` need a BELFEM extension), `SourceFunction` gaps (SIN offset/delay, PULSE, PWL), the node-0-ground remap, the `M`-vs-`MEG` suffix trap, a proposed file/class structure (`NetlistParser`, `NgspiceCircuitFactory`, `spice_number_to_si`), and a 6-phase checklist. Companion to the new usage guide `src/circuit/doc/circuit_usage_guide.md`. **Status:** OPEN — **promoted out of `deferred/` 2026-08-27** (Christian requested it, closing DR-39's "parser that does not exist" leg by building it); Phase-0 decisions taken (6b linkage in `input.conf`, Option A comment directives, no `.subckt` in v1, v1 = parser + `examples/circuit` conversion); audit round 1 complete 2026-08-27 (three P0s adopted as §12 O-items O5–O9; stale §3/§7/§9.2 text struck; by-catch DR-115 diagnosis + DR-117/DR-118 filed); Phases 1–2 may start on approval, Phase 3b gates on O5–O8 sign-off and the DR-115 fix. *(Original proposal Codex+Grok cross-checked 2026-06-28.)*

### Maxwell Visualization: Interface Node Duplication

- **[interface_node_duplication_coil_ferro.md](closed/interface_node_duplication_coil_ferro.md)** - **(2026-07-03)** Extend the `InterfaceProcessor` duplicate-node scheme to ferro-air and air-coil interfaces for ParaView output: ferro-air duplicates hang on their originals with weight 1 (solve unchanged, postproc recovers each side independently → visible B/H jump at the iron surface); air-coil duplicates are **fully decoupled** (coil excluded from computation, nodal fields stay NaN). Key hazards identified: the `MaxwellFactory` ferro source-strip (`cl_MaxwellFactory.cpp:1892-1898`) severs weight-1-hung ferro duplicates, and `unify_duplicates()` derives the original↔duplicate viz pairing from the sources list, which the decouple path removes. Steps R1–R8 (ground truth → explicit classification → decouple path → strip replacement → side-local recovery/NaN → validation). Builds on Codex trace `devlog/dl20260703_interface_duplicate_node_trace.md`. **Codex audit round 1 (2026-07-03) folded in:** recovery is NOT side-local (`recover_fields()` canonicalizes to `original()` and writes to all duplicates — the split needs its own mechanism, O3); the ferro strip is the *previous deliberate decouple design* (`582131de`); master side resolved (`fix_facet_masters()`, higher `DomainType` wins → ferro/coil side is duplicated). **Status:** IMPLEMENTED + SERIAL-VALIDATED (2026-07-03) — R3 (coil decouple, unregistered sourceless duplicates) and R5 (side-local recovery + same-PP sibling merge: cuts continuous, interfaces crisp, coils dead) landed; per-side ownership fix in `Mesh::update_ownerships()` after Codex audit round 2. Open: R6 guard, R7/R8 parallel validation, R4 re-scoped (strip active yet run physical), O4/O5. See `devlog/dl20260703_interface_decouple_sidelocal_recovery.md`.

### Test Coverage

- **[test_normals_and_pipette_coverage.md](closed/test_normals_and_pipette_coverage.md)** - **(2026-08-05)** OPEN, opened by the `tests/old` triage: the eight `fn_normal_*` tests and the pipette checks were deleted rather than revived because they were written against the retired `FEM_geometry.hpp` / `geometry::` namespace; the functionality now lives on `Calculator::normal_tri_straight/_curved`, `normal_quad_straight/_curved`, `normal_hex` (`cl_FEM_Calculator.hpp:1249-1284`) and `cl_Pipette`, all of which have **zero** coverage. The reference fixture survived as `tests/fem/support/test_database.hdf5` and is still valid ground truth. Steps R1-R5 (confirm the outward-normal convention first — physics call → `tests/fem/test_Normals.cpp` → surface/volume via `Pipette` → CMake wiring without the old `/tmp` copy + global `gDatabase` pointer). Supersedes the "keep" verdict for `test/fem/fn_normal_*.cpp` in `closed/tests/existing_tests_triage.md`, which predates the header removal. `cl_IntegrationData_Interface.cpp` deliberately not restored — superseded by `tests/fem/test_FacetIntegrationPoints.cpp`; `tests/old/maxwell/` deliberately not restored — `MaxwellJob::initialize_test`/`run_test` and the `MAXWELL_HPHI_*` enum entries are all gone.

*(Caching Correctness and First-Run / Public-Facing both closed out on 2026-08-11 — see
[`closed/bfm_stale_cache_detection.md`](closed/bfm_stale_cache_detection.md) and
[`closed/example_deck_and_material_db_repair.md`](closed/example_deck_and_material_db_repair.md).)*


### Documentation & Instruction Set

- **[instruction_doc_currency.md](closed/instruction_doc_currency.md)** - **(2026-08-11)** IN PROGRESS — **R1–R4 all done; DR-60, DR-61, DR-62 and DR-63 all closed**; the only thing left is O3, a ruling on whether `BELFEM_ERROR`-on-non-convergence is correct in `ElementMapper::evaluate_general` (2026-08-11 sweep: R4 was ticked after verifying the amendment in the tree). `debt_register.md` DR-60…DR-63. The `CLAUDE.md` sweep that removed the proprietary `literature/` inventory from the open-source tree (and retired the `paperN` aliases, converting all 14 `src/**.cpp/.hpp` comment sites) found **nine claims false against the tree, four of which would have produced wrong code**: "no `share` overload exists for `Cell`" when `commtools.hpp:1648` has one — the same false row sat in `coding_philosophy.md:525`, from which `CLAUDE.md` had copied it; `aligned_alloc(64,…)` taught as the SIMD pattern though `src/` contains no aligned allocation and no backend gives 64 bytes; `BELFEM_ERROR` prescribed for convergence failures, which `coding_philosophy.md:621` calls a design error; and a Blaze/Armadillo default inversion **introduced by the sweep itself** and caught by the audit — the sharpest argument in the file for auditing a doc rewrite like a code change. Six further build/architecture facts were wrong (`-O3 -fno-exceptions`, `make tests`, the executable list, `fvm/` and `visualizer/` build status, "static and shared"), and two protocol rules had no presence in the bootstrap file at all — the **nonfree AI-exchange ban** (vendor leak surface) while the checklist told sessions to *use* the exchange, and "reviewed ≠ verified". `ai_collaboration_protocol.md:171` mandated the retired alias format *and* named `CLAUDE.md` as its definition, so retiring it in one file would have been inoperative under the precedence order; that line moved too. **Closed the same day (DR-60, DR-61):** `doc/literature_references.md` gained a Computational Topology and Optimal Cuts section (Gross & Kotiuga, Pellikka, Mrozek, Dey, Chen & Freedman, Dunfield & Hirani, Costantini, Haken, CLRS — all with DOI or ISBN) and a Lecture Notes subsection, closing the gap where the `homology/` module's own references had no public citation; and the alias conversion ran across 17 documents (117 occurrences), taking `src/` to zero and leaving only the decoder table and the dated devlogs. `paper5`/`paper6` are why that could not be a blind delete — both read "Alves et al. 2022" in the prose, so they had to become `2022a`/`2022b` explicitly. **DR-62 closed too:** Christian ruled to build the mechanical check — against Claude's recommendation to wait for the unwritten input-contract checker — and `scripts/check_doc_claims.py` now guards **32** claims, each anchored on a searchable token rather than a line number (re-run 2026-08-11: 32/32 hold). **Validated by running it against the pre-sweep `CLAUDE.md` from `HEAD`, where it independently finds 13 of the defects**; a green run on just-fixed documents would have proved nothing. Building it surfaced two bugs *in the checker*, both the same class: a probe returned an empty CMake branch, so the check reading it passed vacuously while appearing to run — hence the depth-aware `cmake_branches()` parser. Run by hand; no hook, no CI. **R4 / DR-63 closed 2026-08-11** (`7432385a`) by *amending the convention* rather than manufacturing 22 reference sections: `doc/documentation_guidelines.md` now names `doc/literature_references.md` as the single source of truth, says module docs link to it, and records that 1-of-23 is the expected state — an empty References heading being worse than none. Nothing compiled or run. Session record: `devlog/dl20260811_claudemd_stale_content_sweep.md`.

### Cross-Cutting Register

- **[debt_register.md](debt_register.md)** - **(seeded 2026-08-05, correction pass 2026-08-09, commit-state sweep 2026-08-10)** THE one table of open items across all campaigns, with the **blocking-1.0** column as the September release lens (protocol §11). 53 rows. Most still carry `[seeded — confirm]` — they were lifted mechanically from the June+ devlogs and their severities and blocking flags are *proposals* until Christian's pass. The 2026-08-09 sweep verified 13 rows against the tree: DR-03, DR-04, DR-05, DR-10, DR-20, DR-25 and DR-47 **closed**, DR-13 downgraded (uint8 → uint16), DR-02/DR-23/DR-33 reworded, and two new rows added — **DR-52** (Anderson residual = linear-solve roundoff, P0) and **DR-53** (wall-side fusing vs decoupled viz sheets → Δt collapse). **2026-08-10 sweep:** the working tree is clean and the previously-uncommitted bundle is in history — DR-52's Anderson fix in `4f2c11cd`, the PID controller in `1d6ef305`, the material work in `bc578b5e` — so DR-52 keeps only its greg3 run gate, DR-57/DR-58 record their commit, and DR-59's reachability is refreshed (the retired parallel builder is now *unreferenced* dead code, not abort-shadowed). **DR-54 is the next open item and needs a ruling, not code.** Aggregated view: `devlog/campaigns/release_1.0.md`.

### 2-D Thin Shells

- **[2d_thinshell_todo.md](deferred/2d_thinshell_todo.md)** - **(2026-07-01 v1 / 2026-07-10 v2; status corrected 2026-08-11)** The dependency-ordered 2D checklist, 10 of 53 ticked, **PAUSED** — the two-week sprint window from 2026-08-04 expired, so the day-by-day calendar is an ordering, not a schedule. **Christian reports the 2D thin shells are WORKING (2026-08-11)**, on a deck that is not in the repository — which most likely settles B8, E4, C3 and M1, none of it reproducible by anyone else. **T0 (author a 2D reference deck into `examples/`) is therefore the highest-value item on the list**, and would also give DR-02's §4.1 deck matrix its first 2D case. 2D periodic (C1, D5, M4) is deferred on Christian's instruction. Health warnings: citations baselined to 2026-07-10 have drifted, and the tick state has lagged reality (B6 sat done-but-unticked for four weeks) — verify against the tree before working an item.
- **[closed/2d_thinshell_gap_analysis.md](closed/2d_thinshell_gap_analysis.md)** - **(2026-07-01, ARCHIVED and MOVED to `closed/` 2026-08-11)** The procedure-by-procedure MISSING/STUB/WRONG audit behind the 2D checklist. Its reasoning remains the reference for how the 2D path differs from 3D; its currency does not — several gaps it analyses are fixed and every citation is baselined to 2026-07-01. Not a tracker, not maintained.
### Controller Refactor

- **[iterate_refactor_plan.md](cancelled/iterate_refactor_plan.md)** - **(2026-06-22, swept 2026-08-09)** ACTIVE but **substantially overtaken.** Its premise — that `iterate_magnetic()` is a stale copy missing the coupled path's logic — is now only half true: the matfix port brought across the Picard↔Newton handoff, the damped Newton entry, the watchdog re-anchor and the solver soft-fail path, and `impose_voltage_bcs()` (execution step 2) is extracted. **Three bugs survive and are re-anchored in the file:** the unclamped ω in `iterate_magnetic` (`set_omega` at `:1376` *before* the clamp at `:1377`, while coupled and thermal clamp first), the shared `mNumIterationsDiv`, and its cross-timestep leak (now zeroed on retry only). The **backtracking line-search parity gap is the largest remaining drift** — `tBacktracks` exists only in `iterate_coupled`. Recommendation: land the surviving bug fixes as their own commit; hold the helper extraction until the PID timestep work and the DR-52 residual fix have settled, since all three touch the same functions.

### Other Active Items

- **[ferro_undulator_mu_and_interface_bugs.md](closed/ferro_undulator_mu_and_interface_bugs.md)** - **(2026-08-12)** OPEN — two defects reported by Gregory from a ferro undulator benchmark, his changes uncommitted on his machine. **D1: `cl_FEM_Calculator.cpp:96` calls `constant_property( mu )` unconditionally**, which is an assert-then-return — so a ferro with a B-H curve *aborts in debug*, while release compiles the assert out, returns the stored `NaN`, and gets the right answer by accident. **His proposed fix is more dangerous than the bug**: it leaves `tIsConstantMu0` uninitialized on exactly the ferro path, which the bulk branch reads at `:297`, and a truthy garbage byte selects `compute_mu_0` — **the iron yoke behaves as vacuum**, silently, with the Newton tangent zeroed. `-Wno-error=maybe-uninitialized` means the build will not stop it. Fix is a short-circuit one-liner, **applied 2026-08-12 (R1)**; telling Gregory (R2) still owed. **D2: he commented out one of the three "Interfaces must be disabled!" errors** (`cl_IWG_Maxwell.cpp:349,354,359`) — but `mFunMKF` is never reset between groups, so an empty case body assembles that sideset with **the previously linked group's kernel**, quieter than the abort. Those types are not admitted by `MaxwellFactory::select_sidesets` (`:1940-1988`, the only `set_sidesets` caller), so the errors are a backstop and **the admission route is unidentified** — needs his exact message and deck. Also found: `InterfaceFerroAir` has the same fall-through for non-HPhi formulations — **closed 2026-08-12 (R6)** with an `else` + `BELFEM_ERROR`. Two jury rounds done: round 1 on the diagnosis (Codex confirmed; Grok leg failed on the wide brief), round 2 on the applied fix (**both auditors, both edits confirmed, no P0/P1**). R1 + R6 in source, syntax-checked, not yet built or run. **Gregory then answered R3 and refuted the D2 reconstruction (§8): his error is `fn_check_facet_orientation.hpp:78` and his model has no coils** — the real defect is **D4**, `check_facet_orientation` is 3D-only (needs a shared node *pair*; distinct 2D facets share one node) and `fix_facet_masters`' BFS aborts on any 2D cohomology model where a same-type facet touches a cross-type seed. His comment-out silently degrades tape orientation to gmsh element-numbering luck. **Fix applied (2D chain rule, R10) and jury-audited: both auditors accept, no P0**; the round refuted the O4 side-claim — interior-junction contacts can tear a same-type chain in the BFS (latent P1, pre-existing, endpoint contacts safe). **The R8 run then exposed D5** (§9): the 2D cut pipeline silently needs a *uniform master side per tape*, supplied only by gmsh element numbering — R10's selective flips break it and both sides get duplicate nodes ("No surface found on tape/boundary"). Fixed as **R13/option B** (no same-type BFS propagation in 2D; 3D bit-identical). §9.2 corrects the earlier verdict: Gregory's comment-out was harmless in 2D (flip-all preserves the invariant). Open: R11 (revert, corrected rationale), **D3** `InterfaceCondFerro` admitted-but-no-case (R9), O4/O5 fix directions, R12 raw-id comparison in the 3D loop, R8 debug run gating all four fixes.
- **[input_conf_configurator_plan.md](deferred/input_conf_configurator_plan.md)** - Python tool to author/validate `input.conf`: schema-driven validator + mesh introspection, with a wxPython skin. **IN PROGRESS (2026-08-11): the defect-repair half is done and jury-reviewed; the tool itself has not been started.** Landed: 8 code fixes across 6 files, 6 documentation fixes, the two-artifact policy in `CLAUDE.md`, and the new `doc/input_schema.yaml` (all 9 sections, 114/115 anchors resolving). The campaign found **six real BELFEM defects** independent of the tool — chief among them `tape`/`shell`/`curve`/`cut` parsing as valid topology types while building nothing, and `background dirichlet` corrupting a neighboring boundary condition. **The schema is audited but NOT finished** — four gaps (T1–T4) block the validator. Open: **O6** (close T1–T4 first, or start the wx skin against the schema as it stands), and a `corc`/`helix` smoke run — nothing has been built or run.
- **[thin_shell_overhang_bug_analysis.md](cancelled/thin_shell_overhang_bug_analysis.md)** - (anchors re-baselined 2026-08-09; the §9 relink is still unimplemented, so this stays ACTIVE — and `fix_facet_masters`, which its "upstream normalization exists" argument leans on, has since acquired a second known weakness on all-air 2-D meshes, DR-16). Overhang bug partially addressed (3-AI re-audit 2026-06-22): first-facet sampling fixed + `fix_facet_masters()` normalization; the robust §9 consistent-block relink remains unimplemented. Related: the multi-sideset block-id analysis, folded into `src/homology/doc/homology_usage_guide.md` (the standalone `doc/cohomology_block_id_bug.md` was retired 2026-08-14).

---

## Deferred ([`deferred/`](deferred/))

Parked behind the current production path. **Spot-checked 2026-08-09:** all seven are still
genuinely unimplemented (no `EF_QUAD9TS`/`EF_PENTA18TS` in `src/fem/interpolation/nedelec`,
no Nitsche coupling in the Maxwell kernels, no Buffer branch in `cl_CutFactory`), so the
parking is accurate. One entry's premise has changed, though:
**`maxwell_postprocessor_gpu_acceleration.md` says "implement element-field caching first;
GPU may then be unnecessary" — the caching is done** (see
`closed/maxwell_postprocessor_element_caching.md`, verified obsolete 2026-06-22). Its own
stated precondition is met, so it is due a keep-or-drop decision rather than indefinite
parking.

**Moved here 2026-08-11 (currentness sweep) — two files whose own text already said "parked":**

- **[coreduce_complexPellikkaGeneralized_performance_findings.md](deferred/coreduce_complexPellikkaGeneralized_performance_findings.md)** - Simplicial-complex performance ledger for `coreduce_complexPellikkaGeneralized()` (re-audited 2026-06-22, anchors re-baselined 2026-08-09). **Cold since 2026-07-03 and deferred by its own status line**; tracked as `debt_register.md` DR-29, not on the 1.0 path. The one live item — #1/#2a, porting `pCoreduce`/`coreduceOmit` to Algorithm 6.1 worklists — was explicitly sequenced *behind* `thin_cut_nonunit_rectification_implementation.md`, because the port changes which cells survive reduction and therefore the generator representatives, and T5's guards plus the T8/T9 tests are the Betti/cut-count instrumentation needed to validate it safely. T8 is only part done. **Revive when** those tests land, or when a profile puts `coreduce` back on the critical path.
- **[mumps_symmetric_triangle_extraction.md](deferred/mumps_symmetric_triangle_extraction.md)** - **(2026-08-29)** ALGORITHM DRAFT, explicitly not for implementation yet (Christian: "worth thinking about implementing tomorrow, but not today"). What it would take to make MUMPS `SYM != 0` supported instead of refused: extract the lower triangle from the COO triplet view the wrapper already hands MUMPS. Key simplification — `create_coo_indices()` materialises whichever index array the storage order lacks, so **CSR and CSC need no separate path**; and `set_indexing_base()` shifts `mRows`/`mColumns` together, so the `rows[k] >= cols[k]` predicate is base-invariant (both verified 2026-08-29). Algorithm splits by lifetime: a gather map built once per PATTERN, an O(nnz/2) value gather per SOLVE, no allocation in the solve path. **The extraction loop is the easy part** — the risk is that discarding the upper triangle silently SYMMETRIZES a nonsymmetric matrix, the same silent-wrong class the guard was added to remove, so §4 weighs four verification strategies and recommends structural-once-per-pattern plus numeric-on-first-solve. **AUDITED 2026-08-29 and corrected:** MUMPS accepts EITHER triangle (the draft's "requires lower" was false — lower is a BELFEM policy choice; what is fatal is supplying BOTH, since duplicates are summed); the audit also produced a BETTER validation design than the draft's (a transpose-partner map built during pattern setup, so symmetry is checked on every gather for one extra load and compare, instead of a first-solve check that is false comfort for a nonlinear assembly). O2/O3 resolved and struck; memory cost corrected 12 → 20 bytes/entry; **O6 is now blocking** — the proposed pattern-invalidation test is unsound and nothing else can be trusted until it is settled. **Status:** DRAFT, audited — not for implementation until O6 and O1 are decided.
- **[conditioning_petsc_native.md](cancelled/conditioning_petsc_native.md)** - **(2026-08-29)** R7 of the shift-invert plan, unblocked now that shift-invert works and prints every timestep: let a PETSc field report `sigma_max/sigma_min` from `KSPComputeExtremeSingularValues` instead of paying for an eigen round. **The number is NOT kappa_2** — it is the ratio for the PRECONDITIONED operator, which a good preconditioner deliberately makes small ( measured `kappa_2 = 6.44e7` on the thermal system; a healthy ASM-GMRES might report 1e2 ), so the footer must name it as a third quantity. Covers the GMRES case ONLY: PREONLY does zero Krylov iterations and has nothing to estimate from ( `MatMumpsGetRinfog` is the separate route there ), and CG needs `KSPComputeEigenvalues` instead. **Its accuracy runs opposite to the eigen path's** — the better the preconditioner, the fewer iterations, the worse the estimate — so R1 gates on an iteration floor ( O2, undecided ). Also proposes retiring the `!= SolverType::MUMPS` type-switch and the `reinterpret_cast` in `Solver` in favour of a `Wrapper` virtual, since `get_cond0()` already is one — sequenced AFTER the estimate is proven, never bundled. **Status:** PLAN — pending audit.
- **[conditioning_shift_invert_fallback.md](closed/conditioning_shift_invert_fallback.md)** - **(2026-08-28)** Implements Christian's ruling: report the metric the active solver supplies natively, fall back to **shift-invert** when it cannot, RAM is the user's problem, the diagnostic stays opt-in. Rests on a MEASUREMENT of the tapestack3d thermal Jacobian (`sysdump_thermal_3.hdf5`): symmetric to 3.5e-16, SPD, `lambda_min` 6.176e-8, `lambda_max` 1.830, **kappa_2 = 2.96e7**, small-end relative gap **9.5e-10** against 2.8e-2 under shift-invert. That gap is why no polynomial Krylov method reaches this end at any `ncv`/`maxit`/tolerance, and why the spectral fold is dead. Blocker re-verified: **neither wrapper can reuse a factorization** — MUMPS runs JOB 5 (factorize+solve) on every repeated call and has no JOB 3 path (`cl_SolverMUMPS.cpp:291-301`), PETSc re-pushes values and re-sets-up the PC every solve (`cl_SolverPETSC.cpp:103-131`) — so R1/R2 add an opt-in reuse flag that leaves production JOB sequences bit-identical. R4 notes BELFEM has **no symmetric ARPACK driver at all** (`dsaupd`/`dseupd` absent) although this matrix is symmetric. O1 (where the reverse-communication loop lives), O2 (three different quantities, one footer label) and O3 (retire the fold?) are open; O2 needs Christian. **Has a known-good answer to validate against.** **Status:** PLAN v2 ( 2026-08-29 ) — O1 resolved ( C++ master-loop + collective dedicated-MUMPS solves, NO PARPACK ), mode-3 ground rules verified against scipy's vendored driver ( §2.2, incl. the ipntr(3) trap ), C1-C7 folded into R1-R3; v2 audit pending; a SEPARATE Fable session implements from the file.
- **[arpack_small_end_configuration.md](cancelled/arpack_small_end_configuration.md)** - **(2026-08-28)** Make the eigenvalue conditioning fallback produce a number for the thermal system instead of `n/a`, and stop it costing 4747 ms per timestep to fail (measured on `tapestack3d`, n = 97095, PETSc+ASM, 4 ranks). The cause is structural, not tuning: `job = 0` is `which = 'SM'` in regular mode, and IRAM converges from the **exterior** of the spectrum of `OP = A`, with no resolution where a diffusion operator's eigenvalues are densest. Replaced by a **spectral fold** — `sigma = max|lambda| * (1 + 1e-2)`, `OP = sigma*I - A`, `lambda_min = sigma - mu`, both ends now exterior `'LM'` requests — which needs no factorization and therefore does not touch the JOB 5/6 refactorization blocker that deferred [slepc_eigensolver_integration.md](deferred/slepc_eigensolver_integration.md). `sigma` cancels exactly, so only the folded run needs a tight tolerance. Also: `ncv`/`maxit`/`tol` sized automatically (they were hardwired and no caller in the tree reached the setters), a broadcast `EigenOutcome` so ranks cannot diverge into a collective, and a latch after repeated failure. **Codex ran three rounds and changed the design twice** — folding about the signed `lambda_max` returns `lambda_max` on a spectrum straddling zero, and an underestimated `rho` does the same, hence the magnitude shift and the margin. Six defects recorded, including **D1: an OpenMP data race in the serial matvec** (`a`/`b` shared, not private) that predates this work and affects any single-rank multi-threaded run. **Status:** IN PROGRESS — R1–R5 landed, syntax-checked only, **nothing built or run**; R6 open, R7 (PETSc `KSPComputeExtremeSingularValues`, available here because the thermal deck is GMRES+ASM rather than the `KSPPREONLY`+LU the SLEPc plan assumed) waits on Christian's test run. Grok absent from every round — CLI cannot start on this machine.
- **[slepc_eigensolver_integration.md](deferred/slepc_eigensolver_integration.md)** - **(2026-08-10)** SLEPc as a second eigenvalue backend beside the repaired ARPACK-ng driver, reaching the small end of the spectrum through `STSINVERT` where mode-1 Arnoldi structurally cannot. **Deprioritized by measurement on 2026-08-10, before this sweep:** at the size measured the repaired ARPACK path is projected at ~70–150 ms in a release build, which a shift-invert factorization would not beat, so the case now rests on large problems only. The production line for the conditioning number is [conditioning_diagnostic_backends.md](cancelled/conditioning_diagnostic_backends.md) instead. Nothing here was refuted — R1–R7, the gap table, O1–O3 and Appendix A's four scaffold defects (chiefly `EPS_SMALLEST_MAGNITUDE` with no `ST`, the same trap the ARPACK driver was in) all stand. **Revive when** ARPACK is re-measured at CORC scale and stops being competitive.

- **buffer_cut_implementation_plan_v2.md** - Reroute the cohomology cut through the buffer layer's phi nodes (single-cut approach). Buffer domain plumbing is done; the cut rerouting itself is not yet written.
- **edge_function_quadratic_shells.md** - Quadratic thin-shell Nédélec edge functions (`EF_QUAD9TS`, `EF_PENTA18TS`). Linear shells are the active production path.
- **maxwell_postprocessor_gpu_acceleration.md** - Potential GPU acceleration for extreme-scale postprocessing. Implement element-field caching first; GPU may then be unnecessary.
- **slave_side_enrichment.md** - Slave-side bubble enrichment for both sides of h-φ thin-shell interfaces. Not yet implemented.
- **thinshell_automatic_sublayering.md** - Automatic Lobatto-point sub-layering for thick thin-shell layers (follow-up to the resolved corner-edge transition problem).
- **thinshell_hphi_formulation.md** - Per-layer H–φ formulation switch in the thin-shell stack (insulator / high-contact-resistance interlayers). Sections marked ⚠️ need re-derivation before implementation.
- **thinshell_selective_nitsche_coupling.md** - Selective DG (Nitsche) coupling at high-contrast layer interfaces. Long-term design concept — needs literature verification and prototyping.

---

## Closed ([`closed/`](closed/))

Completed, resolved, or done. Retained for the historical record.

**Moved 2026-08-26:**

- **[ferro_cut_flux_island.md](closed/ferro_cut_flux_island.md)** — **✅ CLOSED 2026-08-26, opened and
  closed the same day: overnight investigation → root cause → fix → gate.** Opened as a
  ferromagnetic-cut formulation question (the gantry's upper yoke arm solving to |B| ≡ 0) and
  **re-rooted by measurement to something material-blind**: `Basis::mNumberOfSources` was a `uint8_t`,
  and `CutSet::create_duplicates` hands `set_sources` a list of [N abstract current nodes, then the
  ORIGINAL node]. At the deck's 464 cuts, 465 wrapped to 209 — the original was dropped, the cut
  constraint became the constant 209·I, flux stopped crossing, and the region behind it was severed.
  The proof was a law, not an argument: 32 distinct N values from 425 to 464, every group pinned to
  exactly ((N+1) mod 256)·I. Christian widened the counter and added `BELFEM_ERROR` guards
  (`d0f4db7d`); a probe measured that the widening costs **zero bytes** at any width up to 64 — the
  counter sat in alignment padding, so the byte `uint8_t` was saving had never been saved. Gate G1 ran
  clean: 0 wrapped nodes (was 302), 1371/1371 cut jumps exact, 1366/1366 cut pairs with B continuous,
  no dead islands anywhere, the dead arm carrying 0.46 T. Three claims from the opening brief are
  refuted inside the file (the "inconsistent MMFs", the "mixed jump signs", and the iron-piercing
  framing itself — 205 of the 302 corrupted nodes were AIR). **Residuals:** gates G2 (np = 10) and G3
  (`make check-fast`) survive in the register; D4 became DR-110; Christian's ferro-cut input flag is
  preserved on its own merits but is explicitly *not* the fix for this. DR-108 struck under the
  DR-42/49 exception. Devlog: `dl20260826_dr108_uint8_source_overflow.md`.

**Moved 2026-08-14:**

- **[nedelec_tex_extraction_plan.md](closed/nedelec_tex_extraction_plan.md)** - **✅ COMPLETE 2026-08-14, same-day plan → two audit rounds → extraction → file audit → close.** The still-valid theory from Christian's pre-BELFEM LaTeX notes (`tmp/nedelec/`, Lagrange-multiplier era) now lives in the module docs: the weak-form chapter (fundamental lemma, divergence corollaries, L2 projection, heat warm-up, b-conform a/a-v marked derived-not-implemented, the implemented h-conform) in `src/fem/maxwell/doc/maxwell_weak_forms.md`, and the element derivation (N/E/B/C operators, Lagrange triangles with the Jᵀ pitfall, Whitney TRI3/TRI6/TET4/TET10 incl. the full TET10 set) in `src/fem/interpolation/doc/nedelec_derivation.md`; `tmp/nedelec/drift.md` records the drift, led by the missing cohomology machinery (deferred to Gregory's paper and thesis) and the missing hanging-node concept, plus 18 transcription errata. Both files passed a blind Codex + Grok equation audit (all weak-form signs re-derived clean, TRI3 an exact code match, TRI6 circulation-1/2 convention recorded, TET10 narrowed after Grok refuted the 1/2 overclaim); prose polished under the no-em-dash / keep-the-voice rules. Everything reviewed, not verified. **Residual: D1/D2 edge-function defects → [nedelec_edge_function_defects.md](closed/nedelec_edge_function_defects.md).** Devlog: `dl20260814_nedelec_tex_extraction.md`.

**Moved 2026-08-11 (currentness sweep):**

- **[bfm_stale_cache_detection.md](closed/bfm_stale_cache_detection.md)** - **DONE and VERIFIED END TO END 2026-08-11** (`e81bb411`, `6b5fa441`, `2c0a7382`); `debt_register.md` DR-21 struck. A `.bfm` caches the *enriched* mesh — cuts, thin-shell layers, coating walls, periodicity, hanging entities — while the reuse test compared only a checksum of the *raw* geometry, so editing a layer thickness or an `edge coating width` reloaded the stale cache **silently**. Fixed by a second stamp beside the untouched checksum (`meta/config` + `meta/config_text`); reuse now requires both. Design ruled by Christian after two blind jury rounds: a tuple of the mesh-defining **numbers**, FNV-1a over `%.12g` text, so a mesh mailed to another machine tags the same. **What verified it, and why the distinction matters here:** running it found two defects that compiling could not (an id list `1,3,5,7,9,11` collapsing to `1`; the terminal filter missing corc's `input curves` spelling), and the post-implementation jury found two more (a loaded `.bfm` dropped its stamp on re-save; circuit `terminal pair` ids were not in the tag at all, so `examples/circuit` reused stale cuts). Then `hphirun` closed it end to end: tag round-trips through HDF5, a 1.6 → 1.7 µm edit rebuilds *naming the changed value*, a restored deck reuses silently, and a genuinely pre-feature `.bfm` rebuilds exactly once. R8's schema half landed in a parallel session — `doc/input_schema.yaml` now carries the `mesh_config_tag` vocabulary and the rule that a new enrichment key outside a `wholesale` section must also reach `fn_mesh_config_tag.hpp`. Accepted gaps on record: material *definitions* untracked (Christian's ruling), `1:3` vs `1,2,3` causing a false rebuild.
- **[example_deck_and_material_db_repair.md](closed/example_deck_and_material_db_repair.md)** - **ALL FOUR DEFECTS CLOSED (2026-08-10), definition of done fully met, moved 2026-08-11.** Two defects blocked a fresh checkout's first successful run: a ≥2-rank cache-miss run aborted in METIS before the mesh was read, and a stale rho cache shadowed regeneration before dying with *"Dataset RRR … does not exist"*. The parallel rho builder turned out to have **never worked and never been needed** — its partition call disables `aSetProcOwners`, so it computed everything on rank 0 regardless; the branch is gone and its retired bodies were deleted in `6ce93549` (DR-59). D4 was closed by replacement rather than repair: Christian swapped `examples/corc` for his working 6-tape model, verified from a clean directory (61,418 nodes / 341,712 elements, `create_periodic`, both rho databases built in-run, 12 timesteps to t = 300 ms, warm restart from its own memdump). Along the way the failure's mechanism was corrected twice — it is **not** the gmsh 2.2 format, and the two directories were different models. Register rows DR-54, DR-57, DR-58 and DR-59 are all struck. **One residual carried out rather than closed:** the `buffer`-label `Tape_Quench/BuiltinMat` deck is now `debt_register.md` **DR-65**, and it needs a ruling, not code.

**Moved 2026-08-09 (currentness sweep):**

- **cross_review_tooling.md** - **DONE 2026-08-05, committed `bd73583b`.** R1–R11 all ticked; `scripts/cross_review.sh`, `scripts/install_autoreview_hook.sh`, `scripts/review_status.sh`, `.claude/commands/cross-review.md` and `METHODOLOGY.md` are all in the tree. The last `[◐]` (METHODOLOGY "drafted, NOT committed") was stale. One inherited confirmation for Christian — the `scls_env` g++ guard being existence-only — is parked in `falsification_tooling.md`.
- **cut_pocket_removal_rules.md** - **IMPLEMENTED 2026-07-16.** `Cohomology::remove_cut_pockets( bool aFireTierA )` runs inside `clean_spfa()` with `fire_node_coboundary()`; Rules 1–7 and O1–O4 are now ticked. Its field outcome was the decisive diagnostic of the campaign: the census exposed the SPFA *density* defect and then proved the `match_edges` abort was the Rule-7 seam-emission defect, **not** a removable pocket — exactly as the rules predicted. Tier-B enablement stays parked behind validation cases (a gate, not work in flight).
- **periodic_cap_cut_emission.md** - **FUNCTIONALLY COMPLETE 2026-07-16.** Option **(d)** — duplicate-aware edge rebuild — was chosen and shipped; (a)/(b)/(c) struck as not taken, and the file's own recommendation (try (a) first) is recorded as overtaken by the probe. The facet-mediated `match_edges` rewrite with combo-III ties, the half-cut policy and the trace-twin alias asserts are live in `cl_Mesh_PeriodicityFactory.cpp` and committed. Kept as the mechanism/formulation record for the seam. Residual half-cut-on-conductor case is a guarded **watch item**.
- **handoff_double_corc_d1_session.md** - **SPENT.** Every premise consumed: D1 is fixed (`Calculator::allocate` gates MaxwellData on non-Air blocks present in the peer dof manager), the "nothing is committed" inventory is committed (`ce6a0e8b`), and all seam probes are stripped. Residual §5 items were re-homed to DR-23/DR-26 and the thin-cut plan; its §4 gate numbers survive as reference values.
- **handoff_thermal_iwg_collapse_session.md** - **DONE 2026-07-21.** Both mission items landed the same day; T1–T7 and S1–S7 all ticked. Re-verified: `cl_IWG_MaxwellThermal.cpp` is 105 lines with an algorithm-based `T_h_picard`/`T_h_newton` pick, `mt_thermal_h.cpp` is 129 lines, and all 13 legacy `T_h_*` bodies are gone. The `T_h_newton` wiring that was "pending its Codex gate" is live.

- **sideconnector_bridge_plan.md** - **SUPERSEDED (2026-08-08), never executed.** The edge-wall impedance bridge (interior side edges fused into one `h_in` per station + curve kernel `r′·BᵀB` on [h_in; g]) was overtaken by the HEX8TB wall element campaign: the degenerate wall became a real degenerate volume element with its own edge function and is wired through the active **[hex8tb_phase2_fem_wiring.md](closed/hex8tb_phase2_fem_wiring.md)**. The coated `h_in` variant survives there as open question O2. Still-valid material: the D1 baseline finding (only outermost layers hang on air φ; interior side edges are free unknowns) and the recyclables inventory.
- **side_connector_effective_resistivity.md** - **CLOSED as a plan (2026-08-08); kept as the standing r′ theory reference.** Problem statement and modeling analysis for HTS tape side conductivity (rev. 4, Codex-audited): the wrap-resistance derivation `r′ = ρ_Cu·h_path/t_w + ΣR_ct/t_w` (slit-edge-interface-controlled, 10⁻⁶–10⁻⁴ Ω·m), the resistor-in-h-formulation argument (currents are jumps/circulations of H — no discrete curl, no new DOFs), and the failed-routes record (explicit copper mesh, in-element DOFs, enrichment, HEX8TS wrap). The implementation route is now the HEX8TB wall element (**[hex8tb_phase2_fem_wiring.md](closed/hex8tb_phase2_fem_wiring.md)**, theory home `src/fem/maxwell/doc/side_coating_wall_element.md`); this file remains the source for the r′ calibration numbers behind plumbing step P4.
- **meshfile_refactor_plan.md** - **DONE (2026-07-01, tri-AI audited).** Repaired and refactored the `.bfm` (HDF5) mesh save/load around the ID-keyed `ProtoMesh` reconstruction engine so a Maxwell run can stop, change parameters, and relaunch from the last timestep with a *fully enriched* mesh (skipping cohomology/cut/thin-shell factories on reload). Retired the legacy index-keyed `HDF5Writer`/`HDF5Reader`; added the `BfmFile` backend (base topology + groups + facets/edges/faces, node duplicates, hanging/T-matrix data, periodicity, thin shells + ghosts, curves/terminals, vertices, control points, base checksum) plus a separate `Mesh::save_fields`/`load_fields` restart file (fields + time cursor). The R12 end-to-end restart test passes on CORC (periodic + thin-shell + cohomology) in **serial and on 2/4/8 MPI procs** after the last two blockers — D24 (thin-shell facet aliasing serialized twice) and D25 (`Topology::run()` on the enriched mesh diverging from fresh, fixed via `run_on_enriched_mesh()`). Format documented in `src/mesh/doc/bfm_file_format.md`. **Residual follow-ups (additive, non-blocking):** coupled-circuit dynamic-state restart (R6b, → `restart_circuit_state.md`) and the `/meta/format_version` attribute. (R6a global variables is done.)

- **bdf_jacobian_scaling_bug.md** - **DONE (2026-07-03, Codex-audited).** BDF2-5 were non-functional: `compute_bdf_coefficients()` was never called, so `mAlpha`/`mBeta` stayed NaN → all-NaN system matrix (BDF1 worked only because it never reads them). Fixed in `cl_IWG_Timestep`/`cl_TimestepMatrices`: **B1** lazy coefficient wiring (`mCoeffsDirty`), **B2** startup order ramp (step 1 BDF1 → step p BDF-p, reject-safe via `mStepCount`), **B3** α=1 for order≤1, **B4** α scaling of `dMdX_times_x` in `assemble_dJdx` (Hairer & Wanner II.4), **B6** docs, **C1** (Codex) `mMethod` gate so temporary `MassOnly`/`StiffnessOnly` switches in `EigenValues` aren't hijacked by stale ramp state, **A1** `mOrder` in the no-stiffness branch (`bdfN_nok` had history depth 1), **A2** reset/retry step-size history. Same-day follow-ups (in `devlog/dl20260703_bdf_jacobian_fix.md`): Controller-owned input `method:` hook + two linkage fixes, and the **thermal savepoint** (`make_savepoint`/`restore_savepoint`) resolving the long-standing `initialize_thermal()` todo — a rejected magnetic step now restores all N thermal sub-step shifts instead of a single broken un-shift. Verification (V1) + the B7 history-derivative producer contract spun out to active **[bdf_nonlinear_mass_verification.md](closed/bdf_nonlinear_mass_verification.md)**.

- **restart_circuit_state.md** - **DONE (2026-07-03, tri-AI audited).** Coupled-circuit dynamic-state restart: the `/circuit` group in `memdump.hdf5` stores `{time, delta_time, x, prev_x}` and scatters it back into the rebuilt circuit on load (`ElectricalCircuit::save_state`/`load_state`, `update_components()`, R4 length guard, Controller wiring). Circuit topology is input-defined and never enriched, so the first `shift()` re-seeds the L/C/terminal-pair registers from the scattered solution — matching the FEM cold-start philosophy. The verification/doc tail (R5 IWG cold-start check, R6 end-to-end test, `/circuit` format doc) was spun out to the active **[restart_circuit_verification.md](closed/restart_circuit_verification.md)**.

- **timestep_controller_cleanup.md** - **DONE (2026-07-03).** Resolved the `#HACK` in `adjust_timestep()` (`cl_FEM_Controller.cpp:1098`): removed the invalid error-proxy machinery (`mEpsilonFirst*Δt`, dimensionally incoherent, would relentlessly shrink Δt), kept a clean iteration-count controller `sqrt(mIterationTarget/mIteration)` clamped `[0.5, 1.5]` (paper1 §4). The parked items (save-point restore of an unvalidated step; inert ω line-search — both actionable bugs) and the future-direction ideas (dual accuracy/convergence controller, contraction-rate signal, feed-forward pre-shrink, TR-BDF2) live on in the closed file for the record; revive as their own todo if scheduled.

- **user_material_backend_boundary_refactor.md** - **DONE (2026-07-02, verified in tree 2026-07-03).** Split the user-material API into a backend-light tier: spline ownership/evaluation moved out of `Material` into `material::SplineLookupTable`, with Metal/Alloy/Magnesia reparented to it; `cl_Material.hpp` is now fully backend-decoupled (no `cl_Spline.hpp`/`cl_Vector.hpp`/`cl_SpMatrix.hpp`), and the polynomial API takes `Cell<real>`/`std::vector<real>`. **Deferred, non-blocking:** R7 (harmless duplicate inline collapse), R11/R12 (full build + material-value check — only `-fsyntax-only` done), R13 (standalone user-material CMake/`MatData` must import host BELFEM compile defs + TPL includes), O3 (linalg-import policy). See `devlog/dl20260702_spline_lookup_table_refactor.md`.

- **periodic_volume_edge_seam_fix.md** - **DONE (2026-07-10, tri-AI audited).** Fixed the silently missing periodic ties on bulk-conductor volume edge DOFs at periodic seam planes, first exposed by the helix example (first solid conductor crossing a periodic plane; corc = thin shells only). Three defects fixed: D1 (cap facet wrappers carry no edges → `collect_edges` zero, `match_edges` silent early-out → zero ties; fix: master-element fallback in `collect_edges`), D2 (silent failure → named `BELFEM_ERROR` coverage check), D3 (BFM parallel reload: `unflag_all_edges` not-finalized branch called `tElement->edge(e)` for static per-type count on air elements with no containers; fix: `has_edges()`/`has_faces()` guards). Supporting fixes: PART 1 skip-guard in `create_hanging_edges_and_facets` (periodic wins on slave rim), deferred cascade flattening in `create_dofwise_t_matrices_master`. Verification: cap-pair Jz correlation +0.999 (was −0.96). Residual watch items (non-blocking): R5 corc regression, R6 parallel helix rerun, O3 λ₁/λ₀ = 0.9899 structural gap.

- **periodic_bc_fix_plan.md** - Master plan for periodic BC handling (DOF / MPI / input). The z-periodic CORC reproducer now solves end to end.
- **periodic_thin_cut_continuity_fix.md** - Cohomology thin-cut φ-continuity defect under periodic BCs; reproducer solves (`devlog/dl20260618_periodic_corc_solves.md`).
- **periodic_dofs_as_hanging_dofs.md** - Design note: a periodic slave DOF is a single-source hanging DOF, reusing the static-condensation / `mHangingDOFs` stack.
- **periodic_input_extension_plan.md** - Explicit `topology{periodic{source;target}}` sideset syntax; implemented.
- **periodic_parallel_source_closure_fix.md** - Distributed-mesh source-closure gap that aborted z-periodic CORC in `MUMPS::solve()`; fixed (`devlog/dl20260620_periodic_parallel_source_closure.md`).
- **ghost_penalty_parameter_provenance.md** - Literature trace of `eta=4.0` / `k_reg=1e-3` in `h_ghost()`: both ad-hoc engineering picks (method + harmonic-mean weighting are literature-grounded; the constants are not). Includes a calibration/test plan.
- **beta_angle_normal_sign.md** - RESOLVED — the β-from-tape-normal computation now wraps the dot product in `std::abs` at all 9 sites (`mt_maxwell_h.cpp`, commit `d2d3a41d`), forcing β∈[0,π/2]; the `β→π−β` sign ambiguity is gone.
- **corner_edge_transition_problem_2026-03-13.md** - RESOLVED — not a code bug; insufficient through-thickness layering exposed by the (correct) edge-deduplication fix.
- **maxwell_postprocessor_review.md** - Review complete — no optimization needed.
- **maxwell_postprocessor_element_caching.md** - Obsolete (3-AI verified 2026-06-22): the proposed per-element field cache is effectively implemented. The old node-patch `process_recovery()` was rewritten; base recovery (`cl_FEM_Postprocessor.cpp:927-980`) now visits each element once, and `compute_element_data()` (`cl_MaxwellPostprocessor.cpp:477-528`) integrates physics once per element with Gauss weights.
- **orientation_claude.md** - Consolidated facet-orientation / integration-point test plan; the orientation-table work it scoped is implemented (`test_facets.cpp`, TET20/TET35/HEX64 specializations).
- **tests/** - Test-suite planning campaign (modules 00–13); modules 1–10 implemented (`devlog/dl20260324_test_suite_bugs.md`).

---

## Completed (historical, no file)

- **Thin shell to conductor coupling** - Edge-to-edge DOF coupling (H-H with orientation signs).
- **Edge-to-node hanging DOFs** - H-φ interface (thin shell H to air φ).
- **Cohomology cuts** - Automatic generation via `src/homology/`.
- **STRUMPACK parallel interface** - MPI barriers and BLR tuning fixes.
- **Quaternion module** - Review complete, no open items.
- **Side-connector removal** - HEX8TS wrap removed as unphysical (`devlog/dl20260620_side_connector_removal.md`).
- [restart_timestep_consistency.md](closed/restart_timestep_consistency.md) — **CLOSED 2026-08-30 as obsolete** — warm-restart loader contract: the first post-restart step integrates the DUMPED Δt while the clock advances the CAPPED one (measured, reproducible, module not yet localized), plus DR-112's missing-BDF-triple hole. R1 is localization and blocks every other step. Replaces the withdrawn DR-149 under the extend-don't-branch policy.
- [Gauging context from the quench run](deferred/gauging_context_from_the_quench_run.md) — measured κ, the tolerance result, and the one test that decides whether a gauge term helps here.
