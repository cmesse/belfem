# Devlog — `todo/` Currentness Sweep and Relocation

**Date:** 2026-09-03
**Topic:** triage of all 57 active `todo/` plans against the debt register and the tree;
relocation of the non-current ones; repair of four false Status lines
**Module:** `todo/` (documentation only — no `src/` file touched)
**AIs involved:** Claude (triage + reconciliation + edits), Codex `gpt-5.6-terra`/medium
(batch A, 29 files), Grok `grok-4.6`/high (batch B, 28 files)
**Verification:** **Reviewed, not verified.** Nothing was compiled or run. The tree checks
below are `grep`/`ls`/`sed` reads at named `file:line`, which is the evidence ladder's static
rung. The one thing that would raise it — a build — was neither needed nor attempted.

## What was asked

Christian: for each `todo/*.md`, is it still current, and has the topic been addressed,
postponed or abandoned; then move the non-current ones into `closed/`, `deferred/`, or a new
`cancelled/`.

## Method

Pre-registered in `tmp/ai_exchange/todo_currentness_sweep.md` before either auditor ran, with a
stated falsifier: if the auditors merely echoed each file's own Status line, the round added
nothing and the Status lines are self-certifying. **The falsifier did not fire** — four Status
lines turned out to be factually wrong, and that was the round's yield.

Batches were disjoint (Codex 29, Grok 28) so neither could anchor on the other. Claude's own
verdicts were sealed to `scratchpad/claude_verdicts.md` before the auditors returned.
Agreement was 24/29 with Codex; the disagreements were bucket-boundary
(`SPLIT` vs `CLOSE?`), not factual.

Context that made the triage cheap: the debt register carries only **10 open rows** (DR-29,
DR-39, DR-71, DR-105, DR-113, DR-119, DR-122, DR-131, DR-133, DR-155), so "does this plan
still have a live register link" is a ten-way check.

## Outcome — nine relocated, 48 active

To `closed/`: `dr106_dr144_joint_session`, `dr144_mumps_lifecycle_plan`, `doc_currentness_fixes`,
`handoff_20260831_open_defects`, `weekly_report_20260828`, `mumps_workspace_cap_raise`.
To `deferred/`: `cut_representative_options`, `side_edge_fusing_cut_aware_plan`,
`mumps_symmetric_triangle_extraction`.

**No `cancelled/` was created.** Only one file was a supersession (`mumps_workspace_cap_raise`)
and its own successor plan, `mumps_memory_budget_plan.md:285`, instructs the move to `closed/`.
A third bucket holding a single document costs more in lookup than it saves.

## The four false Status lines

These are the findings that justify the round. Each was verified at its named site.

1. **`nitrogen_eos_completion.md` — comprehensively false** (Grok found it). Status said "six of
   nine functions correct on disk; three outstanding. Not wired; `EoS_Nitrogen` is unreachable
   from `Gas` today", and a re-verification paragraph named three specific syntax defects.
   Tree: the file is **13691 bytes** (Status said 7050), all three bodies are implemented
   (`cl_GM_EoS_Nitrogen.cpp:246,304,338`), the stray `* *` is gone (zero matches), `kk` is
   declared at `:232`, the file **is** in `SOURCES` (`src/physics/gasmodels/CMakeLists.txt:14` —
   the Status asserted the opposite as a warning to future readers), and `Gas` constructs it at
   `cl_Gas.cpp:845`. The work was evidently completed during the gasmodels migration and never
   recorded. Its "reconstruction insurance" framing is now inverted: the bodies in the plan are
   the *older* text and must not be pasted over the source. What is genuinely owed — and never
   existed — is a numerical gate against Span et al. 2000.

2. **`mumps_error_analysis_opt_in.md` — self-contradictory within one file.** Header said "only
   R8b remains, blocked on O3"; R8b is `[x]` landed at `:343` and O3 is RESOLVED at `:463`,
   explicitly "R8b UNBLOCKED". All 16 steps are done; what remains is the G3–G5 gate matrix.
   This matters beyond tidiness: it is the difference between a one-gate `CLOSE?` and six owed
   gates. The definition-of-done line at `:524` is stale the same way and is flagged in place.

3. **`powerlaw_jc_n_field_derivatives.md` — dangling gate.** The β channel was "gated on DR-69's
   signed-angle work"; **DR-69 is struck** (zero occurrences in `debt_register.md`, struck in
   `debt_register_closed.md`). β is unblocked but unscheduled — the one state it should not
   stay in.

4. **`gauge_newton_tangent_plan.md` — records a reversed ruling.** O3 records Christian's
   2026-08-27 decision to default `chi = 1e-4`. That was **reversed to opt-in on 2026-09-01**:
   the Controller's absent-block branch now writes `set_penalty( 0.0, 2 )`
   (`cl_FEM_Controller.cpp:4454-4461`, whose comment states the reversal and its reason). The
   IWG constructor default is still 1e-4 (`cl_IWG_Maxwell.cpp:54`), which is why the Controller
   writes the 0 explicitly. G3's premise flips back with it.

## Cross-reference repair

Nine stale `todo/` paths repaired across the live plans, `README.md`, `debt_register.md` and
`debt_register_closed.md`; four of those were created by this session's own moves and caught by
re-running the check afterwards. Every `README.md` link now resolves. Files **inside** `closed/`
that cite other closed files were left as written — they are historical records.

## Owed from this sweep

- **Seven DR-42/DR-49 rulings** (finished but for one gate, so not ours to close):
  `current_sign_2d_fix`, `dr111_quad4ts_winding`, `nedelec_edge_function_defects`,
  `pid_timestep_controller_plan`, `falsification_tooling`, `gauge_newton_tangent_plan`,
  `make_doc_repair`.
- **One live contradiction:** `rho_lambda_argument_convention` (T-first doc pass) against
  `rho_lambda_b_first_reorder` (Christian's B-first reversal). Executing the doc pass while
  B-first stands would be wasted work in the wrong direction. One must be retired first.
- **`nedelec_edge_function_defects` Status understates its own file** — it says "OPEN — awaiting
  Christian's ruling. No source modified" while R1/R2 are `[x]` fixed by Christian and
  probe-verified. Left for his ruling rather than rewritten, since the correction and the
  close are the same decision.
- **By-catch, outside this sweep's scope:** two live documentation files cite a `todo/` path
  (`src/fem/maxwell/doc/thin_shell_virtual_domains.md:46` and
  `src/mesh/doc/thin_shell_geometry_and_periodicity.md:61`, both pointing at the now-closed
  `periodic_thin_cut_continuity_fix.md`). Per `CLAUDE.md`, a `doc/` or `src/*/doc/` file must
  not cite a todo file at all, so these are policy violations rather than stale links, and the
  fix is to state the fact inline — not to re-point the path. Not done here.


---

# Round 2 — the method above was wrong, and this is how

**Christian rejected round 1's result:** "I think there are far more that are obsolete, such as
`sch04_history_purge.md`, `matvec_csc_omp_stack_overflow.md`, `hex8tb_phase2_fem_wiring.md` and
many more", then added `ngspice_parser` and `thermal_matrices_cleanup_and_newton_plan`.

He was right, and the error was **methodological, not clerical**. Round 1 asked each plan *what
does your text say is left?* A plan is written when work is planned and is rarely updated when the
work lands, so both of its signals — unchecked `- [ ]` boxes and an OPEN/IN PROGRESS Status line —
are biased in one direction: they **overstate** remaining work. Round 1 found the four cases where
the bias ran the other way (Status understating completion) precisely because those are the ones a
text read *can* catch. It was structurally blind to the common case.

The tell was in round 1's own output and I did not read it as one: `matvec_csc_omp_stack_overflow`
has **17 unchecked boxes** and I called it ACTIVE, when the fix is fully landed
(`splinalg.f90:13,48,63,70` gated on `BELFEM_OMP`, `CMakeLists.txt:88` defaults
`USE_BELFEM_OPENMP` OFF) and most of those boxes describe a design that was **rejected and never
built**, kept deliberately as a record. Box-counting inverted the answer.

## Round 2 method

One question per plan: **is the work already in the tree?** Verdicts required a grep result, never
a quotation. Both auditors were briefed with the `matvec_csc` failure as a worked example and told
that a `LIVE` verdict needs a *negative* tree result naming what was searched for. Disjoint batches
of 24 (Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high).

**Disposition rule, now recorded in `todo/README.md`:** a plan stays in `todo/` only while code or
tests are genuinely missing. Once the substance has landed and the sole remnant is a run gate, a
build or a sign-off, it moves to `closed/` with that gate named in its Status — an owed gate is
tracked in the debt register, not by keeping a finished plan on the active list.

## Result: 19 more moved, 57 → 29 active

To `closed/` (17, remnant = gate or sign-off only): `matvec_csc_omp_stack_overflow`,
`nitrogen_eos_completion`, `current_sign_2d_fix`, `dr111_quad4ts_winding`,
`nedelec_edge_function_defects`, `pid_timestep_controller_plan`, `mumps_error_analysis_opt_in`,
`thermal_matrices_cleanup_and_newton_plan`, `ngspice_parser_plan`, `maxwell_kernel_collapse_plan`,
`hex8tb_phase2_fem_wiring`, `falsification_tooling`, `make_doc_repair`, `gauge_newton_tangent_plan`,
`ferro_undulator_mu_and_interface_bugs`, `dr128_spap_low_field_axis`,
`rho_lambda_argument_convention`.

To the new `cancelled/` (2): `arpack_small_end_configuration` — its spectral fold was **replaced**
by shift-invert, and the source says so itself at `cl_FEM_DofMgr_EigenValues.hpp:129-136`
("REPLACED by shift-invert"), with the implementation at `.cpp:1160`. And `sch04_history_purge` —
**Christian's ruling, asked and answered in-session**: the history rewrite is abandoned. Worth
recording that the tree disagreed with the "obsolete" reading here and the ruling settled it:
`share/material/fysc.hdf5` is still tracked in HEAD and `git rev-list --objects --all` still
returns 6 sch04 objects, so R5–R8 were genuinely undone. Obsolete by decision, not by completion —
which is exactly the distinction `cancelled/` exists to hold.

## Plans the tree flatly contradicted

The yield box-counting could not reach. Each verified at the cited site:

| Plan | What it claimed | What the tree holds |
|---|---|---|
| `nedelec_edge_function_defects` | "OPEN … **No source modified.** no numeric probe has run" | Both fixes in `cl_EF_TET4.cpp` / TET10, plus a 22/22 circulation battery `tests/fem/test_EdgeFunctions.cpp` wired into CMake — and its own R1/R2 boxes are `[x]` |
| `dr128_spap_low_field_axis` | framed as unimplemented | Transition at `cl_JcFunction_Database.hpp:440-480`, applied by `eval`/`deval_dB` at `:580-671` |
| `nitrogen_eos_completion` | three unwritten bodies, file absent from `SOURCES`, unreachable from `Gas` | All nine written; `CMakeLists.txt:14`; `cl_Gas.cpp:845`. Grok found this understated even by round 1's correction — `Nitrogen_Caloric` already FD-checks `phir_dd` |
| `thermal_matrices_cleanup_and_newton_plan` | awaiting `T_h_newton` | 77 lines at `mt_thermal_h.cpp:57`, dispatched from `cl_IWG_MaxwellThermal.cpp` |
| `maxwell_kernel_collapse_plan` | "no doc mentions the kernels" | False — Grok found the documentation |
| `ngspice_parser_plan` | Phase 4 `circuitrun` + shared `CircuitSolver` pending | 0 hits under `src/`; design absent and rejected, residue owned by open DR-39 |

## The 29 that stay are LIVE on negative tree evidence

Each names a thing searched for and not found: `belfem_probe_link` (linalg verification),
`tests/fem/test_Normals.cpp` (0 hits), callers of `parmetis_nd` (none outside itself),
`drho_powerlaw_dbeta` (0 hits — the β channel is genuinely missing code, not a stale box),
`belfemConfig.cmake`, Phase-4 min-cost L1, the designed `stall slope` key (0 hits in schema and
`src/`), a PETSc-native conditioning metric (`KSPSetComputeSingularValues`, 0 hits), and the
still-built `src/executables/electricalCircuit.cpp` that DR-39 exists to delete.

## Honest limits of round 2

- **Nothing was compiled or run.** Every verdict is a static tree read. For the 17 moved to
  `closed/`, "the code is present" is verified; "the code is correct" is not, and each owed gate
  is named in its Status line rather than waived.
- **The gate-vs-code boundary is a judgement.** `mumps_error_analysis_opt_in` (16/16 steps, gates
  only) was moved while `mumps_memory_budget_plan` (R12 owed, two days old, actively being worked)
  was kept, though both are "code landed, gates owed". Age and campaign activity broke the tie, and
  that is a defensible call rather than a derived one.
- **Round 1's four in-place Status corrections stand** and were not re-litigated.

---

# Round 3 — purpose-served and relevance (Fable, after Christian's review of the 29)

Christian, on the 29 left by round 2: "most of them are stale", with two calibrating examples.
`ac_loss_postprocessing` should close because AC losses **are** computed and stored in the exodus
— and they are: `Controller::collect_dotQ()` (`cl_FEM_Controller.cpp:2333`) sums the global every
step and the mesh writes it out. Round 2's auditor had called it LIVE because "no in-tree code
integrates `dotQ` to energy" — true, and beside the point. The plan's *purpose* (a reportable
loss) is served; its *refined scope* (energy, decomposition, honesty checks) is a nicety. That
is the standard round 3 applied throughout: **purpose served in the tree → closed.**

The second lens was relevance: a plan can be genuinely undone and still not deserve to be
carried. Christian proposed cancelling those of least relevance; nine went that way, each with a
one-line reason stamped under its title.

**29 → 6.** Dispositions and reasons are in each file's banner and in `todo/README.md`; the
one ruling asked in-session was `rho_lambda_b_first_reorder` — cancelled, keep T-first.

Two by-catch items from the checks:
- `config/globals.cmake:20` still listed `hphirun hphiTrun` in `BELFEM_INSTALL_EXECUTABLES`,
  two executables that are not built. Removed on Christian's approval; `check_doc_claims.py`
  37/37 afterwards.
- `iterate_refactor_plan`'s three small controller bugs (ω clamped after `set_omega` in the
  magnetic path; the divergence counter shared across three adapts; the streak leaking across
  accepted steps) are recorded in its cancellation banner as debt-register candidates. They were
  re-verified live on 2026-08-11 and not re-verified today.

**Doxygen sweep.** Christian asked for the ~245 unapplied corrections to be implemented today.
Four reports had no `_applied.md` at all (`core_containers_comm`, `homology_math`,
`module_docs_fem`, `module_docs_nonfem`) and `materials` had 6 findings applied of ~50. Five
applier agents were dispatched in parallel under `apply_brief.md` (documentation-only edits, the
six closed homology units untouchable, one log line per finding). Their outcome is appended
below when they return.

## Doxygen applier wave — outcome

| Report | applied | already in place | closed-unit | code-side |
|---|---|---|---|---|
| `core_containers_comm` | 0 | 20 | — | — |
| `homology_math` | 31 | 7 | 4 (Gregory) | — |
| `module_docs_fem` | 30 | 6 | — | — |
| `module_docs_nonfem` | 112 | 3 | 1 (Gregory) | — |
| `materials` (remainder) | 31 | — | — | 1 (D15) |

**204 applied, 36 found already in place.** The "~245 unapplied" figure in the plan was
overstated: `core_containers_comm` had been fully applied on 2026-09-02 and only its log was
lost, and a scatter of findings across the other reports had been fixed on disk by concurrent
sessions. New code-side items D15–D18 and the Gregory list are recorded in the plan's Status.

**Verification of the edits, and its limit.** For all 39 touched source files, the working tree
and `HEAD` were compared after stripping C/C++ (and Fortran `!`) comments: identical for every
applier-touched file — the edits are comment-only by construction. `make doc` (R3) has not run;
Christian runs builds. **Reviewed, not verified.**

**Concurrent-session note.** The same comparison flagged six `src/numerics/integration/` files
(the pyramid Gauss tables and `fn_intpoints.cpp`'s tet/pyra order switch) as carrying real code
changes stamped "Regenerated 2026-09-03". None of the five appliers had numerics in its partition
and none of their logs touches those files; the change is a concurrent session's work on the
gauss-table items (D11/D12) and is deliberately left untouched here, per the shared-checkout rule.
