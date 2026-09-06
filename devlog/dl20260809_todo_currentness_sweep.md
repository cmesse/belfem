# Todo Directory Currentness Sweep — Technical-Debt Reassessment

**Date:** 2026-08-09 (overnight; follows `dl20260809_sideconnector_viz_regression_round.md`)
**Purpose:** Re-read every active `./todo/*.md` not authored that day, check each claim against
the tree / devlogs / `git log`, and correct what had gone stale — so the debt picture going
into the September cut is real rather than inherited.
**Module:** meta (`./todo/`)
**AIs involved:** Claude (sweep + corrections). No Codex/Grok round — this was verification
against the tree, not adjudication of a disputed claim.
**Claude Confidence:** high for every item marked "verified in tree" below (each rests on a
direct file read at a named anchor); medium where a run gate is involved, since no runs were
performed.
**Verification:** source-trace tier throughout. Nothing was built or executed; every "fixed"
verdict is a code read at a cited `file:line`, and every "still open" verdict is a failed
search plus a positive check that the old anchor moved rather than the feature landing.
**Scope:** documentation only. No source file was modified.

## Summary

33 files reviewed. **Five plans had finished without being closed** and were moved to
`closed/`; **nine had status lines that actively misdescribed the tree**; four analysis
documents kept their findings but had every line anchor drift out from under them. The debt
register got its first correction pass since it was seeded: 13 rows verified, 7 closed, one
downgraded, two added.

The pattern behind the staleness is worth naming: **the plans that went stale are the ones
that were overtaken sideways.** Not one of them was abandoned. `maxwell_kernel_collapse_plan`
went stale because Christian leapfrogged R6–R10 in a single 2026-07-21 rewrite, so the
per-family boxes were never ticked even though the whole family tree is gone.
`thermal_matrices_cleanup_and_newton_plan` went stale because its cleanup half was delivered
by a *different* plan. `2d_thinshell_todo` went stale because Gregory's alternation report
pulled E5 forward out of sprint order. `anderson_picard_acceleration_plan` went stale because
it was committed the same night it was written and nobody went back to change the word
"Uncommitted". Sequential plans survive; leapfrogged ones rot.

## Closed (moved to `todo/closed/`)

| file | why |
|---|---|
| `cross_review_tooling.md` | R1–R11 done; everything committed in `bd73583b`. The one `[◐]` said METHODOLOGY.md was "NOT committed" — it is. |
| `cut_pocket_removal_rules.md` | `Cohomology::remove_cut_pockets` in tree since 2026-07-16; Rules 1–7 and O1–O4 ticked. Tier-B enablement is a validation gate, not work in flight. |
| `periodic_cap_cut_emission.md` | Option (d) shipped; (a)/(b)/(c) struck. The file's own recommendation (try (a) first) is recorded as overtaken by the probe rather than quietly deleted. |
| `handoff_double_corc_d1_session.md` | Every premise consumed: D1 fixed, the "nothing is committed" inventory committed (`ce6a0e8b`), all seam probes stripped (zero hits). Residual items re-homed. |
| `handoff_thermal_iwg_collapse_session.md` | All of T1–T7 and S1–S7 already ticked; only the status line still said HANDOFF. |

## Corrected against the tree

Verified fixes that no todo file had recorded:

- **D13 (kernel collapse) dissolved.** Conductive layers are FEM Blocks on the magnetic side
  too, and `Calculator::allocate` builds `MaxwellData` for every non-Air block present in the
  peer dof manager (`cl_FEM_Calculator.cpp:1111-1137`). The "sideset MaxwellData design" R9
  was waiting on was never needed. → DR-04 closed.
- **The thermal Newton tangent is written and wired.** `T_h_newton` (`mt_thermal_h.cpp:54-127`)
  fills all four blocks — including the unsymmetric (Bᵀ∇T)⊗(dλ/dT·N) conductivity block in a
  preallocated `Btg` workspace — and `IWG_MaxwellThermal::link_to_group` selects it on
  `algorithm()`. O1 was decided as option (c) (analytic `Material` API), not the
  finite-difference option that was recommended. → D8 closed, DR-03 closed.
- **`T_phi`'s placeholder properties are gone** — it reads density/cp/lambda from the group
  material (`mt_thermal_phi.cpp:14-40`). → DR-05 closed.
- **The Anderson bundle is committed** (`6b2a2b98`, 2026-07-30 23:31), with its five-case
  regression test. → DR-10 closed. `devlog/campaigns/controller_anderson.md` still says
  UNCOMMITTED and needs the same correction.
- **`mConnectorWidth` is no longer hardcoded** — derived per shell from
  `edge_coating_width()` (`cl_ThinShellFactory.cpp:306`). → DR-20 closed.
- **DR-13 downgraded**: `Vertex` facet counters are `uint16_t` now, so the wrap moved
  256 → 65536. Residual: `allocate_facet_container` still narrows a `uint` into it unguarded.

Two new register rows, both lifted from findings that existed only in devlogs:

- **DR-52 (P0, already fixed)** — the residual reported under Anderson mixing was the linear
  solver's own roundoff, by algebraic identity on the bootstrap iterate. **This row is also
  the sweep's own correction story, and it is the methodological point worth keeping:** it was
  first written up as an open P0 straight from the jury devlog's headline, then re-checked
  against the tree and found already fixed — the same devlog's *addendum 1* records Christian
  executing the fix that same day. Verified in tree: the `mFieldValues` refresh is gone with
  the reasoning inline (`cl_FEM_DofMgr_SolverData.cpp:2776-2787`), `fixed_point_residual()`
  exists as a diagnostic, and the defaults are rolled back to BDF1 + Anderson depth 0/0, i.e.
  opt-**in** again. The fix **reverses D5 of the Anderson campaign** and is **uncommitted and
  not compiled**. Lesson: a devlog headline and its same-day addendum can tell opposite
  stories; the addendum is the one that matches the code.
- **DR-53** — wall-side fusing versus the decoupled viz sheets (today's Δt collapse). Recorded
  with the observation that it lands on **the same 90/348 rim stacks** the 2026-08-06
  authority census identified, so "ship = free rims" now rests on two independent failure
  modes rather than one.

## Escalations and one withdrawal

- **DR-33 — escalation drafted, then withdrawn.** The 2026-08-07 opt-out change made BDF5 the
  default, which would have put every key-less deck on the multi-step history path the
  never-written V1 test guards. It was rolled back to BDF1 the same day
  (`cl_FEM_Controller.hpp:189`), so the row stays P2. The B7 half is done and was ticked; V1
  remains the only guard on a path decks must now select explicitly.
- **DR-23 sharpened.** `clean_spfa`, `rectify_greedy_sweeps`, `fire_node_coboundary` and
  `remove_cut_pockets` are committed production code with **zero** checked-in coverage, and
  the only test harness lives in a session scratchpad.

## Refuted

`nonlinear_iteration_strategy_near_quench.md` (2026-02-03) recommended **against** line search
and Anderson as tools for material nonlinearity. Both shipped since and are load-bearing —
the backtracking line search in `iterate_coupled`, and Anderson(m) mixing (opt-in again after the 2026-08-07 rollback, but implemented, tested and wired). The file now
carries that correction explicitly so the advice is not cited as settled. Its Δt-control half
belongs to `pid_timestep_controller_plan.md`; the genuinely unclaimed idea left in it is the
quench-margin trigger.

## Missed on the first pass, corrected on Christian's question

`maxwell_kernel_collapse_plan.md` R12 was first written up here as having one thing left (the
§4.1 gate). It has two: the gate, **and the `src/fem/maxwell/doc/` kernel architecture note,
which was never written.** No file under `src/fem/maxwell/doc/` or `src/fem/kernel/doc/`
mentions `h_picard` / `h_newton_mu0` / `h_newton_mu` / `T_h_picard` / `T_h_newton` — the only
prose describing the collapsed kernel set lives in devlogs, which do not count as
documentation. A reader of `src/` has no account of why 38 variants became five, nor of the
§6.0 two-way algorithm pick that replaced the material tree. Worth noting *what* the sweep
missed: an undone sub-item inside a box already marked `[◐]`. Partial boxes hide residue —
checking the box state is not the same as reading the sub-items under it.

## Stalled, not failing

`src/fvm` has had no code change since 2026-07-09 and is excluded from the build
(`add_subdirectory( fvm )` commented out, `src/CMakeLists.txt:12`; the `fvmtest` executable
block likewise). Both FVM files now say so. This is a scheduling fact — DR-41 already assumes
FVM is out of 1.0 — but it deserved to be visible in the plan rather than inferred from commit
dates.

## The single highest-leverage next action

Three plans now gate on **one uncommitted, uncompiled working tree**: the Anderson residual
fix (DR-52), the Picard line-search retirement, and the PID timestep controller all live
there together. Anderson R8, matfix R9 and PID R5 each need it to build and run before their
A/B numbers mean anything. Nothing else in `./todo/` unblocks as much per unit of effort.

## Open questions routed to Christian

1. **FVM in or out of 1.0** — a month uncompiled; the longer it sits, the more of the PENTA6
   extrusion phase will have to be re-derived against a moved `cl_ThinShellFactory`.
2. **`ngspice_parser_plan.md`** — a feature proposal with no consumer pressure and no bearing
   on the cut. Candidate for `deferred/`.
3. **`maxwell_postprocessor_gpu_acceleration.md`** (deferred) — its own stated precondition
   ("implement element-field caching first, GPU may then be unnecessary") is now met. Keep or
   drop.
4. **The `[seeded — confirm]` pass** on `debt_register.md` and the six `devlog/campaigns/`
   pages. This sweep verified 13 register rows; the rest are still mechanical lifts, and at
   least one campaign page is now factually wrong.

## Hygiene finds (small, unrelated to any plan)

- `scripts/fluidprop/__pycache__/*.pyc` are git-tracked — remove before the open-source cut.
- `tests/physics/gasmodels` has **five** `DISABLED_` tests, not the four the migration plan
  records.
- `Hex8TbUnitCirculation` is still marked DISABLED although its blocker (DR-47, the missing
  HEX8TB Lagrange factory case) landed on 2026-08-08 — a free test re-enable.

## Files touched

33 under `./todo/` (5 moved to `closed/`, 27 edited, `README.md` rewritten with a sweep
triage table and entries added for the four active files that had never been indexed:
`debt_register.md`, `2d_thinshell_todo.md`, `2d_thinshell_gap_analysis.md`,
`iterate_refactor_plan.md`). No source modified.

## Addendum (same night): three-AI triage round + first fixes

**Thread:** `tmp/ai_exchange/triage_debt_register.md` (pre-registration frozen before
dispatch; blind Codex + Grok legs on an identical brief via the `ask_*.sh` wrappers —
deviation from `/cross-review` noted openly: its hard-coded defect-review priming does not
fit a triage, so the wrappers were driven directly with the protocol shape kept intact).

**Verification yield:** Grok 13/14 citations confirmed, 1 refuted (its DR-32 "latent qold
bug" — `set_timestepping_method` rebuilds the tables at `cl_IWG_Timestep.cpp:222-248`);
Codex all load-bearing citations confirmed, 1 circular (DR-01 = HIGH citing only the
register row back at itself). Each auditor caught something the other missed, and Codex
caught an error in this very sweep (the BDF5-default claim left in
`restart_circuit_verification.md` after the same-day rollback — corrected).

**Register motion:** DR-01 closed as stale (sign-split moot per the R12 note — Grok's
find), DR-32 closed as already-fixed (Codex's find), DR-33's CN/Galerkin clause struck
(hard-error before any tangent), DR-19 and DR-49 stale text refreshed (DR-49's spline
test gate is doubly wrong: production `Spline` uses SuperLU unconditionally,
`cl_Spline.cpp:427-429`).

**Fixed in source (Christian's explicit go-ahead), all pending his build:**
- **DR-24** (Christian): `mPrescribedCurrents` fully removed — member, accessor, read loop.
  Residual: the `set_currents` prefix-order assert is still owed.
- **DR-14** (Claude): `case 1` identity in the edge overload of `to_master_orientation`
  (`fn_to_master_orientation.cpp`); the ±1 sign correctly stays with the caller's
  node-index comparison. Default-branch format string repaired en route.
- **DR-44** (Claude): `EntityType::CELL` fall-through to the ELEMENT body in
  `Distributor::select_sources`, matching the pairing already used at `:596-597`.
- **DR-12** (Claude, after Christian's ruling that a .bfm is correctly oriented by
  construction): MeshChecker ctor now reports its flip count (warning if nonzero — the
  counter existed with zero consumers); `BfmFile::load()` is the single authority for the
  trust flag, which also fixes the direct-`.bfm` deck branch that never set it (latent
  edges-exist abort). Factory's manual flag set removed. NOT done: moving the checker to
  `src/mesh` (blocked by `Pipette`'s home in `fem/kernel`). Run gate: a CW-surface deck.

**Consensus triage** (proposal, not a vote result — full table in the thread): HIGH =
DR-52, DR-02, DR-23, DR-21, DR-22, DR-42, DR-34 + the run-gate block DR-06/08/09/15/38
(both auditors ranked run gates HIGH; Claude's pre-registration had them MEDIUM-as-debt
and defers). DR-18/19 conditional on the side-connector scope ruling. Structural answer,
3/3 independent: one ranking plus a `kind` column (DEFECT / VERIFICATION / DECISION) —
"the kinds need different verbs, not different severities."

## Addendum 2 (overnight, autonomous): DR-02 Gate A attempted, abandoned, and what it turned up

**Mandate:** Christian offered 2 of his cores and handed over the session, asking me to
build the DR-02 Gate A harness. Late in the night he redirected: work only from the current
`sideconnectors` commit. Gate A — matrix-level equivalence against the pre-collapse
kernels — was therefore **dropped**, because an A/B needs the legacy kernels and those exist
only in history. What DR-02 gets instead is physics-level verification (Gate B).

**Operational note:** all work used separate build trees (`belfem-build-gateA`) at
`nice -n 19`, `-j2`. `cmake-build-debug` and the running sidecoatings job were never touched —
a concurrent `make` in a shared tree produces partial-`.o` failures that masquerade as code
bugs.

### Delivered

- **The working tree COMPILES CLEAN** — `hphirun` links, 0 errors, 0 warnings (Blaze, Debug,
  MPI). That covers the whole uncommitted DR-52 bundle (Anderson residual fix, Picard
  line-search retirement, PID controller), Christian's DR-24 removal, and DR-12/14/44. This
  retires the compile half of DR-52's gate, which three plans were waiting on → DR-55.
- **Gate B, helix at HEAD: serial vs 2 MPI ranks give IDENTICAL residual sequences**
  (1.000000 / −156.54 dB; 0.001037 / −29.84 dB; 0.001498 / −28.25 dB), 2-rank reaching
  t = 5.8 s. Parallel consistency verified for the bulk-conductor + periodic families.
- **A reusable comparator** (`compare_system.py`) plus an env-gated assembly dump
  (`BELFEM_DUMP_SYSTEM`, after `compute_jacobian_and_rhs` and before the BC patch/solve, at
  both the coupled and magnetic sites). Self-validated in both directions; its self-test
  caught a genuine bug in itself (zero-size `LHS` crashed the reduction). Kept for whenever a
  two-commit comparison is wanted.

### Found (the more valuable half)

- **DR-57 (P1, public-facing, NEW):** a run on ≥2 ranks **aborts if a metal rho database must
  be built from scratch** — `populate_rho_database_parallel()` partitions its work graph
  through the `Mesh::partition` overload defaulting `aForceContinuousPartitions = true`, the
  graph is non-contiguous, METIS returns −4, and the abort lands in `create_materials()`
  *before the mesh is read*. Proved both ways: helix/2 ranks aborts without a database, and
  the identical run proceeds to t = 5.8 s once the serially-built `Copper_RRR50.hdf5` is
  present. Invisible in practice because every working run dir already carries its cached
  databases — so it hits first-time users. **Must not be fixed by changing those overload
  defaults; the `false` siblings are deliberate (Christian).**
- **DR-54 (P1, public-facing):** `examples/corc` cannot run — `corc.msh` is gmsh **2.2**
  (no `$Entities`, hence no vertices) while the deck declares `periodic` in vertex IDs;
  it aborts in `set_master_plane`. `examples/helix` ships no mesh at all. `examples/` is the
  first thing a new user runs.
- **DR-02 baseline was wrong:** the plan named `ce6a0e8b` (07-17), but the R6 alloy flip
  landed at `abb71c2e` (07-15). Dispatch kernel-assignment counts: `93353369` = 28 (fully
  legacy), `abb71c2e`/`ce6a0e8b` = 26, HEAD = 6. Recorded for any future comparison.
- **DR-56 (historical, closed):** PARDISO could not compile at the 2026-07-14 baseline —
  `int_t` declared inside the first function instead of at module scope, no `use omp_lib`,
  and a `real*`/`const int_t*` signature mismatch. All fixed in `6b2a2b98` (07-30);
  **verified working at HEAD**. An earlier version of that row wrongly said `USE_PARDISO` is
  ON by default — it is OFF at both commits.

### Judgment calls, for the record

- Reordered Christian's request: compiled the current tree *before* the baseline archaeology,
  because it was half the work and gated three plans.
- Patched the baseline's Fortran twice rather than disabling PARDISO, to keep both arms on
  identical flags; fell back to disabling only after a third independent breakage, having
  first checked that `SolverType` is an unconditional enum so no value could shift.
- Both moot once the baseline arm was dropped — the worktree and its build dir are removed.

### Honest accounting

Roughly three-quarters of the night went to build archaeology rather than the comparison, and
the comparison was ultimately cancelled. The findings that fell out of the detour (DR-57 above
all) are worth more than the Gate A verdict would have been — but that was luck, not planning.
The generalisable lesson is the one DR-56 and DR-57 share: **configurations nobody exercises
rot silently**, and both were found only by doing something slightly unusual — building an old
commit, and running a deck in a directory without cached artifacts.

### Addendum 2b: the §4.1 deck inventory — 1 of 5 cases actually runs

Attempting thermal `T_h` coverage after the helix result turned up the real obstacle to any
future §4.1 gate. Tested at HEAD:

| §4.1 case | verdict at HEAD |
|---|---|
| helix | **runs** serially — but ships no mesh (undocumented `gmsh -3 helix.geo`), and needs a pre-built rho database on ≥2 ranks (DR-57) |
| corc | **dead** — gmsh-2.2 mesh has no vertices, deck declares `periodic` in vertex IDs (DR-54) |
| Tape_Quench/BuiltinMat | **dead since 2026-04-08** — its custom material is labelled `buffer`, reserved by `8eae5f00` ("add MgO"); also carries a pre-format-change rho cache (DR-58) |
| Tape_Quench/CustomMat | untested; under ephemeral `tmp/` |
| Validation set | untested; under ephemeral `tmp/` |

`maxwell_kernel_collapse_plan.md` §4.1 records these as "all present in-tree, checked
2026-07-13". They are present. Three of them do not run, one of them for four months.
**The critical path for DR-02 is repairing the decks, not scheduling the run** — and that
should be settled before anyone budgets time for the gate.

Three of tonight's findings (DR-54, DR-57, DR-58) share one shape: **an artifact that stopped
matching the code, with nothing checking.** A gmsh-2.2 mesh against a vertex-ID deck; a rho
cache with no version stamp; a deck label that became reserved. None is hard to fix; none
would have surfaced without running something in a directory that lacked the usual cached
state. That is the argument for a from-scratch smoke over the examples before the 1.0 cut.

## Addendum 3 (2026-08-10, Fable session takeover): DR-57 fixed three layers deep

Christian handed the session to Fable mid-campaign ("Opus is struggling") and later granted
implement-then-jury autonomy ("I trust you, and we can always call the jury after you are
done"). Full record: `todo/example_deck_and_material_db_repair.md` §7; thread
`tmp/ai_exchange/review_example_deck_repair.md` (plan round + post-implementation round).

**The fix, in one sentence:** the material rho-database cache-miss path now runs the former
"serial" build on EVERY rank in lockstep — master owns all nodes and evaluates, workers hold
empty tensor-mesh copies and receive the projected table inside the `Database` ctor's
existing collective pair — because the layer beneath (`Projector::project`, the Database
ctor, `populate_tensor_mesh`) turned out to be symmetric by design all along; the entire
`populate_rho_database_parallel` apparatus (with its Distributor crash, stuck gather
counter, and dense-packing mismatch) was never necessary and is retired unreachable.

**The three attempts are the story:** (1) un-forcing METIS contiguity was necessary but only
exposed the Distributor dying on empty worker sets; (2) master-only-build + all-ranks-load
DEADLOCKED — the "serial" builder ends in a collective ctor, and the pre-flight
"no hidden collectives" check had read function bodies but not callees (the same
shallow-check failure mode, one level deeper, that the plan jury had just caught in the R1
rationale); (3) the lockstep design, enabled by one semantic guard —
`TensorMeshFactory::create_bsplines` early-returns on an element-less mesh, since its loops
iterate the config-sized grid rather than the container.

**Acceptance:** fresh dir + 2 ranks + no cache = the exact first-contact scenario — database
built in-run (0.7 s projection), residuals digit-identical to serial, zero aborts. R2
verified against the genuine old-format cache (warn + rebuild). Determinism characterized
along the way: single-threaded database builds are bit-identical; multithreaded builds carry
~1e-7 relative MKL/MUMPS solver fuzz that predates all of this (serial-vs-serial with the
old binary shows the same scale) — the plan's "bit-identical h5diff" DoD was unachievable
from the day it was written, and is amended to the single-thread statement.

Post-implementation jury dispatched on the focused diff (`tmp/ai_exchange/matdb_fix.diff`)
with the blast-radius question named as the primary target: the new empty-mesh guard fires
for every tensor mesh on a non-master rank, not just this path.
