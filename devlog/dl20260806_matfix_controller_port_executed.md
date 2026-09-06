# Matfix Controller Hardening Port — Executed

**Date:** 2026-08-06
**Purpose:** Execution record for `../todo/closed/matfix_controller_port_plan.md` (PLAN v2,
Codex-audited): harvest the Jul 27–28 matfix controller/solver hardening
(`8c52161f`, `7e5df34c`, `2371339c`) onto `sideconnectors`, reconciled with the
ts17/Anderson work that exists only there. Root cause of Gregory's convergence
regression — see `dl20260806_greg_corc_convergence_regression.md`.
**Module:** `src/fem/kernel` (Controller, DofManager, SolverData), `src/sparse`
(SolverWrapper, MUMPS, STRUMPACK), `src/fem/iwg` (Timestep guard)
**Executor:** fresh Fable session (Christian's routing rule for
correctness-critical refactors); plan + Codex audit were prior sessions.
**Working tree only — nothing committed.** matfix untouched (Gregory's reference).

## What landed (R1–R8, all mechanisms H-A…H-N except the H-J second wave)

- **R1 / H-G — solver soft-fail contract.** `cl_SolverWrapper.hpp`,
  `cl_SolverMUMPS.cpp`, `cl_SolverSTRUMPACK.cpp` taken verbatim from matfix
  (branch diff contained only this work): armed wrappers record failures
  instead of `BELFEM_ERROR`; MUMPS relies on rank-uniform error propagation
  (INFO(1) sign) + forces JOB=6 re-analysis on retry; STRUMPACK parallel path
  reduces the verdict with `MPI_Allreduce` (codes are NOT reduced by
  STRUMPACK). DofManager gained `solve_failed()` / `absolute_residual()`;
  SolverData got `clear_failure()` at solve entry and failed-LHS write gates on
  all four update paths (Direct, Newton, Picard, RHS-matrix) — merged
  surgically around the Anderson methods and the Newton `mFieldValues` refresh,
  which are untouched. Controller: constructor + `set_thermal_kernel` arming,
  four `solve_failed()` consumers (coupled magnetic/thermal,
  `iterate_magnetic`, `iterate_thermal` with its Δt₂/20 vs `reset_thermal`
  mirror), `mSolverFailCount` (abort at 8 consecutive) cleared unconditionally
  in `finalize()` when `!mReset`.
- **R2 / H-A + H-F — acceptance clause + moved-baseline detector** in the
  coupled magnetic line search. Absolute-accept clause now bounded to +1.0
  decade of the reference (kills the multi-decade Newton-kick spiral, ts34).
  Moved-baseline: two consecutive flat reject pairs (< 0.05 decade) under
  genuine ω decrease (< 0.6×) accept the trial (quench baseline drift, ts1727).
  Per O1 ruling (recorded lean, no override): the accepted trial's staged
  Anderson pair is **discarded and the window cleared** — its fixed-point
  residual mixes two thermal baselines. Per C2: the accept path restarts the
  magnetic watchdog best-tracker (same idiom as the escalation restart), else
  the watchdog cuts the step the detector just rescued. Every reject still runs
  `anderson_discard + anderson_clear` (C3 preserved).
- **R3 / H-B — damped first Newton entry** at all four promotion-copy sites
  (coupled magnetic/thermal, `iterate_magnetic`, `iterate_thermal`). Per C1:
  `ω_N = clamp( tDamp · ω_P )` — the sideconnectors **copy carry** is kept, the
  matfix min-merge is NOT reintroduced (except the escalated-stall case in R5).
  O2: factor stays 0.5; R9 traces decide any rescale.
- **R4 / H-C + H-K + H-L.** Trust growth ×2 for a first-trial-accepted Newton
  step that halves the residual (coupled path only; `tRanNewton` captured
  before the loop so a same-call escalation cannot misroute the ω store or the
  trust test). Escalation guard `mEpsilon < 1E1` in `try_escalate_to_newton`
  (ts9: escalating a blown-up iterate disables the +10 dB cut via mForceNewton
  → 75-iteration flail). Print unclamp: `residual_string()` (scientific above
  10) replaces the `min(ε, 9.0)` clamp at all four print sites + the circuit
  print. AIMD growth stays ACTIVE under Anderson (O3 stale premise, plan v2).
- **R5 / H-D + H-E.** Stagnation-latch ω branch: promotion-Newton stall resumes
  Picard from **its own** ω at half throttle; escalated-Newton stall keeps the
  conservative min-merge (Codex DQ2). mForceNewton is read BEFORE it is
  cleared. Latch re-arm union: `reset_timestep` re-arms `mFirstFlip/-2`
  (self-contained reset, Codex DQ1), `reset_thermal` re-arms `mFirstFlip2`;
  the drivers' initializers already re-armed on both branches.
- **R6 / H-H + H-N + H-M.** Absolute tolerance defaults flipped 1e-12 → 0.0
  (disabled; raw ‖Ax−b‖ is dimensional, Codex+Grok RQ3). SolverData stores +
  broadcasts `mAbsoluteResidual`; `mEpsilonAbs2` consumed in `run_coupled` /
  `run_thermal`. The MAGNETIC `mAbsoluteEpsilonTarget` stays parsed-but-unused
  by design — documented as such in the theory doc's input table. Thermal
  flat-stall exit: 5 bit-flat coupled iterations (< 0.001 decade) with magnetic
  converged trip `mThermalStalled` → timestep accepted with
  `print_thermal_stall_warning()` (max-T + material `T_max` clamp diagnostic);
  flat evidence reset on freeze (Codex RQ2) and re-armed at
  initialize/reset. `save()` verifies the converged time lies on the
  `save every` grid instead of trusting the sticky flag (off-grid frame bug).
- **R7 / H-I — CN/Galerkin hard-error** in `set_timestepping_method`, both the
  stiffness branch and the no-stiffness aliasing branch (no silent
  `bdf1_nok` alias). Sideconnectors' BDF5 `mHDropped` shift/restore
  (`cc90ed90`) preserved — the matfix-side removal hunks were NOT ported.
- **R8.** Incremental builds green after R1 and after R2–R7.
  `src/fem/kernel/doc/nonlinear_controller_theory.md` updated to the merged
  behavior (damped entry, decade-bounded acceptance, moved-baseline, trust
  growth, half-throttle stagnation resume, soft-fail net, thermal flat-stall
  exit, save-grid verification, opt-in absolute tolerance incl. the dead
  magnetic knob). Test suite: the shared build tree had `USE_TEST=OFF`;
  temporarily reconfigured ON, ran the suite, restored OFF afterwards.
  (CLAUDE.md's `make tests` is a stale target name — the target is `check`.)
  **Port-relevant suites all pass:** containers, linalg, comm, math, sparse,
  mesh, fem (incl. all 5 AndersonMixing tests, run explicitly), ode, physics,
  core, gastables (needs `BELFEM_DATA=<repo>/share` — `data_path()` appends
  `/fluidprop`, and the relative fallback misses ctest's working directory by
  one level). Not port-related: io/kepler/manta don't compile (see below +
  nonfree `typedefs.hpp` include paths), gasmodels has 6 failing Cubic/Methane
  tests — consistent with the staged gas migration (share/ data excluded,
  `f6babfb0`; Methane critical-region tables live only in
  `nonfree/physics/tables`), and the module didn't compile at all before the
  `aResult` fix below, so these tests have never run on this tree.

## Out-of-scope findings (pre-existing, exposed by enabling tests)

- `src/physics/gasmodels/cl_GM_EoS_Cubic.cpp:431` returned undeclared `aResult`
  (function computes `tResult`) — compile error whenever tests are enabled (the
  gasmodels lib is not in the default target). **Fixed in place** (1 line) to
  unblock `make check`.
- `tests/io/test_HDF5.cpp` (`78c1692e` "more tests") calls a 2-arg
  `HDF5::create_group( label, parent )` that exists on neither branch — test_io
  does not compile; the tests were committed against an API that never landed.
  **Fixed later the same session on Christian's call:** the overload was added
  to `cl_HDF5.hpp/.cpp` — the group-navigation state is a stack, so the parent
  handle must be the active group (BELFEM_ERROR-verified; every test use
  satisfies this, the argument makes the intent explicit). test_io compiles and
  all 43 HDF5 tests pass, including the three nested-group tests that had never
  run. Register: DR-50 (closed). D7 of the fusing campaign is DR-51 (closed).

## Kept intact on sideconnectors (verified, not regressed)

Progress watchdogs, `try_escalate_to_newton` semantics + `mNewtonEscalated`,
Anderson stage/commit/discard/flush handshake with growth ACTIVE, thermal
flip-count latch, retry hygiene, Newton `mFieldValues` refresh + Anderson
SolverData methods, `set_thermal_kernel` Anderson-depth forwarding, BDF5
`mHDropped` restore.

## Open / next

- **R9 (Christian):** Greg CORC deck A/B vs matfix (`algorithm : Picard` AND
  `: Newton`, `target iterations : 20`); ts34-class coupled tape trace (Newton
  entry kick ≤ ~3 dB, no ω-floor pinning, no multi-decade accepts). C4 retune
  decision and the O2 damping factor hang on this.
- **R10 (H-J thermal/material tangent wave):** deferred per O5, own session;
  FVM/QUAD4TS pseudoinverse in/out per O6 (reconcile with the dl20260803
  Blaze-defect line, don't double-port).
- O1's discard+clear ruling was executed from the recorded lean; flip to
  commit-on-accept only if R9 traces argue for it.

**Plan:** `../todo/closed/matfix_controller_port_plan.md` (Status updated, R1–R7 ticked).

## Addendum 2 (same session): R9 trace findings — H-F gate + two probes

Christian's first R9 Garber trace (magnetic-only, `mFuseEdges` re-enabled,
oscillations present) delivered two findings:

1. **"Newton does nothing", characterized:** accepted Newton iterates sit
   bit-flat (±0.01 dB) while ω sweeps 0.6→0.07 (ts11 its 4–8 at −57.96 dB,
   ts14 its 21–25 at −41.73 dB) — the correction has (near-)zero effect on the
   residual at any step fraction. From the same iterate, Picard descends
   decades and both timesteps reach the −156.5 dB machine floor, so the fused
   constraint system is CONSISTENT and the ported controller routes around the
   dead tangent correctly (stall guard → half-throttle Picard resume). The
   defect is in the Newton tangent's interaction with the fused configuration.
2. **Moved-baseline detector mis-fire (controller bug, FIXED):** the H-F accept
   fired on the magnetic-only run (ts11 it32: accepted −5.26 dB from a
   −49.85 dB reference after three flat ω halvings; ts14 it9 similar) although
   its premise — the staggered thermal update moving the baseline — cannot
   hold without a thermal kernel. Gate added on Christian's instruction:
   `mKernel2 != nullptr && ! mThermalFrozen` in front of the accept branch
   (`iterate_coupled`); matfix's motivating coupled-quench case (ts1727) is
   preserved. Theory doc §3 updated. Plan amendment A1.

**Two opt-in probes installed** (env-gated, DR-25 hygiene — strip later).
UPDATE later same day: `BELFEM_PROBE_NEWTON_DX` and the `#bearing` setup print
were REMOVED again after the diagnosis closed (near-null φ gauge mode of the
point bearing — not a bug; see `src/fem/kernel/doc/bearing_gauge_eigenmode.md`;
remedy = pure Picard for net-current decks). `BELFEM_PROBE_FUSED_ROWS` stays
for the oscillation campaign, upgraded to v2 ( per-source
`basis-id:dof-type:weight` triplets, decidable sign coherence and λ content ):

- `BELFEM_PROBE_NEWTON_DX` — prints `|dx|`, `|r|`, ω per Newton solve
  (`cl_FEM_DofMgr_SolverData.cpp`, Newton branch). Splits the paralysis:
  |dx| ≈ 0 = rigid tangent (assembly/scaling on fused dofs); |dx| large with
  flat residual = Newton direction in the lagged operator's near-nullspace.
- `BELFEM_PROBE_FUSED_ROWS` — dumps every LINE2 edge-on-node constraint row at
  dof-manager conversion (`cl_FEM_DofMgr_DofData.cpp`) to
  `fused_edge_rows_rank<r>.txt`: edge id, own node order (the ±1 sign
  convention), original ids, periodic/duplicate flags, hang-source node ids,
  final weights, midpoint coordinates. Setup-time only — usable from a run
  killed after the first residual line, with all BCs (periodicity, cuts)
  exactly as in production; sign coherence incl. the periodic map is checked
  offline and correlates with blob positions via the coordinates.

## Addendum (same session): fusing-campaign defect D7 fixed

On Christian's go-ahead after reviewing the latest Christian↔Gregory exchange
(assessment against `tmp/ai_exchange/side_edge_fusing_handoff.md`, two verified
jury rounds): `ThinShellFactory::connect_side_nodes` now indexes the layer node
container with `tOrg->original()->index()` — the D1 pattern, matching every
other `Layer::Nodes` consumer. Verification sharpened the finding: cut-duplicate
side-curve stations carried the `gNoIndex` sentinel (all master nodes reset in
`create_node_container`, only sorted originals re-indexed), so the old code was
an out-of-bounds access, not merely a wrong-node hang. Dormant while
`mFuseEdges = false`; unblocks the no-cut fusing experiment (gate 2). Open
before any cut/periodic fused run: the §6 orientation gate and the
`reset_source_container()` clobbering audit. Ticked as D7 in
`todo/hex8tb_phase2_fem_wiring.md` (Status refreshed).

---

**Addendum 2026-08-25.** The "kepler/manta don't compile" line above described the
symptom correctly but not the cause, and ctest never reported them as *failed* —
the executables did not exist, so the status was "Not Run". Cause:
`config/scripts/Add_Test.cmake` took its `core`/`comm`/`containers`/… include paths
from `${SSF_SRC_DIR}`, a pre-BELFEM name never defined in this repository, and the
open-source `tests/` tree only compiled because the root-level `banner` executable
leaks `Add_Executable.cmake`'s `include_directories` into it. `nonfree/tests` is
added before that point and saw only `-I/core`. Fixed in `Add_Test.cmake`; see
`devlog/dl20260825_add_test_include_paths.md`.
