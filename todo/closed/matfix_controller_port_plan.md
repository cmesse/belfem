# Matfix Controller Hardening → Sideconnectors Port

**Date:** 2026-08-06
**Purpose:** Harvest the Jul 27–28 controller/solver hardening that exists only on
`matfix` (commits `8c52161f`, `7e5df34c`, `2371339c`) onto `sideconnectors`, reconciled
with the ts17/Anderson work that exists only there. Root cause of Gregory's convergence
regression: the two branches carry divergent controller fix lines that never merged
(devlog `dl20260806_greg_corc_convergence_regression.md`; threads
`../../tmp/ai_exchange/review_greg_corc_regression.md`).
**Module:** `../../src/fem/kernel` (Controller, DofManager, SolverData), `../../src/sparse`
(SolverWrapper, MUMPS, STRUMPACK), `../../src/fem/iwg` (Timestep guard); second wave:
`../../src/fem/thermal`, `../../src/fem/kernel/cl_FEM_Calculator.hpp`
**AIs involved:** Claude (archaeology + plan), Codex (plan audit + prose), executor =
fresh Fable session (Christian's routing rule for correctness-critical refactors)
**Status:** EXECUTED 2026-08-06 (Fable session, Christian's "execute") — R1–R7 landed on
the sideconnectors working tree, R8 build green + theory doc updated; test suite run via
`make check` (USE_TEST temporarily enabled, restored to OFF afterwards — the tree had it
disabled; note: CLAUDE.md's `make tests` is a stale target name). R9 traces = Christian.
R10 deferred per O5. Devlog: `dl20260806_matfix_controller_port_executed.md`.

**Re-verified 2026-08-09 (currentness sweep): accurate; R9/R10/O5/O6 still open.** Spot
checks confirm the port is in the working tree — `iterate_magnetic` now carries the
Picard↔Newton handoff, the damped first Newton entry, the watchdog re-anchor and the
solver soft-fail → `reset_timestep` path (`cl_FEM_Controller.cpp:1327-1400`), and
`impose_voltage_bcs()` is a shared member (`:425`).
**Caveat that affects when R9 can run, not how:** the greg3 jury
(`../../devlog/dl20260807_greg3_false_convergence_jury.md`) showed the residual reported under
Anderson was the linear solver's roundoff. That was **fixed the same day** — the refresh is
removed, ε is the pre-update force residual for both paths, and Anderson is opt-**in** again
(default depth 0). So R9's ε is meaningful again. **Commit state corrected 2026-08-10:** the fix is
**committed** (`4f2c11cd`), as are the Picard line-search retirement and the PID timestep
controller (`1d6ef305`) — all three change what an R9 trace measures, so run R9 at
`bc578b5e` or later and record which commit the trace came from. See `anderson_picard_acceleration_plan.md` §4.4 and `../debt_register.md`
DR-52.
Follow-on note for R9's ω observations: the unclamped-ω ordering bug in `iterate_magnetic`
(`:1376` before `:1377`) is still live and is tracked in `../iterate_refactor_plan.md` — worth
knowing before attributing any ω excursion to the port.

> **Scope guards:**
> - Port + reconcile ONLY. No new controller features, no constant retuning beyond what
>   the merge forces (R9 gate decides if a retune session is needed — separately).
> - Anderson mixing semantics (`fn_FEM_anderson_mixing.hpp`, SolverData handshake) are
>   NOT touched; only the line-search restructure must keep its hooks correct.
> - Side-edge fusing / ThinShellFactory is a separate campaign — out of scope here.
> - Target branch: `sideconnectors`. matfix stays untouched (Gregory's reference).

---

## 1. Harvest inventory (by mechanism, not by commit)

Citations: `matfix:` = `git show matfix:src/fem/kernel/cl_FEM_Controller.cpp`;
bare `:NNN` = sideconnectors working tree.

| ID | Mechanism | matfix source | sideconnectors state |
|----|-----------|---------------|----------------------|
| H-A | Acceptance clause A/R0: absolute accept only within 1.0 decade of reference | `matfix:686-696` | hole at `:750` (`\|\| tLogEpsilon < 0.8` unbounded) |
| H-B | Damped first Newton entry ×0.5 (B/R1, `tDamp = mFirstFlip ? 0.5 : 1.0`) | `matfix:584-600` | full-ω entry `:683-687` |
| H-C | Trust growth ×2 (C/R2: `tRanNewton && tBacktracks==0 && ε<0.5ε₀`) | `matfix:960-980` | ×1.3 crawl only |
| H-D | Stagnation-latch ω branch (R3: promotion-stall → Picard's own ω; escalated-stall → min-merge) + half-throttle resume H2 | `matfix:~500-530` (`:524`) | unconditional copy `:623` |
| H-E | Latch re-arm H1 (`reset_timestep`/`reset_thermal` re-arm first-flip latches) | matfix reset fns | partial overlap with `initialize_timestep` hygiene `:112-138` — port the UNION |
| H-F | Moved-baseline detector (2× consecutive flat reject pairs, 0.05 decade, genuine ω decrease → accept; quench limit-cycle fix) | `matfix:640-744` | absent |
| H-G | Solver soft-fail contract: `set_soft_fail` wiring, `solve_failed()` consumers, wrapper soft-return, MUMPS rank-uniform error path, STRUMPACK MPI_Allreduce verdict, failed-LHS write gates, `mSolverFailCount` finalize clear, late-attached thermal arming | `matfix:59-69, :673-680, :847-855, :1073-1081, :1198-1213, :1588-1595, :2208-2212`; `2371339c` in `cl_SolverWrapper.hpp`, `cl_SolverMUMPS.cpp`, `cl_SolverSTRUMPACK.cpp`, `cl_FEM_DofManager.*`, `cl_FEM_DofMgr_SolverData.*` | absent — solver failure = hard BELFEM_ERROR abort. Merge SURGICALLY: keep Anderson SolverData methods (`cl_FEM_DofMgr_SolverData.cpp:2607-2741`), the Newton `mFieldValues` refresh (`:2190-2200`), and `set_thermal_kernel`'s Anderson-depth forwarding (`cl_FEM_Controller.cpp:2201-2206`) |
| H-H | Absolute THERMAL residual semantics: defaults 0.0 (disabled — raw ‖Ax−b‖ is dimensional, Codex+Grok RQ3), raw residual stored/broadcast by SolverData, `mEpsilonAbs2` vs `mAbsoluteEpsilonTarget2` consumed in `run_coupled`/`run_thermal`. The MAGNETIC `mAbsoluteEpsilonTarget` is parsed-but-unused on BOTH branches — do not claim consumption unless a new behavior is deliberately designed | `matfix:cl_FEM_Controller.hpp:80-88,157-160`; `matfix:859, :1218, :1811-1814, :1838-1840`; SolverData absolute-residual plumbing rides with H-G | parses both keys (`:1910-1913`, `:2011-2013`), consumes neither |
| H-I | CN/Galerkin hard-error in `set_timestepping_method` (both branches incl. no-stiffness aliasing) | `8c52161f` in `cl_IWG_Timestep.cpp` | verify current state (36-line branch delta); port ONLY the guard hunks — preserve sideconnectors' BDF5 `mHDropped` restore |
| H-J (2nd wave) | D11–D18/N3 thermal-tangent wave: `compute_drhodT` clamp (D11), `T_phi_*` on MaxwellData (D12/D13), ThermalFactory Air removal (D15), table-derivative windows (D17), Bloch-Grüneisen derivative wiring (D18) — INCLUDING the required `src/physics/materials/*` APIs (`drho_i_dT`, `ddebyedT`, BG wiring, window clamps: `matfix:cl_Material.hpp:796-906,1978-2034`, `matfix:cl_Material_Metal.cpp:288-306`) | `8c52161f` in `cl_FEM_Calculator.*`, `src/fem/thermal/*`, `src/physics/materials/*` | partially propagated in Calculator; thermal/material pieces absent. Harvest by DIFF not commit. FVM/QUAD4TS pseudoinverse fixes: decide in/out explicitly (O6) |
| H-K | Escalation guard: never hand a blown-up Picard iterate to Newton (`mEpsilon < 1E1`; ts9 trace: 75 iterations flailing at +11 dB with the +10 dB cut disabled by mForceNewton) | `matfix:445-456` | absent at `:507-510`; preserve sideconnectors' Anderson clear + watchdog reset in the same function |
| H-L | Residual print unclamp: `residual_string()` — scientific notation above 10 instead of the 9.000000 clamp | `matfix:422-430, :1640-1695` | absent (`:1674-1723`) |
| H-M | Save-grid verification: `save()` verifies the converged time lies on the `save every` grid | `matfix:2252-2272` | absent (`:2245-2248`) |
| H-N | Thermal flat-stall exit: `mThermalFlatCount`/`mThermalStalled`, absolute-residual re-arm, warning box | `matfix:cl_FEM_Controller.hpp:162-174,356-360`; `matfix:229-231, :874-900, :1381-1384, :1699-1757` | absent |

Excluded (inspected, not part of this port): `src/executables/hphirun.cpp:56-58` /
`hphiTrun.cpp:47-49` on matfix set Nitsche `gRhoMin`/`gRhoMax` only — a different line of
work, not soft-fail plumbing.

**Keep intact on sideconnectors (do NOT regress):** progress watchdogs (`:546-593`),
`try_escalate_to_newton` semantics + `mNewtonEscalated`, Anderson
stage/commit/discard/flush handshake with growth ACTIVE (O2 amendment, `:992-998`),
thermal flip-count latch (`:868-871`), retry hygiene (`:112-138`), Newton `mFieldValues`
refresh (`cl_FEM_DofMgr_SolverData.cpp:2190-2200`) + Anderson SolverData methods
(`:2607-2741`), `set_thermal_kernel` Anderson-depth forwarding
(`cl_FEM_Controller.cpp:2201-2206`), BDF5 `mHDropped` restore (`cc90ed90`).

## 2. Known merge conflicts (resolve as specified, else escalate to O-item)

| # | Conflict | Resolution |
|---|----------|------------|
| C1 | ω-carry semantics: matfix damped **min-merge** vs sideconnectors **active-store copy** (nine sites: `:623, :685, :694, :842, :851, :1054, :1063, :1192, :1201`) | Keep the copy carry everywhere; apply B/R1's `tDamp` ONLY to the four Picard→Newton promotion copies (`:685`, `:842`, `:1054`, `:1192`): `mOmegaNewton = clamp( tDamp * mOmegaPicard, ... )`. Demotion/latch copies stay undamped. The min-merge is NOT reintroduced except where R5 deliberately keeps it (escalated-Newton stall). |
| C2 | Stall guard (flat band) + moved-baseline detector vs progress watchdog — three Δt-cut/exit mechanisms on one loop | All three coexist: band catches oscillating-flat, watchdog catches no-new-best creep, moved-baseline rescues coupled-quench false rejects. Verify the moved-baseline ACCEPT path also resets the watchdog's best-tracker (else the watchdog cuts a step the detector just rescued). |
| C3 | Line-search loop body: matfix restructures accept/reject (A/R0 + flat-pair tracking) while sideconnectors threads Anderson commit/discard/flush through the same branches | Interleave: EVERY reject still runs `anderson_discard + anderson_clear`; the flat-pair bookkeeping adds no new accept path except moved-baseline (see O1). Flush triggers per plan §3.1 of `anderson_picard_acceleration_plan.md` stay complete. |
| C4 | Residual semantics: matfix thresholds (0.3 / 0.8 / +1.0 / 0.05 dec) were tuned on the LAGGED residual; sideconnectors measures the refreshed residual (`cc90ed90`) | Port constants unchanged; R9 trace gate decides whether a retune session is warranted. Do NOT tune blind in this port. |

## 3. Ordered steps

- [x] **R1** — Solver soft-fail contract (H-G): wrappers first (`cl_SolverWrapper.hpp`,
  MUMPS, STRUMPACK incl. Allreduce), then DofManager/SolverData plumbing, then the four
  controller consumers + `set_soft_fail` wiring + executable diffs. Independent of the
  line-search work; land + build first.
- [x] **R2** — Acceptance clause A/R0 (H-A) + moved-baseline detector (H-F) into
  `iterate_coupled`, honoring C2/C3. (after: R1 for clean testing, not a hard dep)
- [x] **R3** — Damped first Newton entry (H-B) per C1 at all FOUR Picard→Newton
  promotion-copy sites: magnetic coupled (`:685`), thermal coupled (`:842`),
  `iterate_magnetic` (`:1054`), `iterate_thermal` (`:1192`).
- [x] **R4** — Trust growth (H-C) in the coupled magnetic ω-adaptation path, plus the
  escalation guard (H-K) and print unclamp (H-L). Keep AIMD growth ACTIVE under
  Anderson exactly as sideconnectors has it (`:992-998`, O2 amendment / ts16 log) —
  do NOT reintroduce the stale growth hold.
- [x] **R5** — Stagnation-latch ω branch + half-throttle resume (H-D); latch re-arm
  union (H-E).
- [x] **R6** — Absolute THERMAL residual + flat-stall semantics (H-H/H-N): defaults 0.0,
  `mEpsilonAbs2` plumbing, consumption in `run_coupled`/`run_thermal`, re-arm points,
  warning box; save-grid verification (H-M). The magnetic absolute knob stays
  parsed-unused unless deliberately designed otherwise (document as such in the theory
  doc's input table).
- [x] **R7** — CN/Galerkin guard (H-I): port only the disabled-scheme guard hunks;
  preserve sideconnectors' BDF5 `mHDropped` restore in `cl_IWG_Timestep.cpp`.
- [x] **R8** — DONE 2026-08-06 (target is `make check`, not the stale `make tests`;
  `USE_TEST` toggled ON → suite → restored OFF). Port-relevant suites ALL PASS:
  containers, linalg, comm, math, sparse, mesh, **fem (incl. all 5 AndersonMixing
  tests)**, ode, physics, core, gastables (needs `BELFEM_DATA=<repo>/share`).
  Not port-related: io/kepler/manta don't compile (D2 + nonfree include paths;
  kepler/manta resolved 2026-08-25 — `Add_Test.cmake` used the undefined
  `SSF_SRC_DIR`, so ctest reported them "Not Run", never failed),
  gasmodels has 6 failing Cubic/Methane tests (staged gas migration, share/ data
  excluded per f6babfb0 — the module didn't even compile before D1). Theory doc
  updated. Build + `make tests` (incl. `test_AndersonMixing`); update
  `nonlinear_controller_theory.md` to describe the merged behavior (acceptance rule,
  damped entry, trust growth, moved-baseline, soft-fail — the doc currently documents
  the unhardened rules as current behavior).
- [ ] **R9** — Trace gates (Christian runs): (a) Greg CORC deck A/B vs matfix,
  `algorithm : Picard` AND `: Newton`, `target iterations : 20`; (b) a ts34-class
  coupled tape trace — Newton entry kick ≤ ~3 dB, no ω-floor pinning, no acceptance of
  multi-decade kicks. C4 retune decision hangs on this.
- [ ] **R10** (2nd wave, separable) — Thermal/material tangent harvest (H-J) by diff,
  including the `src/physics/materials/*` derivative APIs; only after R8 is green. May
  be its own session. FVM/QUAD4TS pseudoinverse in/out per O6.

## 4. Defects

Found during R8 (both pre-existing, exposed only because the shared tree keeps
`USE_TEST=OFF` so these units never build in the default target):

- [x] **D1** — `src/physics/gasmodels/cl_GM_EoS_Cubic.cpp:431` returned undeclared
  `aResult` (function computes `tResult`); compile error with tests enabled. FIXED in
  place (1 line) to unblock `make check`.
- [x] **D2** — RESOLVED 2026-08-06 (Christian's call: fix now): added the missing
  `HDF5::create_group( label, parent )` overload — the navigation state is a stack, so
  the parent handle must be the active group (BELFEM_ERROR-verified; in every test use
  the parent IS the active group, the argument makes the intent explicit). test_io now
  compiles and all 43 HDF5 tests pass, incl. the three nested-group tests that could
  never run. Original finding: `../../tests/io/test_HDF5.cpp` (commit `78c1692e` "more
  tests") called the 2-arg form which existed on NEITHER branch — committed against an
  API that never landed.

## 4b. Post-port amendment (R9 evidence, 2026-08-06)

- [x] **A1 — H-F thermal gate (deviation from the verbatim matfix port, Christian's
  "fix it now"):** the first R9 Garber trace (magnetic-only, fusing ON) showed the
  moved-baseline detector accepting a 44 dB kick (ts11 it32: ref −49.85 dB → −5.26 dB)
  and a +19 dB kick (ts14 it9) — its premise (staggered thermal update moved the
  baseline) cannot hold without a thermal kernel. Gate added: the accept branch now
  requires `mKernel2 != nullptr && ! mThermalFrozen` (thermal exists AND updated last
  iteration). matfix's motivating ts1727 coupled-quench case is preserved exactly.
  Theory doc §3 updated. Same trace also characterized "Newton does nothing": accepted
  Newton iterates bit-flat (±0.01 dB) across ω sweeps 0.6→0.07 — a tangent problem in
  the fused configuration, not a controller problem (stall guard exits correctly,
  Picard reaches −156.5 dB machine floor both timesteps).

## 5. Open questions

- [x] **O1** — EXECUTED 2026-08-06 as the recorded lean (no counter-ruling at dispatch):
  DISCARD+CLEAR — the pair's fixed-point residual mixes two thermal baselines. The
  moved-baseline accept also restarts the magnetic watchdog best-tracker (C2). Flip to
  commit-on-accept only if the R9 traces argue for it.
- [x] **O2** — EXECUTED 2026-08-06 as the lean: 0.5 kept at all four promotion sites,
  applied to the copy carry (`clamp( 0.5 · ωPicard )`). R9 decides any rescale.
- [x] **O3** — RESOLVED 2026-08-06 (Codex audit, verified): stale premise. The
  AIMD-growth hold no longer exists on sideconnectors — the O2 amendment (ts16 log,
  `:992-998`, theory doc `:70-71`) keeps growth ACTIVE under Anderson because holding ω
  was a one-way ratchet to the floor. Trust growth therefore also stays active; the
  flush-on-reject path is the overshoot guard.
- [x] **O4** — RESOLVED 2026-08-06 (Codex audit, verified `matfix:986-990, :1137-1142,
  :1262-1267`): matfix did NOT disable the divergence-counter ω reset; the "(turned
  off)" text at `matfix:cl_FEM_Controller.hpp:74` is a stale comment. Both branches keep
  the reset; port nothing, optionally fix the comment.
- [ ] **O5** — Does the R10 thermal wave ship in this port or as its own plan?
  (Christian; Greg's deck is magnetic-only, so it does not gate his regression.)
- [ ] **O6** — FVM/QUAD4TS pseudoinverse fixes sitting in the same matfix commits:
  include in R10 or exclude explicitly? (Christian. Note: sideconnectors already has its
  own eager-inverse/pseudoinverse line from the Blaze-defect session, dl20260803 —
  reconcile, don't double-port.)

## 6. Definition of done

R1–R8 landed on `sideconnectors`, `make tests` green on both backends, R9 traces show:
no full-ω Newton kick, no acceptance-hole spirals, solver failure ⇒ Δt cut (not abort),
Greg deck convergence ≥ matfix baseline with `algorithm : Newton`. Theory doc updated.
R10 dispositioned per O5. Plan moved to `` with DONE summary.

## 7. Audit trail

- `../../tmp/ai_exchange/review_greg_corc_regression.md` — round-1 jury (branch archaeology,
  controller findings F1–F8, reconciliation)
- `../../tmp/ai_exchange/review_side_edge_fusing_physics.md` — round-2 jury (fusing physics;
  separate campaign, listed for context)
- Devlog: `dl20260806_greg_corc_convergence_regression.md`
- matfix-side design record: `matfix:devlog/dl20260727_controller_line_search.md`
