# Refactor plan: Controller::iterate_coupled / iterate_magnetic / iterate_thermal

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): the helper extraction was judged not worth its diff risk; the three small controller bugs it re-verified (ω clamped after `set_omega` in the magnetic path, the shared divergence counter, the streak leaking across accepted steps) are candidates for a debt-register row, not a refactor plan. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-06-22
**Purpose:** De-duplicate the three `iterate*()` functions and fix the drift between
them, without changing convergence behavior (except where explicitly flagged).
**Status:** ACTIVE, but **substantially overtaken** — re-verified against the tree
2026-08-09. Tri-AI reviewed 2026-06-22 (Claude proposal, Codex + Grok audit, both high
confidence, both recommend Option A over a full unification).

> **2026-08-09 currentness sweep.** The premise of the Verdict below — that
> `iterate_magnetic()` is "a stale single-shot copy that never received" the coupled
> path's logic — is now only **half** true. The matfix controller port
> (`closed/matfix_controller_port_plan.md`, executed 2026-08-06) plus the intervening Anderson /
> ts17 work brought most of the drift across. Current state, all re-checked in
> `src/fem/kernel/cl_FEM_Controller.cpp`:
>
> | Item | Status |
> |---|---|
> | Execution step 2 — extract `impose_voltage_bcs()` | **DONE.** Now a real member (`:447`), called from `iterate_coupled` (`:935`) and `iterate_magnetic` (`:1589`). |
> | `iterate_magnetic` Picard↔Newton handoff, damped first Newton entry, watchdog re-anchor, `mForceNewton`, solver soft-fail → `reset_timestep` | **All present** — `iterate_magnetic` starts at `:1517`; no longer stale relative to `iterate_coupled` (`:772`). |
> | Bug: unclamped ω in `iterate_magnetic` | **STILL LIVE.** `set_omega( tOmega )` at `:1569` precedes `std::clamp` at `:1570`. Coupled clamps at `:875` before `set_omega` at `:914`; coupled-thermal clamps at `:1756` before `set_omega` at `:1759`. The segregated path can still exceed `max relaxation` for one iteration. |
> | Bug: shared `mNumIterationsDiv` | **STILL LIVE.** One member (`cl_FEM_Controller.hpp:99`) written by the magnetic, thermal-standalone and coupled-thermal adapts alike (`cl_FEM_Controller.cpp:1439,:1474`, `:1678,:1694`, `:1855,:1867`). |
> | Bug: `mNumIterationsDiv` never reset per timestep | **PARTLY FIXED.** `reset_timestep` (`:1943`) now zeroes it (`:2000`) along with the Anderson history and the divergence bookkeeping — so a *retry* is clean. A streak still leaks from one accepted timestep into the next; only an improving iterate zeroes it otherwise. |
> | Parity decision — backtracking line search in `iterate_magnetic` | **STILL OPEN, and now the largest remaining gap.** `tBacktracks` and the whole accept/reject loop exist only in `iterate_coupled` (`:895` onward). |
> | Extractions 3–6 (`adapt_relaxation`, `select_algorithm`, `diverged`, `stagnated`) | **Not done.** Still copy-pasted; the copies have grown since 2026-06-22, so the de-duplication payoff is larger now, and so is the diff risk. |
>
> **All `file:line` citations in the body below are from 2026-06-22 and are wrong now**
> (e.g. `iterate_magnetic` is at `:1517`, not `:625`). The table above carries the current
> anchors; re-locate anything else by symbol.
>
> **2026-08-11 currentness sweep.** All three surviving bugs re-checked directly in the tree
> and **all three are still live**; only the anchors moved, and the table above is re-baselined
> to `bc578b5e` + the working tree. The ω-ordering asymmetry was read out of the source in this
> sweep rather than carried over: `iterate_magnetic` sets then clamps, both other paths clamp
> then set. The recommendation below is unchanged — and the reason to hold the extraction is
> now weaker, not stronger, since the PID work and the D4 residual fix are both committed
> (`1d6ef305`, `4f2c11cd`); what still argues for holding is the unrun controller campaign
> (DR-52 greg3 A/B, matfix R9), which an extraction diff would confound.
>
> **Recommendation (unchanged in shape, sharpened):** do the three surviving bug fixes as
> their own commit — the ω-clamp ordering one-liner is trivial and the divergence-counter
> split is small — and treat the helper extraction as a separate, later decision. Do NOT
> start the extraction while the controller is still absorbing the PID timestep work
> (`pid_timestep_controller_plan.md`) and the D4 residual-semantics fix
> (`anderson_picard_acceleration_plan.md` §4.4); both touch the same functions.

## Verdict

Worth doing — the motivation is **drift, not aesthetics**. `iterate_magnetic()`
(`cl_FEM_Controller.cpp:625`) is a stale single-shot copy that never received the
backtracking line-search / grace / accept logic added to `iterate_coupled()`, and the
segregated executable path runs it (`hphiTrun.cpp:148-150`). The Picard↔Newton handoff,
voltage-BC block, ω-adaptation, divergence-reset and stagnation guard are copy-pasted
2–4× and have already diverged in non-obvious ways. **Option A** (extract narrow helpers,
keep the three entry points) captures most of the value at low risk; **Option B** (a single
`iterate(KernelContext&)`) is rejected for now — the thermal path differs in backtracking
applicability, solve ordering, reset action, stagnation test and a forced Picard-only mode,
so a single context would be a bag-of-booleans leaky abstraction. Revisit B only after A
plus characterization tests.

## Bugs to fix (verified at file:line; fold into the relevant helper extraction)

- [x] **Dead coupled thermal divergence check.** `iterate_coupled` increments `mIteration`
  (`:526`) but never `mIteration2`, while the reset check reads
  `mEpsilon2 > 1E1 && mIteration2 > mMinNumIterations2` (`:519-520`) → `0 > 2` is always
  false, so a diverging thermal residual is never caught in the coupled path. Decide:
  advance `mIteration2`, or check against `mIteration`.
  *Fixed 2026-07-09: `iterate_coupled` now advances `mIteration2` inside the
  `mKernel2 != nullptr` block; counter is zeroed per timestep in
  `initialize_timestep` / `reset_timestep`, so the guard fires as intended.*
- [ ] **Unclamped ω in `iterate_magnetic`.** `set_omega(tOmega)` is called before the
  `std::clamp` (`:662` before `:663`); coupled/thermal clamp first (`:394`, `:817`). The
  segregated path can exceed the configured `max relaxation`. Compounded by `mOmega0`
  being set at construction (`:41`) before `mOmegaMax` is parsed from input (`:1534`).
- [ ] **Shared `mNumIterationsDiv`.** One member (`cl_FEM_Controller.hpp:74`) written by
  both the magnetic (`:765/770`) and thermal (`:869/874`) ω-adapt; in the **segregated**
  path these cross-talk. (NOT a coupled-path bug — coupled thermal adapt at `:608-620`
  never touches it.) Give thermal its own counter.
- [◐] **`mNumIterationsDiv` never reset per timestep.** A divergence streak leaks across
  timesteps (only an improving step zeroes it). Reset it in `initialize_*`.
  *Partly fixed: `reset_timestep` zeroes it (`:1774`), so retries are clean; the leak from
  one accepted timestep to the next remains.*
- [x] **Coupled thermal residual uses `mIteration`** (`:523`) not `mIteration2` (standalone
  uses `mIteration2` at `:830`). Reconcile when fixing the dead check above.
  *Fixed 2026-07-09 together with the dead check: coupled path now passes
  `mIteration2` to `residual()`.*

## Helpers to extract (signatures per Codex's completeness notes)

- [ ] `select_algorithm(IWG_Timestep*, real& omegaPicard, real& omegaNewton, real eps,
  real epsSwitch, SolverAlgorithm terminal, uint iteration, real omegaMin,
  bool& firstFlip, bool justPicard, bool& algorithmChanged) -> SolverAlgorithm` — the
  handoff (`:365`/`:632`/`:794`); must report the algorithm change so the caller clears
  `mResidualHistory` (`:388`, `:654`).
- [x] `impose_voltage_bcs()` — the byte-identical rank-0 RHS patch. **DONE** — extracted as
  `Controller::impose_voltage_bcs()` (`:425`), called by `iterate_coupled` (`:831`) and
  `iterate_magnetic` (`:1382`); the `mCommRank == 0` guard sits inside the helper as
  specified.
- [ ] `adapt_relaxation(real& omega, real eps, real eps0, real omegaMin, real omegaMax,
  uint& numDiv)` — the `mBeta + mGamma*atan` grow / `mAlpha` shrink / streak-reset / clamp
  (`:540`/`:608`/`:761`/`:865`). Pass the divergence counter in so magnetic and thermal use
  separate ones.
- [ ] `bool diverged(real eps, uint iter, uint maxIter, uint minIter)` — the
  `eps > 1E1 && iter > minIter` + nan + max-iter trigger.
- [ ] `bool stagnated()` — the ShiftRegister deviation `< mStallBand` test (magnetic only;
  thermal keeps its cruder single-step latch, or gets a documented variant).

Then reduce each `iterate_*()` to a short composition of these.

## Decisions for Christian

- [ ] **Parity:** should `iterate_magnetic` gain the backtracking line search (so the
  segregated path stops running stale logic)? This is a deliberate behavior change for
  `hphiTrun`. (Recommend yes, but separate commit.)
- [ ] **Bug fixes:** apply the verified bug fixes as part of this refactor, or split into a
  preceding bug-fix commit so the refactor stays strictly behavior-preserving? (Recommend
  the latter: fix bugs first, then refactor.)

## Conventions / constraints (from the audits)

- `KernelContext` binding `mEpsilon` vs `mEpsilon2`: use a POD of non-owning pointers (or
  per-call refs in small helpers), not a fat reference-struct — same zero cost, clearer
  intent, and several fields are policy-dependent.
- Hot path: helpers must be inlinable, no hidden allocations (these run every nonlinear
  iteration). Backups are once-per-timestep so `Cell`/`Vector` there is fine.
- MPI: reset paths must stay all-rank decisions; preserve `residual()`'s barrier+broadcast
  ordering (`cl_FEM_DofMgr_SolverData.cpp:2350-2353`); never add an unmatched collective in
  a helper. (Pre-existing latent hazard: the rank-0 voltage-BC RHS patch updates only rank
  0 before the collective solve — out of scope here, but note it.)

## Execution order (each step behavior-preserving + `mpicxx -fsyntax-only` clean)

- [ ] 1. (optional) Preceding bug-fix commit for the 5 verified bugs.
- [x] 2. Extract `impose_voltage_bcs()` (simplest, byte-identical). *(done)*
- [ ] 3. Extract `adapt_relaxation()` (with per-kernel `numDiv`).
- [ ] 4. Extract `select_algorithm()` (+ `algorithmChanged` out-param).
- [ ] 5. Extract `diverged()` and `stagnated()`.
- [ ] 6. Reduce the three entry points to compositions; diff behavior carefully.
- [ ] 7. (decision) Bring `iterate_magnetic` to backtracking parity.
