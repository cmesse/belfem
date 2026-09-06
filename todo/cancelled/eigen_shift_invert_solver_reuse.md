# Cut the Eigen Diagnostic's Second Resident MUMPS Factorization

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): low relevance — per-call release of the dedicated eigen solver is a memory refinement nobody is blocked on. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-30
**Purpose:** The conditioning diagnostic keeps a second MUMPS instance, with factors, resident for
the whole run. Remove that cost. **Two options; the audit reversed which one to try first.**
**Module:** `src/fem/kernel` (`cl_FEM_DofMgr_EigenValues.cpp`), `src/sparse`
**AIs involved:** Claude (draft), Codex `terra`/`xhigh` + Grok `grok-4.6`/`xhigh` (plan audit r1)
**Status:** OPEN — **round-1 audited 2026-08-30; both auditors said REVISE BEFORE CODE, and the
draft's central premise was wrong.** Plan rewritten the same day: the recommendation is now
**Option A (free per call)**, with the borrow demoted to Option B behind a measurement. No code.

> **Scope guards:**
> - **OUT of scope:** changing *what* is computed. If a number moves, the change is wrong.
> - **OUT of scope:** the non-MUMPS-production case — it keeps its dedicated instance either way.
> - **Sequencing:** `todo/closed/mumps_workspace_cap_raise.md` wants `gMaxMemoryRelaxation` 240→960 on this
>   same wrapper. One constant, two plans. Sequence them (Grok).

---

## 0. What the audit changed

Four of the draft's load-bearing claims were false or stale. Recorded so they are not re-derived.

| Draft claim | Verdict | Truth |
|---|---|---|
| "Same matrix — confirmed" (`mK = mParent->jacobian()`) | **FALSE under Newton** (both auditors) | `SolverData::jacobian()` returns **`mSystemMatrix`** (`SolverData.hpp:679-681`), while a Newton body solves **`mJacobianMatrix`** (`SolverData.cpp:2428`) — separate allocations (`:396`, `:405`), with `dJdx` added only to the Jacobian (`assemble_newton`, `:1683-1709`). Claude verified independently. The draft trusted an accessor named `jacobian()` to return the Jacobian. |
| "`compute_conditioning()` reassembles, so production factors are of a previous assembly" | **FALSE on the live path** (Grok) | `mMatrixFlag` is set **only** in `compute_lambda_max()` (`:1650`) and cleared by `EigenValues::reset()` on every assembly (`:1805`). `compute_conditioning()` never sets it, so `:941-943` does not run at finalize. The "previous assembly" fact is still true — but because of the DR-127 head assembly and the Newton/system split, not this flag. |
| DR-106's "created and freed per call" (quoted into H1) | **STALE** (both) | The dedicated instance is deliberately retained for the `EigenValues` lifetime (`:1511-1520`). Today's isolation is "separate object", not "per-call free". |
| Motivation: "observed laddering to 240 in a STRUMPACK/PETSc run" | **Misattributed** (Grok) | On a STRUMPACK/PETSc field the borrow **never fires** — the saving there is exactly zero. The 227k-dof / `ICNTL(14)=240` evidence comes from a run whose **magnetic production solve was MUMPS**. The observation that prompted this work is not the case it would fix. |

**And the finding that reverses the recommendation:** the comment defending the retained instance
(`:1511-1520`) justifies it as *"a per-call free/re-init burns a fresh slot in the mumpstools
registry each time ( … observed as the tapestack3d step-8 abort )"*. **DR-143 fixed exactly that on
2026-08-29** — slot recycling, gated green with `SolverLifecycle.MUMPSPoolRecyclesFreedSlots` and a
tapestack3d rerun reaching *"43 conditioning calls with zero pool-exhaustion aborts where the old
allocator died on the 8th"*. The comment's primary rationale is one day out of date.

## 1. Current Behaviour

`run_shift_invert` builds a dedicated MUMPS solver (`:1206-1208`) and **keeps it, with its factors,
until `~EigenValues`** (`:1511-1520`). Production MUMPS also holds factors until `JOB -2`, and
`Solver::free()` has no FEM caller. So after the first diagnostic call, **two factor sets are
resident for the rest of the run** — confirmed by both auditors, and by the DR-106 log in which
production *and* the diagnostic each walked the `ICNTL(14)` ladder to 240 on one step.

**Bottom line:** the cost is persistent double residency, not a transient peak — and the reason the
instance is kept is a defect that has since been fixed.

## 2. Option A — free the dedicated instance per call **(recommended)**

`free()` (JOB `-2`) the diagnostic's own solver at the end of each `run_shift_invert`.

- **Keeps the isolation the current design states outright** (`:1237-1242`, *"A DIAGNOSTIC MUST NOT
  KILL THE RUN"*). No borrow, no shared state, none of H1–H5 below.
- Removes the persistent second factor set. The transient peak during the ARPACK loop remains.
- **Cost: a fresh JOB `-1` + JOB 6 per timestep** — the analysis reuse the comment's *second*
  justification names (`:1516-1517`) is genuinely lost. That is a time-for-memory trade, and it is
  the whole risk surface.
- **Unblocked by DR-143.** The slot-exhaustion objection is fixed and gate-proven.

- [ ] **A1** — `free()` the dedicated solver at the end of `run_shift_invert`, on every rank
      (`free()` is collective). Keep the object; free the MUMPS instance.
- [ ] **A2** — Rewrite the `:1511-1520` comment. It currently defends the opposite decision on a
      stale premise; leaving it would mislead the next session as it misled this plan's draft.
- [ ] **A3** — Measure. Peak **and** steady-state RSS, and per-timestep diagnostic wall time,
      before/after, on a **MUMPS-production** deck. If the re-analysis cost is unacceptable, that
      result is what justifies looking at Option B.

## 3. Option B — borrow the production solver (only if A3 says the peak still hurts)

Audited in detail and **not recommended as a first move**. Both auditors: revise before code.

### 3.1 The hazards, corrected and extended

| # | State | Status after audit |
|---|---|---|
| H1 | `ICNTL(14)` escalation (30→60→120→240) | Real. **But restore the *saved pre-diagnostic* value, not a hard 30** — production may legitimately sit at 120/240 under DR-106 persistence, and writing 30 would undo its ladder. The C++ slot is sufficient for the *parameter* (the shim writes `ICNTL(14)` every call, `mumpstools.f90:421`); it does **not** release workspace already allocated at 240. |
| H2 | `set_soft_fail(true)`, `clear_failure()` | **Overstated in the draft** (both auditors). The controller already sets production soft-fail at construction (`Controller.cpp:115-122`), `finalize()` only runs on accepted steps, and the next production solve clears the latch anyway (`SolverData.cpp:2412`). Keep as hygiene for non-controller callers; it is not a controller-visible latch bug. |
| H3 | freeze / unfreeze across the ARPACK loop | **The load-bearing hazard** (Grok). A leaked freeze on a *Picard* body is a **silent wrong production solve**: identity is pointer/shape/storage, so changed values are undetectable (`cl_SolverMUMPS.cpp:620-624`). On Newton it aborts instead. |
| **H4** | `mMatrix` / JOB selection | **MISSED by the draft.** The diagnostic's solve sets `mMatrix` to `mSystemMatrix`; a failure nulls it. The next production Newton body then takes JOB 6 instead of JOB 5 — extra analysis, not a wrong factor, *unless* freeze leaks. Must be a documented policy, not a scalar restore. |
| **H5** | Fortran-side workspace after a `-9` ladder | **MISSED.** A C++ `ICNTL(14)` restore does not `JOB -2`, so peak RSS on the borrowed instance can stay high — which attacks the plan's own justification. |

Also: the borrow **inherits production's parameters** (BLR, reordering, refinement), where a
dedicated `new Solver(MUMPS)` gets defaults with BLR off (`cl_SolverMUMPS.cpp:144-148`). If
production runs BLR, the diagnostic's numbers can move — which would violate the scope guard.

### 3.2 Corrections to the mechanics

- [ ] **B1** — Two pointers (owned + non-owning), never an `mOwns` flag. **Refresh the borrowed
      pointer every call**: `SolverData::set_solver()` deletes and replaces its solver
      (`:2044-2050`), so a cached borrow dangles. Update the "owns four raw pointers and deletes all
      four" comment (`cl_FEM_DofMgr_EigenValues.hpp:292-294`) or it will lie.
- [ ] **B2** — **Broadcast the borrow decision.** `SolverParameters::synchronize()` does *not*
      synchronize solver type (`cl_SolverParameters.cpp:323-370`) and `Solver::mType` is fixed at
      construction, so `type() == MUMPS` is a *local* predicate. Ranks disagreeing enter `DMUMPS` on
      different slots. Allreduce the bit; on disagreement all ranks take the dedicated path.
- [ ] **B3** — RAII guard for H1 (saved value) + H3 (unfreeze) + H2 (hygiene). No `clear_failure()`
      on a borrowed solver.
- [ ] **B4** — Uninitialized production wrapper: let `Solver::solve` initialize it
      (`cl_Solver.cpp:159-162`). **Do not fall back to a dedicated instance** — that re-creates the
      second slot the plan exists to remove. (Draft had this backwards.)
- [ ] **B5** — Pick and write down the H1 policy: restore-to-pre-diagnostic (isolation) *or* leave
      the escalation (DR-106 persistence). These conflict; it is a policy choice, not a detail.

## 4. Gates

- [ ] **G1** — Diagnostic value unchanged, fresh process, same ranks, **MUMPS-production deck**.
      Invalid if production BLR is on (Option B inherits it). No in-tree test calls
      `compute_conditioning()` at all — `tests/` search is empty.
- [ ] **G2** — ~~force a `-9` inside the diagnostic and assert production returns to 30~~
      **NOT RUNNABLE as drafted** (both auditors): in the observed log production had *already*
      laddered to 240 before the diagnostic ran, so a borrowed diagnostic starts at 240 and can
      never exercise 30→240. Needs **deterministic fault injection** (force `MemoryRelaxation` low,
      or a test hook in `escalate_workspace`). Option B must not land without it.
- [ ] **G3** — ~~a deck where production fails after a diagnostic call~~ **Insufficient**: production
      already soft-fails and re-clears the latch. Replace with a **unit test of the RAII guard**.
- [ ] **G4** — Peak **and** post-diagnostic RSS, per rank, MUMPS deck, fresh process. This is the
      gate that decides A-vs-B and it is mandatory for either.
- [ ] **G5** — `make check` — a build regression only; it does **not** exercise this path.

## 5. Open Questions

- **O1 — RESOLVED 2026-08-30 → not a blocker.** Both auditors traced it: after an accepted timestep
  the production instance is idle, synchronous, initialized, and holds factors of the *last body*
  solve. `compute_conditioning()` runs only from `finalize()`, which drivers skip on a reset
  (`hphiTrun.cpp:121-137`, `belfem.cpp:202-217`), so a failed last solve never reaches it. The one
  edge — no body solve has ever run — is handled by `Solver::solve`'s lazy initialize, **not** by a
  dedicated fallback.
- **O2 — REDIRECTED.** DR-144's slot invariants are *not* the risky neighbor (a borrow reduces pool
  pressure). The real coupling is **DR-106's `ICNTL(14)` persistence** on a shared instance (H1/B5),
  and the `gMaxMemoryRelaxation` collision with `todo/closed/mumps_workspace_cap_raise.md`.
- **O3 — ANSWERED, and it is a requirement, not a question:** the borrow decision must be
  broadcast (B2).
- **O4 — RESOLVED 2026-08-30, Christian → not a defect, and not to be "fixed".** Under Newton the
  eigen diagnostic factors `mSystemMatrix` while production factors `mJacobianMatrix`, so
  `|λ|max/|λ|min` and `MUMPS ADD COND1` are computed from **different matrices**, not merely
  different quantities of one. **Christian's ruling: that is the nature of the thing.**
  `|λ|max/|λ|min` of the system matrix is a good estimator of *problem difficulty*, which is what
  that diagnostic is for; COND1 bounds the forward error of the solve production actually
  performed. Both are doing their job. No row is filed, and a future session must not "align" the
  two operators — recorded here precisely because the discrepancy looks like a bug to anyone who
  finds it cold. This also **answers O3 of `todo/closed/mumps_error_analysis_opt_in.md`**, which was the
  last thing blocking that plan's R8b.

## 6. Definition of Done

- [ ] Option A measured (A3/G4) before Option B is considered at all.
- [ ] If B: H1–H5 each have a written disposition, B5 policy chosen, G2 injectable, G4 green.
- [x] ~~O4 filed as its own row rather than carried here.~~ **Resolved by ruling instead
      (2026-08-30, Christian): expected behaviour, no row.**
- [ ] Code round after implementation — this touches solver lifetime and MPI collectives, the
      §9.1 safety-boundary row.

## 7. Audit Trail

`tmp/ai_exchange/eigen_solver_reuse.md` — round 1, 2026-08-30, Codex `gpt-5.6-terra`/`xhigh` and
Grok `grok-4.6`/`xhigh`, parallel, on an identical brief that named O1 as the blocking question and
invited a "this kills it" answer. **Both returned revise-before-code.** They agreed independently
on the Newton/system matrix split, the incompleteness of H1–H3, and that G2 is not runnable as
written. Codex additionally traced the `mMatrix`/JOB-selection consequence (H4) and the parameter
inheritance; Grok additionally refuted the `mMatrixFlag` reassembly premise, caught the misattributed
STRUMPACK/PETSc motivation, found the `synchronize()` solver-type gap, and **proposed Option A**,
which is now the recommendation. Every claim above was re-verified against the code before being
written in: `SolverData.hpp:679`, `SolverData.cpp:2428`/`:396`/`:405`/`:1683`,
`EigenValues.cpp:941`/`:1650`/`:1805`/`:1511-1520`, `EigenValues.hpp:292-294`, and the DR-143 strike
record in `todo/debt_register.md`.
