# DR-106 + DR-144: One Round Over the MUMPS Wrapper

**Date:** 2026-08-30
**Purpose:** Close DR-106 (a MUMPS `-9` workspace soft-fail is answered by a timestep cut instead
of a workspace retry) together with DR-144's two deferred steps R12/R13, because DR-106's retry and
DR-144's R13 edit **the same two functions**. Doing them as separate rounds means two plan+audit
cycles over the same forty lines of `cl_SolverMUMPS.cpp`.
**Module:** `src/sparse` (wrapper + Fortran shim), `src/fem/kernel` (controller abort text, eigen diagnostic)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** **COMPLETE — both rows STRUCK 2026-08-30 ( Christian's ruling ) after both gates ran
the same afternoon.** Christian rebuilt ( shim + wrapper one pass ) and restarted
`tape_quench_usermat` on the new binary: two `-9` events at step 791, both recovered by the ladder
( 30→60→120→240, step converged, zero soft-fail boxes; the second ladder was the conditioning
probe's own per-call instance — per-instance persistence as designed ), and `check-fast` ran 15/15
suites green incl. `sparse`. Archived rows in `todo/debt_register_closed.md` carry the citations.
Earlier state of this line, kept for the record: R1–R6 implemented
in one session under the frozen protocol: plan jury (Codex terra/high + Grok 4.6/high) → amendments
→ implementation → code-diff jury (both at xhigh, MPI-collective/ABI safety boundary) → four
auditor-prescribed fixes applied and re-probed. Evidence level: `g++/gfortran -fsyntax-only` clean
on all touched TUs under Armadillo AND Blaze defines and a no-MUMPS arm — **below the compile/link
rung; nothing is verified until R7/R8 run.** Audit trail: `tmp/ai_exchange/dr106_dr144_joint.md`
(distilled into `devlog/dl20260830_mumps_minus9_retry.md`). Nothing here may be built while
Christian's `tape_quench_usermat` run holds the tree.

> **Scope guards:**
> - **OUT:** raising the hardcoded `ICNTL(14) = 30`. Framework-wide, and both auditors ranked it
>   last. If a stopgap is wanted before this round lands, the deck-only lever is
>   `library : strumpack` on the magnetic block — zero rebuild, and DR-105/106 record zero `-9`
>   on the identical deck.
> - **OUT:** DR-140's symmetric triangle extraction (`todo/mumps_symmetric_triangle_extraction.md`).
> - **OUT:** changing the soft-fail contract, the frozen-factorization scope, or `select_job()`.
> - **IN:** the `-9` retry, the three `soft_fail()` arms, R12, R13, and the controller's abort text.

---

## 0. Why this is a joint session, and why it is unblocked now

DR-144's R12 and R13 were marked **"DEFERRED — blocked on ownership, not on difficulty: the block
is a different session's uncommitted work."** That session (the one that fixed DR-143 and filed
DR-144) **closed on 2026-08-30 at ~04:09** — its socket is gone and its PID is dead — after
committing everything it held. `git status` is clean at `aecea892`. **The ownership block no longer
exists.** Confidence: high (verified directly: socket absent, process absent, working tree clean).

The overlap that makes this one round rather than two:

| | DR-106's retry | DR-144 R13 (G11) |
|---|---|---|
| file | `cl_SolverMUMPS.cpp` | `cl_SolverMUMPS.cpp` |
| site | inside the soft-fail block, **before** `flag_failure()` | the `mSolverID <= 0` bail, same functions |
| overloads touched | both (single-RHS and matrix-RHS) | both |

Two rounds over the same two functions is the waste. One round, one audit pair, one gate.

---

## 1. Where each row actually stands

**DR-106 — open, fix shape fully specified, no code.** Diagnosis is not in dispute: it was filed
2026-08-24 from the tape_quench four-run comparison (BDF5/mumps hit `-9` 41 times) and independently
re-derived 2026-08-30. The jury round added the constraints in §2 that the row's original text did
not name.

**DR-144 — GATE-PASSED 2026-08-29, `make check` green (Christian), 9 of 11 steps done.** G1–G9
landed. The green build was confirmed to have actually exercised the new test rather than inferred
from the suite result. **Read the plan's own caveat before trusting the gate:** five steps (R3, R4,
R5, R6, R11) land *reviewed-not-verified* — they have no reachable in-tree trigger, so a green build
says nothing about them. Full detail: `todo/dr144_mumps_lifecycle_plan.md`.

**Neither row is live in the running calculation.** Verified 2026-08-30:
- DR-144 G10's trigger is `set_equation()` → `SolverData::reset()` deleting `mSystemMatrix`. The
  only caller of `set_equation()` is `Kernel::create_field()` (`cl_FEM_Kernel.cpp:631`), which is
  setup-only, so it cannot fire inside the timestep loop. Confidence: high.
- DR-106 has produced **zero** `-9` events in the current run. **This is not evidence the exposure
  closed** — the mesh was refined between runs (227 468 free magnetic dofs against 148 400), so it
  is a different matrix with different pivoting. The previous mesh threw 24 soft fails on the same
  deck.

---

## 2. Gap table

| # | Gap | Class | Citation / rationale |
|---|---|---|---|
| **G-A** | `-9` retry cannot be keyed on a rank-local code. The shim exports `INFO`, not `INFOG` — the integer copy-out is `forall( k=1:40 ) aInfo( k ) = tMUMPS%INFO( k )` while only the REAL array is `RINFOG`. MUMPS sets `-9` on the failing rank and `-1` elsewhere. A retry written `if ( mInfo( 0 ) == -9 )` fires on ONE rank and strands the others outside the collective call | (c) | `mumpstools.f90` copy-out at the close of `mumpstools_solve`; both auditors, independently; re-verified against the shim |
| **G-B** | `print_soft_fail`'s `-9` arm is **dead code** on a multi-rank job: it is rank-0 only, and rank 0 was never the failing rank in any observed run (failing ranks 7, 4, 1, 6). Every box read "error on rank N"; the `-9` was legible only from raw MUMPS stdout | (c) | `cl_SolverMUMPS.cpp` `print_soft_fail`; run log `tape_quench_usermat/out.txt` |
| **G-C** | Re-analysis already happens and is **not** the missing ingredient. Post-fail `mMatrix = nullptr` forces JOB 6 every retry but never raises `MemoryRelaxation`, so each retry re-analyses and still meets `ICNTL(14) = 30`. MUMPS prescribes raising `ICNTL(14)` and repeating the FACTORIZATION for `-9`, reserving re-analysis for `-5`/`-7` | (c) | `cl_SolverMUMPS.cpp` soft-fail arms; `select_job()`; MUMPS 5.9.0 error `-9` |
| **G-D** | The retry must sit **inside the wrapper, before `flag_failure()`**. The controller aborts after eight consecutive solver failures (`mSolverFailCount`, reset on an accepted step), so a retry above the wrapper still loses a run needing nine. A warm-restart session died to exactly that cap after 8 fails from Δt 0.5 ms to 0.0039 ms | (c) | `cl_FEM_Controller.cpp` `mSolverFailCount` |
| **G-E** | That abort names the wrong failure — it reports a recoverable workspace event as "persistently singular or not converging". Third-category miss per `doc/coding_philosophy.md` | (c) | same site as G-D |
| **G-F** | `ICNTL(23)` is a documented **fourth** option, written nowhere in `src/`. MUMPS 5.9.0 recommends a per-rank working-memory cap in PARALLEL because parallel analysis estimates can be inaccurate. A sizing lever, **not** a substitute for the retry | (b) → **O1** | MUMPS 5.9.0 §5.12; grep of `src/` returns zero writes |
| **G-G** | The `.gt. 0` guard means a future `memory relaxation : 0` deck key would **skip** the `ICNTL(14)` write, inheriting a recycled instance's previous value instead of MUMPS's default. Latent today (C++ forces 30) | (c) | `mumpstools.f90` guarded write; recorded on DR-106 under extend-don't-branch |
| **G10** | `EigenValues::mK` captured once, can dangle on a WORKER; the block's "link_matrix() only runs on the master" comment is FALSE | (c) | `cl_FEM_DofMgr_EigenValues.cpp:370`, `:456`, `:605` — full text in the DR-144 plan |
| **G11** | `mSolverID <= 0` bail in BOTH solve overloads returns silently even when `soft_fail()` is false; the `BELFEM_ASSERT( is_initialized() )` below is unreachable for that case | (c) | `cl_SolverMUMPS.cpp`, both solve overloads |

**§2.1 Cross-cutting.** G-A is the one that turns a working fix into a hang. It was absent from
DR-106's text for six days and both auditors found it independently; it is now on the row. Treat any
patch that does not establish a rank-uniform code before retrying as incorrect regardless of how it
tests on one rank.

---

## 3. Ordered steps

- [x] **R1 — Make the error code rank-uniform.** *(DONE 2026-08-30: `aInfoG` export, five call sites, create/free untouched.)* Export `INFOG(1:40)` from `mumpstools_solve`
      alongside `INFO`, or allreduce the failing rank's `INFO(1:2)`. Decide which in **O2**.
      *(closes G-A; prerequisite for everything below)*
- [x] **R2 — Revive the `-9` reporting arm.** *(DONE: `print_soft_fail`/`error_message`/error arms on `mInfoG`; the `-9` box names the exhausted ICNTL(14); the shortage count deliberately dropped — INFOG(2) is in millions when negative and does not fit the field honestly.)* With a rank-uniform code, `print_soft_fail` can name
      the workspace instead of only the rank. *(closes G-B; after R1)*
- [x] **R3 — The retry itself.** *(DONE: shared `escalate_workspace()` helper, both overloads; cap + zero-hang guard + overshoot clamp per the code jury.)* Inside the wrapper, before `flag_failure()`: on a uniform `-9`,
      double `MemoryRelaxation`, retry JOB 2/5 **without** nulling `mMatrix`, bounded (see **O3**),
      and only then soft-fail. Must cover the single-RHS arm and the create arm; the multi-RHS arm
      is not live for this deck family (`num_rhs_cols == 1`) but should not be left inconsistent.
      *(closes G-C, G-D; after R1)*
- [x] **R4 — Correct the controller's abort text** *(DONE, 4 sites; wording names no wrapper-specific retry — STRUMPACK/PETSc take the same path.)* so a workspace event is not reported as
      "persistently singular". *(closes G-E; independent of R1-R3)*
- [x] **R5 (= DR-144 R13, G11).** *(DONE, both overloads, else-shape.)* `if( ! this->soft_fail() ) BELFEM_ERROR( false, … ) ; else { this->flag_failure() ; return ; }`
      in both overloads, so the silent path exists only for the caller that asked for it.
      **Same functions as R3 — land them together.** *(closes G11)*
- [x] **R6 (= DR-144 R12, G10).** *(DONE: unconditional recapture + BELFEM_ERROR, comment corrected.)* Assign `mK = mParent->jacobian()` unconditionally on every entry
      with a `BELFEM_ERROR( mK != nullptr, … )` beside it, and correct the false master-only comment.
      *(closes G10; independent of the rest)*
- [x] **R7 — GATE.** *(RAN 2026-08-30 14:50, Christian: `check-fast`, 15/15 suites green incl.
      `sparse` on the rebuilt tree; the full suite rides the nightly CI once committed.)* `make check` green, including `tests/sparse/test_Solver.cpp`. Note the existing
      `SolverLifecycle.MUMPSPoolRecyclesFreedSlots` is a **ceiling gate at pool size 8**, not an
      ID-recycle assertion — it does not cover this round's changes.
- [x] **R8 — GATE, the one that matters.** *(RAN GREEN 2026-08-30 14:51, Christian's restart:
      `tape_quench_usermat/out.txt:34758-34812` — step-791 magnetic solve laddered 30→60→120→240
      and CONVERGED; conditioning-probe instance laddered independently; zero soft-fail boxes.
      End-to-end reproducer, the top rung.)* Rerun a deck that reproduces `-9`. The red/green
      reproducer on record is the tape_quench BDF5/mumps configuration. **A green `make check` does
      not test the retry** — the retry has no unit-level trigger. Without R8 this round is
      reviewed-not-verified, exactly like DR-144's R3/R4/R5/R6/R11.

---


### 3.1 Accepted design (2026-08-30, this session — Christian resolved O1/O3/O4)

**Decisions:** O1 = ICNTL(23) **out of scope** (follow-up row if ever wanted). O3 = bound expressed
as a **maximum relaxation value 240** (30 -> 60 -> 120 -> 240, i.e. at most three doublings from the
default), **persist-bounded**: an escalated value stays in `mIParameters` for later solves, so a hard
phase does not repay the failed factorizations every step; `free()` restores the default. O4 = **no
deck key** this round; G-G stays recorded as latent. O2 (Claude's recommendation, auditors to
challenge): **export `INFOG`**, not allreduce — it mirrors what the shim already does for `RINFOG`
and keeps the rank-local `INFO` untouched for `check_warnings()`.

**Concrete edits, per file:**

**E1 (R1) — `mumpstools.f90` + `mumpstools.hpp` + wrapper members.**
- `mumpstools_solve` gains an `aInfoG` output (40 x `int_t`), placed between `aInfo` and `aRInfoG`
  in both the Fortran argument list and the C prototype; copy-out
  `forall( k=1:40 ) aInfoG( k ) = tMUMPS%INFOG( k )` beside the existing `INFO`/`RINFOG` lines.
- ~~`mumpstools_create_solver` returns `INFOG( 1 )` instead of `INFO( 1 )`~~ **struck in the
  plan-audit round (Grok, must-fix):** whether `INFOG(1)` is populated after a *failed* `JOB = -1`
  is unproven, and a zero there would skip the occupancy rollback — reopening the slot leak
  DR-144 just closed. Create and `mumpstools_free_solver` stay on `INFO`. The export is from
  `mumpstools_solve` **only**.
- `cl_SolverMUMPS.hpp`: new member `Vector< int_t > mInfoG ;` sized 40 in the ctor; all **five**
  `mumpstools_solve` call sites pass `mInfoG.data()` — :611, :639 (vector master/slave), :764
  (matrix master, Armadillo), **:786 (matrix master, Blaze `#else` — the arm a skim misses)**,
  :814 (matrix slave). Container adjudication (Codex wanted `Cell`, Grok accepted `Vector`):
  `Vector< int_t >` to match the siblings `mInfo`/`mIParameters` — extend-don't-branch; a lone
  `Cell` beside them would be the odd one out.
- Both solve overloads key the **error arm** on `mInfoG( 0 ) < 0` and decode from `mInfoG` (the
  sign was already rank-uniform via MUMPS error propagation; the *code* now is too). The warning
  arm stays on the rank-local `mInfo( 0 ) > 0` — warnings are not propagated, per the existing
  comment, and that behaviour is deliberately unchanged. The **hard-fail decode moves too**:
  `error_message( mInfoG.data(), … )` at both call sites, or its `-1`/`-9` arms keep printing the
  propagated payload (Grok #7). `mInfoG.fill( 0 )` joins `mInfo.fill( 0 )` at both solve entries.

**E2 (R2) — `print_soft_fail`** reads `mInfoG` instead of `mInfo`, so the `-9`/`-8` arms decode the
real code instead of the propagated `-1` ("error on rank N"). The create soft-fail arm fills
`mInfoG( 0 ) = tInfo` **alongside** the existing `mInfo` fill (the create/retry split, per the
audit: the create arm gets E2's reporting fill, never a relaxation ladder — `JOB = -1` cannot
produce a `-9`). The `-9` arm
additionally names the relaxation the attempt died at: `ICNTL(14) = %i exhausted` read from
`mIParameters`.

**E3 (R3) — the retry, both solve overloads.** File-local constants in the `.cpp`:
`tDefaultMemoryRelaxation = 30` (also used by `init_defaults()` and restored by `free()`) and
`tMaxMemoryRelaxation = 240` (the bound, stated in code per the definition of done). Structure: the
master's pre-work (`create_coo_indices`, base switch, `select_job`, `broadcast( tJob )`) stays
before a retry loop; the loop wraps only the parameter stores and the `mumpstools_solve` call
(master and slave arms); the master's base restore stays after the loop. Loop policy, decided
identically on every rank from rank-uniform data:

```
call the shim
if mInfoG( 0 ) != -9            -> break        ( success, or a different failure )
if tJob is not 5 or 6           -> break        ( a frozen JOB 3 never retries )
if relaxation already >= 240    -> break        ( bound reached -> existing soft-fail arms )
relaxation *= 2 ( persisted in mIParameters )   ( rank-uniform: every rank runs this line )
tJob = 5                                        ( analysis from the failed attempt is valid;
                                                  MUMPS prescribes repeating factorization only )
rank-0 notice line at InfoLevel::Minimal naming the new relaxation
```

On fall-through the existing arms run unchanged: a `-9` that survives the ladder soft-fails
exactly as today (`mMatrix = nullptr`, invalidate, print, return). A retry that succeeds reaches
the warning check with `mInfoG( 0 ) >= 0` and never touches `flag_failure()`.
**Retry job adjudicated in the plan-audit round — JOB 5 stands, on these MUMPS 5.9.0 quotes**
(`~/transfer/userguide_5.9.0.pdf`; both auditors objected *conditionally on the absence of a
quote*, and the condition is met):
1. *"-9 … the user should increase the value of ICNTL(14) before calling the factorization
   (JOB= 2) again"* — factorization repeat is the prescription;
2. *"JOB= 5 combines the actions of JOB=2 and JOB= 3. It must be preceded by a call to MUMPS with
   JOB= 1 on the same instance"* — satisfied: the failed JOB 6's analysis sub-phase completed,
   and `-3` fires only if analysis *"was not performed (or failed)"*;
3. *"ICNTL(14) … Phase: accessed by the host both during the analysis and the factorization
   phases"* — the raise takes effect at factorization time, no re-analysis needed;
4. JOB 2 **alone** is wrong through this shim (Grok): the shim copies `aY` into `aX` before the
   call, so a factorize-only JOB returns the RHS copy as "solution". JOB 5 is the only
   single-call shape that is both the documented recovery and correct here.
Failure containment: if MUMPS rejects the JOB 5 (`-3`), `INFOG( 1 ) != -9` breaks the ladder on
every rank → the existing soft-fail arms → status quo. No hang, no wrong answer.

**Audit-mandated additions:** (a) the doubling line runs on **every rank** — a master-only
escalation leaves workers at 30 while the host sends 240, and ICNTL(14) then disagrees across the
communicator (only the rank-0 *notice line* is rank-guarded); (b) `free()` gains a new
**unconditional** restore `mIParameters( MemoryRelaxation ) = default`, placed next to
`mMatrix = nullptr` (NOT inside `if( mInitialized )`) — today the value is written only in
`init_defaults()`, which only the ctor calls, so a torn-down wrapper would hand 240 to its next
`initialize()`; (c) the constants live in an anonymous namespace in the `.cpp`
(no `t` prefix at file scope); (d) the master's indexing-base restore stays **after** the loop —
restoring C++ COO between attempts would hand 0-based indices to the retry; (e) the stale
`mumpstools.f90:307` reference in the `init_defaults()` comment is corrected while there.

**E4 (R4) — `cl_FEM_Controller.cpp`, all four `mSolverFailCount` aborts** (lines 1237, 1655, 2193,
2409: linear + thermal, two code paths each): replace "persistently singular or not converging"
with wording that does not name a cause the controller cannot know:
`"linear solver failed on %u consecutive timestep attempts - the cause is in the solver messages above ( e.g. a singular system, non-convergence, or workspace exhaustion that survived the in-solver retry )"`
(and `thermal solver` in the two thermal sites).

**E5 (R5 = DR-144 R13, G11) — both `mSolverID <= 0` bails:**
```cpp
if( mSolverID <= 0 )
{
    if( ! this->soft_fail() )
    {
        BELFEM_ERROR( false, "MUMPS::solve() called without a solver instance ( initialize failed or was never called )" );
    }
    else
    {
        this->flag_failure() ;
        return ;
    }
}
```
(`else`-shape per the audit: `BELFEM_ERROR` never returns, but the `else` keeps a later reader
from thinking `flag_failure()` runs on the hard path.)

**E6 (R6 = DR-144 R12, G10) — `cl_FEM_DofMgr_EigenValues.cpp` ~1278:** unconditional
`mK = mParent->jacobian() ;` with `BELFEM_ERROR( mK != nullptr, ... )` beside it, and the comment
rewritten to drop the false "link_matrix() only ever runs on the master" claim (it runs at lines
370, 456 and 605) — the point of the unconditional assign is that a captured `mK` can be stale
after `SolverData::reset()`, on any rank.

**Deliberately untouched:** `mumpstools_free_solver` still exports rank-local `INFO` (teardown
path, same G-B class, out of scope — noted for a future row); the `.gt. 0` `ICNTL(14)` guard in
the shim stays (G-G latent, O4 = no deck key); `-8` is *not* on the retry ladder (same ICNTL(14)
remedy per MUMPS, but the scope guard names `-9` and no `-8` was observed — logged as **O5** below).

**Build atomicity note for Christian (Grok R-2):** the widened `mumpstools_solve` signature means
the Fortran shim and `cl_SolverMUMPS.cpp` must be **one rebuild** — a stale `.o` on either side of
the `bind(c)` boundary is an argument-shift smash, not a link error.

**Gates this session:** R7/R8 are **deferred** — the `tape_quench_usermat` run holds the shared
build tree (verified live, 8 ranks, started 03:43). Everything lands *reviewed-not-verified* until
Christian can build; the todo boxes for R7/R8 stay open and say so.

## 4. Open design questions (logged, not decided)

- [x] ~~resolved 2026-08-30: **(a) out of scope**~~ **O1 — Does `ICNTL(23)` belong in this round?** MUMPS recommends it in parallel. It is a sizing
  lever, not a retry. Options: (a) out of scope, (b) set a conservative per-rank cap alongside the
  retry, (c) expose both as deck keys. **Christian's call.**
- [x] ~~resolved 2026-08-30 (plan jury, both auditors concur): **export from `mumpstools_solve` only**; allreduce would not carry INFOG(2); create/free stay on `INFO` (slot-leak regression risk)~~ **O2 — `INFOG` export vs allreduce for R1.** Export is the smaller diff and matches what the shim
  already does for `RINFOG`; allreduce avoids widening the C interface. No strong recommendation
  yet — decide with the auditors.
- [x] ~~resolved 2026-08-30: **cap at ICNTL(14) = 240, persist-bounded, free() restores the default**~~ **O3 — Retry bound.** How many doublings before giving up, and does the escalated relaxation
  persist for later solves or reset per call? Persisting risks a ratchet; resetting risks paying
  the retry on every step of a hard phase. Both auditors asked for a bound; neither proposed a number.
- [x] ~~resolved 2026-08-30: **no key this round; G-G stays latent**~~ **O4 — Does `memory relaxation` become a deck key?** If yes, the two-artifact rule fires
  (`doc/input_file_reference.md` **and** `doc/input_schema.yaml`, same session) and G-G must be
  designed out: `0` must force an explicit write, not a skip.

---

- **O5 — Does `-8` (integer workspace) join the retry ladder later?** Same ICNTL(14) remedy per
  MUMPS 5.9.0, but never observed in the runs on record and outside this round's scope guard.
  Logged, not decided.

## 5. Definition of done

- [x] Every gap-table row mapped to a step or an open question.
- [x] R1 landed **before** any retry code — a rank-local retry is a hang, not a bug.
- [x] Retry bounded, and the bound stated in the code, not just the plan (`gMaxMemoryRelaxation = 240`, anonymous namespace, `cl_SolverMUMPS.cpp`).
- [x] `make check` green (R7) — `check-fast` 15/15 @ 2026-08-30 14:50; full suite via nightly CI post-commit.
- [x] **Deck rerun reproducing `-9` (R8)** — the only gate that exercises the retry. RAN GREEN, see R8.
- [x] DR-106 and DR-144 rows updated in the same session the code lands, then **STRUCK and
      archived** on the gate evidence; `check_doc_claims.py` 37/37 after both the update and the
      strike.
- [x] Devlog written (`dl20260830_mumps_minus9_retry.md`, strike addendum included); `devlog/README.md` index line added.
- [x] `src/sparse/doc/solver_memory_and_compression.md` updated (audit find, Codex): its two
      "`INFOG` is not copied back" claims (lines ~113, ~251) and the "hard-coded 30 %" row
      (~249) become false when this lands — fixed in the same session as the code (Codex language
      sweep dispatched on the touched cells).

---

## 6. Audit trail

- `tmp/ai_exchange/review_mumps_icntl14_ceiling.md` — the 2026-08-30 jury round (Codex + Grok,
  blind/parallel) that produced §2's constraints, plus a verification pass in which every auditor
  citation was re-checked against the tree and a reconciliation table with severities.
- `todo/dr144_mumps_lifecycle_plan.md` — DR-144's full plan, two jury rounds, D-tracker including
  D20/D21 (= G10/G11).
- `devlog/dl20260829_mumps_pool_recycling.md` — DR-143, the pool fix that removed G10/G11's trigger.
- Register rows: `todo/debt_register.md`, DR-106 and DR-144.

**Findings that were REFUTED in the jury round — recorded so they are not re-derived:** there is no
Δt threshold (failures are not monotone in Δt; one session failed at every Δt down to 0.0039 ms
while another succeeded at 0.0420 ms); "raise the hardcoded 30 to 100–150" is the weakest option;
and DR-144's instance-lifecycle class is **OUT** for the `-9` symptom, because `ICNTL(14)` is written
inside `mumpstools_solve` on every call and a recycled instance cannot carry a stale relaxation.

---

## 7. One thing the next session should check first

The 148 400 vs 86 215 free-magnetic-dof discrepancy on the **identical** 11 566-node mesh
(`tape_quench_usermat/out.txt`, warm-restart sessions against the cold start) is **unexplained**.
It is not DR-106 and not DR-144. It is either what `tape.bfm` restores versus a cold homology build
by design, or a restart-path defect. It surfaced during this jury round and was deliberately not
diagnosed. It may deserve its own row.
