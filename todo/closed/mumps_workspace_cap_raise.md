# MUMPS ICNTL(14) Retry Ladder: Raise the Workspace Cap

**Date:** 2026-08-30
**Purpose:** The `-9` (out-of-workspace) retry ladder tops out at ICNTL(14) = 240 after three
doublings. On the `tape_quench_usermat` deck the ladder is exhausted while the missing workspace is
still only ~2.4 MB, the step is soft-failed back to the controller, and Δt is halved — repeatedly.
Worse, because the escalation persists on the solver instance, once the cap is reached *every*
later solve gets a single attempt with no ladder at all. Raise the cap so the ladder has headroom
and, critically, so the retry mechanism is not permanently disarmed after the first hard step.
**Module:** `src/sparse`
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** **SUPERSEDED 2026-09-01 by `todo/mumps_memory_budget_plan.md`** — O1 resolved there at 480
(one rung), D2 (no ICNTL(23) ever written) became that plan's spine, D4 (the conditioning instance
walks the same ladder) its gap row 6. Kept for the evidence in §1. Previously: PLAN — **audited
2026-08-30 (jury: Codex gpt-5.6-terra/xhigh, Grok grok-4.6/xhigh).**
Verdict: revise before implementing. Mechanism sound; plan over-claimed. No P0. Four P1s (D1-D4)
logged below and folded into the steps. **Blocked on O1 (the ceiling value) — Christian's call.**
No code written. Exchange: `tmp/ai_exchange/review_mumps_workspace_cap.md`.

> **Scope guards:**
> - `gDefaultMemoryRelaxation` (30) is **OUT of scope** for this plan. Raising the floor costs
>   memory on every solve including easy ones and interacts with the persist ratchet; it is a
>   separate decision (see O2).
> - The underlying conditioning (the spectral ratio (5.5e13 - 2.4e14)) is **OUT of scope**. This plan buys margin, it does
>   not fix the physics.
> - No `input.conf` key is involved — ICNTL(14) is not a deck key
>   (`src/sparse/doc/solver_memory_and_compression.md:249`), so the Input Contract rule does not
>   trigger.

---

## 1. Evidence

Observed in `cmake-build-debug/tape_quench_usermat/out.txt`, step 47, after a warm restart at
t = 6.0 ms. `INFOG(2)` is the missing workspace in real entries (it is in *millions* only when
negative — `cl_SolverMUMPS.cpp:1758-1762`; all values below are positive, so they are direct counts).

| Δt (ms) | ICNTL(14) attempts | INFOG(2) deficit | outcome |
|---|---|---|---|
| 0.5 | 30 → 60 → 120 → 240 | 1,905,777 → 606,895 → 3,425,627 → 4,243,613 | ladder exhausted, soft fail |
| 0.25 | 240 only | 2,327,741 | soft fail (no rungs left) |
| 0.125 | 240 only | **302,097** (~2.4 MB) | soft fail (no rungs left) |
| 0.0625 | 240 | — | succeeded |

Two facts drive this plan:

- **F1.** ~~At Δt = 0.125 ms the step died 302,097 entries (~2.4 MB) short.~~
  **CORRECTED 2026-08-30 after audit (Grok, D5).** `out.txt:1947-1951` shows `Magnetic Picard 1,
  residual 0.917532` printing *before* the `-9`: that step's FIRST factorization succeeded at 240
  and a later iterate failed. The clean case is instead **step 48** (`out.txt:1990-1993`):
  `INFOG(2) = 89764` (~0.7 MB) with no preceding Picard line — a genuine first-solve failure,
  single-shot at 240 with no rungs left. Confidence: high for the corrected reading; **medium**
  that 480 would have carried it (the pivoting path is not monotone).
- **F2.** The escalation is a persist-bounded ratchet by design (`cl_SolverMUMPS.cpp:556-560`), so
  after step 47 the instance sits at 240 and `escalate_workspace()` returns false immediately on the
  `>= gMaxMemoryRelaxation` gate (`:566-569`). Every later step therefore gets **one attempt and no
  retry**. Confidence: high (read directly from the gate).

**Counter-evidence, recorded deliberately.** At fixed Δt the deficit is **non-monotone** — it fell
to 606,895 then rose to 3,425,627 and 4,243,613 as workspace grew. That is not the signature of a
merely undersized workspace; it says the pivoting takes a different, worse path on each attempt,
which is what the spectral ratio (5.5e13 - 2.4e14) does. **A higher cap will rescue marginal steps like F1; it will not make
genuinely hard steps easy.** Anyone reading this later should not expect a cure.

## 2. Change

- [ ] **R1.** `src/sparse/cl_SolverMUMPS.cpp:34` — `gMaxMemoryRelaxation` 240 → 960.
      Ladder becomes 30 → 60 → 120 → 240 → 480 → 960 (five doublings; 960 = 30 · 2⁵).
- [ ] **R2.** Comment sweep — **THREE sites, not one** (D1). `grep -n "240"` returns exactly
      `:30`, `:34`, `:437`, `:558`:
      - `:28-31` "30 -> 60 -> 120 -> 240, at most three doublings"
      - `:437` "( up to 240 ) to whatever initialize() comes next"  **(missed in v1)**
      - `:558` "the bound: 30 -> 60 -> 120 -> 240"  **(missed in v1)**
- [ ] **R3.** `src/sparse/doc/solver_memory_and_compression.md:249` — states the "240 % ceiling
      (30 → 60 → 120 → 240)". Update to match R1. **NB this file already has uncommitted edits.**
- [ ] **R4.** ~~Confirm the ladder reports rungs above 240 and steps that soft-failed now
      factorize.~~ **REWRITTEN after audit (D4): as written this would be misread.** The
      conditioning instance also walks the ladder (E2), so "a rung above 240 appeared" does not
      show the *production* solver was rescued. Revised gate:
      1. Disable conditioning for the gate run, or label the instance in the log.
      2. Require the first failing **magnetic** solve after a fresh instance to report 480/960.
      3. Score on the *marginal* single-shot failures (step 48's 89,764; step 54's 305k) now
         factorizing — **do NOT require Δt = 0.5 ms to succeed** (its deficit GREW to 4.2M).
- [ ] **R4a.** **NEW (D9).** The gate run MUST use `-v 4`. As of the 2026-08-30 log-hygiene change
      (`devlog/dl20260830_mumps_minus9_retry.md`, Addendum 2) `ICNTL(1)/(2)/(3)` follow the info
      level, so the `On return from DMUMPS, INFOG(2)=` lines — the only source of the deficit
      numbers R4 scores on — are **suppressed below `-v 4`**. Without this the gate cannot be read.
- [ ] **R5.** **NEW (D2).** Record per-rank peak RSS and MUMPS's post-analysis memory estimate at
      480 and 960 before accepting a framework-wide cap. `mumpstools.f90` writes ICNTL 1-7, 14 and
      others but **never ICNTL(23)**, so per-rank allocation is unbounded from BELFEM's side.
- [ ] **R6.** **NEW (D6).** Register this file in `todo/README.md`.

~~No other site references the cap.~~ **WRONG — see R2 (D1).** `cl_SolverMUMPS.hpp:151-152`
mentions `gMaxMemoryRelaxation` without a number (no edit needed; this cite was challenged by Grok
as `:152-153` and my original is correct — `:153` is `*/`). No test references it: the search path
in v1 said `test/`, which **does not exist** — the tree is `tests/` (D7). Conclusion survives;
both auditors grepped `tests/` independently and found none.

## 3. Gaps / risks

| # | Risk | Assessment |
|---|---|---|
| ~~G1~~ | ~~**MUMPS-internal 32-bit workspace sizing.**~~ **REFUTED 2026-08-30 (D3).** MUMPS sizes its work arrays with its own `MUMPS_INT`, 32-bit in a default build, independent of BELFEM's `BELFEM_INT64` (`src/core/typedefs.hpp:48`). At 960 % the workspace is ~10.6× the analysis estimate. | For this deck the estimate is ~10⁶ entries, so ~10⁷ — nowhere near the 2.1e9 int32 ceiling. But this constant is **framework-wide**: a much larger deck could plausibly reach it, and MUMPS's failure mode on internal overflow is not obviously a clean `-9`. **This is the main thing the audit must weigh.** Confidence: low — I have not established MUMPS's actual build-time int width, only that the default is 32-bit. |
| G2 | Persist ratchet now ratchets 4× higher. Once a run touches 960 it holds 960 until `free()`. | Memory measured: 6.85 GB across 8 ranks at the current 240, against 39 GB available. Headroom is ample *for this deck*. Not established for larger decks. |
| G3 | More rungs = more failed factorizations before giving up, so a genuinely hopeless step now costs ~2 extra full factorizations of wall clock before the controller gets it back. | Accepted: a factorization here is seconds, a Δt halving costs far more. But it is a real cost and should be stated, not hidden. |
| G4 | 960 is a judgement call, not a derived number. | Alternatives: 480 (one rung, minimal), 1920 (six doublings). See O1. |

## 4. Open questions

- [ ] **O1.** Is 960 the right ceiling, or is 480 sufficient / 1920 safer? F1 only demonstrates that
      *one* more rung would have helped one step. Nothing in the evidence justifies five doublings
      specifically over four.
- [ ] **O2.** Should `gDefaultMemoryRelaxation` rise off 30 as well? The analysis estimate is
      clearly under-predicting for this matrix class, so 30 guarantees a climb on every hard step.
      Deliberately out of scope here; logged so it is not silently dropped.
- [ ] **O3.** Should the cap be a solver parameter rather than a compile-time constant, so a hard
      deck can raise it without a rebuild? Larger change; not proposed now.

## 5. Pre-registration (before dispatch, per the frozen protocol)

What I expect the auditors to say, recorded so that agreement cannot be mistaken for confirmation:

- I **expect both** to accept R1–R3 as mechanically correct and low-risk for this deck.
- I **expect at least one** to raise G1 (int32 workspace sizing) as the substantive concern, and I
  regard that as the finding most likely to change the number.
- I **expect** challenge on O1 — that 960 is unjustified by the evidence, which I agree with; I
  would accept 480 on argument.
- I would treat **unanimous approval with no G1 discussion as a weak result**, not a strong one:
  it would mean neither auditor checked the one thing I could not establish myself.


---

## 6. Audit findings (jury round, 2026-08-30)

Full record and verification pass: `tmp/ai_exchange/review_mumps_workspace_cap.md`.
No P0. Verdict: **revise before implementing.**

- [ ] **D1 (P1, Codex + Grok).** R2 incomplete — stale `240` literals at `cl_SolverMUMPS.cpp:437`
      and `:558`. Verified by Claude: `grep -n "240"` returns exactly `:30 :34 :437 :558`.
      Folded into R2.
- [ ] **D2 (P1, Codex).** No ICNTL(23) memory limit is ever written
      (`grep "ICNTL( 23 )" src/sparse/mumpstools.f90` → no match, verified). Per-rank allocation is
      unbounded from BELFEM's side, and the cap raise applies to **both** MUMPS instances. R4 had
      no memory gate. Folded into new R5. **Single raiser — needs human adjudication.**
- [x] **D3 (P1, Codex + Grok).** **G1 REFUTED — my own claim was wrong.**
      `tmp/MUMPS_5.9.1/src/dfac_mem_dynamic.F:23-33`: `INTEGER(8) :: MAXS` / `ALLOCATE(S(MAXS))`;
      `dfac_driver.F:218-220`: `INTEGER(8) :: LWK, LIWK`. The main real workspace is **64-bit**
      sized, so the "2.1e9 MUMPS_INT ceiling on the -9 path" does not exist. Verified by Claude
      against MUMPS source. Closed: the risk is replaced by D2, not carried forward.
- [ ] **D4 (P1, Codex + Grok).** R4 could not distinguish the production solver from the
      conditioning instance reaching a high rung. Matches Claude's independent E2. Folded into R4.
- [x] **D5 (P2, Grok).** F1 overstated — the 302,097 failure was a later iterate, not the first
      factorization. Grok supplied a **better** case (step 48, 89,764, first-solve). Corrected in §1.
- [x] **D6 (P3, Codex).** Not registered in `todo/README.md`. **Done 2026-08-30.**
- [ ] **D9 (P1, Claude, post-audit).** The R4 gate reads `INFOG(2)` deficits, which the same-day
      log-hygiene change suppresses below `-v 4` (`mumpstools.f90`, ICNTL(1)/(2)/(3) now follow the
      info level). Found while checking whether the 16:14 edit invalidated the field evidence — it
      does **not** (log hygiene only, no `-9` semantics touched), but it does break the gate's
      readability. New R4a. Confidence: high.
- [x] **D7 (P3, Codex + Grok).** v1 cited a `test/` directory; the tree is `tests/`. Corrected.
- [x] **D8 (P3, Grok).** κ quoted as a constant 5.5e13 and conflated with `MUMPS ADD COND1`.
      Corrected in §1 / §G-table.
- **FALSE POSITIVE (Grok, retracted by Claude's verification).** Grok corrected the header cite to
  `cl_SolverMUMPS.hpp:152-153`; the original `:151-152` is right (`:153` is `*/`). No change.

### Open question elevated by the round

- [ ] **O1 (unchanged, now explicit).** Neither auditor would derive 960 from the evidence, and
      Claude agrees. Corrected F1 supports "one more rung would have rescued one first-solve
      failure" and nothing stronger. **480 is the defensible minimum; 960 is headroom.**
      **Design call — routed to Christian, not settled by the round.**
