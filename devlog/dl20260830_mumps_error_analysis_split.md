# `mumps error analysis`: splitting the MUMPS ADD diagnostic out of `compute conditioning`

**Date:** 2026-08-30
**Purpose:** Session record — a new deck key that makes the MUMPS ICNTL(11) error analysis
opt-in and independent of the eigenvalue conditioning estimate; the naming decision behind it;
the two-round plan audit that blocked and then cleared the design; and the code-audit round.
**Module:** `src/fem/kernel` (+ input contract, `src/fem/doc`, `src/fem/maxwell/doc`)

## Why

Christian: the `MUMPS ADD COND1` footer value is noise for his purposes — κ is actionable,
COND1 is not — and it should be the user's choice. Investigation found the coupling ran both
ways: one flag pair gated the ARPACK eigen estimate *and* the ICNTL(11) numbers, so neither
diagnostic could be requested alone. (The cost premise inverted along the way: the ADD pair is
the *cheap* half — triangular solves against the existing factorization — while the eigen
estimate is the expensive one. The noise argument, not the cost argument, carried the split.)

A Gemini-provided explanation of the MUMPS-vs-ARPACK discrepancy was checked against the MUMPS
5.9.1 source and **rejected for the documentation**: its headline claim (COND1 is computed on
the scaled matrix; "ICNTL(7)" as scaling) is the exact error this repo refuted on 2026-08-29
(`dsol_driver.F:6101` — "Notice that D is always the identity"). Its one useful residue is O3
in the plan: whether the operator ARPACK iterates carries the same Dirichlet/penalty
contributions as the one MUMPS factors. Open; blocks only R8b (the causal "why" paragraph).

## The name

Blind jury round (`tmp/ai_exchange/mumps_cond_key_name.md`, Codex terra/high + Grok
grok-4.6/high; Claude's pre-registration sealed in the scratchpad so the auditors could not
read it). Split verdict: both ranked `compute cond1` < `compute mumps conditioning` <
`mumps error analysis` — the middle option judged *harmful* because it contains the word the
split exists to separate — but Grok proposed dropping the vendor (`error analysis`) on the
capability-in-the-key convention. Rejected on a ground neither auditor weighed: the flag is a
**no-op without MUMPS** (the arming functions return early on solver type), so a vendor-free
name would sit silently inert in every STRUMPACK deck. Christian pinned **`mumps error
analysis`** — permanent, since renaming a deck key after release breaks user decks.

Second ruling, same session: the κ₂/`|λ|max/|λ|min` footer label branch is **deleted** — both
fields always print `|λ|max/|λ|min`. Correct for both matrix classes, never over-claims, no
glossary needed; matches INC-213, which had already recorded the κ₂ label as part of a
misreading.

## Plan: blocked, amended, cleared

`todo/mumps_error_analysis_opt_in.md`, per the template. Round 1 (terra/high + grok/high,
parallel, same brief): **both blocked, both independently found D1** — `arm_conditioning_thermal`
did two jobs behind one early return, and re-gating it whole would have silently dropped
`set_symmetric(true)`, sending every existing deck's thermal field through the nonsymmetric
ARPACK driver (reachable: the shipped `examples/3D_tapestack` deck is exactly that shape).
Plus: hoisting the footer tail was insufficient (D2), the thermal warning was structurally dead
at parse time — a *pre-existing* hole (D3, repaired this session on Christian's O4 ruling),
G1 was vacuous (D4), R8 asserted what O3 forbids (D5).

Christian's label ruling turned D1 from "accept or work around" into "design out":
`set_symmetric` moved to a `setup_thermal_eigen()` helper called from **both** attach paths,
and the arming functions became pure ICNTL(11) toggles.

Round 2 (both at xhigh): no blocking defect. O5 (does the thermal EigenValues survive attach?)
resolved safe by two independent lifetime traces — the object is built once in the DofManager
ctor and nothing rebuilds it. Round 2 also caught the plan contradicting itself: a value pin
cannot detect a lost `set_symmetric` (both drivers converge to the same ratio), so the real
tripwire is R12's `BELFEM_ASSERT( is_symmetric() )` — which is why `is_symmetric()` survives
R3c as a test-facing accessor.

## What landed (all in `cl_FEM_Controller.{hpp,cpp}` + docs)

- `mMumpsErrorAnalysis`/`2` + presence bools; three-site parse mirroring `compute conditioning`
- arm/capture re-gated onto the new flags in one edit (INC-294 discipline: `Full` only in the
  two arms, `None` only in the two captures)
- `setup_thermal_eigen()` — the one `set_symmetric(true)` write, both attach paths
- footer rebuilt: per-row gating, no outer if/else, single unconditional tail, κ₂ ternaries gone
- `check_magnetic_diagnostics()` / `check_thermal_diagnostics()` — no-MUMPS warnings for both
  keys plus the migration notice (fires on key *absence* for a MUMPS field; explicit false is
  informed silence); thermal checks moved to kernel-attach time, repairing the never-fired
  pre-existing warning
- R12 tripwire assert; input contract updated in both artifacts same session; timestepping
  strategy §4 rewritten with the post-split footer; coulomb_gauge_penalty_theory's
  "already-scaled matrix" claim corrected against the MUMPS source; example deck annotated with
  an explicit `false`

## Code audit

Both PASS (terra/high + grok/high), scoped **by symbol** — the shared checkout carried another
session's uncommitted work in the same file, so `git diff` would have audited the wrong code.
Post-round fixes: **C1** (both found it) — the thermal-diagnostics one-shot was
first-caller-wins, so `set_thermal_kernel` before `set_params` would have consumed it against
default flags; replaced with a both-prerequisites latch (`mParamsSet`), any call order
converges. **C2** (Grok) — three stale comments, one of which was the text most likely to talk
a future session into re-opening INC-294. **C3** — the R12 assert is debug-only, so gate G3b
needs `USE_DEBUG=ON`.

## O3 answered by a second plan, then ruled on

The last open question — whether the eigen path iterates the operator MUMPS factors — was
answered from an unexpected direction. Christian proposed reusing the production MUMPS solver for
the shift-invert diagnostic instead of building a second instance
(`todo/eigen_shift_invert_solver_reuse.md`). Its audit round found that **the matrices differ
under Newton**: `SolverData::jacobian()` returns `mSystemMatrix` (`SolverData.hpp:679`) while a
Newton body solves `mJacobianMatrix` (`:2428`) — separate allocations, `dJdx` added only to the
Jacobian. My draft had claimed "same matrix — confirmed"; I had trusted an accessor named
`jacobian()` to return the Jacobian and verified one hop short.

**Christian ruled it not a defect:** the ratio estimates problem difficulty, which is what that
diagnostic is for, while COND1 bounds the forward error of the solve that ran. That unblocked
**R8b**, which now explains the gap on two axes — different quantity, and different matrix under
Newton (coinciding under Picard) — a better answer than the material that raised the question.

That reuse plan itself came back **revise-before-code** from both auditors, with three further
draft claims refuted, and the recommendation reversed to Grok's alternative: `free()` the
diagnostic's own instance per call, which keeps the isolation and removes the *persistent* double
residency. Decisive find: the comment defending the retained instance cites a slot-exhaustion abort
that **DR-143 fixed on 2026-08-29** — a one-day-stale rationale. Left OPEN for Christian's call.

## Status

**All 16 steps landed.** Executed this session: two syntax checks (debug + release flags), the
scoped G1 grep, two Codex prose sweeps (each of which flagged a real technical error, not just
prose — the COND2 omission predicate understated `cond2_row_wanted`'s omega2==0 case, and
"well-conditioned right-hand side" is nonstandard), and Christian's own `tape_quench_usermat` run,
which verified three changes in real execution: R5a's single-field warning firing, R3c's unified
`|λ|max/|λ|min` label with no `κ₂` anywhere, and R4's box structure measured at 76 display columns.
Everything else is reviewed, not verified. Owed: Christian's deck runs only — G3 (four-combination matrix, successful first
iterate), G3b (pre-R10 fixture, **debug build**, since the R12 assert compiles out in release),
G3c (two MUMPS fields, per-field override), G3d (multi-iterate ICNTL(11) disarm), G4 (both
no-MUMPS warning paths), G5 (magnetic-only MUMPS). Ten fixtures, a runner that refuses a stale
binary, and pre-registered expectations are staged in `build/gates/`.

Nothing was committed: the tree carries at least three sessions' uncommitted work, and
`todo/mumps_workspace_cap_raise.md` touches the same `gMaxMemoryRelaxation` constant the reuse
plan would.
