# Split the MUMPS ADD Diagnostic Out of `compute conditioning` and Make It Opt-In

**Date:** 2026-08-30
**Purpose:** The `MUMPS ADD COND1 / COND2` footer rows are printed whenever `compute conditioning`
is on, so a user who wants the eigenvalue κ must also take the MUMPS Arioli–Demmel–Duff numbers —
and, symmetrically, cannot get the MUMPS numbers without paying for an ARPACK run. Split them:
a new `mumps error analysis` boolean, default false, gates the MUMPS half independently, and
`compute conditioning` keeps meaning exactly the eigenvalue estimate.
**Module:** `src/fem/kernel` (+ `doc/` input contract, `src/fem/doc/`)
**AIs involved:** Claude (exploration + plan + verification), Codex (naming audit, plan audit), Grok (naming audit + dissent, plan audit)
**Status:** PLAN — round-1 audited 2026-08-30 (Codex `gpt-5.6-terra`/high + Grok `grok-4.6`/high).
**BOTH auditors blocked implementation on the same P0**, found independently: D1 below. Amendments
applied the same day — R3 and R4 rewritten, R5 split, R8 split, gates rebuilt, five gap rows added.
Christian then ruled on D1 (2026-08-30): resolve it by **removal** — always print
`|λ|max/|λ|min` (R3c) and move `set_symmetric` to `set_thermal_kernel` (R3b) — rather than by
accepting the regression. No source modified. Key name pinned (Appendix A). O4 ruled 2026-08-30 → the pre-existing dead
thermal warning is repaired here, in R5b. **Round 2 (both auditors, `xhigh`, 2026-08-30) returned
NO blocking defect: "implementable for R1–R4", "I would start R1".** O5 resolved safe by two
independent lifetime traces; twelve gate and step refinements folded in, including one contradiction
of this plan's own reasoning (G3b). Approved 2026-08-30. **IN PROGRESS — R1, R2, R3, R3b, R3c, R4,
R5a, R5b, R9, R12 landed (10 of 12), D1/D2/D3 closed. Code-audit round 2026-08-30: BOTH auditors
PASS, no blocking defect; post-round fixes C1 (call-order latch, `mParamsSet`) and C2 (three stale
comments) applied and re-syntax-checked debug+release.** Docs landed same day: R6, R7, R8a, R10, R11 — 15 of 16
steps done. **CORRECTED 2026-09-03 (currentness sweep, Grok — verified in this file): this line
used to read "only R8b remains, blocked on O3", and both halves of that are false. R8b is `[x]`
landed 2026-08-30 (`:343`) and O3 is `[RESOLVED]` the same day, explicitly "R8b UNBLOCKED"
(`:463`). All 16 steps are done. What actually remains is the gate matrix — G3, G3b, G3c, G3d,
G4, G5 (`:408+`) — which is why this file is NOT a single-run close.** The definition-of-done
line at `:524` still says "O1 and O3 remain open" and is stale for O3 in the same way.
G1 run and green (scoped grep). G2 run and applied same day —
17 edits landed, 2 declined against the code, and the sweep caught a real claim error (the COND2
omission predicate understated `cond2_row_wanted`'s omega2==0 case; both docs fixed). Remaining:
gates G3–G5 and R8b (blocked on O3).

**Gate fixtures are prepared and pre-registered** (2026-08-30): ten run directories under
`build/gates/` derived from the tapestack3d deck (cold-start, 12 ms ≈ 2–3 steps; g3d 30 ms), a
runner `build/gates/run_gates.sh`, and `build/gates/EXPECTATIONS.md` written **before any run** —
per-gate PASS criteria plus explicit anti-expectations. The debug tree is `./build` on SCLS
`/opt/scls/debug`, so R12's assert is live and G3b's `USE_DEBUG=ON` requirement is met. **The
binary must be rebuilt first**: the 21:35 build predates the C1 latch fix (source 21:50), and the
runner refuses to execute against a binary older than `cl_FEM_Controller.cpp`.

A genuine **pre-change baseline** was found for G3b in
`cmake-build-debug/tapestack3d/out2.txt` (Aug 29): `κ₂ of Thermal System : 6.46e7 / 5.92e7` and
`|λ|max/|λ|min, Magnetic : 7.2e10 / 6.3e10`. Different run conditions, so it pins the **magnitude
band**, not digits — the driver identity itself is checked by R12's assert, not by the value. Syntax-checked clean; no gate has run, so nothing here is verified.
Next: the documentation steps. Steps R1–R12, gates G1–G5.

> **Scope guards:**
> - The **name is settled** and is not reopened: `mumps error analysis` (decided 2026-08-30,
>   Christian; see Appendix A). Renaming after the end-of-August-2026 release breaks user decks.
> - **OUT of scope:** adding PETSc/SuperLU backends for the ADD numbers — that is
>   `todo/conditioning_diagnostic_backends.md`, and it stays open.
> - **OUT of scope:** changing what either diagnostic *computes*. This is a gating and
>   documentation change only. No numerical behaviour moves.
> - **IN scope:** the new key, its parse + warning, the `print_footer` re-gating, both
>   input-contract artifacts, and the two documents that currently conflate the two quantities.
> - Compatibility: an existing deck with `compute conditioning : true` keeps its eigen κ rows and
>   **loses** its `MUMPS ADD` rows. That is the intended, announced behaviour change.

---

## 1. Current Behaviour and How It Fails

One flag pair (`mComputeConditioning`, `mComputeConditioning2`;
`cl_FEM_Controller.hpp:134-135`) gates two unrelated diagnostics:

| Diagnostic | Produced by | Cost per timestep per field | Slot |
|---|---|---|---|
| eigen κ | `EigenValues::compute_conditioning()`, `cl_FEM_DofMgr_EigenValues.cpp:929` | an ARPACK/PARPACK run — the expensive one, timed at `cl_FEM_Controller.cpp:3712-3721` | 0 |
| MUMPS ADD COND1/COND2/omega2 | `ICNTL(11)=1` armed around the first solve, `mumpstools.f90:404` | omega statistics + Hager estimator (ITMAX 5, twice if the second row category is non-empty); triangular solves against the existing factorization, no refactorization | 1, 2, 3 |

| Failure | Mechanism | Evidence |
|---|---|---|
| Unwanted output — the ADD rows cannot be suppressed without also losing κ | both are gated on the same flag | `cl_FEM_Controller.cpp:3184`, `:3259`, `:3294`, `:3748-3785` |
| Unwanted cost — the ADD numbers cannot be had without an ARPACK run | same flag; `compute_conditioning()` runs unconditionally when the flag is on | `cl_FEM_Controller.cpp:3025-3027`, `:3184-3200` |
| Documentation conflates the two quantities | the reference row describes MUMPS error analysis and the eigen estimate as two routes to one κ | `doc/input_file_reference.md:219`; `src/fem/doc/timestepping_strategy.md:154-177` |
| Stale footer diagram | the ASCII footer sketch predates the current row set | `src/fem/doc/timestepping_strategy.md:154-177` vs `cl_FEM_Controller.cpp:3736-3818` |

**Bottom line:** two independent diagnostics share one switch, so neither can be requested alone,
and the documentation has grown to describe them as one thing.

## 2. Architecture: Why a Second Independent Boolean

A second boolean parsed in the same three places, with the same override chain, is the smallest
change that makes both diagnostics separately requestable. Rejected alternatives:

- **Enum-valued `compute conditioning : none|eigen|mumps|both`.** More expressive, but breaks the
  existing `get_bool` parse and silently changes the meaning of every deck that already sets it.
- **Subordinate flag** (only effective when `compute conditioning` is also on). Smaller edit, no
  footer re-gating — but preserves exactly the coupling this plan exists to remove.
- **Delete the ADD rows.** The user asked for opt-in, not removal; and the numbers are the correct
  input to a forward-error bound for anyone who wants one.

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation |
|---|---|---|---|---|---|
| 1 | `mMumpsErrorAnalysis` / `…2` members | gating the MUMPS half | no — does not exist | (c) | new, beside `cl_FEM_Controller.hpp:134-135` |
| 2 | parse of `mumps error analysis` (solver shorthand + 2 per-field overrides) | deck contract | no | (c) | mirror `cl_FEM_Controller.cpp:4459-4477` |
| 3 | `arm_conditioning_magnetic` gate | stop paying ICNTL(11) when unwanted | keyed on the eigen flag | (c) | `:3259` |
| 3b | `arm_conditioning_thermal` — **two jobs behind one early return** | `set_symmetric(true)` belongs to the EIGEN flag; only `set_mumps_error_analysis` belongs to the new one | keyed on the eigen flag; a naive predicate swap breaks the eigen path | (c) | `:3270-3288`; see **D1** |
| 4 | `capture_conditioning_*` gate | must match the arm gate exactly, or slots stay NaN | keyed on the eigen flag | (c) | `:3294`, `:3322` |
| 5 | `compute_conditioning()` eigen calls | must NOT run when only the MUMPS half is asked | keyed on the eigen flag (correct already) | (a) | `:3184`, `:3193` — no change |
| 6 | `:3025` guard that calls `compute_conditioning()` | must not fire for a MUMPS-only request | keyed on eigen flags (correct already) | (a) | `:3025-3027` — no change |
| 7 | eigen-timing footer row | must not print when no eigen run happened | keyed on eigen flags (correct already) | (a) | `:3712-3721` — no change |
| 8 | κ footer rows | print on the eigen flag | keyed correctly, but nested inside the shared block | (c) | `:3733-3760`, `:3786-3797` |
| 9 | ADD footer rows | print on the new flag | nested inside the eigen-flag block | (c) | `:3762-3785`, `:3799-3822` |
| 10 | outer `if/else` around timestep-time + wallclock rows | must print in all four flag combinations | two branches with an identical tail | (c) | `:3733` vs `:3834-3846` |
| 11 | eigen-failure banner | must not fire for a MUMPS-only request | keyed on eigen flags (correct already) | (a) | `:3683-3708` — no change |
| 12 | no-MUMPS warning for the new key | the flag is a silent no-op without MUMPS | no | (c) | mirror `:4493-4516` |
| 13 | `doc/input_schema.yaml` entry | machine contract | no | (c) | beside `:248-251`, `:282` |
| 14 | `doc/input_file_reference.md` row | human contract | conflated | (c) | `:219`, `:183` |
| 15 | `src/fem/doc/timestepping_strategy.md` §4 | user-facing explanation of the two quantities | conflated + stale diagram | (c) | `:154-177` |
| 16 | text of the existing no-MUMPS warning | its closing paragraph presents MUMPS error analysis as what `compute conditioning` gives you — the conflation this plan removes | wrong after the split | (c) | `:4508-4512` |
| 17 | `doc/input_file_reference.md` §4 preamble | names `compute conditioning` as the *sole* key read in both places; after R2 there are two | wrong after the split | (c) | `:181-184` |
| 18 | `src/fem/maxwell/doc/coulomb_gauge_penalty_theory.md` §5.3 | attributes κ under `compute conditioning` to "MUMPS from its own error analysis" — third conflating user-facing doc | wrong today, wronger after | (c) | `:262-265` |
| 19 | stale code comments naming the wrong flag | send the next session to the wrong gate | wrong after R3 | (c) | `hpp:134`, `hpp:667-668`, `cpp:1252`, `cpp:3122-3127`, `cpp:4448-4451` |
| 19b | κ₂ label branch in `print_footer` | user-facing row label | branches on `is_symmetric()` | (c) | `:3753-3760`, `:3790-3794`; R3c |
| 19c | docs naming the footer row `κM` / `κT` | must match what the footer prints after R3c | stale after R3c | (c) | `doc/input_file_reference.md:219`, `src/fem/doc/timestepping_strategy.md:165` |
| 20 | `examples/3D_tapestack/input.conf` | shipped deck that silently loses its ADD rows | untouched | (c) | `:19` (magnetic MUMPS), `:48` (thermal PETSc), both `compute conditioning : true` |

### 3.1 Cross-cutting findings

- **Rows 3/3b and 4 must move together, and BOTH mismatch directions are permanent.** The disarm
  exists only in the capture path (`:3315-3316`, `:3333-3334`); the only `Full` writes are the two
  arm sites (`:3263-3264`, `:3287-3288`). Grok tabulated both diagonals and the plan originally
  described only one:

  | Arm gated on | Capture gated on | `eigen=off, mumps=on` | `eigen=on, mumps=off` |
  |---|---|---|---|
  | new flag | old flag | arms, capture returns → **stays Full** | arm skips; capture reads NaN |
  | old flag | new flag | arm skips; capture reads NaN | arms, capture returns → **stays Full** |

  "Stays Full" is not hypothetical: it is **INC-294** (`doc/lessons_learned_evidence.md:386`),
  where `MumpsErrorAnalysis::Full` was set once and never cleared and MUMPS ran error analysis on
  **every one of 30+ solves per timestep** while one value was read. An asymmetric edit re-opens a
  closed incident.
- **Row 3b is the P0 (D1).** `set_symmetric(true)` sits *between* the eigen guard and the MUMPS
  guard, and is the only write of `mSymmetric` in the tree.
- **Row 10 is the only structural edit, but hoisting alone is NOT sufficient** — see D2 and the
  rewritten R4.
- **Rows 5, 6, 7, 11 are genuinely class (a)** — confirmed independently by both auditors. They
  read no MUMPS state. Listed so the audit can see they were considered, not to schedule work.
- **MPI:** arm and capture are deliberately unguarded by rank (`cl_FEM_Controller.hpp:669-671`).
  The new bools are parse-time and rank-identical because `set_params` runs on every rank. **Do not
  rank-guard the new flag** — it would desynchronize ICNTL(11) across ranks.

### 3.2 Defects found in round-1 audit (2026-08-30, Codex + Grok)

- [x] **D1 — CRITICAL. A naive R3 silently breaks the thermal eigen diagnostic.**
  Found independently by **both** auditors. `arm_conditioning_thermal` (`:3268-3289`) does two
  jobs behind one early return: it sets the thermal `EigenValues` object symmetric, *then* arms
  MUMPS. `set_symmetric(true)` (`:3282`) is the only write of `mSymmetric`, whose default is
  `false` (`cl_FEM_DofMgr_EigenValues.hpp:180`); it selects `dsaupd`/`dseupd` over
  `dnaupd`/`dneupd` and is what earns the footer label `κ₂ of Thermal System` rather than
  `|λ|max/|λ|min, Thermal` (`:3790-3794`). Re-gating the whole function on the new flag makes that
  write a no-op whenever `mumps error analysis` is off — which, after this split, is **the default
  for every existing deck**. Root cause: the gap table treated one function as one job.
  **RULED 2026-08-30, Christian — resolved by removal, not by acceptance.** The label half is moot:
  R3c drops the κ₂ branch and always prints `|λ|max/|λ|min`, which is correct for both matrix
  classes and does not ask the reader to know κ₂. The driver half is fixed structurally by R3b,
  which moves `set_symmetric(true)` out of the arming path into `set_thermal_kernel`, so
  `arm_conditioning_thermal` no longer has a second job to forget. Net effect: no existing deck
  changes its thermal number, the label under-claims nowhere because it no longer claims κ₂
  anywhere, and R3 reduces to the straight predicate swap the draft wanted.
  **Fixed 2026-08-30** by R3b + R3c (Claude). `set_symmetric( true )` now lives in the new
  `Controller::setup_thermal_eigen()`, called from `set_thermal_kernel` *and* the two-argument
  constructor's `aKernel2 != nullptr` branch; `arm_conditioning_thermal` is a pure ICNTL(11)
  function carrying a comment that names this defect so it is not re-made. Verified: clean
  `g++ -fsyntax-only` against the project's real defines/includes. NOT verified: no gate has run.
  Gate owed: G3b.

  **Proof it is reachable, not theoretical:** `examples/3D_tapestack/input.conf` ships with
  magnetic MUMPS `compute conditioning : true` (`:19`) and thermal **PETSc**
  `compute conditioning : true` (`:48`). Thermal is not MUMPS, so today it returns at the MUMPS
  check *after* `set_symmetric` has run, and earns κ₂. Naive R3 changes the driver, the label and
  the number, with no gap row to catch it and no G3 combination that would see it.
  *Fixed by: the rewritten R3 + gate G3b.*
- [x] **D2 — HIGH. R4 as originally written does not produce the four combinations.**
  Found by Grok. The ADD rows are nested *three* deep: outer `if (eigen_m || eigen_t…)` at
  `:3733-3734`, inner `if (mComputeConditioning)` / `if (mComputeConditioning2 && mKernel2)` at
  `:3748` / `:3786`, then the MUMPS type check at `:3763` / `:3798`. Hoisting the duplicated tail
  and leaving those predicates yields `(eigen=off, mumps=on)` → outer false → **no ADD rows at
  all**, and `(eigen=on, mumps=off)` on a MUMPS field → **ADD rows still print**. Both are exactly
  the behaviours the split exists to create and remove. *Fixed by: the rewritten R4.*
- [x] **D3 — MEDIUM. R5 could not have worked as specified.** Found by Codex, corroborated by
  Grok. `set_params` runs with `mKernel2 == nullptr` (`cl_MaxwellFactory.cpp:983`; `belfem.cpp:144`
  then `:166`), and the existing thermal warning branch requires `mKernel2 != nullptr`
  (`:4497-4500`). So the *existing* thermal no-MUMPS warning is already dead on the coupled path,
  and mirroring it would copy the hole. `set_thermal_kernel` (`:4595-4623`) forwards headroom and
  compression but not this warning. **Fixed 2026-08-30** by R5a/R5b (Claude): `check_magnetic_diagnostics()` runs at `set_params`
  where the magnetic kernel exists; `check_thermal_diagnostics()` runs from **both** attach paths
  (`set_thermal_kernel`, and `set_params` when the two-argument constructor already supplied the
  kernel), one-shot so a caller reaching both sites warns once. The pre-existing dead eigen warning
  is repaired in the same move — its thermal branch is stripped from the parse-time block. Reviewed,
  not verified — gates owed: G4.
- [x] **D4 — LOW. G1 was not a gate.** Found by both. `scripts/check_doc_claims.py` reads
  `CLAUDE.md`, `doc/coding_philosophy.md` and a debt-register count (`:41-47`) — never
  `doc/input_schema.yaml` or `doc/input_file_reference.md`. A green G1 after R6/R7 would have
  proven nothing. **Fixed 2026-08-30**: G1 replaced by the scoped grep; run same day — token present in schema (4), reference (3), controller (14), and both table rows state the split.
- [x] **D5 — LOW. R8 asserted what O3 forbids.** Found by Grok. R8 said "explaining why the two
  numbers differ"; O3 says that must not be asserted until the operator question is settled.
  **Fixed 2026-08-30**: R8a landed without a causal paragraph (the measured pair is quoted as a measurement only); R8b remains blocked on O3.

**Retracted from the draft:** the claim that the `if`/`else` tails are "verbatim" copies. Both
auditors checked: they differ in comments and in `"┘"<<` vs `"┘" <<`. **Semantically** identical,
so the hoist stands; the word was wrong.

## 4. Ordered Steps

> **AMENDED 2026-08-30** after the round-1 audit. Superseded step text is struck through rather
> than deleted, per the plan conventions.

- [x] **R1** — Add `mMumpsErrorAnalysis` / `mMumpsErrorAnalysis2` (default `false`) beside
      `cl_FEM_Controller.hpp:134-135`. The comment must state they gate `ICNTL(11)` **only** and are
      independent of `mComputeConditioning`; it must not repeat `hpp:667-668`'s claim that all four
      arm/capture functions key off the eigen flags, which R3 makes false.
- [x] **R2** — Parse `mumps error analysis` in `set_params`, mirroring `:4459-4477` exactly:
      solver-section shorthand sets both, `linear magnetic` / `linear thermal` override per field.
      *(after: R1)*
- [x] **R3** — Re-gate arming and capture. **~~Re-gate `arm_conditioning_magnetic` /
      `arm_conditioning_thermal` and `capture_conditioning_*` onto the new flags.~~**
      *(struck 2026-08-30: a straight predicate swap on the thermal arm is D1.)*
      **~~Split the two jobs inside `arm_conditioning_thermal`, keeping `set_symmetric(true)` under
      `mComputeConditioning2`.~~** *(struck 2026-08-30, superseded by R3b — Christian's ruling on
      the label makes the cleaner move available.)* The required shape:

      - **R3b first** (below) moves `set_symmetric(true)` out of the arming path entirely.
      - `arm_conditioning_magnetic` and `arm_conditioning_thermal` then become **pure ICNTL(11)
        functions** and both are a straight predicate swap onto `mMumpsErrorAnalysis` /
        `mMumpsErrorAnalysis2`. No second job survives in either.
      - `capture_conditioning_magnetic` / `capture_conditioning_thermal` — same swap. **Must land in
        the same edit as the arm sites** (§3.1: either mismatch direction leaves `ICNTL(11)` armed
        for the rest of the run — INC-294).
      - Do **not** add a rank guard.

      *(after: R1, R2, R3b)*
- [x] **R3b** — Move `set_symmetric( true )` out of `arm_conditioning_thermal` (`:3282`) into
      `set_thermal_kernel` (`:4595-4623`), beside the existing late-attach block that already
      repeats the headroom and compression gates `set_params` could not run against a null kernel.
      Unconditional — it selects the ARPACK driver and is inert unless the eigen diagnostic runs, so
      it needs no flag. This is what makes R3 safe by construction rather than by comment: after it,
      no MUMPS-arming function carries a side effect on the eigen path.

      **VERIFIED round 2 (O5): safe.** The object is created once in the `DofManager` constructor
      and never rebuilt — full trace in O5. Confidence upgraded medium → high.

      Delete `:3282` **and** its "set on EVERY call" comment (`:3279-3281`) together: that comment
      defends a per-call write against a one-shot *guard bool*, which R3b does not introduce, so
      leaving it would mis-explain the new code.

      **Also write it in the two-argument constructor's `aKernel2 != nullptr` branch**
      (`cl_FEM_Controller.cpp:71-91`). `Controller( Kernel *, Kernel * = nullptr )`
      (`hpp:466`) can attach a thermal kernel without ever reaching `set_thermal_kernel`. It has no
      production caller today — `MaxwellFactory::create_controller` passes the magnetic kernel alone
      (`cl_MaxwellFactory.cpp:981-983`) — but it is public API, and both auditors flagged that a
      setter-only write reintroduces D1 there. The constructor already dual-paths `set_soft_fail`
      for exactly this reason (`:111-114`); mirror that. The write is idempotent, so a deck that
      takes both paths is harmless.

      **MPI:** `set_thermal_kernel` runs on every rank and `mSymmetric` needs no broadcast
      (`cl_FEM_DofMgr_EigenValues.cpp:770-772`). **Do not rank-guard this write.** *(after: R1)*
- [x] **R3c** — Drop the κ₂ label branch. `print_footer` currently selects between
      `"κ₂ of Magnetic System"` and `"|λ|max/|λ|min, Magnetic"` on `is_symmetric()` (`:3753-3760`,
      and the thermal twin `:3790-3794`). **Always print `|λ|max/|λ|min`** (decided 2026-08-30,
      Christian): it is correct for symmetric and nonsymmetric matrices alike, it never over-claims,
      and it does not ask a reader to know what κ₂ means. This also matches **INC-213**
      (`doc/lessons_learned_evidence.md:305`), which recorded the κ₂ label as itself part of a
      misreading — `|λ_max|/|λ_min|` equals κ₂ only for a normal matrix, and the h-φ Jacobian is not
      even symmetric. Removes the footer's only dependency on `is_symmetric()`. Both format strings
      already exist at the correct column width, so this is a deletion, not a re-layout.

      **Do NOT also delete `is_symmetric()` or `set_symmetric()`.** Both auditors confirmed the
      footer is the accessor's only caller, so it becomes unused by production code — but
      `mSymmetric` still selects the ARPACK driver, serial (`:487`, `:527`) and PARPACK (`:773`,
      `:816`) alike, so `set_symmetric` stays load-bearing and R3b is not dead code. The getter
      survives because **G3b needs it**: after R3c the printed value can no longer distinguish the
      two drivers, so the accessor is the only handle a gate has on whether R3b worked. It moves
      from display-facing to test-facing. *(after: R1)*
- [x] **R4** — Re-gate `print_footer`. **~~Hoist the shared tail out of the `if/else` and gate the
      κ rows on the eigen flags and the ADD rows on the new flags.~~** *(struck 2026-08-30: the
      hoist alone leaves both wrong behaviours — D2.)* Four edits, all required:

      1. **Drop** the outer `if` at `:3733-3734` entirely — the hoisted tail already covers the
         all-off case, so no widened predicate is needed.
      2. **Split** the inner `if (mComputeConditioning)` (`:3748`) and
         `if (mComputeConditioning2 && mKernel2)` (`:3786`) so each wraps **only** its κ row.
      3. **Gate** the ADD rows on the new flags, *keeping every existing guard*:
         `mKernel2 != nullptr` (thermal), `solver() != nullptr`, `type() == SolverType::MUMPS`,
         and `cond2_row_wanted` for the COND2 row.
      4. **Hoist** the tail once and delete both copies.

      **Preserve the current row order** — physics stats → failure banner → eigen timing →
      postprocessing → magnetic κ/ADD → thermal κ/ADD → tail. Emitting ADD after the box close, or
      leaving it nested under the eigen `if`, is D2 again.

      *(after: R1, R2, R3, R3c)* — **R3 is required**: if R4 lands first, `(0,1)` prints ADD rows as
      `n/a`, because nothing has armed or captured them yet. R3c is required so R4 does not restore
      the ternary it deletes.
- [x] **R5a** — Magnetic no-MUMPS warning at parse time, symmetric to `:4493-4516`: warn when
      `mumps error analysis` is on for a magnetic field whose solver is not MUMPS. *(after: R2)*
- [x] **R5b** — Thermal warnings, **both of them**, from `set_thermal_kernel` (`:4595-4623`)
      rather than `set_params` — the thermal kernel does not exist at parse time (D3).

      1. **New key:** warn when `mumps error analysis` is on for a thermal field whose solver is not
         MUMPS.
      2. **Repair the pre-existing hole** (O4, ruled in scope 2026-08-30, Christian): the *existing*
         eigen warning's thermal branch requires `mKernel2 != nullptr` (`:4497-4500`) inside a block
         that runs while `mKernel2` is still null on the coupled path, so **it has never fired for a
         thermal field**. Move that branch to `set_thermal_kernel` too.

      **Emit two independent per-field messages.** The existing message composes both fields into
      one `sprint` — `"magnetic and thermal"` / `"magnetic"` / `"thermal"` with a matching plural
      (`:4505-4515`) — and that cannot survive the split. ~~Or defer the whole message until both
      kernels are known.~~ *(struck 2026-08-30, round 2: the deferred option is unimplementable —
      there is no call site where both kernels are known. `hphirun` never calls
      `set_thermal_kernel` at all, and `set_params` cannot know whether one will be attached
      later.)* Two messages is also semantically the more accurate shape: the fields genuinely can
      use different libraries, which is the reason the existing comment gives for warning per field
      in the first place.

      **Strip the thermal half out of the `set_params` block** (`:4497-4500`) rather than leaving it
      in place. It is harmless on production paths, where `mKernel2` is null, but it would
      double-fire on the two-argument constructor path.

      **Ordering trap:** the two-argument constructor attaches *before* `set_params` runs, so a
      warning fired only at attach would read default flags there. Put the check in a small helper
      and call it from whichever of the two operations completes second.

      Note for the record: this repair means some decks will start seeing a warning they have never
      seen, because it was structurally unreachable. That is the point, but it should be expected
      rather than treated as a regression. *(after: R2)*
- [x] **R6** — Update `doc/input_schema.yaml`: new key beside `:248-251`, and the
      "may also appear inside linear magnetic / linear thermal" note at `:282`. Anchor by
      searchable token (`'"mumps error analysis"'`), never by line number. *(after: R2)*
- [x] **R7** — `doc/input_file_reference.md`: rewrite `:219` so it stops describing the two
      diagnostics as two routes to one κ, add the sibling row, **and amend the §4 preamble at
      `:181-184`**, which names `compute conditioning` as the sole key read in both places — after
      R2 there are two. Per Grok, the new row's first sentence must contain the footer strings
      `MUMPS ADD COND1` / `MUMPS ADD COND2` verbatim, so a user reading the footer can grep to the
      key. **Also drop `κM` / `κT` from the row text** — after R3c the footer never prints those
      (gap row 19c). *(after: R6, R3c)*
- [x] **R8a** — `src/fem/doc/timestepping_strategy.md:154-177`: un-conflate the two quantities and
      refresh the footer diagram **against the post-R4 code**, not against `:3736-3818`, which R4
      moves. The diagram at `:165` still shows the old two-cell `κM … κT` layout and must be redrawn
      with the current one-row-per-quantity shape and the `|λ|max/|λ|min` labels (gap row 19c).
      Unblocked. *(after: R4, R7, R3c)*
- [x] **R8b** — Add the causal explanation of why the two numbers differ by orders of magnitude.
      ~~**BLOCKED on O3.**~~ *(unblocked 2026-08-30: O3 resolved by ruling.)* **Landed
      2026-08-30** — two axes stated (different quantity; different matrix under Newton, coinciding
      under Picard), no causal claim beyond what O3 settled. The measured pair (κ₂ = 2.96e7 vs COND1 in 1e4–1e5) may be quoted as a
      measurement — it already is, at `cl_FEM_Controller.cpp:3181-3183` — but no *why* until O3 is
      settled. *(after: R8a, O3)*
- [x] **R9** — Migration notice, **per field, at two call sites**. Rank-0 `InfoLevel::Minimal`,
      same channel as `:4502-4515`, one line naming the new key. Fires when `compute conditioning`
      is on **and** `mumps error analysis` is absent **and** that field's solver is MUMPS. An
      explicit `false` is informed silence; an absent key is not.

      **~~A single check in `set_params`.~~** *(struck 2026-08-30, round 2 — two independent
      reasons.)*

      1. **Presence is lost by R2.** `Section::key_exists` does distinguish absence from an explicit
         `false` — `create_key` inserts either way (`cl_Input_Section.cpp:104-124`) and `key_exists`
         tests membership (`:225-227`); the house already uses this shape for `restart`
         (`cl_FEM_Controller.cpp:4434-4436`). But R2's member bools throw the provenance away
         (default `false` and explicit `false` are the same bit), and `Controller` does not retain
         `aSection`. **R1/R2 must therefore record two presence bools**, one per field, each true
         iff *either* the solver-section shorthand *or* that field's `linear` override carried the
         key. Testing only the shorthand false-positives on
         `linear magnetic { mumps error analysis : false ; }`.
      2. **The thermal solver does not exist at `set_params`** — D3 again. So: magnetic R9 in
         `set_params` against the live magnetic solver; thermal R9 in `set_thermal_kernel` against
         the live thermal solver, using the stored presence bit.

      Do **not** implement the predicate as `!mMumpsErrorAnalysis` — that warns on an informed
      `false`. Do **not** read the `library` string as a substitute for `solver()->type()`.
      *(after: R2, R5b — shares R5b's helper and its call sites)*
- [x] **R10** — Fix the shipped example. `examples/3D_tapestack/input.conf:19,48` is the exact
      shape that regresses; add the key (value Christian's call) and a comment so the shipped
      example is not a trap. *(after: R2)*
- [x] **R11** — Correct the stale comments in gap row 19 (`hpp:134`, `hpp:667-668`, `cpp:1252`,
      `cpp:3122-3127`, `cpp:4448-4451`) and the warning text in gap row 16 (`:4508-4512`), plus the
      third conflating doc at `coulomb_gauge_penalty_theory.md:262-265` (gap row 18).
      *(after: R3, R4, R5b)* — `hpp:667-668` is the arm/capture comment R3 invalidates, and
      `:4508-4512` sits inside the block R5b splits, so both must land first or R11 rewrites text
      that is about to move.

- [x] **R12** — Add the R3b tripwire in the thermal branch of `compute_conditioning`
      (`:3193-3202`):

      ```cpp
      BELFEM_ASSERT( mKernel2->dofmgr()->eigen_values()->is_symmetric(),
          "thermal EigenValues lost mSymmetric after attach" );
      ```

      `BELFEM_ASSERT`, not `BELFEM_ERROR`: this is an invariant check on a logic bug, and the eigen
      branch is per-timestep rather than setup code. It is the only cheap thing that can observe a
      failed R3b — see G3b. *(after: R3b, R3c)*

**Gates.** None are satisfied by the build alone, and none are optional.

- [x] **G1** — ~~`scripts/check_doc_claims.py` clean.~~ *(struck 2026-08-30 — D4: that script does
      not read either input-contract artifact.)* Replaced by a grep assertion: the token
      `mumps error analysis` appears in `doc/input_schema.yaml`, in `doc/input_file_reference.md`
      and in `cl_FEM_Controller.cpp`, and `compute conditioning` is no longer described as the
      switch for the MUMPS ADD rows **in the files this change owns** — the schema,
      `input_file_reference.md`, `timestepping_strategy.md`, `coulomb_gauge_penalty_theory.md` and
      the controller comments. Scoped deliberately: an unqualified "anywhere" would fail on this
      plan file and on INC-213, which are records of what was true when written and are not
      rewritten (round 2).
- [x] **G2** — Codex language sweep over R7/R8a prose (`gpt-5.6-terra`, `medium` — dense document).
      `timestepping_strategy.md` is user-facing and IS swept; this plan file is not.
- [ ] **G3** — Deck matrix, all four `(compute conditioning, mumps error analysis)` combinations on
      a **MUMPS** deck. Require **finite COND1** on `(0,1)` and `(1,1)` — not merely that the row
      label exists. A label check would pass on `n/a` while `ICNTL(11)` stayed armed; finite COND1
      proves capture ran, and capture is what disarms. Require the *Time for eigenvalue analysis*
      row **absent** on `(0,1)`, not merely small. State the expected presence/absence of **every**
      row for all four cases, not just those two. **Use a timestep whose first solve succeeds**: the
      failed-first-solve path bypasses capture entirely (`:1228-1241`, thermal twins `:1653-1659`,
      `:2418-2421`), so a step that fails its first solve would leave `ICNTL(11)` armed for reasons
      unrelated to this change and mis-attribute the result.
- [ ] **G3b** — **The shipped shape, and the gate for D1/R3b.** Run a **fixture copy** of
      `examples/3D_tapestack` in its **pre-R10** form: magnetic MUMPS + thermal PETSc,
      `compute conditioning : true` on both, new key **absent**. Not the shipped file itself —
      R10 adds the key to it, and after R10 either value breaks this gate (`false` → R9 must not
      fire; `true` → magnetic ADD still prints). Caught in round 2.

      Assert: the magnetic ADD rows are gone; R9's magnetic notice fired; the thermal value is
      unchanged against a pre-change run.

      **The value pin is a "nothing moved numerically" check, NOT proof that R3b worked.**
      ~~the gate that proves the `set_symmetric` move did not silently drop the Lanczos driver~~
      *(struck 2026-08-30, round 2 — and it contradicted this plan's own finding that the number
      does not change: both ARPACK drivers converge to the same spectral ratio on a symmetric
      matrix, so a lost `set_symmetric` prints an identical value.)* The actual driver tripwire is
      **R12**'s `is_symmetric()` assertion, which is why R3c keeps the accessor. **G3b must run a
      `USE_DEBUG=ON` build** — `BELFEM_ASSERT` compiles out in release, so a release G3b cannot
      observe a lost `mSymmetric` (code audit, C3).
- [ ] **G3c** — Per-field override, on a deck with **two MUMPS fields** (a non-MUMPS thermal field
      cannot prove isolation, because its ADD rows are suppressed by the solver-type guard whatever
      the flag says — round 2). Solver-section shorthand on, with
      `linear thermal { mumps error analysis : false ; }`. Assert magnetic ADD present and thermal
      ADD absent. Proves R2 precedence and that overrides do not leak between fields.
- [ ] **G3d** — Multi-iterate: a MUMPS deck at `(0,1)` or `(1,1)` on a timestep with more than one
      Newton/Picard iterate. The footer cannot show whether `ICNTL(11)` stayed armed after iterate
      0 — this is INC-294's actual cost. Confirm later iterates are not running error analysis
      (MUMPS `ICNTL(4)` chatter, or first-vs-later solve time).
- [ ] **G4** — Non-MUMPS decks for the warnings: a **magnetic-only** STRUMPACK deck with the new
      key true (proves R5a at parse time), **and** a coupled deck with non-MUMPS thermal, with
      **both** flags on, asserting **both** distinct thermal warnings appear — the new key's and the
      repaired eigen one (proves R5b from `set_thermal_kernel`; the parse-time path structurally
      cannot cover either, D3).
- [ ] **G5** — Magnetic-only MUMPS deck with the new key on and `mKernel2 == nullptr`. Cheap, and
      it is the null-kernel guard on the thermal ADD rows.

## 5. Open Design Questions

- **O1 — Does the new key belong in the `solver` section at all, or only per-field?**
  `compute conditioning` is read in both, and R2 mirrors that. Logged because the mirror is being
  copied on the strength of consistency rather than a fresh argument.
- **O2 — RESOLVED 2026-08-30 → unreachable; no defensive disarm needed.** Both auditors
  independently traced the write set: the only writes of `mComputeConditioning*` are the member
  defaults (`cl_FEM_Controller.hpp:134-135`) and the `set_params` parse (`:4461-4476`), and
  `set_params` is called once, from `MaxwellFactory::create_controller`
  (`cl_MaxwellFactory.cpp:983`). `load_memdump` (`:4874-4950`) restores Δt, BDF state and fields —
  not these bools. Confidence upgraded **medium → high**. Recorded as an implementation invariant,
  not a language guarantee: the members are mutable, so a future mid-run write would re-open it.
- **O3 — RESOLVED 2026-08-30 → answered, and ruled not an issue (Christian). R8b UNBLOCKED.**
  The question was whether the eigen path iterates the same operator MUMPS factors. **It does not,
  under Newton**: `EigenValues` takes `mParent->jacobian()`, which returns `mSystemMatrix`
  (`cl_FEM_DofMgr_SolverData.hpp:679`), while a Newton body solves `mJacobianMatrix`
  (`cl_FEM_DofMgr_SolverData.cpp:2428`) — separate allocations (`:396`, `:405`) with the `dJdx`
  tangent term added only to the Jacobian (`assemble_newton`, `:1683-1709`). Found during the
  `eigen_shift_invert_solver_reuse` audit by both auditors and verified independently.
  **Christian's ruling: expected, not a defect** — `|λ|max/|λ|min` on the system matrix is a good
  estimator of problem difficulty, which is what it is for, while COND1 bounds the forward error of
  the solve that actually ran. So the two footer numbers differ on **two** axes: different quantity,
  and (under Newton) different matrix. That is now a statable fact rather than an open question,
  which is exactly what R8b needed.

- **O5 — RESOLVED 2026-08-30 → yes; R3b is safe. Both auditors traced it independently.**
  `DofManager`'s constructor is the only assignment of `mEigenValues` (`cl_FEM_DofManager.cpp:60`),
  its destructor the only delete (`:96`). It is constructed inside `Kernel::create_field`
  (`cl_FEM_Kernel.cpp:616`), which `ThermalFactory::create_thermal_kernel` calls before returning
  (`cl_ThermalFactory.cpp:161`) — so the object is live by the time any executable calls
  `set_thermal_kernel` (`belfem.cpp:166`, `hphiTrun.cpp:89`). Nothing rebuilds it afterwards:
  `EigenValues::reset()` clears only `mMatrixFlag` (`:1803-1806`), `DofManager::reset()` never
  touches `mEigenValues` (`:157-168`), `link_matrix()` rebinds the matrix rather than the object
  (`:74-126`), a solver swap replaces only `SolverData::mSolver`, and `load_memdump`'s
  `initialize(true)` no-ops on an initialized manager. **The existing code already assumed this
  model**: `cl_FEM_DofMgr_EigenValues.cpp:770-772` says `mSymmetric` "is set once at setup, never
  derived from matrix values, so it needs no broadcast". The per-timestep write was belt-and-braces
  against a rebuild that does not exist.

- **O4 — RESOLVED 2026-08-30, Christian → repair in this plan.** The pre-existing eigen thermal
  no-MUMPS warning has never fired on the coupled path (D3): its branch requires
  `mKernel2 != nullptr` (`:4497-4500`) but sits in a block that runs while `mKernel2` is still null.
  Folded into **R5b**, which now moves both thermal warnings to `set_thermal_kernel`. Expect decks
  to start emitting a warning that was previously unreachable.

## 6. The Key Contract

```
solver
{
    compute conditioning : true ;    // eigen κ₂ / |λ|max/|λ|min   -> slot 0, ARPACK cost
    mumps error analysis : true ;    // MUMPS ADD COND1/COND2      -> slots 1-3, ICNTL(11) cost

    linear magnetic { mumps error analysis : false ; }   // per-field override wins
}
```

| Property | Value |
|---|---|
| Type | bool |
| Default | `false` |
| Sections | `solver` (shorthand for both fields), `solver.linear magnetic`, `solver.linear thermal` |
| Precedence | per-field key overrides the solver-section shorthand |
| Schema anchor | `'"mumps error analysis"'` |
| Effect without MUMPS | none; warns — magnetic at parse time (R5a), thermal at kernel attach (R5b) |
| Independent of | `compute conditioning` — either may be set alone |

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row (1–20) mapped to a step, or marked (a) "already correct, no edit".
- [ ] Each claimed gap backed by a citation, not assumption.
- [ ] Ordered steps with dependencies — note R3/R4 depend on **R2**, not R1 alone.
- [ ] D1–D5 each closed with a "Fixed YYYY-MM-DD" line naming what changed and who verified it.
- [ ] Open questions logged, not decided — **O1 and O3 remain open**; O2, O4 and O5 resolved 2026-08-30.
- [ ] G1–G5 all green, with the deck names and configurations recorded here.
- [ ] Round-2 audit of the amended steps returned before R1.

## 8. Audit Trail

- **Plan-audit round 1, 2026-08-30** — `tmp/ai_exchange/mumps_error_analysis_plan.md`. Codex
  `gpt-5.6-terra`/`high`, Grok `grok-4.6`/`high`, dispatched in parallel on an identical brief that
  told both the key name was closed. **Both blocked implementation, and both found D1
  independently** — the strongest signal available, since the brief did not point at
  `arm_conditioning_thermal`; it asked them to attack the class-(a) rows, and the defect was
  adjacent to one. Codex additionally found D3 and D4; Grok additionally found D2 and D5, produced
  the both-directions arm/capture table, tied it to INC-294, and found the shipped
  `examples/3D_tapestack` deck that makes D1 reachable. Everything above was re-verified against
  the cited code before amendment: `arm_conditioning_thermal:3270-3288`,
  `mSymmetric = false` at `cl_FEM_DofMgr_EigenValues.hpp:180`, `cl_MaxwellFactory.cpp:983`,
  `check_doc_claims.py:41-47`, `examples/3D_tapestack/input.conf:19,48`,
  `lessons_learned_evidence.md:386` (INC-294), `input_file_reference.md:181-184`,
  `coulomb_gauge_penalty_theory.md:262-265`. No auditor claim was accepted unverified.
- **Code-audit round 1, 2026-08-30** — `tmp/ai_exchange/mumps_error_analysis_code.md`. Codex
  terra/high + Grok grok-4.6/high, parallel, scoped **by symbol** because the shared checkout
  carries another session's uncommitted work in the same file. **Both PASS: no blocking defect.**
  Both confirmed: INC-294 not re-opened (`Full` only in the two arms, `None` only in the two
  captures, all four on the new flags); no third `mKernel2` attach path; footer gating correct in
  all four combinations with every pre-existing guard intact; migration presence bits correct per
  deck shape; `message()` calls safe. Post-round fixes applied and syntax-checked debug+release:
  **C1** (found by both) — the thermal diagnostics one-shot was first-caller-wins, so
  `set_thermal_kernel` before `set_params` would have consumed it against default flags; replaced
  with a both-prerequisites latch (`mParamsSet`), any call order now converges. **C2** (Grok) —
  three stale comments fixed, including the magnetic capture call-site text that could have talked
  a later session into re-gating capture onto the eigen flag and re-opening INC-294. **C3**
  (Grok) — recorded on G3b: the R12 assert is debug-only, so that gate needs `USE_DEBUG=ON`.
  Grok's residual open risks 3–5 are pre-existing shapes, LOW, recorded here rather than patched:
  a null solver at check time is a silent no-op; thermal κ prints without a solver guard (matches
  the old inner `if`); COND2 prints n/a on a failed first solve (G3 already requires a successful
  first iterate).
- Naming round: `tmp/ai_exchange/mumps_cond_key_name.md` — jury (parallel blind), 2026-08-30.
  Codex `gpt-5.6-terra`/`high`, Grok `grok-4.6`/`high`. Split 1–1; see Appendix A.
- Claude's pre-registration was sealed outside the exchange file (the wrappers instruct auditors to
  read the thread first, so writing it there would have biased the round) and was folded into the
  thread only after both verdicts landed.
- All auditor citations were re-verified against the cited code before inclusion. Two corrections:
  Codex's `mumpstools.f90:400` for ICNTL(11) is actually `:404`; and Codex's portability rebuttal
  cited SuperLU `dgscon` (norm-based) while the same section
  (`todo/conditioning_diagnostic_backends.md:136-140`) says `dgsrfs` returns the ADD quantities
  proper — so the portability concern was stronger than Codex allowed.
- **Plan-audit round 2, 2026-08-30** — same thread. Codex `gpt-5.6-terra`/`xhigh`, Grok
  `grok-4.6`/`xhigh`, parallel, on a brief that named the closed items (key, label ruling, O2, O4)
  and put **O5** first. **Neither returned a blocking defect.** Codex: "block only for three precise
  plan amendments"; Grok: "implementable for R1–R4… I would start R1". Both traced O5 to the same
  answer independently, and both found the two-argument `Controller` constructor hole, the R9
  provenance problem, and that G3b could not prove what it claimed.

  Round-2 findings unique to one auditor: **Codex** — R4 needs `after: R3` or an intermediate build
  prints unarmed ADD rows; the constructor attaches *before* `set_params`, so R5b's helper must fire
  from whichever operation completes second; G3c needs two MUMPS fields. **Grok** — G3b as written
  *fights R10*, since R10 adds the key to the very file G3b requires it absent from, so G3b needs a
  pre-R10 fixture; R5b's deferred-message option has no call site and is unimplementable; G1's
  "anywhere" wording would fail on this plan file and INC-213 and must be scoped; the
  failed-first-solve path bypasses capture, so G3 must use a successful first iterate.

  **One round-2 finding contradicted this plan's own reasoning and was accepted:** G3b claimed to
  prove R3b by pinning the thermal value, while the plan elsewhere establishes that the value does
  *not* change when the driver does. Both auditors caught it. The driver tripwire is now R12's
  `is_symmetric()` assertion, and R3c keeps the accessor for it.

  Re-verified before folding in: `Controller( Kernel *, Kernel * = nullptr )` at
  `cl_FEM_Controller.hpp:466` and its body `:71-91` (attaches `mKernel2` without
  `set_thermal_kernel`); `is_symmetric()` has exactly two callers, both the footer branches R3c
  deletes; `mEigenValues` assigned only at `cl_FEM_DofManager.cpp:60`.

---

## Appendix A — Decision: the key is named `mumps error analysis`

**Question.** What to call the boolean that gates `ICNTL(11)`, pinned permanently because the
release is end of August 2026 and renaming a deck key breaks every user's `input.conf`.

**Candidates and verdicts.** Both auditors independently ranked the three offered names identically,
worst to best: `compute cond1` < `compute mumps conditioning` < `mumps error analysis`.

- `compute cond1` — names one of two printed numbers and none of the switch; freezes the current
  print set into the deck contract while `get_forward_error()` (`cl_SolverMUMPS.cpp:1940-1948`) is
  already wrapped and unused. No numbered-jargon precedent in the deck.
- `compute mumps conditioning` — contains the very word the split exists to separate. Both auditors
  independently judged it the *harmful* middle option for that reason, not the safe compromise.
  Both also found `compute ` is a **singleton**, not a house prefix: `compute conditioning` is the
  only key in the deck that uses it.
- `mumps error analysis` — names the switch, matches the internal vocabulary exactly
  (`MumpsErrorAnalysis` in `en_SolverEnums.hpp`, carrying `/** MUMPS only */`;
  `set_mumps_error_analysis` in `cl_Solver.hpp:160-161`), and covers COND1 + COND2 + omega2 because
  those are outputs of error analysis rather than "a condition number".

**The dissent, and why it did not win.** Grok proposed a fourth name, `error analysis`, dropping the
vendor, on the ground that the house convention is capability-in-the-key and vendor-in-the-note —
verified: `matching` is STRUMPACK MC64 and is not called `strumpack matching`; `preconditioner`,
`krylov method` and `initial guess` are PETSc-only and carry no vendor; `metis nodendp`
(`doc/input_schema.yaml:353`) is the lone vendor-named key and exists because NodeND-vs-NodeNDP *is*
the distinction. Grok argued that `todo/conditioning_diagnostic_backends.md` being open (R3/R4
unchecked, O1 requiring one quantity under one name across backends) means a `mumps `-prefixed key
decides by accident that SuperLU and PETSc get a second key.

Rejected on a ground **neither auditor weighed**: the flag is a *no-op without MUMPS today*.
`arm_conditioning_magnetic` and `arm_conditioning_thermal` both return early on
`solver()->type() != SolverType::MUMPS` (`cl_FEM_Controller.cpp:3259-3262`, `:3272-3274`). A key named
plain `error analysis` would therefore sit inert in any STRUMPACK deck — the project's first-choice
solver — with nothing in its name to say why. Grok traded a real, present silent no-op against a
hypothetical future rename. If SuperLU later gains `dgsrfs`, a second key can be added and this one
stays accurate forever, whereas a generic key pinned today would have to be retrofitted with
per-backend semantics.

| | Verdict |
|---|---|
| Pinned name | `mumps error analysis` |
| Decided | 2026-08-30, Christian |
| Jury | Codex: same name, high confidence. Grok: dissent (`error analysis`), medium-high |
| Deciding argument | the flag is a no-op without MUMPS; the name must say so |
| Kept from the dissent | the R7 wording constraint (footer strings verbatim in the first sentence) and the R5 no-MUMPS warning |
