# DR-39 closing plan made implementation-ready, then audited

**Date:** 2026-08-30
**Purpose:** Record the second DR-39 ruling (close by executing (a) and re-homing (b)), the eight
defects found in the harvest plan's own text before dispatch, and the jury round that found nine
more — three of them P0.
**Module:** `todo/circuit_demo_harvest_plan.md` (+ `tests/circuit`, `src/executables`, docs — none
touched yet)

## The ruling

Christian, asked what "close DR-39" covers: **execute (a), re-home (b), strike the row.** (b)
`Pulse`/`PWL` is a feature, not debt — it has a home at `todo/ngspice_parser_plan.md` Phase 6 and a
user-visible behaviour today (the parser names and refuses both). (c) was discharged 2026-08-27.

Christian also pinned R4 at Δt = 50 µs, and asked whether PWL is even in BELFEM's scope. The answer
went into the plan as a recommendation to Phase 6: **`PULSE` earns its place, `PWL` does not on
current evidence.** `PULSE(v1 v2 td tr tf pw per)` is a trapezoid — ramp, plateau, ramp, repeat —
which is the standard HTS magnet duty cycle and the standard AC-loss cycle, and BELFEM can only
approximate it today with `Ramp`/`Sigmoid`. `PWL` is the general breakpoint table whose only magnet
use is replaying a measured current profile, for which no BELFEM deck syntax exists either — and it
is the leg that costs most (variable-length payload ⇒ `share`/`receive`, `Cell< real >` only under
the plugin-ABI ban). Recorded as a recommendation, not a decision.

## Eight defects found in the plan before any auditor saw it

The two that would have cost most were both in R4's switching gate:

- **The pre-switch amplitude window was half a period.** At 50 Hz the period is 20 ms = 400 steps at
  Δt = 50 µs, not 200. An amplitude gate over a half period can miss the peak depending on phase.
- **The window included a post-switch sample.** `shift()` fires the latch *before*
  `compute_MNA_matrix()` stamps, so the solution recorded at the end of step 500 already belongs to
  the closed-switch regime. The window must stop at 499.
- **`t_switch` at exactly 25 ms is a floating-point tie, not determinism.** `mTime` accumulated over
  500 additions of `50e-6` differs from `1.25/50` at ~1e-18 with unpredictable sign, so the latch
  fires at step 500 or 501 depending on the build — and possibly *differently* in the hand-built and
  netlist runs, which would break R6's equivalence gate for a reason having nothing to do with the
  parser. Moved to 24.975 ms, half a step early; the physical firing instant is unchanged at 25.0 ms.

Also: R1′ was mandatory on a premise that does not survive reading the code (the demo's own ω update
is a **no-op on the first iteration** — enters at `BELFEM_REAL_MAX`, multiplies by ≈1.3, clamps to 1
— so it cannot protect the iterate most at risk), and R1 would have created a second test TU,
duplicating `solve_attempt`/`take_step` on a *de-duplication* row.

## The jury round: three P0s

Codex `gpt-5.6-terra`/high and Grok `grok-4.6`/high, blind, both read-only (verified against a
pre-round snapshot).

1. **R8's documentation inventory was wrong by nine sites — it said four, there are thirteen.**
   Codex found one, Grok found four, and re-running the grep independently found four more. The one
   that matters most was found by neither auditor: **`scripts/update_doc_index.py:124` is the
   generator** that emits the `executables` blurb producing `doc/doxygen_nav.dox:52` and
   `doc/groups.dox:266`. Patching the generated files without it means the dead binary name
   *regenerates* on the next doc-index refresh. Also newly found: `Doxyfile.in:1166` and
   `examples/scripts/Allrun:120-125` (a shipped script). Grok's own best catch,
   `examples/scripts/README.md:63-72`, is user-facing and **already false** — it says the binary
   *reads* `circuitAnalysis.txt`; it writes it.
2. **R6 would not compile on a `USE_HDF5=OFF` tree.** Both auditors caught it independently:
   `FileGuard` is defined inside the `#ifdef BELFEM_HDF5` block, and R6 writes a `.cir`, not HDF5.
   The right pattern already exists as `write_scratch` in `test_HybridFactory.cpp:43-50`.
3. **Consequently P6 and the Definition-of-Done grep were false as written**, and R7's compatibility
   premise is true as *invocation* but was false as *documentation*.

Six P1s, all confirmed at `file:line`: the R2/R6 1e-12 gate is undefined at the sine's zeros, which
a 400-step window contains exactly (k = 200, 400) — it needs an absolute floor; Appendix A's ε₀
anchors are wrong (`:2143` is `dump_system_if_requested` in the *magnetic* Newton, `:155` is a print
line — the live assignments are `:2547` and `:162`) and R9 would have regressed the register by
copying them; the netlist parser case-folds every label, so an R6 compare against `"V1"` fails;
R5's `−1e-9 V` rectification floor is tighter than the Newton *relative*-residual tolerance can
guarantee; `solve_attempt` returns `bool`, so R1′'s decision evidence is unrecordable as specified;
and **R1′ cannot repair a first-iterate overflow at all** — `solve()` mutates `mX` with no rollback,
so ω would scale an update against a poisoned state.

The last one interacts with Grok's best technical finding: **the failure mode named in the plan was
the wrong one.** A 10 V source cannot drive any bridge diode to the 18.45 V `exp` overflow bound; the
real risk is **commutation** at the source zero crossings, where the conducting pair swaps and the
newly forward diodes start from reverse bias — which the 0.31 V-per-step warm-start argument does not
cover. If R5 flakes, that is where to look, and the mitigation is a smaller Δt, not relaxation.

Grok's O2 ruling was adopted: keep ±1 % on R4/R5, do not tighten to 0.1 %, because tightening
converts a topology gate into a discretisation gate that §2 explicitly refused.

## Four stale doc claims fixed out of band (Christian approved)

The audit surfaced documentation that was false **independently of DR-39** — not caused by the
pending deletion, and therefore wrong to park behind R7's gates. Christian approved fixing them
separately, and they are the only tree edits this session made:

- `examples/scripts/README.md:65-67` and `examples/scripts/Allrun:120-121` both claimed the binary
  **reads** `circuitAnalysis.txt`. It reads no input file at all — the only stream in the file is an
  `ofstream` (`electricalCircuit.cpp:41-42`). Both now say it *writes* the file. `Allrun` re-checked
  with `bash -n`.
- `src/circuit/doc/circuit_usage_guide.md` §11's runner bullet was false in **three** ways, not the
  one the audit named: the MNA matrix is recomputed on *every* attempt rather than only when Δt
  changes (`:120`); `shift()` runs at the **top** of each attempt (`:119`), not the end; and the
  solution is written **after acceptance** (`:182-183`), so each output line is that step's own
  result, not the previous step's. All three had been true before the DR-118 fix of 2026-08-29.
- The same section's ngspice bullet still said "deferred, nothing written" and pointed at
  `todo/deferred/ngspice_parser_plan.md` — a dead path and a false status, since v1 shipped
  2026-08-27. It now records the importer as shipped, names `NetlistParser` /
  `NgspiceCircuitFactory` and the hybrid `circuit → file` key, and lists `circuitrun` and PULSE/PWL
  as the remainder.

No language sweep was run on these: they replace wrong statements with right ones, which is
currentness work rather than prose. The sweep stays owed on whatever R8 rewrites. R8's site table
was updated so the implementer does not redo them — what remains at those sites is only the
deletion churn.

## State

Plan is audited and implementation-ready. **No C++ touched and nothing built or run — reviewed, not
verified.** The only tree edits are the four documentation corrections above. Exchange: `tmp/ai_exchange/review_dr39_close.md` (pre-registration frozen
before dispatch, per-citation verification, reconciliation table). Next: implementation dispatch,
then `make check` as the gate, then R9's strike.
