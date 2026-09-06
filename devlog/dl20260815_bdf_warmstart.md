# BDF Warm Start: the Memdump Learns the Integrator State

**Date:** 2026-08-15
**Purpose:** Record the campaign that made a warm restart resume at its
earned BDF order instead of re-anchoring the ramp at BDF1 — and the
last-completed-step subtlety both reviewers caught independently
**Module:** fem/iwg, fem/kernel, tests/fem
**Round:** `tmp/ai_exchange/bdf_warmstart_plan.md` + `bdf_warmstart_diff.patch`

## The trigger

A same-evening restart of the tapestack3d quench study re-entered at
BDF1/50 ms into an actively heating state. The BDF1 tangent is a different
operator from the BDF5 one the dump was written under; the magnetic Newton
promotion detonated (residual 8e-6 → 5.6e-3 on the promotion iterate), the
thermal source followed it to 18.9, and the step burned ~40 minutes and one
rejection before a Δt halving cured it. The memdump already carried the full
field history (15 vectors, hundreds of MB) — the run was discarding it for
want of the integrator bookkeeping.

## The design, and the trap at its center

The missing state is three members of `IWG_Timestep`, per equation:

| member | size | role |
|---|---|---|
| `mH` | 5 reals | previous step sizes (slot 4 is a never-written spare) |
| `mStepCount` | 1 uint | ramp position; `mOrderActive` derives from it lazily |
| `mDeltaTime` | 1 real | **the last COMPLETED step h_n — the trap** |

Both reviewers independently found the P0 in the original two-member plan:
`shift_fields` pushes the *IWG's own* `mDeltaTime` into `mH(0)`, and the
controller assigns the upcoming step only *after* the shift — so at dump
time the IWG holds h_n while the already-dumped `/meta/delta_time` holds
h_n+1, two different quantities. A fresh IWG initializes `mDeltaTime = 1.0`,
so restoring only `mH` + `mStepCount` would have the resumed first shift
inject a fictitious one-second step into the history — the same failure
class as the cliff, quieter. The thermal twin's own comment states the
invariant: "shift BEFORE assigning the new delta_time, so that mH(0)
records the size of the completed step."

## What shipped

- `save_history_state` / `restore_history_state` on `IWG_Timestep` — a
  deliberate public API, not a friend leak. Restore is a closed validation
  gate (no member written before the last check): step count zero rejects
  (it would silently ramp to BDF1 behind a full history), the last dt must
  be positive and finite and is *not* clamped to the deck window (it is
  history), entries must be finite and non-negative, and the slots a
  resumed step and its retry read — `aH(0..r-1)`,
  `r = min(aStepCount−1, mOrder−1, 4)` — must be positive; `aH(4) == 0` is
  the normal spare. On success the counter caps at `mOrder` (the ramp
  saturates there). `mHDropped` stays NaN: it is unreachable before the
  first shift refills it (both reviewers; seeding it from `mH(4)` — always
  zero — would have poisoned the first retry).
- `/meta` keys `bdf_h`, `bdf_step_count`, `bdf_last_dt` (+`2`-suffixed for
  the thermal equation), written on rank 0 beside the existing
  `delta_time`. Presence means all three: a partial triple is a named
  hard error ("incomplete BDF history ( magnetic )"), not a generic HDF5
  failure. Old dumps lack the keys and cold-start the ramp — announced via
  raw cout inside the WARM RESTART banner's guard, because
  `message(Default)` is silent at Minimal info level and a cold-start must
  never be silent (Grok). A coupled run resuming a magnetic-only dump
  announces the thermal re-anchor separately.
- Load side: keys read while `/meta` is open, found-flags + payload
  broadcast after the field synch and Anderson clear, restore on every
  rank. Collectivity verified for all four combinations (old dump, fresh
  run, magnetic-only dump into hphiTrun, coupled dump into hphirun).
- Tests (`test_BdfTimestepMethod.cpp` §5): the probe gained the manual
  five-line shift rotation — `shift_fields` needs a DofManager, and a
  restore-then-compute test without the rotation would have passed with
  the P0 alive. The round-trip compares a resumed object against an
  uninterrupted twin **bitwise** (α, all β, all mH); both auditors
  confirmed it goes red if `mDeltaTime` is dropped from the design. Full
  validation matrix including ±Inf, negatives, cap-at-order, BDF1
  harmlessness, and per-gate reject-state-unchanged.
- Docs: four passages taught the old contract ("BDF spacing deliberately
  not persisted") — `input_file_reference.md`, `timestepping_strategy.md`
  twice, `bdf_timestepping_theory.md` — all rewritten, not appended-to.
  The friend-seam comment now names the fourth history writer.

## By-catch

**DR-74:** segregated hphiTrun calls `load_memdump` but never
`save_memdump` — its restart path can only consume dumps written by other
execution modes. Found independently by both reviewers; needs a decision,
not silently absorbed.

## Process

Full two-phase round. Phase 1: both reviewers approved-as-amended with the
convergent P0 plus Grok's validation spec, savepoint/spill analysis
(magnetic rejects use `reset_fields`, not savepoints — my Q3 premise was
wrong), and the docs-contradiction sweep. Phase 3: cores PASS from both;
peripheral FAILs (partial-key policy, message level, test-matrix
completeness, one more stale doc, comment rot) all fixed and re-gated.

Status: **reviewed, not verified** — all syntax gates green; executable
gates are `make check-fast` (the new §5 cases) and a deliberate stop/resume
A/B after the next rebuild. First real beneficiary: any restart of the
running quench study once the binary carries this — the next cliff simply
does not happen.
