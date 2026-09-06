# Night session: the reset seed hole, chi default-on, and the BDF1 reversal

**Date:** 2026-08-28 (00:00–01:45, continuous with dl20260827_evening_quenchfront_campaign)
**Purpose:** Record the diagnosis chain that ended the restart ladder — the incomplete soft
reset — plus the chi default ruling, the watchdog calibration finding, the BDF-reset jury, and
the BDF1-over-BDF5 empirical reversal.
**Module:** fem/kernel (Controller), fem/iwg (Timestep), doc/input contract

## 1. The seed hole (jury round `review_reset_vs_warmstart.md`)

Christian's thesis — a soft timestep reset should equal "load a memdump of the last accepted
step" (n = 1) — sent a jury to diff `reset_timestep()` against `load_memdump()`. Grok found
the mechanism; verification confirmed all four legs on source: **the reset restores the mesh
FIELDS but never re-seeds the DOFS, and the Picard residual reads the dofs** (`Calculator::q()`
assembles from fields; `mFieldValues`/`Dof::value()` carry the failed attempt's last iterate;
`load_memdump` calls `seed_dof_values()` with a comment naming this exact hazard — the reset
path never did). A non-finite iterate survives `x -= ωΔ` unchanged → 10–15-rejection cascades
non-responsive to Δt, cleared only by restarts (which seed). The out3/out4 same-dump
divergence showed FP-ordering nondeterminism modulates WHICH attempt poisons, while the
15-in-a-row in-run failures vs 1–2-try restart successes showed the state damage dominates.
Claude's pre-registered prime suspect (BDF un-shift field slot) was REFUTED by both auditors
(the retry's own re-shift heals it before assembly) — the jury structure caught its
dispatcher's error, twice in one night.

**Fix (Christian: "OK, let's fix it!"):** `seed_dof_values()` on both kernels at the end of
`reset_timestep()`, and the same call in `reset_thermal()` (covers the new `belfem`
executable's segregated loop; verified against `belfem.cpp`'s two loops). Auditor-authored
shape; step 2 (magnetic savepoint) deliberately withheld so the experiment discriminates.

**Validated in production the same hour:** out5/out6 crossed the 4300–4320 ms band — which had
consumed three runs and a manual restart ladder — with every rejection recovering in 1–2
same-order retries. No cascade has formed since the fix, across four runs.

## 2. Build-system interlude

The libbelfem refactor's reconfigure had flipped `USE_MKL`/`USE_PARDISO` OFF → phantom
`-lcblas -lblas` link failures (netlib fallback names absent on SCLS). Restored via
`-DUSE_MKL=ON -DUSE_PARDISO=ON`; binaries at 23:09/23:16 carry the seed fix + gauge tangent.
`test_core` link under the refactor still fails the same way — refactor's list.

## 3. chi default-on (Christian's ruling, input contract updated in-session)

"Let's set gauge chi to 1e-4 by default. I feel better when it's on." Landed: ctor
`mPenalty(2) = 1e-4`, Controller absent-block branch 1e-4, opt-out `chi : 0`;
`input_file_reference.md` + `input_schema.yaml` updated same session (DEFAULT CHANGED stamps);
Codex sweep applied — its by-catch caught a third doc site still claiming absent = off.
Gauge-campaign memory re-scoped (vacuous = sub-critical only).

## 4. Watchdog calibration finding

With the seed fix in, out6's remaining ratchet traced to the thermal watchdog
(`watchdog window : 8` vs `max iterations : 15`): it cut attempts at −67.6 dB — 3 dB from
target, relaxation recovered to 1.0, one iterate after best-ever — because the coupled
alternation wobbles the thermal residual ±15 dB/iterate. With max-iterations as backstop, the
window-8 watchdog only pre-empts 7 iterates of a stuck attempt while demonstrably executing
converging ones. Third live exhibit for the floor plan's D2/M1 (trend test, not scatter test).

## 5. BDF-reset jury (`review_bdf_reset_on_reject.md`) and the BDF1 reversal

Christian's proposal (reset to BDF1 on every rejection) was jury-rejected unanimously:
in-house counter-evidence (the 2026-08-15 restart cliff — BDF1 re-entry = different tangent,
detonated a Newton promotion), production-code precedent (order→1 is the discontinuity-reinit
/ repeated-failure rung, never the first response to a convergence failure), and the post-fix
distribution (modal rejection recovers on one same-order retry). Survivor: the N≥2 variant
(= M10) with Grok's implementation-trap catalog (mStepCount write point, mH must not be
zeroed — div-by-zero at order ≥ 3, both equations, PID-hold interplay). Claude's brief
mis-cited DR-112 as precedent (it is counter-evidence) and overclaimed down-switch safety —
both corrected in reconciliation.

**The deck-level lever the jury pointed at was then measured:** `method : bdf1` (out7) crossed
the 4314–4331 band holding Δt at 0.7–1.9 ms with 3-iterate convergence — where BDF5+Newton
(out6) had ratcheted to 0.1 ms — reaching the campaign's furthest point. Two mechanisms:
L-stable damping + history-free tangent (physical), and the 1.5× growth clamp at order 1 vs
1.2 at order ≥ 4 (administrative; noted honestly). This matches the tree's own standing
guidance ("prefer bdf1/bdf2; treat bdf5 as experimental") and Messe 2023 §4's original design.

## 6. The night's synthesis

Much of the accumulated stabilization apparatus (aggressive watchdogs, scheme experiments,
restart rituals) had been compensating for the undiagnosed seed hole. With the root fixed,
the paper's simple configuration (BDF1, halve-on-failure) outperforms the elaborate one, and
each remaining knob can be re-evaluated on honest ground. Register rows owed (morning): seed
hole (fixed, validated), circuit rollback gaps (shift_back misses Switch/sources; memdump
saves no component histories), stale reset ω comment, watchdog calibration, ledger entries.
