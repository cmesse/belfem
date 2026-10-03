# DR-105 weekend campaign: the "parallel wall" is a floor retry that replays a stalled Picard attempt twenty times; the hybrid completes 0.3 s twice; the 1770 K is the model outside its tables; plan v2 approved

**Date:** 2026-10-03 (campaign 2026-10-02 18:02 → 2026-10-03 12:56; survey and plan rounds the same afternoon)
**Purpose:** Christian offered the weekend for DR-105 as an unattended run. The campaign discharged the
row's long-rerun residue twice, replaced its "unexplained parallel wall" with a measured controller
mechanism that is not parallel at all, turned the "physics ruling" into a material-table validity
envelope, and produced an approved plan for the controller change. Nothing in `src/` was edited.
**Module:** `src/fem/kernel` (`Controller`: the stagnation guard, `try_escalate_to_newton`, the floor
retry in `reset_timestep`, `initialize_timestep`, memdump save/load) — read only; `tmp/dr105_weekend/`
(the recovered T5 deck, driver, parser, ten run directories); `todo/dr105_picard_stall_plan.md` (new);
`todo/debt_register.md` (DR-105)
**AIs involved:** Claude (campaign design and driver, monitoring, the brief, the plan); Codex
`gpt-6-astra`/xhigh + Grok `grok-4.7`/xhigh on two blind rounds (the survey on the brief, the plan
jury). Monitoring and the mid-campaign additions ran on Opus 5.5 after a `/model` switch; the brief,
plan and rounds on Fable 5.1
**Status:** campaign COMPLETE; plan v2 APPROVED with rulings O1–O3 (O4 open); R0 diagnostics launched
at close (`r0a_stall0`, `r0b_anderson0`, `r0c_np4`); implementation is a later round with its own code
jury; DR-105's row restated
**Literature References:** Messe et al. 2023 §4, Eqs. 12–14 and the 10⁻¹¹ terminal criterion (the brief's
§2.7 pointer was wrong — §2.7 is the kernel layout; Grok); Walker & Ni 2011 (Anderson type II as
implemented); Anderson 1965
**Verification:** executable throughout — ten runs on a private HEAD `5da57f7` debug binary (MKL, Blaze,
Open MPI, pinned copy), logs parsed by `analyze.py` (calibrated on the August logs: np=8 263 accepted /
213 floor warnings / 81 %; serial 1.1 %), fields compared with `scipy.io.netcdf_file`. Both auditors
re-read the logs and the source; their citations were re-read by me (all CONFIRMED). No build of a
change, no `make check` (no change exists)
**Exchange:** `tmp/ai_exchange/review_dr105_picard_stall.md` — pre-registration, two survey entries,
verification + reconciliation, plan pre-registration, two plan-jury entries, verification +
reconciliation. Brief: `tmp/dr105_weekend/STALL_BRIEF.md`; pre-registered run plan and its three
amendments: `tmp/dr105_weekend/ANALYSIS_PLAN.md`; driver-written `STATUS.md`

## Summary

The reconstructed T5 deck was not lost: it survived in `~/codes/belfem_transfer` with both August
logs, and its two old-API plugins build unchanged against HEAD. A four-minute smoke test showed the
serial debug binary at 24 ms/h where the 2026-08-27 binary did 1 ms/h, so the 0.3 s rerun had become
a 13-hour job. Christian approved an unattended campaign under `setsid`: six floor discriminators,
then the serial long run beside the first clean parallel configuration, both on the debug binary.

**The August np = 8 wall at 0.375 ms does not reproduce**: np = 8 ran 28 ms clean in its hour. What
every pure-Picard configuration does instead is abort on the floor-retry limit somewhere between
10.6 and 33 ms into the quench — np = 4, np = 8 with MUMPS, np = 8 with STRUMPACK on the second try,
and the serial long run at 15.41 ms. The two hybrid runs (`algorithm : Newton`, serial and np = 8)
completed the full 0.3 s with no floor warning and agree with each other to 1e-5 on temperature.

**The mechanism, established by the logs and confirmed by both auditors:** in the quench a Picard
attempt starts from a predictor already within 5–30× of tolerance; the first body lifts the
residual by 9–30 dB; a flat cluster follows. The stagnation guard (5-sample ring, 0.2 dB mean
absolute deviation) ends any attempt whose fifth sample is still flat — at iterate 6, every time. A
pure-Picard deck cannot escalate to Newton, so Δt is cut; at the floor the cut becomes a retry that
raises `max iterations` and the watchdog ×2, ×4 — budgets the guard never lets the attempt reach —
and the next attempt restarts at ω = 1, so the twenty retries replay one attempt to printed
precision (one distinct trace among 21 attempts in three of the four aborts). The hybrid takes the
same guard verdict and escalates instead (414 times per run). Whether the cluster would resolve if
the guard let it is not established: every surviving attempt at a larger Δt did resolve (8–27
iterates), but the guard censors any plateau longer than four. That is what R0a measures.

**The 1770 K.** Christian's prior was a Newton overswing. The logs show none (largest per-step
T_max change 0.4 K; the hybrid's T_max equals pure Picard's to 0.1 K at every common time to 33 ms;
monotone Newton descents). Both auditors found the real issue in the plugin: the tables end at
199 K (Hastelloy conductivity), 282 K (HTS heat capacity), 300 K (copper), 303 K (Hastelloy heat
capacity), 390 K (Hastelloy resistivity), the endpoint is held above, and copper's resistivity is
clamped at 1358 K. The 0.3 s result is the implemented model outside its tables from ≈ 47 ms on.
Christian ruled: report against the 199 K envelope; extending the tables is a separate task.

**Plan.** My first design — a patience factor on the guard that grows with the retry count — was
refuted by both plan auditors on arithmetic: with overlapping windows it delays the cut by at most
three iterates and never reaches the raised budgets. Plan v2 suspends the guard's Picard verdict
while a floor retry is in force, passes the cut reason to `reset_timestep` (20 call sites, so a
missed one is a compile error), and lets the raised watchdog and iteration budget bound the
attempt. Christian approved it, unconditional, with the five pure-Picard shipped decks left as they
are.

## What the auditors corrected in my brief

Recorded so the next reader does not inherit them: the plateau level is not Δt-independent (it
drops ≈ 6 dB per halving; the floor is where the next halving does not exist); the footer's
`|λ|max/|λ|min` is not κ₂ and the residual sits three decades below that ratio times ε, not five
above; the "natural 4-iterate plateau" is censored by the 5-window guard; thermal enters a stalled
attempt unconverged and converges in one body; memdumps DO persist the BDF history (my "re-anchors
at BDF1" came from the stale INC-339 record); the operative ω reset on the coupled driver is
`initialize_timestep`; Messe §4, not §2.7.

## The campaign in numbers

| run | ranks | algorithm | outcome |
|---|---|---|---|
| p1_np1 / np2 | 1 / 2 | Picard | 9.7 / 12.6 ms in the 60-min cap, clean |
| p1_np4 | 4 | Picard | floor abort 11.745 ms (93.6 K) |
| p1_np8 / np8_cond | 8 | Picard | 28.3 / 26.8 ms in the cap, clean; κ-ratio 8.5e12–1.7e14 |
| p1_np8_mumps | 8 | Picard, MUMPS | floor abort 33.15 ms (146.6 K) |
| p2_np8_strumpack | 8 | Picard | floor abort 10.63 ms (same deck/binary as p1_np8) |
| p2_np1 | 1 | Picard | floor abort 15.41 ms (100.1 K), 2.1 h |
| p2_np8_newton | 8 | hybrid | 0.3 s, 15826 steps, 211 rejections, 3.1 h, peak 1770 K |
| p2_np1_newton | 1 | hybrid | 0.3 s, 15824 steps, 208 rejections, 13 h |

Serial vs np = 8 hybrid: I_1 to 2e-13, U_1 to 1.5e-5, T to 1e-5 relative at every 1 ms save.

## Open

- R0a/R0b/R0c running (`tmp/dr105_weekend/r0*`); their result sizes nothing in R1 any more but
  decides whether the plateau resolves at all (R0a) and whether Anderson shapes it (R0b).
- Implementation of R1 (approved plan, own code jury, `make check`, the reproducer gate), R2, R3a–b.
- O4: `target iterations : 100` in the T5 deck (Christian).
- The table extension (O3's separate task).

## Files Updated

- `todo/dr105_picard_stall_plan.md` (new, v2 approved), `todo/README.md`, `todo/debt_register.md` (DR-105)
- `devlog/README.md`, `devlog/campaigns/release_1.0.md`
