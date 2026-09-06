# Devlog 2026-08-30 — DR-29 coreduce analysis: three-way jury, fix proposal, no code touched

**Date:** 2026-08-30
**Topic:** DR-29 (homology coreduce worklist port + stale-adjacency cleanup #10) — analysis-only round; problem characterization, fix proposal, impact assessment
**AIs involved:** Claude (Fable), Codex (gpt-5.6-sol, xhigh), Grok (grok-4.6, xhigh) — blind jury, per-vendor exchange slugs, reconciled by Claude with citation re-checks
**Claude Confidence:** high on the confirmed defects; the one medium claim (F6) was refuted in-round
**Codex Audit Confidence:** high
**Grok Audit Confidence:** high on inventory, medium on abort reachability
**Literature References:** Mrozek & Batko 2009 §4/§5/Alg. 6.1 (Thm 4.1, Thm 6.2); Pellikka et al. 2013 §2.2; Giard et al. 2026 draft (tmp/cohomologyPaper.pdf), Eqs. (21)/(22), Algorithms 2-4
**Verification:** reviewed, not verified — static source trace + literature consistency + reconciled AI audits only. Nothing was built or executed. Christian's standing instruction for this session: *we are not touching it.*

## Summary

Christian asked for a three-way analysis of DR-29 before anyone touches the cohomology
algorithm again, with the explicit history that nearly every past change broke it. Claude
pre-registered eight findings (F1-F8) in the exchange, then dispatched Codex-sol and Grok
at xhigh with a neutral blind brief (identical five-question docket, per-vendor slugs).
Both returned full audits; every load-bearing Grok citation was re-checked against the
tree and held. Deliverable: `tmp/dr29_cohomology_report.md`, written for Gregory (module
author), emailed by Christian.

## Key Findings

- **Production path correction (3-way):** the register's "live path = `coreduce_complexCCR()`"
  was wrong. Production is `coreduce_complexPellikkaGeneralized()` (`cl_MaxwellFactory.hpp:48`
  default; `mSuggestHomologies` hardcoded true at `cl_CutFactory.hpp:91` so `reduce_*` never
  runs). Register row corrected in place.
- **#10 confirmed (3-way):** zero-valid-neighbor case in `pGeneralizedCocombine` is an
  elementary coreduction pair; the two cleanup blocks (`cl_SimplicialComplex.cpp:1175-1196`)
  are skipped while `remove_kcochainFromMap` still runs (`:1223-1224`) → dead IDs in live
  boundary maps. Reachable mid-pass (created by earlier neighbor-merges); vacuous in 2-D,
  so the annulus battery structurally cannot see the worst of it.
- **The fragility mechanism, named:** no cochain-side cleanup ever removes a dead ID from a
  survivor's coboundary map (the chain side does — `cpp:730`, `cpp:909`; asymmetry is real).
  Several `Map::operator()` reads abort on a dead key in every build type
  (`cl_Map.hpp:230-246`); they are safe today ONLY by phase ordering, which is checked and
  documented nowhere. Both auditors independently constructed the same production-reachable
  abort cascade (3-D, p=1 zero-neighbor skip → tet boundary shrinks to dead raw singleton →
  `pCoreduce(2)` picks it → abort at `hpp:369`). Order-dependent via unordered_map iteration —
  hence "works until you change something."
- **F6 REFUTED — the trap for the next person:** Claude and Grok both initially concluded the
  Smith-stage matrix `mD(1)` misses the phase-2 merge projection (row copies not updated).
  Codex refuted with a Schur derivation: BA=0 forces the pivot-row identity αd+λD=0, so the
  transformed lower operator is [0; D] — surviving rows unchanged, and the live-restriction
  `createMatrixFromCoboundaryMap` computes is exactly the retracted operator. Claude
  re-verified the derivation independently; refutation stands. **Do not** mirror merges into
  (p−1) coboundary rows; **do not** switch the matrix source to boundary maps (Grok's
  alternative fix — wrong, discarded at reconciliation). Residual true gap: no reciprocal
  unit-pivot validation (`b.coboundary[a]` vs `a.boundary[b]`).
- **Latent, not live (both auditors):** `coreduceOmit` could omit a second vertex from the
  same component only if a stale entry stalls `pCoreduce(0)` first — not producible from the
  fresh production sequence. Handled by component-aware seeding in the eventual port.
- **`pCocombine` is not a worklist template (3-way):** `Cell::pop` is LIFO, `unique()` sorts
  the queue, unguarded pop dereference, raw size-2 test, no ±1 check. The 2026-07-03 "we
  already have a worklist in-module" precedent is a category error.
- **`pCoreduce` lacks Mrozek §4's ±1 pivot check** (`hpp:364-367`) — narrow reachability
  today, mandatory before any reordering (connects to the deferred #11).

## Changes Made / Proposed

- Proposed fix (design only, strictly ordered; full text in the report): (0) debug-level
  invariant — live↔live row/column incidence agreement + mD(k+1)·mD(k)=0, the check that
  would have caught everything; (1) #10 cleanup unconditional AND guarded (incl. the
  `(p+2)` walk `:1190-1194`, today shielded by the very skip being removed) + reciprocal
  unit-pivot check; (2) `pCoreduce` pair test live-and-unit; (3) Alg. 6.1 worklist port
  last — FIFO + membership set, dequeue revalidation, IDs not pointers, component-aware
  omit — still sequenced behind rectification T5/T8/T9, in separate commits from the
  correctness steps so cut-shape changes stay attributable.
- Impact: Betti numbers and cut counts must be invariant at every step (gate); cut shapes
  may legitimately move at steps 1-3; step 0 and the guards are inert on green runs.
- **No source files were modified.** Register row DR-29 updated (live-path clause struck in
  place, findings + fix plan recorded); `check_doc_claims.py` 37/37 green after the edit.

## Outcome — handed off, not queued

Christian sent the report to **Gregory Giard**, who develops the homology/cohomology module
(first author of the algorithm paper draft this round read). **DR-29 stays deferred and the
decision is his.** The ordered fix plan is a proposal for the module author, not work queued
for a BELFEM session: no future session should self-start it, and the analysis does not need
repeating. The gates listed below/in the register apply only if he adopts the plan.

## Open Questions

- Does a `p=2` neighbor-merge actually occur on production 3-D decks? Irrelevant for
  correctness after the F6 refutation, but it decides how often #10's cascade can fire.
  A one-line counter probe would answer it (probe policy: piggyback on the fix round).
- Whether historical breakages actually present as bare "Key not found in map" — worth a
  grep of old run logs before the fix round, as the pre-registered falsifier of the
  fragility mechanism.
- DR-29 #11 (unit-coefficient check) stays deferred but must be re-decided at port time —
  the port re-orders which size-1 states arise, invalidating the 2026-07-03 unreachability
  argument.

## Files Updated

- tmp/dr29_cohomology_report.md (new — the deliverable for Gregory)
- todo/debt_register.md (DR-29 row: live-path correction + analysis outcome)
- devlog/dl20260830_dr29_coreduce_analysis.md (this file), devlog/README.md
- tmp/ai_exchange/dr29_coreduce_stale_adjacency.md (pre-registration + reconciliation),
  dr29_brief.md, dr29_coreduce_codex.md, dr29_coreduce_grok.md (ephemeral, will be swept)
