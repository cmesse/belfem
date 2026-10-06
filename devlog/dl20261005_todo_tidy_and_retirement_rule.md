# Todo tidy, DR-105 struck, and the plan-retirement rule at session close

**Date:** 2026-10-05
**Purpose:** Christian's guess that most of `./todo/` was already done or obsolete, checked against the tree and the registers; and a standing rule so that finished plans stop piling up on the active list
**Module:** meta (`todo/`, `doc/ai_collaboration_protocol.md` §6, `.claude/skills/session-close`)
**AIs involved:** Claude Opus 5.5 (primary); no auditor seat (docs and registers only; protocol edits are text, not code)
**Status:** COMPLETE — active list 17 → 6, index rebuilt, rule written into three places
**Verification:** static only. `scripts/check_doc_claims.py` 38/38; `scripts/check_practice_drift.py` reports only its three pre-existing items (two `grok-4.5`/`grok-4.6` allowlist rows, four retired devlog titles), nothing from this diff

## Summary

The guess was half right. `todo/` already had 106 plans in `closed/`, 19 in `deferred/`, 11 in
`cancelled/`; the active list held 17 files, of which 8 were finished by their own Status line
and the register, 1 was low-priority tooling, and 6 still describe code or tests not in the tree.
The README index had grown to 1646 lines of dated "Added" entries and four currentness sweeps.

## What moved

To `closed/` (substance landed; at most a gate owed, named in the plan's Status line):
`dr190_dr191_small_fixes`, `dr113_dr190_dr191_residues` (both DRs struck 2026-09-30),
`picard_line_head_residual`, `release093_precheck_fixes`, `release093_ultra_followups`,
`materials_tangent_cluster` (all five DRs struck; the debug-assert leg of its gates owed),
`ghost_swip_average` (DR-168 struck; its R6 private-tree and `make check-fast` gates owed), and
`dr105_picard_stall_plan` — Christian: "we also dealt with dr105". Its R3b (tolerance sensitivity)
and R4a are struck; R3c (the tables) re-filed.

To `deferred/`: `practice_drift_mutation_suite` — D-m2 and D-m10..D-m14, all low, per the tooling
finish-line rule (Christian, quiz).

Kept active: `test_gap_backlog`, `test_gap_followups` (R2-R6, R12, R13 ruled and open),
`test_gap_tier2_execution` (B4 mesh and B7 comm held for a peer's WIP that the clean checkout
suggests is gone; T-04 / T-43 (b) on a ruling), `periodic_postprocessor_fields` (audited plan, no
source), `pancake_inner_leadout` (awaits Christian's S1/S2 ruling), `mkl_num_threads_allrun`
(Allrun still exports no `MKL_NUM_THREADS`). The backlog's and Tier 2's Status lines were stale
("no test written", "seven of eight batches planned") and are rewritten to the landed batches.
Christian asked about the three test-gap files at close; whether to fold the Tier 2 residue and the
followups' open rows into the backlog is left to him.

## Register

DR-105 struck to `debt_register_closed.md` on Christian's ruling (quiz: strike, re-file the
tables), rebuilt cell by cell. Residue ( 3 ) is DR-193 [MIXED] [P]: the usermat plugin's tables
end at 199 K (`k_Hast`) … 390 K (`rho_Hast`), so the 1770 K hot spot both algorithms reproduce is
the model outside its tables from ≈ 47 ms; ruled 2026-10-03 to report against the envelope, table
extension a separate task. Live 6, archive 182; the preamble's `[P]` 3 / `[W]` 3 count is unchanged.
Links in the registers and the active plans now point into `closed/`; devlogs keep their paths.

## The retirement rule

Christian, mid-session: a session close should remove the plan it worked from, or at least strike
its pending steps. Written as protocol §6 "Todo-file retirement at session close", in the
`session-close` companions table, and in `todo/README.md`: the close that lands a plan's
substance `git mv`s it to `closed/` (or `cancelled/` if abandoned by ruling) in the same commit,
names the owed gate in its Status line and the register, and repoints the index line; a plan that
cannot close gets its finished steps ticked, obsolete ones struck, and its Status rewritten.
"Remove" is implemented as the move, not deletion — `closed/` is where the README already keeps
finished plans, and the eight moved today were exactly the backlog this rule prevents.

## Index rebuilt

`todo/README.md` 1646 → ~220 lines: header and layout rules kept, then one line per plan in
Active / Deferred / Cancelled / Closed (+ `closed/tests/`). Descriptions taken from the old
index entries, falling back to each plan's `**Purpose:**` or title where the old entry carried
sweep rationale instead. The chronological index is `git show 18e57fe:todo/README.md`.

## Files Updated

`doc/ai_collaboration_protocol.md`, `.claude/skills/session-close/SKILL.md`, `todo/README.md`,
`todo/debt_register.md`, `todo/debt_register_closed.md`, `todo/test_gap_backlog.md`,
`todo/test_gap_followups.md`, `todo/test_gap_tier2_execution.md`, nine moved plans,
`devlog/campaigns/release_1.0.md`, `devlog/README.md`.
