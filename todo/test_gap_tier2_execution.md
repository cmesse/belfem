# Test-Gap Tier 2 Execution: 48 Triage Rows Plus the Six Reopened Exclusions

**Date:** 2026-09-22
**Purpose:** Run Tier 2 of `tmp/devlog_sweep/triage.md` (§2, 48 rows) plus N-49 … N-53 reopened by
the exclusion re-read (`todo/test_gap_backlog.md` §3b; N-54 stays Tier 3) through the standing
loop, batched by suite: one plan per batch (a section per row), one plan jury per batch, Opus
implementation, gates in both trees, one code jury per batch diff, one commit per row. Consumer
inventory before any src fix that changes a value or refuses an input; rulings batched for
Christian at the end (the quiz).
**Module:** every suite under `tests/`; src fixes only where a row's defect is live and the jury confirms
**AIs involved:** Claude Fable (orchestration, planning agents, juries, commits), Codex + Grok (juries; terra/high, astra/high on the critical list), Opus (implementation)
**Status:** IN PROGRESS — started 2026-09-22 (afternoon wall clock); planning agents launched for the seven batches

> **Standing approvals (carried from Tier 1, `memory: test-gap-overnight-run`):** tests and jury-confirmed
> local src fixes land without a per-row stop; Claude runs make/ctest; one commit per row, own files only,
> peer WIP never staged (`git diff --cached` empty before every add; shared CMake lines staged apart);
> a red not fixed in one attempt, a vendor quota failure, or an open ruling → register row and move on.
> One `make` per tree. No G4 mutation chain while a jury reads the tree (P5).

## 1. Batches

| batch | suite(s) | rows | plan slug |
|---|---|---|---|
| B1 | math, linalg, containers, core | T-05, T-14, T-19, T-28, T-53 (linalg half), T-54, N-02, N-05, N-23, N-52 | `review_tgap2_b1` |
| B2 | io, homology | T-04, T-06 (O9), T-07 (O1), T-43, T-56, T-77 | `review_tgap2_b2` |
| B3 | sparse | T-09, T-17, T-18, T-20, T-41, T-42, T-57, T-76 (MKL lane), N-22, N-28, N-29 | `review_tgap2_b3` |
| B4 | mesh | T-03, T-39, T-45, T-46 (mesh half), T-48, T-55, N-03, N-08 | `review_tgap2_b4` |
| B5 | fem, circuit, interpolation | T-15, T-21, T-31, T-35, T-37, T-46 (fem half), T-50, T-59, N-01, N-04, N-06, N-07, N-20, N-21, N-25, N-26, N-27, N-30, N-49, N-50, N-51, N-53 | `review_tgap2_b5a` / `_b5b` |
| B6 | physics, gasmodels | T-16, T-22, T-27, T-29, T-30 (O3), T-58, N-24 | `review_tgap2_b6` |
| B7 | comm | T-40, N-28 (comm half), T-53 (commmpi half) | `review_tgap2_b7` |

## 2. Ledger

| id | batch | re-read on `devel` | gate | commit | residue |
|---|---|---|---|---|---|
| (filled per row as it lands) | | | | | |

## 3. Rulings owed (the quiz)

- O1 (T-07) `Ascii` relative path: cwd or `$PWD`? — carried from the triage.
- O3 (T-30) Ni/Fe λ(300 K) within 15 %: land `DISABLED_` until ruled.
- O5 (T-53) logical size vs capacity after shrink.
- O9 (T-06) `.bfm` foreign git hash: refuse or flag.

## 4. Incidents

(none yet)

## 5. Audit trail

Exchange threads `tmp/ai_exchange/review_tgap2_<batch>.md`; plans `..._plan.md`; gate logs
`tmp/ai_exchange/tgap2_*.log`.
