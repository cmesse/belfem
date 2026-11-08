# Jury wrappers: one private work directory per round, receipts that match what the wrappers write, a fail-closed merge lock, and bash 3.2

**Date:** 2026-11-08
**Purpose:** Fix the `cross_review.sh` temp-file collision between concurrent rounds, found
2026-11-07 when the DR-202 and DR-203 rounds ran at the same time, and the empty-array abort that
`ask_copilot.sh` hits under macOS `/bin/bash` 3.2. Christian approved taking this "tooling round"
under the 2026-09-22 ruling, which says tooling scripts owe a jury.
**Module:** meta (`scripts/cross_review.sh`, `.claude/scripts/ask_copilot.sh`, `.claude/scripts/ask_grok.sh`)
**AIs involved:** Claude Opus 5.5 (primary). Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high
sat on the plan jury and code round 1; Codex `gpt-5.6-terra`/xhigh and Grok `grok-4.7`/xhigh on
code round 2.
**Status:** COMPLETE. Plan `todo/closed/jury_wrapper_portability.md` v2 plus amendments R7 and R8,
all approved by Christian. DR-204 filed.
**Verification:** executable gates, run by Claude in the scratchpad with stub CLIs (`CODEX_BIN`,
`GROK_BIN`, `COPILOT_BIN`); no build and no real auditor in the gates. Real tree:
- G1 syntax under bash 5.3 and 3.2
- G2 two concurrent same-slug jury rounds
- G3 `--quick` under `/bin/bash`
- G4/G5 empty and non-empty argument arrays under `/bin/bash`
- G6 drift scan and its self-test

A case-sensitive APFS sparse image carried the rest: same-slug concurrency, a non-canonical slug, a
stale lock, and SIGTERM during the lock wait. Positive controls showed that HEAD's script collides
and that the round-1 receipt names lose audits on a case-sensitive volume.
**Exchange:** `tmp/ai_exchange/review_jury_wrapper_portability.md` (+ `.diff`)

## Summary

The script had four defects. One was latent until this round's own fix exposed it.

- **Literal temp name.** `mktemp ".../crossrev_subject_XXXXXX.md"` is a literal name to BSD
  `mktemp`, which randomizes only trailing X's. Two concurrent rounds collided on it, and a round
  killed with SIGKILL left it behind for the next one. Each round now gets its own
  `mktemp -d .../cross_review.XXXXXX` directory, with the cleanup trap installed only once the
  directory exists.
- **Fixed receipt names.** The jury's per-seat receipts were named after the slug alone, so two
  rounds on one slug deleted each other's receipts (Codex, plan jury). They now carry the work
  directory's suffix, `RUNTAG`.
- **bash-4 syntax.** `${AUDITOR^^}` on the `--quick` path needs bash 4: `/bin/bash` stops with
  `bad substitution` when the line runs, although `bash -n` passes it. It is now a `tr`.
- **Empty arrays.** Under `set -u`, bash before 4.4 treats an empty array expanded as
  `"${arr[@]}"` as unbound. Every array that can be empty is now expanded as `${arr[@]+"${arr[@]}"}`
  (`ask_copilot.sh` `EFFORT_ARGS`; `ask_grok.sh` `PERM_ARGS` and `tDisallow`). `ask_grok.sh`'s
  comment had stated the wrong rule. `ask_copilot.sh` now starts with `#!/usr/bin/env bash`.

## What the code juries changed

**Code round 1, Grok's P1, introduced by my own receipt names.** The random suffix contains upper
case. `mktemp` produced `LpMP0K`, `itHWfB` and `OnvBG1`. All three wrappers lowercase the slug
before using it as a filename (`sanitize_slug`). On a case-sensitive filesystem `collect` therefore
looked for a file that did not exist and recorded `Audit FAILED` for a seat that had succeeded. In
`--quick` the result would be a P0 that never reaches `AUTOREVIEW_P0.flag`.

This Mac's volume is case-insensitive (`touch CaseT; ls caset` succeeds), so the first G2 could not
see the defect. A case-sensitive sparse image reproduced it:
- old receipt names: 0 entries, 4 `Audit FAILED`, 4 orphaned receipts
- fixed: 4 entries, no failures

`RUNTAG` is now lowercased. **Lesson:** a gate that compares paths has to run on a filesystem that
can tell them apart.

**Code round 1, Codex's P2.** Now that the second round no longer dies in `mktemp`, two same-slug
rounds both reach the unlocked final `tee -a`. My answer was a `mkdir` lock, which is atomic
everywhere; `flock(1)` is not on stock macOS.

**Code round 2, Grok's P1, introduced by that lock.** `collect` deletes each receipt as it reads
it, and I called it before the lock wait. For up to 120 s the only copy of both audits was
`$WORKDIR/merged.md`, which the EXIT trap deletes. An interrupt during the wait therefore lost both
finished audits.

**Code round 2, Codex and Grok on the timeout path.** On timeout the script appended without the
lock, and it never cleared a stale lock.

**The fix for both.** The lock is taken first. `collect … | tee -a` streams as it did before, then
the lock is released. A timeout exits 3, keeps the receipts and names the lock. The gates on the
case-sensitive image:

| Case | Result |
|---|---|
| two concurrent same-slug rounds | entries in clean `CODEX GROK CODEX GROK` order |
| stale lock | rc 3, nothing appended, receipts kept with their text |
| SIGTERM 3 s into the wait | rc 143, receipts kept, the foreign lock intact |

**Pre-existing, raised by all three voices.** `cross_review.sh` used the raw slug while the wrappers
sanitize theirs, so a slug such as `review-jury` produced `Audit FAILED` plus an orphan. `SLUG` is
now canonicalized with the wrappers' own expression; `Review-Jury` gives `review_jury.md`, with
2 entries and no failures.

**Refuted along the way:**
- My plan's inference that one round's trap could delete another round's file (Grok:
  `mkstemp` is `O_EXCL`, and the trap fires once).
- Grok's claim that `${SHA: -1}` needs bash 4.2. It ran on 3.2; only a negative *length* needs 4.2.

## Process notes

- **Plan jury merged by hand.** The plan round was dispatched with `nohup … &`, and stopping that
  untracked parent skipped `cross_review.sh`'s merge step. Both seat files were complete, and they
  were merged by hand with a note in the thread.
- **Positive control recognized.** The first code round failed on the very collision under repair
  (DR-202 against DR-203). After the fix, HEAD's script was run under the same conditions as a
  positive control: `mkstemp failed … File exists`.
- **Shim for G3.** `--quick` needs `flock(1)`, and stock macOS has none. G3 ran with a shim on PATH,
  and the gap is filed as DR-204.
- **Gate side effect.** G3 created `tmp/ai_exchange/autoreview.md`, which did not exist before; it
  was removed afterwards, and no P0 flag was written.
- **Removal blocked.** A safety check refused one `rm -rf $R` in the test-bed setup. The volume was
  fresh, so the setup was rewritten without the removal.

## Open Questions

- DR-204: give `--quick` the same `mkdir` lock.
- bash 4.0–4.3 is not an executed endpoint. The idiom is the standard workaround for that range.

## Files Updated

`scripts/cross_review.sh`, `.claude/scripts/ask_copilot.sh`, `.claude/scripts/ask_grok.sh`,
`todo/closed/jury_wrapper_portability.md` (new), `todo/README.md`, `todo/debt_register.md` (DR-204),
`devlog/README.md`.
