# Mutation suite for the practice-drift checker: follow-ups

**Date:** 2026-09-23
**Purpose:** Residue of the review of `scripts/tests/mutate_check_practice_drift.py`, split out of the
closed plan `todo/closed/practice_drift_checker.md` as the plan template prescribes. The suite proves that
`scripts/check_practice_drift.py --self-test` can go red: each mutant edits the checker, and the self-test
must then fail.
**Module:** meta (`scripts/tests`)
**AIs involved:** Grok `grok-4.7`/high (blind review). **Codex produced no audit** — usage limit, resetting
at 19:10 — so this was a one-voice round.
**Status:** OPEN — D-m1 and D-m3 fixed before the commit on Christian's ruling (2026-09-23); D-m2, D-m4,
D-m5 and D-m6 open. A narrow re-review of the fixed suite is owed once Codex is back from its usage limit.

> **What was fixed before the commit:** the round found the suite could report a false green, the exact
> failure it exists to prevent. Christian ruled to fix the baseline and the encoding and commit, with the
> rest filed and a narrow re-review once Codex returns.

## Findings — all confirmed

- [x] **D-m1 (false green — executed) — FIXED**: the unmutated `--self-test` runs first; if it is red the
      suite stops with exit 2. Proven on the same always-red checker: was `68/68, exit 0`, now exit 2.
      Original finding: no baseline run. Against a checker whose self-test always fails,
      the suite reports `68/68 mutants killed` and exits 0. Fix: run the unmutated `--self-test` first and
      exit 2 with a message if it is red.
- [ ] **D-m2 (weak attribution — executed)** — a mutant counts as killed on any non-zero exit, so a
      traceback in an unrelated fixture kills it. `D-c46 closers never close` dies of an `ExtractionError`
      in the *clean* fixture, with zero FAIL lines, before its own fixture runs. Fix: require the expected
      `FAIL <label>` line in the child's stdout, per mutant, rather than the exit status alone.
- [x] **D-m3 — FIXED**: both opens pass `encoding="utf-8"`. Original finding: `open()` without `encoding=` on the checker source, which contains U+2014; under a
      non-UTF-8 locale the read can raise before the first mutant. Add `encoding="utf-8"` to both opens.
- [ ] **D-m4** — the anchor check is an `assert`: count 0 and count 2 both abort the whole campaign with
      "not unique", and `python -O` strips it, turning a missing anchor into a silent no-op. Replace with
      an explicit check that reports the entry and continues or fails loudly.
- [ ] **D-m5** — fixed defects with no mutant: the commit-body side of "one count per record"
      (`commits.update(set(MODEL.findall(body)))`), and `MODEL`'s trailing-period and seat-name rules.
      Their fixtures could be deleted without the suite going red.

- [ ] **D-m6 (surfaced by the D-m1 fix)** — the *checker's* `--self-test` is locale-dependent: it
      writes and reads its fixture files with `write_text()` / `read_text()` and no encoding (25 call
      sites), so under a C locale the `—` in a fixture raises `UnicodeEncodeError`. The live run is
      unaffected (its reads pass UTF-8; exit 0 under a C locale). The new baseline guard reports this as
      exit 2 instead of a meaningless result. Fix in the checker, which takes its own narrow jury.

Also noted, not filed: the temporary work directory is never removed; there is no `__main__` guard; a
hang counts as killed, which the docstring states.

## Audit trail

- Thread: `tmp/ai_exchange/review_mutation_suite.md` (Grok's audit, the Codex failure, verification).
