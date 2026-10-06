# Mutation suite for the practice-drift checker: follow-ups

**Date:** 2026-09-23
**Purpose:** Residue of the review of `scripts/tests/mutate_check_practice_drift.py`, split out of the
closed plan `todo/closed/practice_drift_checker.md` as the plan template prescribes. The suite proves that
`scripts/check_practice_drift.py --self-test` can go red: each mutant edits the checker, and the self-test
must then fail.
**Module:** meta (`scripts/tests`)
**AIs involved:** Round 1: Grok `grok-4.7`/high only — Codex hit its usage limit. Round 2, after the reset:
Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high, both blind.
**Status:** **v11 committed** 2026-09-23 — D-m1, D-m3..D-m9 fixed; the v11 jury (Codex and Grok at `high`,
blind) found no regression and no defect in the fixes. Open: **D-m2** (per-mutant attribution — a design
task) and five low-severity items from the v11 round, D-m10..D-m14. None makes the checker report a
wrong result.

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
      `FAIL <label>` line in the child's stdout, per mutant, rather than the exit status alone. Round 2
      (Grok) refines this: a FAIL line alone is too narrow — `D-c20 OSError uncaught` is correctly killed
      by a traceback escaping the missing-input fixture, so the predicate is "the expected FAIL label, or
      an escaping traceback from that fixture".
- [x] **D-m3 — FIXED**: both opens pass `encoding="utf-8"`. Original finding: `open()` without `encoding=` on the checker source, which contains U+2014; under a
      non-UTF-8 locale the read can raise before the first mutant. Add `encoding="utf-8"` to both opens.
- [x] **D-m4 — FIXED v11**: an explicit per-mutant check; a broken anchor is reported as BROKEN by name,
      the campaign continues, and the run exits 1 (proven with an anchor edited out). Original finding: the anchor check is an `assert`: count 0 and count 2 both abort the whole campaign with
      "not unique", and `python -O` strips it, turning a missing anchor into a silent no-op. Replace with
      an explicit check that reports the entry and continues or fails loudly.
- [x] **D-m5 — FIXED v11**: five mutants added for the four named behaviours (commit-body count, MODEL's
      trailing period, MODEL's seat names, the window date order, the `--strict` expression); 73/73. Original: fixed defects with no mutant: the commit-body side of "one count per record"
      (`commits.update(set(MODEL.findall(body)))`), and `MODEL`'s trailing-period and seat-name rules.
      Round 2 (Grok) adds two: the window's date sort `sorted(by_name, key=date, reverse=True)` — D-c17's
      mutant only rewrites the name tiebreak — and the `--strict` exit expression
      `1 if (strict and drift) else 0`. Their fixtures could all be deleted without the suite going red.

- [x] **D-m6 — FIXED v11**: all 25 self-test `write_text`/`read_text` calls pass `encoding="utf-8"`,
      inserted at each call's closing parenthesis by an `ast` walk; the self-test and the suite pass under a C
      locale. Original finding: the *checker's* `--self-test` is locale-dependent: it
      writes and reads its fixture files with `write_text()` / `read_text()` and no encoding (25 call
      sites), so under a C locale the `—` in a fixture raises `UnicodeEncodeError`. The live run is
      unaffected (its reads pass UTF-8; exit 0 under a C locale). The new baseline guard reports this as
      exit 2 instead of a meaningless result. Fix in the checker, which takes its own narrow jury.

- [x] **D-m7 — FIXED v11**: a work directory the suite created is removed on every exit path. Original: the default `tempfile.mkdtemp()` work directory is never removed. Clean up
      the directory the suite created; leave an explicitly supplied one alone.
- [x] **D-m8 — FIXED v11**: a red baseline prints its exit status or timeout and the last 20 lines of the
      child's output (proven with a crashing baseline). Original finding: a red baseline is undiagnosable: the suite prints one sentence and exits 2,
      never showing the child's stdout, stderr or return code, and a timeout reads the same as a FAIL.
      Print the child's output tail and say which of the two it was.
- [x] **D-m9 — FIXED v11**: both pipes decode as UTF-8 with `errors="replace"`. Original finding: both `subprocess.run` calls use `text=True` without
      `encoding="utf-8"`; D-m3 fixed only the file opens. A fixture that printed a non-ASCII character
      would make the parent raise while decoding the pipe and abort the campaign.

- [ ] **D-m10 (two raisers)** — the work directory leaks on exception paths: `mkdtemp` runs before the
      checker is read, and cleanup happens only in `finish()`. Use `TemporaryDirectory` or `try/finally`.
- [ ] **D-m11** — a red baseline's last-20-lines tail need not contain the failing `FAIL` line, since the
      self-test prints each where it occurs. Print the `FAIL` lines themselves as well as the tail. (The
      docstring's "shows why" was corrected in the v11 commit to state what is actually printed.)
- [ ] **D-m12** — a baseline timeout returns empty output, although `TimeoutExpired` carries the partial
      stdout and stderr.
- [ ] **D-m13** — two self-test launches in the checker (the D-c16 CLI fixtures) still decode with
      `text=True`; the expected stderr is ASCII, so this is latent.
- [ ] **D-m14** — `D-m5 seat names admitted` changes two things (adds `Opus`, widens the separator to
      `[- ]`); the fixture kills the first, so a separator-only regression would survive. Split it.

Also noted, not filed: there is no `__main__` guard; a hang counts as killed, which the docstring states.

## Audit trail

- Thread: `tmp/ai_exchange/review_mutation_suite.md` (Grok's audit, the Codex failure, verification).
