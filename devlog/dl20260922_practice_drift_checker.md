# A checker for documents that have drifted from practice, and the drift it caused in itself

**Date:** 2026-09-22
**Purpose:** Three convention drifts found by hand on 2026-09-21/22 were invisible to every gate
because each side was internally consistent. Build the mechanical pass in the practice → document
direction, which `scripts/check_doc_claims.py` does not cover.
**Module:** meta (`scripts`, `doc/ai_collaboration_protocol.md`, `.claude/skills`)
**AIs involved:** Claude Opus 5 (plan, three revisions, implementation, verification). Plan round 1:
Grok `grok-4.7`/high blind; **Codex produced no audit**, two usage-limit failures. Plan round 2:
Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high, both blind, fresh slug. Code jury: same depth,
dispatched at the close.
**Status:** **v2 shipped.** The first version was not shippable: a two-seat code jury found its
gate could lose two DRIFT probes and stay green, which mutation testing then confirmed executably.
v2 fixes all ten and kills 16/16 mutants. A second code jury found eleven more, five of them silent
passes. On Christian's ruling v2 ships and those eleven are v3's work, with a narrow jury on the
diff.
**Verification:** executable gate — `python3 scripts/check_practice_drift.py --self-test` passes
12/12 and exits 0; `--strict` exits 1 on the live tree's DRIFT findings. `check_doc_claims.py`
38/38. The probes are verified; the *design* rests on two audit rounds and is reviewed, not verified.
**Exchange:** `tmp/ai_exchange/review_practice_drift_checker.md` (round 1),
`review_practice_drift_checker_r2.md` (round 2), `review_practice_drift_code.md` (code jury).

## Summary

`scripts/check_practice_drift.py` compares the convention documents against what the devlogs, the
commit bodies and the auditor wrappers actually do. Seven probes across four subjects, tiered by
whether the relationship they check is **normative**: DRIFT is a stated rule the tree breaks, INFO
is a convention moving where the document declines to impose one. Only DRIFT gates under `--strict`.

On the live tree it reports a lagging Grok wrapper default and six devlogs missing from their own
index, plus four allowlisted-but-unruled models and four retired-form titles as INFO.

## What the session actually demonstrated

The plan's own history became the argument for the tool, three times over.

**Round 1 found the plan describing a world that no longer existed.** Its §1 cited `**Module:**` as
absent from the §6 template and the template's field as `**Topic:**` — both true on 2026-09-21 and
both fixed hours earlier in the same session, by me. `**Topic:**` now occurs zero times in the
protocol.

**Round 2 found that documenting a defect had disabled the probe that detects it.** Revision 2 read
§9.1's *prose* for the list of deliberately-unnamed models, and called the design validated on live
data. At ~02:45 the latest-grok policy paragraph was added to §9.1 with a sentence naming
`grok-4.6` as a live defect. From that moment the prose contained the token, so the probe would have
suppressed the one drift it was built to catch. The measurement was taken at 02:37 and copied into
the plan at ~03:20 — true when taken, false when used. Round 1 could not have found it: the defect
was introduced *between* the rounds, by the fix to a round-1 finding.

The line-number errors that came with it were the symptom. The premise was the defect: §9.1 prose
carries escape-hatch declarations *and* defect descriptions, and no token search separates them.
Exceptions are now declared structurally, on one anchored line the protocol owns, with prose
explicitly freed to describe defects.

**The corpus moved under the work twice.** A field census taken at 02:37 disagreed with one taken
at 01:40 in seven of its figures — nothing miscounted, the window simply slid as concurrent
sessions added entries. And the live run at 04:15 reported six unindexed devlogs where 02:37 had
measured seven, because a peer indexed one meanwhile. Every count the checker prints now carries
its window.

## Rulings taken

| Question | Ruling |
|---|---|
| How to declare a deliberately-unnamed model | one anchored line, `**Allowlisted but named by no row (by design):**`, read by the probe and nothing else |
| Can template drift gate? | split by direction — a declared field nobody writes is DRIFT, a practice field the template omits is INFO |
| The live gate | fixtures decide, the live run only reports, so repairing a drift cannot break the build |
| `ask_copilot.sh` | not a jury member; the reserve seat for a temporarily unavailable vendor; low priority |
| Tooling scripts under `scripts/` | take the full plan+audit → code+audit round |

## Protocol changes

- §9.1 gained the anchored by-design exception line, with the rule that an exception is never
  inferred from a sentence.
- §9.1 gained the Grok-seat policy: always the latest allowlisted `grok`; a wrapper default that
  lags it is a defect, not a tier. This is why most rows carry an effort and no Grok model.
- **§9.2 The Reserve Seat** is new: Copilot substitutes for an unavailable seat, a substituted round
  says so in the record, and the reserve does not make an outage silent.

## Open

- Four of six findings I pre-registered against revision 2 were addressed by neither seat and are
  carried as O5: commit-body corroboration may lose a real case; `--recent 40` is underived; P2c may
  be too strict for a deliberate pin; and both live drifts the plan cites are this team's own
  housekeeping rather than errors discipline cannot catch. That last one is the honest weak spot.
- **The code jury found the gate itself was the weak point, twice over.** `--self-test` does not
  require P2a, P2b, P2d or P4 to fire — deleting the P2a or P2b append leaves all twelve assertions
  green, because the precedence check asserts only that a subject appears once, not which probe
  reported it. And `index_targets()` captures the markdown *label* while its docstring claims the
  link target; it is masked only because `devlog/README.md` repeats the filename in both positions,
  and because the fixture I wrote used the same label and destination, so the very test named "P3
  matches link targets, not substrings" could not separate them. A green gate proved less than it
  claimed. Ten defects recorded as D-c1..D-c10.

  They are fixed in v2, and the evidence for that is **mutation testing**, not the green gate: each
  probe is disabled in turn, and each fixed defect is put back, and all sixteen mutants turn the
  self-test red. The same P2a and P2b mutations applied to v1 leave its gate green, which confirms
  round 1's P0 executably rather than on the jury's word.

- **Code jury round 2** (both seats, `xhigh`) confirmed no regression and found eleven more, all
  reproduced. Five are silent passes on a stated rule: an exception id no wrapper accepts, a Grok
  row lagging the policy, a Codex default satisfied by a Grok-column token, an index link into
  another directory, and a recent devlog without `**Date:**`. Two needed rulings and got them:
  a lagging Grok row is DRIFT, and a missing Date is DRIFT. Both are now stated in the protocol,
  so v3 will enforce written rules. Carried as D-c11..D-c21.
- `ask_grok.sh:134` still defaults to `grok-4.6`. One line, and the checker's own P2c reports it.
  Not fixed here: the file was modified by a concurrent session throughout.
- `.claude/scripts/ask_codex.sh:13` still claims only terra and luna are named by the depth table.
  O1 — verifying that class of claim needs prose understood, so it stays out of the checker's scope.

## Files

- `scripts/check_practice_drift.py` (new, 576 lines)
- `todo/practice_drift_checker.md` (new; three revisions, Dn/Rn/On, R1–R8 ticked)
- `doc/ai_collaboration_protocol.md` — §9.1 exception anchor and Grok-seat policy, §9.2
- `CLAUDE.md`, `.claude/skills/session-close/SKILL.md` — the checker registered
- `.claude/skills/code-with-jury/SKILL.md` — tooling-script row in the exemption table
- `todo/README.md` — plan registered

## After the ship: the wrapper fix, D-c22, the index repair, and a file I destroyed

**The wrapper fix unmasked a false positive (D-c22).** Christian moved `ask_grok.sh`'s default to
`grok-4.7`, and P2c went quiet — but P2c's precedence had been absorbing `grok-4.6` for P2a, which
now reported nine devlog headers and two commit bodies recording rounds run on it *while it was the
latest*. History, not drift. Ruled: for the Grok seat a model on the wrapper's arm counts as
declared for P2a, since §9.1 keeps older groks allowlisted to reproduce rounds; the default and the
rows must still be the latest. Filed for v3 (`505075d`).

**Six unindexed devlogs indexed (`37cca0f`).** Each line written from the devlog itself and placed
among its date neighbours: the contact-impedance theory rewrite, the debt-register split, the DR-108
fix round, the CurveFactory walk fix, the greedy-cleaner trace, the Smith evidence pack for Gregory.

**I wiped `devlog/README.md` while doing it.** A one-liner
`open(dst,"w").write(insert(open(src).read()))` with `src == dst` truncates the file *before* the
argument reads it; `insert()` then raised on the empty text and left the working tree at 0 bytes.
The file carried other sessions' uncommitted index lines. Nothing committed was at risk, and no
backup had been taken — against a standing rule to back up endangered files first.

Recovery: the README rebuilt from HEAD plus my six lines, built to a separate path and verified
before being copied into place. Three lost lines (`dl20260918_rutherford_cable_generator.md`,
`dl20260919_smith_inverse_branch.md`, `dl20260921_test_gap_sweep_and_dr113.md`) came back verbatim
from a scratchpad backup taken at 01:28. Claude Code's file-history held no snapshot, because the
edits had been made through Bash. The session "Devlog sweep test coverage index" re-added its two
(`dl20260922_test_gap_tier1_overnight.md`, `dl20260922_homology_report_for_gregory.md`) on request,
and corrected my attribution twice: the `dl20260922_ask_codex_cli_audit.md` line was not theirs, and
they had made no edit to an existing line.

**Still lost:** the index line for `dl20260922_ask_codex_cli_audit.md` — the devlog is untracked and
was written by a Codex run (`gpt-6-astra`/high), so no Claude session can be asked for it — and one
modified existing line whose author and content are unknown. The five restored lines stay
uncommitted: three of them index devlogs that are themselves untracked, and committing those would
put broken links into HEAD.

The lesson is in the backup memory: back up a shared file before any scripted rewrite of it, and
never read and write one path in a single expression.

## v3

D-c11..D-c22 fixed in one pass: by-design ids must be on their seat's arm, a Grok row must name the
latest grok, the Codex default must pass the Codex arm, index destinations are resolved against the
README's directory, P1c enforces the required `**Date:**`, CLI values are validated, a missing input
is exit 2, the window runs date-descending and name-ascending and names its cutoff, and historical
use of an older allowlisted grok is no longer P2a (D-c22, ruled).

Gates: `--self-test` 56/56; mutation testing 30/30 killed. The first mutation run of v3 left five
survivors, and all five were tests passing for the wrong reason, not code defects — the CLI tests
accepted exit 2, which a mutant run outside the repository also returns when it cannot find its
inputs; two fixtures had a default failing both the arm check and the column check, so disabling
either left a P2c standing; and P2c precedence had become untestable through Grok once D-c22
declared every allowlisted grok, so it is now proven on the Codex side. All five were closed by
tightening fixtures.

A narrow code jury on the v2→v3 diff (Codex `gpt-5.6-terra`/xhigh and Grok `grok-4.7`/xhigh, blind)
requested changes, and was right: **D-c13 was only half done.** "Codex-column row" was inferred from
the `gpt-` prefix and never read from the table's cells, so a Codex default named only in the Grok
cell passes silently — executed, both seats. My own code comment claimed the column semantics; it
was false, and is corrected in the commit to say what the code does. Both seats also caught that
`resolve()` follows symlinks, giving a false P3 on a correctly indexed symlinked devlog — executed.
Five smaller findings from Grok and Codex. On Christian's instruction v3 is committed with all seven
filed as D-c23..D-c29 and D-c13 reopened as partial. None fires on today's tree; live, the checker
reports one DRIFT — the index line owed for `dl20260922_ask_codex_cli_audit.md`.
