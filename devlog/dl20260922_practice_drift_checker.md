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

## v4

Run on autopilot at Christian's request. D-c23..D-c27 and D-c29 fixed, D-c13 completed: the depth
table is now read by its own column headers (`Codex`, `Grok`) and each column is checked against its
own seat's wrapper arm, so a Codex default named only in the Grok cell is P2c and a model in the
wrong column is P2b; index destinations are normalised lexically, so a symlinked devlog is no
longer a false P3 and a malformed destination cannot end the run; an all-headless window reports P1c
instead of failing the guard; two tests that could not fail were rebuilt. `--self-test` 64/64 and
38/38 mutants killed on the first run. Live output identical to v3 — the real columns agree with
v3's prefix inference on today's table.

The narrow jury (Codex `gpt-5.6-terra`/xhigh, Grok `grok-4.7`/xhigh, blind, on the v3→v4 diff)
requested changes, and the load-bearing one was a **regression I introduced**: the D-c28 fence fix
was a boolean backtick toggle. Executed against v3 and v4 — it fixed a simple closed fence, missed
a nested one, and turned an **unclosed** fence from a correct P1c into a silent pass by swallowing
the body. For a checker whose rule is that a silent pass is worse than a loud one, that is the wrong
trade, so D-c28 was **reverted to v3's scanner** (juried and shipped) and reopened with the full
specification. Everything else in v4 was traced clean and stays. After the revert: 63/63, 37/37.

Filed for v5: D-c28 in full (CommonMark fences, unfenced header lines only, unclosed-fence
fallback), D-c30 (strict column parsing — two raisers), D-c31..D-c33 (a still-dead test conjunct,
the older-grok check reading only the Grok column, wording and coverage), and O6 — whether an
all-headless window should also report P1a on every template field. That last one needs Christian.

The pattern across four code rounds is worth stating: every round found real defects, and three of
the four found at least one in the *previous round's fix*. Mutation testing caught weak tests, but
it can only kill mutants I think to write — none of mine modelled a nested or unclosed fence. The
jury found what the mutants could not imagine.

Process slip: the v4 jury was dispatched with `&` inside a compound command, not as its own
background task. It survived; a watcher was attached rather than re-dispatching a paid round.

## v5

Two rulings from Christian first: P1a is suppressed when no entry in the window has a header (O6), and
a malformed §9.1 table fails closed (D-c30). Then autopilot.

D-c28 done properly: CommonMark fences (``` and ~~~, a closer at least as long with nothing after it,
nesting, the backtick info-string rule), only unfenced lines count as header — which also closes the
gap v3 shared, a fenced `**Date:**` standing in for a missing one — and an unclosed fence falls back
to the plain scan. D-c30: the depth table is the first contiguous `|` block outside fences, parsed
with exactly one edge pipe removed per side, one `Codex` and one `Grok` header, uniform width. The
fence logic is one shared helper, so the header scanner and the table reader cannot drift apart.

Gates: `--self-test` 80/80, **49/49 mutants killed**. One survived the first run — the only
unclosed-fence fixture put its Date where the fallback and a bare unfenced scan agree — and was
closed by a discriminating fixture plus a mutant that reproduces v4's body-swallowing regression
exactly. Live output unchanged; the real table passes the strict parser.

The narrow jury (both seats, xhigh) found two more silent passes, both executed: a later unclosed
fence makes the fallback re-admit an earlier *closed* fence's Date, and a closer followed by a
non-breaking space is accepted where CommonMark allows only spaces and tabs. **Neither is a
regression** — v4 fails both identically — so v5 is committed and they are filed as D-c34..D-c38
with three test weaknesses. The fallback's docstring claimed it "errs toward reporting a missing
field rather than hiding one"; the jury disproved that, and the docstring is corrected.

Across five code rounds the useful discipline has been one question asked of every finding before
deciding what to do: *is this a regression against what is committed?* In v4 the answer was yes, and
the fix was reverted. In v5 it is no, and the version ships with the gaps filed.

## v6 — the overnight run, and why it stopped at v6

**For the morning: nothing needs a ruling.** v6 is committed; the run stopped at v6 of a v8 cap
because the last round left no silent pass and no regression. Residue D-c39..D-c42 is low-severity.

Christian set the run going before leaving ("make an overnight run out of this"); the rules written
into the plan were: cap v8, stop early on a round with no silent pass and no regression, file anything
needing a ruling, stop on a vendor quota failure, commit only with the gates green, and revert any fix
that regresses the committed version.

v6 fixed D-c34..D-c38. The header scan now resolves fences with one pass per stray opener: an opener
that never closes is re-read as text and the scan runs again, so fences that closed earlier stay
excluded — v5's fallback restarted a fence-blind scan that re-admitted a closed fence's `**Date:**`.
Closers accept only spaces and tabs after them, per CommonMark. One mutant survived the first run in
an instructive way: with demotion in place, removing the closer-length check was *rescued* — the short
closer ended the outer fence early, the outer closer became a stray opener, demotion made the answer
come out right. A fixture where early closing exposes a fenced Date closed it. 83/83, 50/50.

The round-6 jury (both seats, xhigh) found a **regression against v5**, executed: I had also moved the
§9.1 table reader onto the new resolver — a consolidation no finding asked for — and on a span cut by
a raw heading inside a fenced example, demotion exposed the example table as the depth table. v5 fails
closed there. By the v4 rule the table reader was reverted to v5; the header fix scans the whole file
and keeps its gain. A fixture locks "an example table is never adopted". My first attempt at a mutant
for that regression was ineffective — it filtered lines but fed them into a loop that re-detected the
fence — and the honest mutant is v6's actual reader substituted verbatim, which the fixture kills.
83/83, 51/51.

Codex's other finding, to restore the post-table fenced-row fixture, was refuted by Grok and by
reading: a fence opener after the last row is itself a non-`|` line and ends the block regardless.

Lesson, again and sharper: both of this night's regressions (v4's fence toggle, v6's table reader)
came from changes *beyond* what a finding asked for — a general mechanism applied one place further
than needed. The findings themselves were fixed cleanly both times.

## v7 (2026-09-23)

D-c39..D-c42, each kept to its finding's scope — both earlier regressions came from a fix reaching one
place further than asked. `_depth_span` now bounds §9.1 by unfenced headings, fences resolved over the
whole file, so a fenced example straddling a heading is skipped and the real table read instead of
failing closed; `_table_rows` stays v5's loop, as the v6 revert left it. The fence resolver is one
pass: a suffix table of the longest bare closing run per fence character decides at once whether an
opener has a closer, instead of re-scanning per stray opener. The bare-run rule became one helper,
`_bare_run`, after writing D-c40 had left a second copy of it inside the drift checker itself.

**D-c42's premise was wrong.** Its fixture failed: a closer followed by a form feed never reached the
fence predicate, because `str.splitlines()` breaks lines at form feed, vertical tab and six more
separators. `_lines` now splits on CommonMark line endings only. That touches every extractor's input,
so it was gated before it was kept: none of those characters occurs in 494 live files or the wrappers,
and every live file splits identically.

Gates: `--self-test` 88/88, including the one-pass resolver checked against the old re-scan definition
on 3000 random documents; 57/57 mutants killed on the first run. Two v6 mutants were dropped because
their anchor now lived only in the self-test's reference oracle — they would have mutated the oracle.

The narrow jury (both seats, xhigh) refuted no functional claim. Codex compared the resolver with the
re-scan definition on **1,111,111** documents without a mismatch. Both seats raised one case —
an unclosed fence before a forged `### 9.1` makes the table reader adopt the forged table — and it
was executed on v6 and v7: **identical, so not a regression**; v7 additionally fixes the closed-fence
variant that v6 got wrong. Filed as D-c43 with one test weakness, D-c44. v7 committed.

## v8 (2026-09-23)

D-c43 and D-c44. For the forged-anchor case both seats raised in round 7, v8 took neither suggested
fix: the §9.1 heading must occur exactly once outside fences, or the run fails closed — an application
of Christian's ruling that a malformed table fails closed, not new policy. In the reported case the
stray opener is demoted, the forged and the real heading are both unfenced, and the run stops with
"appears 2 times" instead of picking one. Uniqueness needed the heading matched exactly: the old
`startswith("### 9.1")` would also have counted a future `### 9.10`. The unclosed-example fixture now
passes only on the empty-table failure it is about. 90/90, 59/59 on the first run, live identical to v7.

The narrow jury refuted no functional claim and found no regression. Both seats noted the new anchor
pattern's `\s` accepts a form feed where the file's rule is spaces and tabs — confirmed, but v7's prefix
match accepted that and `### 9.10` besides, so v8 is strictly narrower. One Grok claim was refuted by
execution: the section-end test does not match `#### Note`. Filed D-c45..D-c47, all minor.

## v9 (2026-09-23)

D-c45..D-c47. The depth anchor accepts only a space, tab or end of line after `### 9.1` (a form feed no
longer counts). D-c46 resolved a disagreement inside the checker: the heading search read a stray fence
opener above §9.1 as text while the table reader failed closed on it. Per the malformed-table ruling,
v9 fails closed — a new helper, `_open_fence_at`, scans strictly up to the heading, and a fence still
open there stops the run with the opener's line named. The shared resolver was left alone. The
duplicate-anchor error lists its lines. 92/92, 63/63 on the first run, live identical to v8.

The jury found the logic sound and refuted nothing; Codex judged the fail-closed scope right. Two test
gaps, confirmed: the tab-only and bare heading forms are unpinned (two mutants survive), and the
line-number assertions check phrases, not numbers. Filed D-c48, D-c49. For the first time a round's
findings are about the tests alone.

## v10, and done (2026-09-23)

Test-only — every hunk inside `self_test()`, checked mechanically, so the jury ran at §9.1's code-diff row
(`high`) rather than `xhigh`. D-c48 pinned the tab and bare forms of the depth heading, killing the two
mutants that survived round 9; D-c49 asserts the exact line numbers in the two fail-closed messages,
computed from the fixture text with a plain newline split so the oracle shares no code with the line
reader under test. 94/94, 68/68.

**The first clean round in ten.** Codex approved with no defects; Grok refuted nothing and hand-traced the
same line numbers the self-test asserts. Its three open risks all fail loud — a reshaped fixture goes red,
an empty oracle list aborts — so none was filed.

Over ten versions the checker went from twelve self-test cases and no mutants to 94 cases and 68 mutants,
and every probe has a fixture that must make it fire and one that must keep it silent. Two regressions
were caught by the jury and reverted, and both came from a fix reaching one place further than its finding
asked. The plan is closed.

## The last DRIFT fixed, and the mutation suite given a home (2026-09-23)

The one DRIFT the checker still reported — no index line for `dl20260922_ask_codex_cli_audit.md` — was a
line I had destroyed in the README wipe; the entry's own "Files Updated" lists the README. It was
reconstructed from the devlog and committed with the devlog itself, which had never been committed
(`09d6289`). The checker now reports no DRIFT and `--strict` exits 0.

The mutation suite moved out of sweepable `tmp/` to `scripts/tests/mutate_check_practice_drift.py`, next to
its subject, on Christian's ruling. Its review ran with one voice — Codex hit its usage limit — and Grok
found, confirmed by execution, that **the suite has no baseline**: against a checker whose self-test
always fails, it reports all 68 mutants killed and exits 0. Every earlier run was genuine, because the
real self-test passes; but the suite cannot tell. A crash in an unrelated fixture also counts as a kill.
Christian ruled to fix the baseline and the file encoding before committing. The suite now runs the
unmutated self-test first and stops with exit 2 if it is red — proven on the same always-red checker,
which went from "68/68, exit 0" to exit 2. That guard immediately surfaced something new: the checker's
own `--self-test` is locale-dependent, writing its fixtures without an encoding, so it fails under a C
locale while the live run does not. Filed with the remaining findings as D-m2, D-m4..D-m6 in
`todo/practice_drift_mutation_suite.md`; a narrow re-review is owed when Codex returns. The instrument
built to catch false greens had one of its own — found, as every other one here was, by the jury rather
than by its author.
