# A checker for convention documents that have drifted from practice

**Date:** 2026-09-22
**Purpose:** `scripts/check_doc_claims.py` verifies the convention documents against the **tree**;
nothing verifies them against **practice** — what the devlogs, the commit bodies and the wrapper
defaults actually do. Four drifts of that second class were found by hand on 2026-09-21/22 and two
are still live tonight. This plan adds `scripts/check_practice_drift.py`, the mechanical pass in the
practice → document direction.
**Module:** meta (`scripts`, `doc/ai_collaboration_protocol.md`)
**AIs involved:** Claude Opus 5 (plan, verification, measurements). Round 1: Grok `grok-4.7`/high
blind, 11 findings; **Codex produced no audit** (two usage-limit failures). Round 2: Codex
`gpt-5.6-terra`/high and Grok `grok-4.7`/high, both blind, on a fresh slug — 16 findings.
**Status:** **v9 committed** 2026-09-23 — D-c45 (the anchor accepts only space, tab or end of line
after `### 9.1`), D-c46 (a lone depth heading inside an unclosed fence fails closed, naming the opener),
D-c47 (the duplicate-anchor error names the lines). `--self-test` 92/92, 63/63 mutants. Round 9: logic
sound, no regression; two test gaps filed, D-c48 and D-c49. Chain: v9, v8 `9e74a6b`, v7 `d74a8ab`, v6
`fb1cc13`, v5 `42b5e42`, v4 `e295ce8`, v3 `cff668b`, v2 `fb74ed5`.

> **Overnight run rules (Christian, 2026-09-22, "make an overnight run out of this"):**
> - hard cap **v8** — at most three jury rounds (v6, v7, v8), each both seats at `xhigh`;
> - stop early when a round finds no silent pass and no regression;
> - anything needing a ruling is filed as an O-question and that line stops — no policy guesses;
> - a vendor quota failure stops the run and is reported; the reserve seat is not used overnight;
> - a version is committed only with the self-test green and every mutant killed; a finding that
>   is a regression against the committed version reverts that fix (the v4 rule);
> - own files only, private-index commits, a HEAD guard on every commit.

> **Scope guards:**
> - **Out of scope:** any probe needing prose parsed *for meaning*. A token search inside a prose
>   span is allowed and is load-bearing (see D2). Exactly one candidate probe dies on this rule; it
>   is O1. (O2 and O3 are ordinary open questions, not casualties of it — corrected from the first
>   draft, which miscounted all three.)
> - **Out of scope:** document ↔ tree claims. That is `check_doc_claims.py`'s class and it works.
> - In scope, named explicitly: §6 and §9.1 of `doc/ai_collaboration_protocol.md`;
>   `.claude/scripts/ask_codex.sh` and `.claude/scripts/ask_grok.sh` — allowlists **and** defaults;
>   `devlog/dl*.md` with `devlog/README.md`; recent commit bodies.
> - **`.claude/scripts/ask_copilot.sh` is out of scope.** It is the reserve seat (protocol §9.2),
>   not a jury member, and §9.1 defines no row for it, so there is nothing to compare against. Its
>   default `auto` would otherwise fire P2c against a policy the protocol has not adopted. Revisit
>   only if Copilot gains depth rows; Christian's ruling of 2026-09-22 marks it low priority.
> - **Token alphabet: `gpt-*` and `grok-*` only.** Devlog headers carry seat names §9.1 never lists
>   (`Claude Opus 5`, `Claude Fable 5.1`, and `Claude Opus 5` inside the §6 template itself). Every
>   model extractor uses this alphabet, not just the `Allowed:` scrape.
> - **Advisory by default.** `--strict` exits 1 on DRIFT only. INFO never gates — a checker that
>   blocks on a style frequency gets switched off in a week.

---

## 1. What Drifted, and What Is Still Drifting

> **HISTORICAL (reconciled 2026-09-22).** The first draft of this section described drifts A and B
> in the present tense, citing `**Module:**` as absent from the §6 template and the template's field
> as `**Topic:**`. Both were true on 2026-09-21 and **both were fixed hours before this plan was
> written** — by me, in the same session. `**Topic:**` now occurs zero times in
> `doc/ai_collaboration_protocol.md`. Grok caught it (D1). The corrected table below separates what
> is history from what a checker would catch tonight, because a plan justified by an already-fixed
> defect cannot be gated on reproducing it.

| # | Drift | Mechanism | State tonight |
|---|---|---|---|
| A | `gpt-6-astra` ruled the Codex seat for critical framework code 2026-09-18; §9.1 still said `gpt-5.6-terra` | the ruling was spoken, landed in an AI memory file, and the table was not edited. Auditors read the protocol and can see neither a skill nor a memory file | **FIXED** 2026-09-21 (`doc/ai_collaboration_protocol.md:466`) |
| B | §6's devlog header template described fields practice had left behind | each entry is written from the last entry, not from the template | **FIXED** 2026-09-21 (`:218-230`) |
| C | `ask_grok.sh:134` defaults the Grok seat to `grok-4.6`, which §9.1's rows name nowhere | the default was correct when pinned and nothing re-reads it | **LIVE** — named in recent devlog headers and commit bodies; the "14 records" figure sums two windows that may overlap and is not re-measured |
| D | seven `devlog/dl*.md` files have no entry in `devlog/README.md` | the index line is a separate manual step from writing the entry | **LIVE** — five of them weeks old |

**Bottom line:** each of these is a fact with two sources in the tree that no gate compares. A and
B were each individually self-consistent on both sides, which is why no review round would have
caught them. C is weaker as an argument — `ask_codex.sh:13`'s equivalent stale claim is a tracked
line in a reviewed script, and a reader *could* see it (D15). The honest claim is not "a review
would not have caught these" but "**no gate compares the two sides, and review has not**."

## 2. Architecture: The Tier Follows Normativity

The first draft claimed both probes were pure set difference with "no heuristics, no scoring", then
specified a `--threshold` percentage. Grok was right that these are incompatible (D4). There are
**two** kinds and they get different tiers:

| Question | Tier | Why |
|---|---|---|
| Is a **normative** relationship broken? | **DRIFT** | the document states a rule and the tree breaks it: a model used but undeclared, a wrapper default off the policy, a devlog outside its index |
| Is a **non-normative** convention drifting? | **INFO** | §6 states that `**Date:**` is the only universal field and the rest "are the ones that recur, not a schema". A session may add a field; reporting that as a defect contradicts the document the probe enforces |

The first draft said "set difference is always DRIFT, frequency is INFO". Codex refuted it: P1 is a
set difference *and* is correctly INFO (D-r2-11). **The tier follows normativity, not the shape of
the comparison.** One consequence, ruled by Christian 2026-09-22: the two directions of the template
probe get different tiers, because only one of them is normative — see P1a/P1b below.

Rejected: parsing prose for meaning (O1), and storing an expected-state file to diff against — a
second artifact to keep current is the disease, not the cure.

**Every probe anchors on a searchable token, never a line number** — the rule
`check_doc_claims.py` and `doc/input_schema.yaml` already follow.

**Every *required* extract must fail when it finds nothing — and the exception set is not one.** `scripts/check_doc_claims.py:31-35`: *"A probe that
finds nothing FAILS; it never skips… A checker whose failure mode is silently reporting success is
worse than no checker."* Three bugs of that shape were found in that script while it was written.
The first draft of this plan reproduced the defect its sibling exists to forbid (D3): an empty
extract has frequency 0, misses any threshold, and `--strict` exits 0.

Codex and Grok both refuted the first repair (D-r2-3, two raisers): asserting **every** set non-empty
is wrong, because `exception_models()` is an *exception* set. Empty means "no escape hatch", which is
a valid protocol and is exactly the clean fixture R6 needs — a universal guard fails that fixture
before any comparison, or forces a dummy exception into every clean document. `REQUIRED_FACTS`
protects facts the tree **must** contain. So:

| Extract | Guarded non-empty? |
|---|---|
| §9.1 table rows, both wrapper allowlists, both wrapper defaults | **yes** |
| devlog corpus, index entries, §6 template fields | **yes** |
| `exception_models()` — the anchored by-design list | **no**, empty is legitimate |

**Every reported number carries its window.** Tonight's re-count moved seven of the census figures
this plan was originally built on — not because anything was miscounted, but because concurrent
sessions added entries while it was being written:

| field | first draft | measured 02:37 |
|---|---|---|
| `**Purpose:**` | 33 | 34 |
| `**Module:**` | 27 | 28 |
| `**AIs involved:**` | 11 | 16 |
| `**Verification:**` | 5 | 9 |

A hand census of this corpus decays in hours. That is the strongest argument for the checker, and
the reason output states the window (how many files, ordered how, as of which commit) rather than a
bare percentage.

## 3. Probes

| id | Probe | Tier | Precedence |
|---|---|---|---|
| P1a | a §6 template field used by **no** recent entry | **DRIFT** | — |
| P1b | a practice field absent from the §6 template | INFO | — |
| P2a | model used in recent records, absent from §9.1's rows and from the by-design list | **DRIFT** | suppressed for a model P2c already reports |
| P2b | a §9.1 row names a model **its own seat's** wrapper would reject | **DRIFT** | — |
| P2c | a wrapper default that is not the model the policy names | **DRIFT** | takes precedence over P2a |
| P2d | allowlisted model named by no row and not on the by-design list | INFO | — |
| P3 | `devlog/dl*.md` file whose name is no index entry's **link target** | **DRIFT** | — |
| P4 | devlog title form | INFO | — |

**P1a/P1b are split by Christian's ruling of 2026-09-22.** A template field nobody writes is a
*fact* about a document describing something dead — normative, so DRIFT. A practice field missing
from the template is not, because §6 says the list is not a schema and a session may add a field.
The original motivating drift (`**Module:**` used and undeclared) is therefore P1b and advisory;
the reverse direction is what can fail. This was my pre-registered F1 and it is now decided, not
silently kept.

**P2c encodes the policy, not a lookup.** §9.1 says the Grok seat is *the latest allowlisted grok*.
Checking "the default string appears in some row" only coincides with that today (D-r2-12). The
check is: **the Grok default equals the greatest `grok-*` on `ask_grok.sh`'s own allowlist**, and
the Codex default is named by a row of the Codex column.

**P2b binds to a seat.** Set-differencing every row token against every allowlist fires on a
healthy table — `grok-4.7` is a row and is absent from `ask_codex.sh`, and every `gpt-*` row is
absent from `ask_grok.sh` (D-r2-7). Codex rows check against `ask_codex.sh` only; Grok rows against
`ask_grok.sh` only.

**P3 matches the link target, not the filename anywhere.** `devlog/README.md:18` links
`dl20260919_smith_inverse_branch.md` and *mentions* `dl20260919_smith_inverse_report_for_gregory.md`
later on the same line. A substring search silently drops that file (D-r2-8). Capture the bracket
target of `^(?:- |#{2,4} )\[([^\]]+)\]`, and match on the captured name alone. Both index forms
are matched — a `- [..]` line and a `### [..]` heading (`devlog/README.md:96`).

### 3.1 The discriminator — REDESIGNED after it was disabled by its own documentation

> **SUPERSEDED (2026-09-22).** Revision 2 read §9.1's **prose** for the exception list: a model
> named in a non-table line was treated as a declared escape hatch. Revision 2 also recorded that
> this "fires on `grok-4.6` alone, validated on live data". **Both claims died together.** At
> ~02:45, hours before that revision was written, the latest-grok policy paragraph was added to
> §9.1 with a sentence naming `grok-4.6` as a live defect. From that moment the prose contained the
> token, so the probe would have suppressed the one drift it was built to catch. The measurement was
> taken at 02:37 and copied into the plan at ~03:20 — true when taken, false when used, which is the
> decay §2 warns about, landing on this plan's load-bearing claim. Round 2 caught it; round 1 could
> not have, because the defect was introduced *between* the rounds by the fix to a round-1 finding.
>
> The line-number errors that came with it (`:480` was a fence, the sol sentence had moved, the span
> ran past `:499`) were the symptom. **The premise was the defect: §9.1 prose carries escape-hatch
> declarations *and* defect descriptions, and no token search separates them.** Re-measuring would
> have restored the probe only until the next drift was documented.

**The design now, ruled by Christian 2026-09-22: exceptions are declared structurally.**
`doc/ai_collaboration_protocol.md` carries one anchored line, and the probe reads that and nothing
else:

```
**Allowlisted but named by no row (by design):** `gpt-5.6-sol`
```

Anchor: the literal string `**Allowlisted but named by no row (by design):**`. Models are taken
from that line only. Prose is then free to name a stale model, quote a wrapper, or describe a defect
at length without blinding any probe — the protocol says so in place, directly under the line.

This also settles P2d, which in revision 2 reported `gpt-5.6-sol` while the definition of done
forbade reporting it (D-r2-9). P2d now consults the same by-design list, and its remaining INFO set
is the allowlisted-but-unruled tail: `gpt-5.5`, `gpt-5.4`, `gpt-5.4-mini`, `grok-4.5`.

> **RETRACTED (2026-09-22) — the usage-count discriminator.** Before the prose design, I proposed
> separating a defect from an escape hatch by how often a model is *used*. Refuted: a concurrent
> session used `gpt-5.6-sol` in a real round and recorded it in a devlog header
> (`devlog/dl20260922_test_gap_tier1_overnight.md:9-11`). Kept here so it is not re-proposed.

### 3.2 Extraction hazards, each found by measuring or by an auditor rather than by reading

- **The index has two forms**, and matching must capture the link target (above).
- **`Allowed:` appears twice per wrapper.** `ask_codex.sh:66` lists models, `:78` lists *efforts*
  (`low medium high xhigh`); same split in `ask_grok.sh:140` and `:148`. Filter to the `gpt-`/`grok-`
  alphabet, which §T of the scope guard now applies to every model extractor, not only this scrape.
- **The §6 fence is not the only fence,** and the fence guard can reject the good fence. The header
  example's own first line begins with `#` (it is a title inside the block), so a scan treating `^#`
  as a heading fails the assert on the very template it must read (D-r2-15). Take the first fence
  after the `### Required Header` anchor, be fence-aware, and assert it closes before the named
  successor heading `### When to Write a Devlog`.
- **The window is ambiguous.** Ten dates inside the last-40 window carry more than one entry. Order
  by filename date, then filename lexically, and print the resolved window.
- **Use versus mention.** Devlog *headers* only for P2a's devlog side. Commit bodies are
  **corroborating** evidence for a model already seen in a header, never a sole trigger — `cb86f2e`
  names models while describing a rule. `record_models()` therefore takes an explicit repo path or
  injectable commit texts, so the path is fixtureable (D-r2-13).
- **`**Topic:**` is not a "declared spelling".** It occurs in devlogs and zero times in the
  protocol, so Purpose-versus-Topic is not a contest between two declared tokens (D-r2-16). It is
  P1b: a practice field the template does not declare — INFO. Any count printed with it carries its
  window, and the plan states counts for the window measured, not for the corpus.

## 4. Ordered Steps

- [x] **R1** — extractors, each taking explicit path(s), each returning a set or Counter, no
      globals: `template_fields(proto)`, `devlog_fields(paths)`, `table_models(proto)`,
      `exception_models(proto)`, `wrapper_models(sh)`, `wrapper_default(sh)`,
      `record_models(paths, repo, n_commits)`, `index_targets(readme)`, `devlog_files(dir)`,
      `title_forms(paths)`. The wrappers are named in the caller: `ask_codex.sh` and `ask_grok.sh`
      only.
- [x] **R2** — the required-extract guard of §2: non-empty asserted for every extract **except**
      `exception_models()` (after: R1).
- [x] **R3** — P2a–P2d with the seat binding, the policy form of P2c, and P2c-over-P2a precedence
      (after: R2).
- [x] **R4** — P3, index coverage, matching link targets in both index forms (after: R2).
- [x] **R5** — P1a (DRIFT), P1b and P4 (INFO), each printing its window (after: R2).
- [x] **R6** — **the gate: fixtures decide, the live run reports.** Christian's ruling of
      2026-09-22. Synthetic fixtures under a `--self-test` flag are the only pass/fail: each probe
      must fire on a planted drift and stay silent on a clean tree, including a clean tree whose
      by-design exception list is empty. The live run prints findings and asserts **no** exact set,
      so repairing a real drift can never break the gate. Revision 2 required the live run to report
      "exactly `grok-4.6` and the seven unindexed devlogs" — an oracle that names today's defect
      list, lives only in sweepable scratch, and fails the moment someone fixes one of them
      (my pre-registered F2, independently raised by Grok).
- [x] **R7** — report format and exit policy: grouped by probe, DRIFT and INFO labelled, window
      stamped, `--strict` exiting 1 on DRIFT alone, `--recent N` (default 40), `--commits M`
      (default 40) (after: R3–R5).
- [x] **R8** — register in `CLAUDE.md`'s tooling block and in the `session-close` companion table.
      Cadence is O3's to decide, not this step's.

### 4.0 Audit record — rounds 1 and 2

Blind plan audit, Grok `grok-4.7`/high, 2026-09-22 01:45. Codex failed twice on a usage limit, so
every finding has **one** independent raiser; all citations were re-verified against source by me
and all held. D1–D4 were P0 and are addressed above. D5–D13 (P1) and D14–D17 (P2) are addressed in
§3.1, §3.2, §4 and §5. Full thread and reconciliation table:
`tmp/ai_exchange/review_practice_drift_checker.md`.

Findings that changed the **design** rather than the prose: D2 (the sol refutation → the row/prose
discriminator), D3 (the empty-extract guard), D4 (two comparison kinds instead of one).

### 4.1 Code-jury defects (round 3, blind, Codex `gpt-5.6-terra`/high + Grok `grok-4.7`/high)

All ten fixed in one rewrite. The evidence is mutation testing, not the green gate alone — D-c1's
lesson is that a green gate can be green for the wrong reason. New in the fix, beyond the ten:
**P2e** (the `Allowed:` message disagrees with the case arm that enforces it — a direct consequence
of D-c8, since reading the arm makes the message checkable), the policy grok is declared for P2a
and P2d even when no row names it (D-c10), `wrapper_default` is anchored on the seat's variable
name (my pre-registered F4), and the window stamps the HEAD commit. `--verbose` was removed rather
than given a meaning (D-c6).

All confirmed against source; two executed. Thread: `tmp/ai_exchange/review_practice_drift_code.md`.

- [x] **D-c1 (P0)** — `--self-test` does not gate P2a, P2b, P2d or P4. Deleting the P2a or P2b
      append leaves all twelve assertions green: the precedence check asserts only that the subject
      appears once, not which probe reported it. Two raisers.
- [x] **D-c2 (P0)** — `index_targets()` captures the markdown **label**, not the link target, while
      its docstring claims the target. Masked because `devlog/README.md` repeats the filename in
      both positions. `- [short name](file.md)` would be a false positive. Two raisers.
- [x] **D-c3 (P1)** — `record_models` counts occurrences per line; `devlog_fields` builds a per-file
      set. P2a prints occurrence counts as "recent devlog headers".
- [x] **D-c4 (P1)** — a failed `git log` becomes `commit_texts = ""` while the window line still
      says the commits were read: the vacuous-success mode §2 forbids. `UnicodeDecodeError`
      unhandled.
- [x] **D-c5 (P1)** — `_newest_grok` raises `ValueError` on a non-integer version piece; `main`
      catches only `ExtractionError`. Three raisers, including pre-registration.
- [x] **D-c6 (P2)** — `--verbose` is a no-op. Three raisers.
- [x] **D-c7 (P2)** — `_write_tree` never removes `extra_devlogs`, so fixtures leak into later cases.
- [x] **D-c8 (P2)** — the allowlist is read from the `Allowed:` echo, not the `case` arm; an arm
      edited without its echo is invisible.
- [x] **D-c9 (P2)** — `grok-4.7-build-fast` clips to `grok-4.7`, a depth-row model, hiding a name
      the wrapper deliberately excludes.
- [x] **D-c10 (P2)** — P2d is quiet on `grok-4.7` only because the `--quick` cell still names it.

### 4.2 Code-jury round 2 — open, for v3

Both seats at `xhigh`, blind, fresh slug. Thread: `tmp/ai_exchange/review_practice_drift_code_r2.md`.
Every item below was reproduced with a fixture or confirmed by reading; none fires on the live tree.

- [x] **D-c11 (silent pass)** — a model on the by-design exception line that no wrapper arm
      accepts is unioned into `declared` unchecked. Two raisers.
- [x] **D-c12 (silent pass)** — a Grok depth row naming an older grok passes while the default is
      correct. **Ruled 2026-09-22: DRIFT**; the protocol now says so.
- [x] **D-c13 (silent pass — completed by D-c23 in v4)** — the Codex-default P2c tests membership in the whole table, so a
      Grok-column token satisfies it. Must be the Codex column and the Codex arm.
- [x] **D-c14 (silent pass)** — P3 keeps only the basename of a link destination, so
      `../archive/dl….md` or a URL counts as indexing the local file. Resolve relative to the README.
- [x] **D-c15 (silent pass)** — a recent devlog without `**Date:**` is invisible in a 40-file
      window. **Ruled 2026-09-22: DRIFT**; the protocol now says so.
- [x] **D-c16** — `--recent`, `--commits`, `--threshold` are not validated; `--recent -1` is accepted.
- [x] **D-c17** — the window tiebreak is reverse-alphabetical (`reverse=True` reverses the name
      too), and the report names no file, so the window cannot be reconstructed. Print the cutoff.
- [x] **D-c18** — one self-test assertion is `check(..., True)` and cannot fail.
- [x] **D-c19** — P1b divides by headed files while its help text says the window.
- [x] **D-c20** — a missing protocol or devlog file raises uncaught `OSError`, not exit 2.
- [x] **D-c21** — an allowlist of only suffixed groks leaves `policy_grok` unset and P2c skips
      silently; the policy grok should be a required extract.
- [x] **D-c22 (false positive, surfaced 2026-09-22)** — once the `ask_grok.sh` default moved to
      `grok-4.7`, P2c stopped claiming `grok-4.6`, and P2a began reporting it: nine recent devlog
      headers and two commit bodies record rounds run on it *while it was the latest*. That is
      history, not drift. P2c's precedence had been masking it, so v2 could not show it while the
      default lagged. **Ruled 2026-09-22: for the Grok seat, any model on the wrapper's arm counts
      as declared for P2a**, because §9.1 keeps older groks allowlisted to reproduce earlier rounds.
      The default (P2c) and the depth rows (D-c12) must still be the latest, so real drift stays
      caught; P2a's Grok job narrows to "a grok not allowlisted at all". The Codex seat keeps the
      stricter rule — `gpt-6-astra` was allowlisted, used and unnamed, and that was a defect.

### 4.3 Code-jury round 3 (narrow, on the v2→v3 diff) — open, for v4

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v3.md`. v3 committed with
these open on Christian's instruction; none fires on today's tree.

- [x] **D-c23 (silent pass, two raisers)** — D-c13 was only half done. The Codex-default check reads
      "Codex row" from the `gpt-` prefix, not from the table's Codex cell, so a `gpt-` id written only
      in the Grok cell satisfies it — executed. Parse cell provenance: split each row on `|`, cell 2
      is Codex, cell 3 is Grok. Fixtures: a `gpt-` id only in cell 3 (must fire), only in cell 2
      (must stay silent). The same flattening lets a misplaced id evade P2b.
- [x] **D-c24 (false positive, two raisers)** — `index_targets` uses `Path.resolve()`, which follows
      symlinks, so a correctly indexed symlinked devlog is reported unindexed — executed. Normalize
      `.` and `..` lexically instead.
- [x] **D-c25** — the "P1c does not depend on P1a" assertion cannot fail; the clean devlog carries
      a Date, so P1a never fires on it anyway.
- [x] **D-c26** — the "failures are loud" test runs `quiet=True`, which silences the message it is
      named for. Assert stderr, or print extraction errors regardless of `quiet`.
- [x] **D-c27** — an all-headless window fails the `devlog headers` guard (exit 2) instead of
      reporting P1c for each file.
- [x] **D-c28 (REOPENED, then fixed in v5 — see §4.4)** — the header scan ends at the first `## `, including one inside an opening example,
      which can hide a later `**Date:**` and false-fire P1c.
- [x] **D-c29** — a malformed index destination can raise `OSError` in `resolve()` and abort the
      whole run; one bad line should not end the check.

### 4.4 Code-jury round 4 (narrow, on the v3→v4 diff) — open, for v5

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v4r.md`. v4 committed with
these open and D-c28 reverted to the v3 scanner. None fires on today's tree.

- [x] **D-c28, full specification (two raisers) — v5** — the v4 fence toggle fixed the simple closed fence
      only, missed nesting, and turned an unclosed fence into a silent pass (executed). A correct
      header scan: CommonMark fences (``` and ~~~, a closer at least as long as the opener with no
      info string, nesting), **only unfenced lines count as header** (a fenced `**Date:**` must not
      stand in for a missing one — a gap v3 shares), and an unclosed fence must fall back to the
      plain scan rather than swallow the body. Required fixtures: the three executed cases.
- [x] **D-c30 (two raisers) — v5; a malformed table fails closed, ruled 2026-09-22** — the column parser accepts ambiguous tables: duplicate seat headers,
      rows narrower than the header, doubled edge pipes (`strip("|")` shifts cells), and `|` lines
      anywhere in the §9.1 span. Require exactly one `Codex` and one `Grok` header and uniform width.
- [x] **D-c31 — v5, the dead conjunct dropped** — the D-c25 conjunct `not has(f, "P1a", "Date")` is still dead: one headed sibling
      keeps Date in use. The P1c half is live and carries the claim; drop or rebuild the conjunct.
- [x] **D-c32 — resolved by design, no change: a `grok-*` in the Codex column is already P2b against
      the Codex arm; adding P2c would make one defect two findings** — P2c's older-grok check reads only the Grok column; a `grok-4.6` in the Codex column
      is only P2b. Low.
- [x] **D-c33 — v5** — P2c says "depth row" where P2b says "column"; no fixture removes the Codex header;
      the "Codex column" guard message covers both a missing header and an empty set. Low.

### 4.5 Code-jury round 5 (narrow, on the v4→v5 diff) — open, for v6

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v5r.md`. Both silent passes
below were executed and are **shared with v4** — v5 is no worse anywhere tested. v5 committed.

- [x] **D-c34 (silent pass, shared with v4)** — when a fence opened in the header never closes, the
      fallback restarts the plain scan, which re-admits lines from fences that had already closed, so
      a fenced example `**Date:**` stands in for a missing one — executed. Fix per the jury: treat an
      unclosed opener as ordinary text and resume the fence-aware scan, keeping earlier closed fences
      excluded. Fixture: a closed fenced Date followed by an unclosed fence.
- [x] **D-c35 (silent pass, shared with v4)** — `_fence_closes` accepts any trailing characters
      `str.strip()` removes; CommonMark allows only spaces and tabs. A closer followed by NBSP ends
      the fence early and exposes a fenced Date — executed. Use `.strip(" \t")`; NBSP fixture.
- [x] **D-c36** — resolved by removal: a fence *after* the table cannot matter, because its opener
      is itself a non-`|` line that ends the block; the fixture was deleted and the reason stated in
      `_table_rows`. Original finding: the `fenced_row` fixture cannot fail: the blank line after the table ends the
      block before the fence is reached. Put the fenced block directly under the last row.
- [x] **D-c37** — the fixture named "an unclosed fence falls back instead of swallowing the body"
      accepts the swallow; `unclosed_then_date` is the discriminating one. Relabel or rebuild.
- [x] **D-c38** — `double_pipe` asserts only that no error is raised, not the column each model lands
      in; `fenced_before` would fail by traceback rather than a named FAIL.

### 4.6 Code-jury round 6 (narrow, on the v5→v6 diff) — the overnight run stopped here

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v6r.md`. Grok found a
regression against v5 (executed): routing the table reader through the demoting resolver made it adopt
an example table from inside a fence straddling the §9.1 span cut. Only `_table_rows` was reverted to
v5; a fixture now locks that an example table is never adopted. After the revert the round left no
silent pass and no regression, so **the overnight run stopped at v6** rather than the v8 cap. None of
the items below fires on today's tree, fails silently, or regresses v5.

- [x] **D-c39** — cut the §9.1 span at the first *unfenced* `##`/`###`, with fences resolved over the
      whole protocol file, so a fenced example straddling a heading is skipped and the real table is
      read, instead of failing closed as now.
- [x] **D-c40** — `_fenced_flags` is O(n²) worst case (one re-scan per stray opener). Low: Python on
      devlog-sized input.
- [x] **D-c41** — the fixture "a stray opener is text, so a ## after it still ends the header" is
      also satisfied by a scanner that swallows to EOF; relabel it or add a lock that separates them.
      The discriminating pair already exists.
- [x] **D-c42** — the premise was wrong: the form-feed case was NOT handled by the shared predicate,
      because `str.splitlines()` treats a form feed as a line break before the predicate ever sees it.
      `_lines` now splits on CommonMark line endings only (a no-op on all 494 live files, scanned).
      Original finding: untested: a tilde line inside a backtick fence, and a closer trailed by a form feed.
      Both are handled by the shared predicate; neither is locked.

### 4.7 Code-jury round 7 (narrow, on the v6→v7 diff) — open

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v7r.md`. No functional claim
in the diff refuted; Codex confirmed the one-pass resolver against the re-scan definition on 1,111,111
documents. v7 committed.

- [x] **D-c43 (two raisers, shared with v6 — not a regression)** — v8: the depth heading must occur
      exactly once outside fences, or the run fails closed (the existing malformed-table ruling), and it
      is matched exactly so `### 9.10` does not count. Narrower than either option below. Original finding: an *unclosed* fence placed before a
      forged `### 9.1` heading: whole-file demotion leaves the forged heading unfenced, the span starts
      there, the slice drops the opener, and `_table_rows` adopts the forged table — executed on v6 and
      v7 alike. Fix per Codex: keep source-indexed fence state across the slice, or resolve the anchor
      with a strict CommonMark map in which an unclosed fence runs to EOF and extraction fails closed.
      Fixture: an unclosed opener immediately before a forged section, a real section after it.
- [x] **D-c44** — v8. the `unclosed_example` fixture accepts any `ExtractionError`, so an anchor miss or an
      empty table would also pass it. Assert the specific failure.

### 4.8 Code-jury round 8 (narrow, on the v7→v8 diff) — open

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v8r.md`. No functional claim
refuted, no regression (all executed or traced against v7). One Grok claim refuted: the span end does not
match `#### Note`. v8 committed.

- [x] **D-c45 (two raisers, not a regression)** — v9. Original finding: `DEPTH_HEADING` uses `\s`, which accepts a form feed
      after `### 9.1` where the file's own rule is spaces and tabs. Use `(?:[ \t]|$)`; fixture beside
      the `### 9.10` one. v7's `startswith` matched this and `### 9.10`, so v8 is already narrower.
- [x] **D-c46 (shared with v7)** — v9: a lone §9.1 heading that a strict reading places inside an
      unclosed fence fails closed, naming the opener's line (the malformed-table ruling). Original finding: with exactly one exposed `### 9.1`, an unclosed opener just above it
      is outside the slice, and `_table_rows` reads the first `|` block after the heading. Only the
      two-candidate shape fails closed.
- [x] **D-c47** — v9. Original finding: the duplicate-anchor error reports a count, not the line numbers. Cosmetic.

### 4.9 Code-jury round 9 (narrow, on the v8→v9 diff) — open, test gaps only

Both seats at `xhigh`, blind. Thread: `tmp/ai_exchange/review_practice_drift_v9r.md`. Logic sound, no
functional claim refuted, no regression; Codex judged the fail-closed scope of D-c46 correct. v9 committed.

- [ ] **D-c48** — `DEPTH_HEADING` accepts a tab or the end of line after `### 9.1`, but only the space
      form is exercised: mutants "space only" and "no end-of-line" both survive (executed). Add passing
      fixtures for `### 9.1\tDepth Selection` and a bare `### 9.1`, and both mutants to the suite.
- [ ] **D-c49 (two raisers)** — the duplicate-anchor and unclosed-fence assertions check the phrases,
      not the line numbers; assert the fixture's expected lines.

## 5. Open Design Questions

- **O1** — `ask_codex.sh:13` claims "Only terra and luna are named by the depth table"; false since
  the astra row landed. Verifying that class of claim needs prose understood, so it is out of scope
  for the checker. **Recommend deleting the claim from the comment and citing §9.1 instead** — a
  comment restating a table is the defect this session already fixed twice. Not yet done: the file
  is a one-line executable change and the wrapper is currently modified by a concurrent session.
- **O2** — should P1 and P3 read `nonfree/devlog/`? It follows the same §6 template. **Leaning no**:
  that tree is absent from most checkouts and the script must not fail when it is missing. Pin the
  glob to `devlog/dl*.md` explicitly, or the lean is lost in implementation.
- **O3** *(unchanged, still open)* — cadence. Every session close is too often; the corpus moves slowly and a finding would
  repeat until someone acts. **Proposal, not a decision:** run when the session touched
  `doc/ai_collaboration_protocol.md`, `CLAUDE.md`, a wrapper or `devlog/README.md`, plus monthly.
  Nothing in `session-close` runs on a monthly clock, so the monthly half needs a home.
- **O5** — my pre-registered F3, F4, F5 and F6 were not addressed by either seat in round 2 and
  remain unaudited: commit-body corroboration may lose a real case; `--recent 40` is still
  underived; P2c may be too strict for a deliberate pin; and both live drifts the plan cites as
  justification are this team's own housekeeping lapses rather than errors discipline cannot catch.
  Carried, not closed.
- **O6 — RESOLVED 2026-09-22 → suppress P1a when no entry in the window has a header (Christian).**
  An all-headless window now reports P1a DRIFT for every template field ("used by none of
  the last N entries") alongside P1c for each file. P1a applied literally, but it accuses the
  template when the entries simply have no header. Keep, or suppress P1a when the window has no
  headed entry? Christian's call.
- **O4** — should P3 also flag the reverse, an index entry with no file? Measured tonight: zero such
  cases, so the probe would be dead code today. **Leaning yes anyway** — it is one set difference in
  the opposite direction and catches a renamed or deleted entry.

## 6. Interface Constraints

No new file format. One constraint: **every extractor takes its path(s) as arguments** so R6 can
point them at fixtures rather than the live tree. The first draft asserted this in prose while R1's
signatures took none of it (D6).

## 7. Definition of Done

- [ ] Every probe in §3 maps to a step; every §3.2 hazard has a guard in the code.
- [ ] Every probe anchors on a token, not a line number — including the by-design exception line.
- [ ] Every extract except `exception_models()` is asserted non-empty before use.
- [ ] `--self-test` fixtures: each probe fires on a planted drift and is silent on a clean tree,
      and the clean tree includes one with an **empty** by-design list.
- [ ] A fixture proves P2c-over-P2a precedence: one model, one finding, not two.
- [ ] A fixture proves P2b's seat binding: a healthy table produces no P2b finding.
- [ ] A fixture proves P3 matches link targets, not substrings: an index line that *mentions* a
      second filename does not count as indexing it.
- [ ] The live run prints findings and asserts no exact set.
- [ ] `--strict` exits 1 on DRIFT, 0 on INFO alone.
- [ ] Runs clean with `nonfree/` and `literature/` absent.
- [ ] Every reported count prints its window.

## 8. Audit Trail

- Thread: `tmp/ai_exchange/review_practice_drift_checker.md` (pre-registration, Grok audit,
  verification pass with the reconciliation table, measurements, two Codex failures, timestamp
  correction).
- **Reduced independence, recorded:** Codex hit a usage limit at 01:45 and again at ~02:30, quoting
  the same fixed 03:30 reset both times, so the credits added between attempts were most likely
  never spent. One voice audited this plan. Owed before implementation: a **blind two-voice jury on
  this revision**, Codex `gpt-5.6-terra`/high and the latest Grok at high — §9.1's plan-audit row.
  Not the astra row: this is tooling, not framework code.
- **Ruling, 2026-09-22 (Christian):** a Python tooling script under `scripts/` takes the full
  plan+audit → code+audit round. ~~Owed in `.claude/skills/code-with-jury/SKILL.md`'s exemption
  table, which currently says nothing about `scripts/*.py`.~~ **Already done** — the row is at
  `SKILL.md:24`, written the same day, before this sentence claimed it was missing (D-r2-5).
- **Round 2, both seats blind, fresh slug so neither could read round 1:** Codex
  `gpt-5.6-terra`/high and Grok `grok-4.7`/high. 16 findings, four P0, **three with two raisers** —
  the stale §9.1 citations, the invalid universal non-empty guard, and the rotting live gate. Round
  2 found a P0 round 1 could not have: the discriminator was disabled *between* the rounds by the
  fix to a round-1 finding. Thread: `tmp/ai_exchange/review_practice_drift_checker_r2.md`.
- **Reserve seat:** `.claude/scripts/ask_copilot.sh` may fill a seat whose vendor is temporarily
  unavailable (protocol §9.2, ruled 2026-09-22). It did not fill Codex's empty seat in round 1 —
  the wrapper was created at 02:57, after that round ran at 01:45.
