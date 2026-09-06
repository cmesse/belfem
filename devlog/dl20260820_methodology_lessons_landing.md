# Devlog 2026-08-20 — the lessons layer lands: `doc/lessons_learned.md`, METHODOLOGY.md, and a narrowed prose-sweep rule

**Date:** 2026-08-20
**Purpose:** Record the landing of `doc/lessons_learned.md`, the METHODOLOGY.md additions that describe it, the two corrections made to the drafts before placing them, and the narrowing of the Codex prose-sweep rule to user-facing documentation
**Module:** repository conventions

---

## `doc/lessons_learned.md` landed

388 lines: a 141-line trigger-indexed tripwire layer over 18 evidence-backed cards, distilled
from 298 devlog entries and the debt register (539 catalogued incidents, 20 clusters).

Run scaffolding was stripped on the way in, because a `doc/` file states its facts inline rather
than citing a run's working artifacts: the header's pointers to the catalog, cluster and
enforcement files were replaced with a standard document header; two in-text references to
working files were reworded to state the fact directly; and the distillation's `LOCKED`
pre-registration stamp was removed, since it is a property of the run and not of the document.
The full-fidelity copy, including the lock and its correction note, remains the audit trail.

Wired in so it is actually read rather than merely present:

- `CLAUDE.md` Minimum Session Compliance Checklist gains a Layer 1 read as its second item,
  directly after the protocol read.
- `doc/README.md` §Guidelines and Philosophy gains its index entry.

Status on every card is `proposed`. The open adjudications from the distillation — whether
"presence is not execution" is a distinct cluster, and the five rules proposed by the external
auditors — remain open, and landing the file does not settle them; new cards are additive.

## `doc/lessons_learned_evidence.md` landed, and the Doxygen site grew a section for both

The 539-row incident catalog landed as `doc/lessons_learned_evidence.md` (580 KB) so that the
`INC-NNN` identifiers the cards cite resolve to something durable instead of to a swept
scratch directory. Same scaffolding treatment: run-stage title and spec references replaced
with a standard header, working-file pointers reworded to state their facts inline, and the
per-agent id and chunk-file columns dropped from the table of repaired rows, which were an
artifact of how the mining was parallelised and carry no evidential value.

**Landing the two files broke `scripts/update_doc_index.py`, which is how the site's page tree
is kept honest.** `Doxyfile.in:1046` sweeps all of `doc/` into the Doxygen INPUT, and the index
generator refuses to file a page that has no reader-path entry — it errored with
`no PAGE_GROUPS entry for: doc_lessons_learned, doc_lessons_learned_evidence`. That guard is
doing exactly its job: a document cannot land on the site unfiled.

Resolved on Christian's suggestion by giving them their own reader path rather than mixing
them into the user-facing groups:

- `scripts/update_doc_index.py` — new `PAGE_GROUPS` entry, **"Code development"**, alongside the
  existing "Run a simulation", "Extend a module" and "Project process" paths.
- Both files carry explicit Doxygen anchors (`{#doc_lessons_learned}`,
  `{#doc_lessons_learned_evidence}`) matching the ids the generator derives.
- `doc/mainpage.md` regenerated; `--check` reports 0 files out of date.
- `doc/README.md` gains both entries under §Guidelines and Philosophy.

The reader-path split is what answers the "will this confuse a user?" question: the landing
page now separates *how to run a simulation* from *how this code gets developed*, and the
lessons pages sit in the second.

---

## Summary

`METHODOLOGY.md` grew from 41 to 106 lines. Christian drafted the additions; this session
placed them, after checking the drafts against the artifacts rather than accepting them.

Landed, in file order:

- One sentence in the intro naming `doc/lessons_learned.md` as the method's distilled output.
- **Five new terminology rows**: agreement-as-weakest-rung mapped to correlated failure in
  multiversion programming; jury-brief contamination mapped to clinical-trial unblinding; blind
  slice assignment mapped to independent coding; and two entries with **no published analogue
  found** — the blindness trust boundary including the tooling, and recurrence-triggered
  promotion of a rule into a mechanism.
- **A vocabulary note** distinguishing our "verified" (an executable gate ran) from
  verification in the V&V sense and in formal methods. The collision is three-way, not two-way.
- **Deviations 5 through 8**: agreement is bounded evidence; the blind includes the tooling;
  mechanical checks carry their own falsifiers; written rules graduate into mechanisms.
- **A closing section, "The method, measured"**, reporting the distillation's own numbers and
  the shape of its weaknesses.

`scripts/check_doc_claims.py` passes 34/34 after the edit.

## Two corrections made before placing

Both are recorded because the drafts were good and the errors were the interesting part.

**1. A miscount that had already propagated.** The literature cross-reference summarised its own
findings as "5 cards have strong published lineage" while naming four, having counted one card's
two independent lineages as two cards. The downstream METHODOLOGY draft inherited the figure as
"five of the extracted rules". Corrected to four in both places, and the cross-reference now
carries an in-place correction note rather than a silent fix.

**2. An overclaim about where the unmatched rules cluster.** The same summary asserted that the
six rules with no published analogue "cluster tightly" around epistemic hygiene. Enumerated,
they do not: three fit (the tree is the evidence, scope every absence, re-open every literal),
one is adjacent (sole-cause), and two are ordinary engineering lessons that simply fall outside
the three books consulted — restart is not continuation, and naming what breaks before deleting
the redundant. For two of the six, adjacent literature exists and was **not** consulted. The
claim was a universal statement made without enumerating its members, which is the defect the
scope-every-absence rule exists to prevent, committed by the document reporting that rule.

## One conflict between the two drafts, resolved toward the stricter side

The cross-reference's patch list supplies the Knight and Leveson terminology row **with** the
study's figures (versions tested, test count, confidence level). The METHODOLOGY draft's own
placement note says those figures must wait until the primary paper is checked, and its
deviation 5 deliberately omits them.

Applying both as written would have landed exactly the figures the note says to withhold. The
row was placed **without** the figures, stating only the qualitative result — that the
independence assumption was tested and rejected — which is what the secondary sources support.
The figures can be restored in one edit once the primary paper is read.

## Forward references: one resolved, one open

`METHODOLOGY.md` had two statements that were not yet true when it was written.

1. ~~**`doc/lessons_learned.md` does not exist.**~~ **Resolved the same session** on Christian's
   instruction to land it. The four references (lines 10, 77, 87, 106) now resolve.
2. ~~**The debt register still has no closure fields.**~~ **Resolved pre-commit by rewording**:
   deviation 8 now says the firing signal "is designed to come from closure fields in the debt
   register's protocol; those fields are proposed alongside the rules and adopt with them" —
   a true statement of the proposed state that does not preempt the adjudication of the
   register amendment.

## Pre-commit sweep: the two index defects the evidence file names were fixed before shipping it

Committing `doc/lessons_learned_evidence.md` with INC-538 and INC-539 still open would have
published two known, trivially fixable defects as open findings. Both fixed:

- **INC-539** — the dangling `dl20260520_tycho_macos_build_fix.md` link (a file that never
  existed in this repository; tycho is a `nonfree/` module and its devlogs live there) is
  removed from `devlog/README.md`.
- **INC-538** — the eight devlog entries on disk but absent from the index are backfilled, each
  at its chronological position, each marked *(backfilled 2026-08-20; INC-538)* so the index
  does not pretend they were indexed on time. The greg3 jury round got the full-length entry
  its content warrants.
- The verification then caught its own tool being wrong: the first index-vs-disk diff reported
  two further missing entries, both of which were indexed all along — the check's character
  class lacked capitals and missed stems containing `dT` and `lowT`. The check was fixed and
  the reconciliation now reads **0 unindexed, 0 dangling** with a correct pattern. Stage 0 of
  the distillation used the same style of stem pattern, which is worth knowing when reading
  its coverage claims.

The incident rows themselves stay in the evidence file unchanged — they record what was true
when mined; the fixes belong to this entry, not to a rewritten row.

## The prose-sweep rule, narrowed twice in one session

Christian's ruling, in two steps. First: the sweep is for documentation intended for a *user* to
read, so `METHODOLOGY.md`, `doc/lessons_learned.md` and the devlogs are excluded — the sweep
dilutes, and smoothing costs information per line, which is a bad trade for a document whose
value is its density. Second, after the `./todo/` question was raised: todo files are excluded
too, as not worth the token cost.

Applied in two documents, in precedence order:

- `doc/ai_collaboration_protocol.md` §6 — the standing todo-file prose rule is **withdrawn** and
  says so, with the date and the reason. Technical audits of a todo file's content are
  explicitly unaffected. The protocol outranks `CLAUDE.md`, so the withdrawal has to live here
  to be real.
- `CLAUDE.md` — the compliance-checklist line and §"Prose Gets a Language Sweep" both narrowed
  to user-facing documentation, with a three-item exclusion list and the dilution rationale.

An intermediate state existed and was wrong: when only the first ruling had been given, a
pointer was added to `CLAUDE.md` saying `./todo/` files **keep** their sweep. That pointer was
corrected in the same session rather than left to rot. The stored memory asserting the old todo
rule was deleted, since a remembered rule that contradicts the tree is worse than no memory.

## Owed

- ~~A **Codex prose sweep** over the added prose.~~ **Ruled the same session, and the rule
  itself was refined.** Christian's ruling: the sweep is for documentation intended for a *user*
  to read. `METHODOLOGY.md`, `doc/lessons_learned.md` and the devlogs are excluded, because the
  sweep dilutes — smoothing costs information per line, which is a good trade for a guide read
  once to learn a module and a bad trade for a document whose value is its density. Applied to
  `CLAUDE.md` in three places: the compliance checklist, the opening of §"Prose Gets a Language
  Sweep", and its closing exclusion list. `./todo/` files are explicitly **not** swept up in the
  exclusion — they keep their own standing rule in `doc/ai_collaboration_protocol.md` §6, which
  outranks `CLAUDE.md`, and a pointer saying so was added so the new wording cannot be misread
  as cancelling it. `check_doc_claims.py` 34/34 after the edit.
- Verification of the Knight and Leveson figures against the primary paper (item V1 of the
  cross-reference's verification list) before the figures are restored to the terminology row.

## Terminology review (Christian's ruling, same session)

The METHODOLOGY.md terminology table was reviewed term by term:

- **devlog** — kept as-is, by preference; the published mapping (daybook / ADR) stands beside it.
- **jury, pre-registration, relay** — already adopted working vocabulary; no change.
- **blackboard** — adopted as a **spoken synonym** for the `tmp/ai_exchange/` channel; the
  directory is deliberately **not** renamed. `ai_exchange` is descriptive too, and the rename
  would have touched ~10 live files, the delta hook (`.claude/hooks/ai_exchange_delta.sh` +
  its `settings.local.json` registration), and the `AI_EXCHANGE_SLUG` env var for no
  functional gain. Recorded in protocol §2 ("see blackboard" resolves to the current task's
  exchange file) and in the METHODOLOGY table, both with the one honest nuance: during blind
  jury rounds the blackboard property is deliberately suspended — auditors must not read each
  other's postings until reconciliation, which is the opposite of the published pattern's
  opportunistic sharing.
- **read-only diagnosis** — made the official working name for the investigation mode, stated
  in protocol §7 where the practice is defined, with hypothesis-driven debugging as the
  published analog.

`check_doc_claims.py` 34/34 and the doc index clean after these edits.

## Files changed

- `doc/lessons_learned.md` — **new**, 388 lines (141-line Layer 1, 18 cards, annex).
- `doc/lessons_learned_evidence.md` — **new**, 580 KB, the 539 catalogued incidents.
- `doc/README.md` — index entries for both under §Guidelines and Philosophy.
- `doc/mainpage.md`, `scripts/update_doc_index.py` — new "Code development" reader path.
- `METHODOLOGY.md` — 41 → 106 lines.
- `CLAUDE.md` — Layer 1 read added to the compliance checklist; prose-sweep rule narrowed to
  user-facing documentation with a three-item exclusion list and the dilution rationale.
- `doc/ai_collaboration_protocol.md` — §6 todo-file prose rule withdrawn; §2 gains the
  blackboard synonym; §7 names the read-only diagnosis.
- `devlog/README.md` — index entry for this file.
- `tmp/lessons_distill/07_literature_crossref.md` — in-place correction to its summary paragraph.

`check_doc_claims.py` passes 34/34 after every edit above. Nothing is committed.

`todo/debt_register.md` also shows as modified in `git status`; that change is from the
2026-08-19 evening cross-review pass (mtime 21:49) and is not from this session.
