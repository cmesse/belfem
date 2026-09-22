# The workflow skills: code-with-jury, session-close, and two rules promoted out of memory

**Date:** 2026-09-22 (session opened the evening of 2026-09-21)
**Purpose:** Christian asked whether "plan+audit -> code+audit" — the phrase he repeats at the
start of every coding task — could be a skill. It could, and writing it down exposed two rules
that existed only in Claude's memory files and contradicted the documents the auditors read.
**Module:** meta (`.claude/skills`, `CLAUDE.md`, `doc/ai_collaboration_protocol.md`)
**AIs involved:** Claude Opus 5 (1M context), sole voice. No jury round: the whole session is
text, which the 2026-08-28 ruling exempts.
**Status:** COMPLETE. Two skills written and registered, two protocol rules promoted, four memory
files trimmed. Nothing committed yet — three of the touched files carry a peer's hunks (§6).
**Verification:** static only. `scripts/check_doc_claims.py` 38/38 after each documentation edit;
no executable gate applies to text. Reviewed, not verified.
**Exchange:** none — no auditor was dispatched, so no exchange file was opened.

## Summary

Two skills now sit in `.claude/skills/`, beside the `literature` skill a peer is still writing:

- **`code-with-jury`** — the outer loop. Its content was never in the repository: `/cross-review`
  has always specified *one round*, but the sequencing around it (plan, plan jury, implement, code
  jury, close out), the exemption table, and the depth choice lived only in four memory files, so
  the loop fired on recall rather than on the task. The skill's description matches any non-trivial
  change to executable code, so it loads whether or not the phrase is spoken. Two approval stops:
  after the plan jury (plan *and* edit approval, since read-only is the default) and after the code
  jury (fixes are a separate approved step). Dispatches in between run unattended, because the
  depth is table-driven and there is nothing to decide at that moment.
- **`session-close`** — the devlog and the commit, extracted from practice rather than invented:
  protocol §6, 475 index lines in `devlog/README.md`, and the last 25 commit messages with bodies.

Both are registered in `CLAUDE.md`'s Cross-Review Tooling block so a session that never triggers
them still sees they exist.

## The commit-body finding

The single most valuable thing the extraction produced. **Commit bodies in this repository are
prose paragraphs, not bullet lists** — one paragraph per finding, opening with its id, saying what
was wrong, by what mechanism, what the fix does, and what pins it. The `Code findings tier 5` and
`tier 6` messages are the reference. Left to its defaults an AI writes `- Fixed X` / `- Added Y`,
which would have degraded the style of 500 commits one session at a time. The skill also pins the
four things a body records because they are the audit trail: the juries with their models, efforts
and *what each one changed*; the gates with their results and counts; what was filed for later and
who raised it; and the live/archive register counts after a strike.

## Two rules that were memory-only

Writing the loop down surfaced a class of defect neither review tooling would have caught, because
both artifacts were individually self-consistent.

**The astra rule.** Christian ruled on 2026-09-18 that critical framework code — mesh, kernel, dof
manager, sidesets, blocks, calculator, basis, distributor, maxwell factories — takes
`gpt-6-astra` on the Codex seat for both the plan audit and the code audit. That ruling reached a
memory file and nothing else. For three days the depth table in `doc/ai_collaboration_protocol.md`
§9.1 still said `gpt-5.6-terra` for exactly that subject, and §9.1 is the document Codex and Grok
themselves read — skills are Claude-Code-only and no auditor can see one. The row is now at
`doc/ai_collaboration_protocol.md:448`, mirrored in the `/cross-review` excerpt (three rows became
four), with a bullet under the table explaining why it breaks the rung discipline stated two
paragraphs above it: §9.1 says the second rung raises effort and not model, and this row moves the
model, so without the note it reads as drift instead of the deliberate exception it is.

I introduced a fourth copy of that table in the skill before noticing, and the copies immediately
disagreed — the skill said `astra`, the command said `terra`, for the same subject. That is the
argument for the rule the session settled on, below. The skill now carries a pointer, not a table.

**The devlog header template.** §6's template and practice had diverged far enough that the
template described almost nothing anyone writes. Measured over the last 40 devlog entries:
`**Date:**` 40, `**Purpose:**` 33 against `**Topic:**` 7, `**Module:**` 27 — absent from the
template altogether, and the largest gap — `**AIs involved:**` 11 against `**AIs:**` 2 and
`**Voices:**` 1, `**Status:**` 6, `**Verification:**` 5, `**Claude Confidence:**` 4,
`**Codex Audit Confidence:**` 3, `**Exchange:**` 2. On titles, 206 of 475 files use the template's
`# Devlog <date> — <topic>` form, but only one of the last ten.

§6 now specifies practice: a topic-phrase title, `**Purpose:**`, `**Module:**`, `**Status:**` and
`**Exchange:**` added, and `**AIs involved:**` kept under its own more common name but showing each
seat's model, effort and which jury it sat on — `AIs involved: Claude, Codex` does not survive
contact with §9.1, since a finding is attributable only if the record says how hard the auditor was
asked to look. Per-finding confidence is stated as belonging in the body beside the finding it
qualifies, with the header form left available. The framing line changed too: `**Date:**` is the
only universal field and the rest recur rather than compose a schema, because the old "fill only the
applicable fields" read as a closed list while practice invents fields freely. The ~206 older
titles are explicitly not to be churned.

## The rule the session settled on

Three layers now describe overlapping ground, and they are read by different audiences:
`CLAUDE.md` and `doc/` by Claude, Codex, Grok and humans; skills by Claude only; memory files by
Claude on recall. **Skills carry procedure; `CLAUDE.md` and `doc/` carry contract. A skill may cite
a table those files own; it must not restate one.** The `literature` skill already worked this way —
it says so in as many words and routes its first step back to CLAUDE.md's table — and
`code-with-jury` was brought into line by deleting its duplicate depth rows.

The corollary matters more than it looks: the tempting cleanup of deleting CLAUDE.md's literature or
coding-standards material now that skills exist would blind the jury to precisely the material it
needs, since `AGENTS.md:46` sends the auditors to `CLAUDE.md` and the wrapper preamble feeds them
`AGENTS.md` plus the protocol. Overlap across that boundary is two audiences, not duplication.

## Rulings

- **Namespace.** A `belfem-skills:` prefix is plugin-derived; there is no namespace for project
  skills. Packaging as a plugin was declined for now: `.claude/skills/` ships with the git
  checkout, so every worktree and every fresh clone has the skills with no install step, while a
  plugin keeps its install state in `~/.claude`. Revisit when the same toolchain is wanted in a
  second repository — `nonfree/`, or a new project — which is when the prefix stops being cosmetic.
- **Memory trimmed, not deleted.** `three-vendor-audits`, `comments-free-code-jury`,
  `probes-skip-audit` and `astra-for-critical-code` keep the ruling, its date, Christian's words and
  the reason; the mechanics they used to carry now live in the skill. They still earn their place:
  they fire on recall in sessions where the skill never triggers, and they are the authority behind
  the skill's content for anyone who was not in the room for the ruling.

## Open questions

- `CLAUDE.md`'s `### 3. Literature-First Approach` is a third copy of the paper roster that
  `doc/literature_references.md` owns in full with DOIs, and eight of its entries already appear in
  `## Literature and References` above it — roughly 50 lines loaded into every session for a
  duplicate. Collapsing it to a pointer is the obvious win, but it is inside the `literature`
  skill's blast radius while that work is live, so it was flagged and not touched.
- Whether `doc/ai_collaboration_protocol.md` owes a Codex language sweep. CLAUDE.md scopes the
  sweep to user-facing documentation and excludes `METHODOLOGY.md` and `doc/lessons_learned.md` as
  descriptions of the method written for us; the protocol is the same kind of document but is not on
  the exclusion list. Not dispatched. Christian's call.

## The index race at close

The commit did not go as the new skill described, and the skill was corrected from what happened.
`git add` of the six paths was followed by a separate `git commit`; between the two, a concurrent
session committed, and its commit `9a721cf` ("DofMgr BlockData: one predicate selects a thin shell's
facets...", 01:27 PDT) absorbed both skills, the CLAUDE.md hunk, the protocol edits and this devlog
under its own message — the index is shared per worktree and takes whatever is staged. That session
then rewrote history, leaving `9a721cf` no longer an ancestor of HEAD, so the work was uncommitted
again; it survived only in the working tree, and in a scratchpad copy made before anything was
touched. HEAD moved twice more during the close.

`session-close` now says to commit with the pathspec form, `git commit -F <msg> -- <paths>`, which
ignores the index entirely, and to check `git rev-parse HEAD` before and after. This is the second
rule in this session that came from a collision rather than from reading a document.

## Files updated

- `.claude/skills/code-with-jury/SKILL.md` (new)
- `.claude/skills/session-close/SKILL.md` (new; pathspec-commit rule added after the index race)
- `CLAUDE.md` — both skills registered in the Cross-Review Tooling block
- `doc/ai_collaboration_protocol.md` — §9.1 critical-code row and its rationale bullet; §6 header
  template rewritten to practice
- `.claude/commands/cross-review.md` — critical-code row mirrored, row count corrected
- memory: `three-vendor-audits`, `comments-free-code-jury`, `probes-skip-audit`,
  `astra-for-critical-code` trimmed to ruling plus reason, `MEMORY.md` hook updated
