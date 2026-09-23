---
name: session-close
description: Close a BELFEM session the house way — the dated devlog, its devlog/README.md index line, the same-session companions (campaign page, debt register, todo checkboxes, input contract), then the commit in the repo's own message style. Load it when the user says to close, wrap up or finish the session, asks for a devlog or a commit, says "we are done", "log this", "commit this" or "write it up", or when substantial work is finished and nothing else is pending.
---

# Session close

Two artifacts, in this order: the **devlog** (durable record) and the **commit** (the change).
Neither substitutes for the other — the devlog holds the reasoning and the evidence, git holds the
diff and the who-and-when.

There are 475 entries in `devlog/README.md`. Match them; do not invent a format.

## 1 — Survey before writing

```bash
git status --short && git diff --stat HEAD && git log --oneline -5
```

In a shared checkout other sessions have live edits. Separate **your** files from theirs before
anything else — you will stage explicitly, never `git add -A`.

Then answer, for yourself: what question opened the session, what actually landed, what evidence
ran, what is still open. That is the devlog.

## 2 — The devlog

`devlog/dlYYYYMMDD_topic.md`, lowercase_with_underscores, titled with a topic phrase. Work scoped
to `nonfree/` writes to `nonfree/devlog/` instead and does **not** touch `devlog/README.md`.

**The header block is specified in `doc/ai_collaboration_protocol.md` §6 "Required Header" — use
it, do not work from memory.** `**Date:**` is the only always-present field; `**Purpose:**`,
`**Module:**`, `**AIs involved:**`, `**Status:**`, `**Verification:**`,
`**Literature References:**` and `**Exchange:**` are filled where they apply. On the
`**AIs involved:**` line, name each auditor seat's **model and effort** and which jury it sat on —
a finding is attributable only if the record says how hard the auditor was asked to look.

Sections are free-form beyond `## Summary` — `## What prompted it`, `## What landed`,
`## What the checker found`, `## Open Questions`, `## Files Updated` all appear. Use what the
session shaped.

Two things the devlog must do:

- **Lift the conclusion out of `tmp/ai_exchange/<slug>.md` before it is swept.** That channel is
  ephemeral and `.gitignore`d. What was confirmed, refuted or left open — with `file:line` evidence
  and confidence — is lost at the sweep if it lives only there. Name the slug in `**Exchange:**`.
- **Say "reviewed" unless a gate ran.** Name the gate and its result in `**Verification:**`, or say
  plainly that the work is static-only. AI agreement is the weakest rung and is never written up as
  verification.

## 3 — The index line

**Prepend** to `## Entries` in `devlog/README.md` — newest at the top:

```markdown
- [dl20260921_topic.md](dl20260921_topic.md) — **the one claim that matters, in bold.** Then the
  dense detail: numbers, file paths, what was refuted, what is carried. Ends with the evidence
  status ("Reviewed, not verified", or the gate that ran).
```

These lines are long and information-dense on purpose. This index is how a future session finds
the work; a vague line loses it. Do not smooth it — the density is the point, and `devlog/` takes
no prose sweep.

## 4 — Same-session companions

Skip none of these that apply; they are owed in the session that did the work, not later.

| If the session… | Also update |
|---|---|
| opened or closed a tracked item | `devlog/campaigns/<name>.md` **and** `todo/debt_register.md` — rows are exactly 6 cells, rebuilt cell by cell; quote the live/archive counts in the commit body |
| worked from a `todo/` plan | tick that step's checkbox, strike obsolete steps (never delete), refresh the Status line |
| touched any `input.conf` key | `doc/input_file_reference.md` **and** `doc/input_schema.yaml`, then a Codex prose pass over the touched `.md` sections |
| edited `CLAUDE.md`, `doc/coding_philosophy.md`, `CMakeLists.txt` or the compiler config | `python3 scripts/check_doc_claims.py` |
| added a citation | `python3 scripts/check_citations.py` |
| edited `doc/ai_collaboration_protocol.md`, a wrapper, or `devlog/README.md` | `python3 scripts/check_practice_drift.py` — reports where the convention documents and practice have parted |
| added user-facing documentation | the Codex language sweep (not for `devlog/`, `todo/`, `METHODOLOGY.md`, `doc/lessons_learned.md`) |

## 5 — The commit

**Commit with a pathspec, in one command with the add.**

```bash
git add <your new files> && git commit -F <message-file> -- <all your paths>
```

`git commit -- <paths>` commits exactly those paths and ignores everything else in the index. That
is the protection that matters: the index is shared per worktree, so without the pathspec you
commit whatever another session happened to stage.

A new file has to be `git add`ed first — the pathspec form cannot commit an untracked path — so
chain the two in a single command and keep the staged window to seconds. On 2026-09-22 a `git add`
of two new skills, a devlog and two doc edits sat staged while the message was written; a
concurrent session committed, its commit absorbed all of it under a message about `DofMgr
BlockData`, and a later history rewrite in that session orphaned the commit. The work survived
only because it was still in the working tree and copied to the scratchpad first.

`git add -A` is worse again: it commits a peer's half-finished work under your message.

**Check that HEAD has not moved** before and after (`git rev-parse HEAD`). In an active checkout it
does, and a commit written against a stale view of the tree is a merge problem for someone else.

**Subject:** `<Area>: <what the code now does>` — lowercase after the colon, specific, present
tense, describing the new state rather than the act of changing it. Two clauses joined by `;` or
`,` when one will not carry it.

```
CutFactory: 2-D terminal-curve orientation reads the master's edge, checks its preconditions
Debt register: strike DR-119 on its notice run, trip half waived
```

**Body: prose paragraphs, not bullet lists.** This is the house style and it is deliberate — a
paragraph can carry a mechanism, a bullet list cannot. One paragraph per finding or per change,
opening with its id when it has one (`F25.`, `DR-119.`). Say what was wrong, by what mechanism, what
the fix does, and what pins it.

Four things the body records when they happened, because they are the audit trail:

- **The juries** — which models at which effort, and *what each one changed*: "both approved the
  deletion, the plan jury split the identity check into its own test, the code jury cut two
  overclaims from the comments."
- **The gates** — exactly what ran and its result: "make check on both trees (17 of 17), fresh
  Valgrind on dipole and tapestack3d at two ranks, no record from the site."
- **What was filed** for a later session, and who raised it.
- **Register counts** after a strike: "Live register 5, archive 148."

**Trailer** — name the model that did the work, per the session's attribution instruction:

```
Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
```

Commit when closing is asked. **Push only when the user asks** — and never commit onto `main`;
`devel` is the working branch.

## 6 — After the commit

If the post-commit autoreview hook is installed (`BELFEM_AUTOREVIEW=1`), check
`tmp/ai_exchange/AUTOREVIEW_P0.flag` and report any P0 rather than leaving it in the file.
`scripts/review_status.sh` shows commits since the last auto-review.

## 7 — Sign off with the end marker

**The final message of a close ends with this line, and nothing follows it:**

```
--- THE END ---
```

It exists for one reason: several sessions run at once in separate terminal windows, and a
glance at the bottom of a window has to answer "is this one finished?". A marker is enough —
the client does not need to exit, so nothing is killed and the scrollback survives.

Two rules keep it honest:

- **Only on a genuine close.** Not at the end of a long answer, not when work is merely paused.
  It means: the devlog is written, the companions are updated, the commit is made, and nothing
  is pending on my side.
- **Anything unfinished is named directly above it.** Files deliberately left uncommitted, a
  gate that was not run, a ruling still owed. The marker says "I am done", never "everything
  is done" — a close with residue is still a close, but the residue is listed first.

## Never committed

`.gitignore` excludes `nonfree/`, `archive/`, `/literature/`, `tmp/` (so the whole exchange tier),
`build/` and `cmake-build-*/`. `.claude/` **is** tracked, minus `settings.local.json` and the
runtime state files — so a new skill or command belongs in the commit, and a peer's WIP skill
does not.
