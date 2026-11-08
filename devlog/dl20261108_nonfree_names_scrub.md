# Nonfree module names scrubbed from the open-source tree

**Date:** 2026-11-08
**Purpose:** Christian's ruling, after a stray `git add -A` staged a plan that documented the nonfree
test layout: the names of the proprietary modules do not appear in the open-source repository at
all, not in source, docs, devlogs, todo files or the debt registers. This entry records what was
found, what was reworded, and what stays.
**Module:** repository-wide (no code behaviour touched; one source comment)
**AIs involved:** Claude Fable 5.1 (inventory and rewording); no jury (prose only, no executable change)
**Verification:** `git grep -i` over the tree outside `nonfree/` after the edit: the only remaining
match is a university name in `doc/literature_references.md`, which is a thesis citation and not a
module. `scripts/check_doc_claims.py` 38/38 after the `CLAUDE.md` edit; `check_practice_drift.py`
reports only its two pre-existing P1b/P4 items.

## Summary

The inventory found the three module names in 17 committed open-source files, in 24 places:
`CLAUDE.md` and `doc/ai_collaboration_protocol.md` (a staging-tree path in the nonfree exception),
one source comment (`src/visualizer/cl_VTK_Curve.hpp`), the devlog index and nine devlogs, two
evidence rows in `doc/lessons_learned_evidence.md`, two closed todo plans and the live debt register
row DR-201. Every mention is reworded to "the nonfree suites", "a nonfree module", "the proprietary
spacecraft-thermal application" or "the nonfree Python extension", keeping the fact each sentence
carried. The DR-201 plan that triggered the ruling moved to `nonfree/todo/closed/`, and the DR-201
devlog was split by the session that wrote it, its nonfree half going to `nonfree/devlog/`.

Kept: the thesis citation naming a Linz university in `doc/literature_references.md`; it is a place, not a module.

## What remains

- The struck DR-201 row in `todo/debt_register_closed.md` still names the suites; it is the other
  session's uncommitted edit and is being reworded there before its commit.
- **History.** Three commits already on `origin/devel` and `origin/main` (`2e190f7`, `08c306b`,
  `f81dd84`) and three more local ones carry the names in devlog and todo text. Removing them from
  history means rewriting pushed branches; that is Christian's decision, not done here.
- No mechanical guard exists yet. A checker in the style of `scripts/check_practice_drift.py`,
  run from the pre-commit hook, would refuse a commit that introduces the names outside `nonfree/`;
  it owes a tooling jury under the 2026-09-22 ruling.
