# Removing DR- Tracker Tags from the Source

**Date:** 2026-09-01
**Purpose:** Record the sweep that removed DR- task indices and author-attribution
phrases from code comments in `src/` and `tests/`, and the convention behind it.

## The rule

DR numbers are **ephemeral**. They live in `todo/debt_register.md`, in `todo/` plans and
in `devlog/` entries, where a reader has the register at hand and the number resolves to
something. In code they do the opposite of what they promise: a comment reading
`( DR-127 )` tells the reader that an explanation exists elsewhere without giving it, and
after the row is struck and archived the tag points at nothing. Same for
`( Christian's ruling )` — the code does not care who decided, only what was decided and
why.

**Convention, effective 2026-09-01:** no DR- indices, plan/round step IDs, or ruling
attributions in `src/` or `tests/`. The explanation stays, rewritten to be
self-contained: "the historical stray message", "before the fix", "the original defect",
"what this check exists to close". If a comment needed the tracker to be understood, the
comment was underwritten, and this is where it gets fixed.

## What was changed

Comments only — no executable line, signature, or preprocessor gate was touched. The
sweep was applied as exact-match replacements with a per-site count assertion, so a
missed or duplicated match aborted the run rather than editing blind.

| Tree | Sites | Files |
|---|---|---|
| `src/` | 77 | 36 |
| `tests/` | 46 | 18 |

Heaviest concentrations: `cl_FEM_Controller.cpp` (17, nearly all DR-127 certified-exit
labels), `cl_FEM_DofManager.cpp` (6), `cl_SolverMUMPS.cpp` and
`cl_JcFunction_Database.hpp` (4 each).

Three classes of rewrite:

1. **Drop the parenthetical.** `// longer reroute safely ( DR-135 ). The` →
   `// longer reroute safely. The`. The majority.
2. **Drop the label, keep the sentence.** `// DR-113: the solve/distribute exchange ends
   here` → `// the solve/distribute exchange ends here`.
3. **Replace the tag with the thing it named.** Where the DR *was* the noun, the defect
   is now described: `was DR-129` → `was a real defect`; `where DR-100's stray was born`
   → `where the historical stray message was born`; `the DR-98 sign gate` → `the
   curl-sign gate`. In `cl_FEM_DofMgr_EigenValues.hpp` a cross-reference to "Solver
   ( DR-141 )" became "Solver ( `cl_Solver.hpp` )" — a pointer that stays true.

Orphaned sub-indices went with their parent: `( DR-144 G11 )`, `( DR-128, F3 )`,
`( DR-127, plan D4 )` were removed whole, since `G11` or `D4` alone is worse noise than
the tag was. `( DR-106 G-B, observed failing ranks 7, 4, 1, 6 )` kept its observation.

Three comment blocks were reflowed after the deletion left a ragged short line.

## Not touched, deliberately

- **`todo/`, `devlog/`, `todo/debt_register.md`** — the tracker tier. 1302 and 2028
  occurrences respectively; DR numbers belong there.
- **`doc/`** — 24 DR- references and 6 ruling attributions remain, in
  `input_file_reference.md`, `input_schema.yaml`, `parallel_execution.md` and
  `lessons_learned_evidence.md`. `lessons_learned_evidence.md` is an evidence ledger and
  the citations may be load-bearing; the other three are user-facing and probably should
  be swept the same way. Owed decision.
- **Cross-review bookkeeping in code** — `( round-3 C5 )`, `( round-3 R-G )`,
  `( round-5 audit )` and similar survive at 13 sites, mostly in
  `cl_FEM_Controller.cpp`. Same noise class, not in this sweep's scope.
- **Prose dates** — `Until 2026-08-30 every entry here carried the VOLT exponents` is a
  factual statement about the code's history and was kept. Only the stamp form
  (`! DR-155, 2026-09-01:`) was removed.

## Posture

**Reviewed, not verified.** Comment-only edits, confirmed by diffing every changed line
against a comment-marker filter: the only non-comment lines in the working-tree diff
belong to pre-existing uncommitted work (`cl_SpMatrix.cpp`, the `BELFEM_OMP` gating in
`splinalg.f90` / `arpacktools.f90` / `parpacktools.f90`). Nothing was built or run.

The edits are **uncommitted**, in the Darwin working tree at
`/Users/christian/codes/belfem` on `main` at 7323d528 — 57 modified files under
`src/` and `tests/`. They exist on no branch and on no other machine.

The live peer sessions were notified of the convention. One of them, on a Linux
checkout at bc6fb67d, correctly reported 127 DR- hits still present and challenged
the claim that the trees were clean — my first message said "have just been swept
clean" without saying *where*, which reads as landed. Corrected to all three. The
counts reconcile: 129 occurrences at 7323d528 in `src`+`tests` across 50 files, 135
at bc6fb67d (two commits older), against 123 edit *operations* here — several
operations spanned multi-line blocks or lines carrying two tags, and the peer's grep
covered only `.cpp`/`.hpp` where this sweep also reached a `.f90`, a `CMakeLists.txt`
and one `.md` under `src/fem/kernel/doc/`.

The convention is recorded **only here and in the session transcripts**. Nothing was
written to `CLAUDE.md` or `doc/coding_philosophy.md`; until it is, a fresh session
will reintroduce the stamps. That edit is Christian's call.
