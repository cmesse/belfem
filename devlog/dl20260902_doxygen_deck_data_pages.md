# Two stray Doxygen chapters removed: the deck `data/` READMEs

**Date:** 2026-09-02
**Purpose:** Record why `block3d plugin data` and `undulator2d plugin data` stopped being
top-level pages of the generated site, and what now keeps deck-internal markdown out of it.

## Symptom

The generated site carried two chapters under **Run a simulation**, beside Getting Started and
the Input File Reference:

- `doc_examples_block3d_data` — "block3d plugin data"
- `doc_examples_undulator2d_data` — "undulator2d plugin data"

Neither is a manual chapter. Both are a note left in `examples/<deck>/data/` telling whoever is
already editing that deck's plugin how to drop a measured two-column table beside it and read it
with `usermat::Table`. Addressed to someone standing in the directory; useless as a chapter.

## How they got there

Not by a decision, by two mechanisms meeting:

1. `scripts/update_doc_index.py::collect_examples()` walked **all** of `examples/` and treated
   every `.md` it found (other than `CLAUDE.md`/`AGENTS.md`) as a documentation page — assigning
   an anchor, writing that anchor into the file's heading, and requiring a `PAGE_GROUPS` entry.
2. `Doxyfile.in` has `INPUT` reaching `examples/`, `RECURSIVE = YES` and `*.md` in
   `FILE_PATTERNS`, so any markdown under that tree becomes a page whether it carries an anchor
   or not.

The two `data/README.md` files landed on 2026-08-31 with the usermat migration. The generator's
`--check` then failed with `no PAGE_GROUPS entry for: doc_examples_block3d_data,
doc_examples_undulator2d_data`, and the session that met that error **satisfied** it — filing both
under "Run a simulation" — rather than asking whether the pages should exist. That silenced the
check by promoting the files, which is the wrong direction: the check exists to stop a real page
going unfiled, not to conscript every markdown file in the tree.

## Fix

Three edits, because removing only one of them puts the pages back:

| File | Change |
|---|---|
| `scripts/update_doc_index.py` | `NOT_DOCUMENTATION_DIRS = {"data"}`, pruned from the `os.walk` of `examples/`; the two `PAGE_GROUPS` entries removed |
| `Doxyfile.in` | `EXCLUDE_PATTERNS` gains `*/examples/*/data/*` (plus `*/CLAUDE.md`, `*/AGENTS.md`) |
| the two `README.md` files | injected `{#doc_examples_*_data}` anchors stripped |

`doc/mainpage.md` was then regenerated and lost the two `@subpage` lines.

The Doxyfile half is not redundant. Stripping the anchor alone leaves an **auto-labelled** page —
`md_examples_2block3d_2data_2README`, still titled "block3d plugin data", still listed under
Related Pages. The generator half is not redundant either: without it the next
`update_doc_index.py` run re-injects the anchors and re-fails `--check`.

`*/CLAUDE.md` and `*/AGENTS.md` were added in the same edit because they are the identical leak
one directory over — `examples/corc_solder/python/CLAUDE.md` arrived with the solder example on
2026-09-01 and would have rendered as a chapter titled "CLAUDE.md" on the next `make doc`. The
generator has excluded those two filenames since it was written; Doxygen did not.

## Evidence

- `python3 scripts/update_doc_index.py --check` → **0 files out of date**; 27 modules, 99 pages
  (was 101 — the two data pages).
- The exclusion was run against the **real** `examples/` tree in a scratch Doxygen configuration
  (`INPUT = examples`, `RECURSIVE = YES`, the patched `EXCLUDE_PATTERNS`, `EXTRACT_ALL = YES`):
  markdown-derived pages produced are exactly `doc_examples.html` and `doc_examples_scripts.html`,
  while the deck plugin sources still get pages (`block3d_2src_2matlib_8cpp.html`,
  `corc__solder_2src_2current_8cpp.html`, …). So the pattern removes the two chapters without
  taking the example sources with them. `examples/*/data/` currently holds nothing but the README,
  so nothing else falls under it.
- `make doc` was **not** run — the user runs builds. The scratch Doxygen run is the gate.

## Why the rule is directory-scoped rather than depth-scoped

A depth rule ("anything more than one level below `examples/`") is the more general statement of
the intent, and the generator could express it. Doxygen cannot: its `EXCLUDE_PATTERNS` wildcard
matches `/`, so `*/examples/*/*/*` would also swallow `examples/<deck>/src/*.cpp`, which **is**
documented today. Keeping both halves on the same directory-name rule is what lets them stay in
step. If a second deck-internal directory appears, add its name to `NOT_DOCUMENTATION_DIRS` and a
matching `EXCLUDE_PATTERNS` line — both, or the page comes back.

Documentation and tooling only; no `src/` change. Reviewed plus the scratch Doxygen gate above.
