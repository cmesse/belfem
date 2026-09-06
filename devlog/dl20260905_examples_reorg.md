# Examples directory reorganized for release

**Date:** 2026-09-05
**Purpose:** Record the renaming, merging and documentation pass over `examples/` before the 0.9.0 release, and the audit round that shaped it.

---

## Decision

Christian asked for a category-based reorganization (2d, 3d, user materials, user boundary
conditions, defect functions) and for the `*_christian` / `*_gregory` names to go. I
inventoried the 19 decks, proposed a `2d/`/`3d/` split, and put it to Codex (gpt-5.6-terra,
high) and Grok (grok-4.6, high) as a round-1 plan audit (`tmp/ai_exchange/examples_reorg.md`).

Outcome: **flat directory, descriptive names, feature matrix in `examples/README.md`.**
Dimension is the only attribute that partitions the set (6 × 2D, 13 × 3D); every other
category overlaps (`undulator2d` is 2D + plugin, `tape_quench_usermat` is 3D + thermal +
two plugins). A directory level would also have broken the six per-deck `Allrun`/`Allclean`
wrappers (`exec ../scripts/…`), the `examples/*/input.conf` glob in `python/belfem_conf/cli.py`,
possibly the `*/examples/*/data/*` Doxygen exclusion, and every first-run path in three
documents, for a 6/13 split the README table gives for free. Grok argued flat; Codex accepted
the split conditionally; Christian chose flat and picked the names.

## Audit corrections worth keeping

Three claims in my proposal were wrong and both auditors caught them: the two `tapestack3d`
`.geo` files were byte-identical at `numTapes = 8` (the 2-tape file belongs to `circuit`; I had
misread my own diff), both 2D tapestack decks use `builtin : ybco` so "builtin" was a false
contrast, and `tape_quench_usermat` builds two plugins, not three (the current source is
linked into `userdefect.so`). Grok also found that `pancake` shipped `make_mesh.py`, which the
launcher's `./*/main.py` discovery never finds although the README said it did.

## Changes

Renames (`git mv`; all 167 tracked files followed):

| before | after |
|---|---|
| `corc_gregory` | `corc_periodic_bc` |
| `garber` | `corc_garber_effect` |
| `corc_christian` | `corc_twolayer` |
| `tapestack2d_gregory` | `tapestack2d_racetrack` |
| `tapestack2d_christian` | `tapestack2d_layered` |
| `tapestack3d_gregory` + `tapestack3d_christian` | `tapestack3d` (gregory kept: it carries the explicit `mumps error analysis : false`) |
| `circuit` | `tapestack_circuit` (Christian wrote `tapestack-circuit`; underscore used for consistency, flagged) |
| `block3d` | `racetrack_usermat` |
| `pancake/python/make_mesh.py` | `pancake/python/main.py` |

Deleted: `tapestack3d_christian/`, `tape_quench_usermat/defect.cpp` (68-line older copy of
`src/defect.cpp`, not built), `corc_twolayer/python/CLAUDE.md` (said "this repository"),
`pancake/python/make_mesh_box.py` (obsolete variant).

**Every `input.conf` re-indented with tabs** and its `key : value ;` columns aligned per
block (scratch script `retab.py`; one tab per brace depth, trailing comments aligned). Gate:
the `belfem_conf` parser produces an identical tree (headers, key names, values, nesting) for
all 18 decks before and after. The C++ reader strips tabs in `clean_string` (`stringtools.cpp:107`).

**CORC trademark block** added to `corc_garber_effect/input.conf`; the other two CORC decks
already had it.

**18 deck READMEs written** (`examples/<deck>/README.md`: header block, what it shows, files,
run recipe, shared-launcher paragraph advising to copy the `Allrun`/`Allclean` wrappers from
`examples/helix/`). Facts checked against each deck; four of my first drafts were corrected on
the check (`bhdata.hdf5` ships in `share/material/`, not "written on first use"; the Python
mesh tools call the `gmsh` executable and need numpy/scipy/sympy/matplotlib, not the gmsh
module; `dipole` has no materials section; `plot_iv.py` reads `iv_results.csv`). Codex
language sweep (luna, medium): see the exchange for what was applied.

**References updated:** `examples/README.md` §1 and §5 (feature matrix replaces the starter
table), `examples/scripts/README.md`, `doc/getting_started.md`, `doc/input_schema.yaml`,
`doc/input_file_reference.md`, `src/circuit/doc/circuit_usage_guide.md`,
`src/fem/maxwell/fn_mesh_config_tag.hpp` (comments), `tests/circuit/test_NetlistParser.cpp`
(comment), `scripts/update_doc_index.py` (comment), `racetrack_usermat/{data/README.md,src/CMakeLists.txt}`,
`tapestack_circuit/tapestack.cir` (comment). Devlogs, `todo/` and `doc/lessons_learned_evidence.md`
keep the old names as written.

## Not done / owed

- Nothing built or run. Gate: `make check` untouched by this change; the real gate is
  `./Allrun --dry-run` in `pancake` and `corc_twolayer` (discovery of `python/main.py`) and one
  solver start per renamed deck to confirm the reformatted `input.conf` files load.
- Plugin `src/build/` directories, if any exist locally, embed the old absolute
  `BELFEM_USER_DATA_DIR` and must be rebuilt.
- `corc_twolayer/python/mesh/` and `pancake/python/mesh/` are byte-identical vendored copies;
  dedupe deferred (two `sys.path` lines, easy to break, not coupled to the rename).
- Stale comment `examples/corc_solder` in `examples/scripts/Allrun:244` and
  `src/io/cl_Input_Section.cpp:78` refers to a deck that was never shipped; left alone.
