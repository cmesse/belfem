# Release 0.9.0 final sweep before the history re-init

**Date:** 2026-09-04
**Purpose:** Three-voice read-only sweep answering one question — is there a roadblock to
publishing 0.9.0 tonight — before Christian archives the repository, deletes `.git`, and starts
a clean public history (the old history carries the combustion module, now in `./nonfree/`)
**AIs involved:** Claude (Fable, sweep + reconciliation), Codex `gpt-5.6-terra`/high, Grok
`grok-4.6`/high (blind jury, `tmp/ai_exchange/release_090_final_sweep.md`)
**Status:** reviewed, not verified — nothing built or run except the python and doc gates

## Method

Tracked set only (`git ls-files`, `git grep`), so the gitignored `nonfree/`, `archive/`,
`literature/`, `tmp/` trees were excluded by construction. Mechanical gates run:
`check_doc_claims.py` 38/38, `check_wrapper_policy.py` clean, `review_status.sh` no open P0,
`belfem-conf drift` (4 findings) and `roundtrip` (19/19). `h5dump` on every shipped table.
Literature checked on Christian's suggestion: `literature/books/coding` has nothing on
publishing or history rewriting — Oliveira & Stewart 2006 §17.2 is a 2006 survey of RCS, CVS and
Subversion; Rouson et al. 2011 does not cover it.

## The re-init fact that mattered most

After `rm -rf .git`, the publish filter is `.gitignore` alone. Two things were held only by
`.git/info/exclude`, which dies with `.git`: the stale 51 MB worktree
`.claude/worktrees/debt-register-sweep/` (its own `.git`; `.gitignore` re-includes `.claude/`
wholesale, so a fresh `git add .` would have picked it up), and the other Claude state paths.
Conversely, 30 tracked `more/gmsh/*.msh` element-type fixtures match the `*.msh` ignore rule and
would have vanished silently. Grok found the first independently and rated it the single most
important fact of the round; Codex missed it.

**Fixed (Christian approved `.gitignore` edits, nothing else):** the eleven `.claude/` state
paths moved into `.gitignore`; `!more/gmsh/*.msh` keeps the fixtures; `perf.data`, `*.rej`,
`*.orig` and the hook bookmark `.claude/ai_exchange_pos.txt` are ignored so the tracked junk
(a 5.9 MB `perf record` dump from `b81cc436`, `tests/linalg/test_Matrix.cpp.rej`, two `.pyc`) is
not re-added. Simulated with `git ls-files -i -c --exclude-per-directory=.gitignore` and the
`--others` twin: exactly those five files drop, nothing new is added, the 30 meshes stay.

## Open, cheap, not blocking

- `src/visualizer/visualize.cpp:51` defaults to `/home/christian/codes/belfem/nonfree/share/spacecraft.msh`
  (only built under `USE_VTK`); `scripts/scls_env.sh:17` and `.claude/hooks/ai_exchange_delta.sh:16-17`
  carry personal absolute paths.
- `belfem-conf drift`: `custom` parsed at `cl_MaterialFactory.cpp:265` but not in the schema;
  anchors `"eta"`, `"has_block_global"`, `"reorient_generators"` resolve nowhere.
- `CLAUDE.md:179` and `doc/coding_philosophy.md:63` say there is no `.gitlab-ci.yml`; the file has
  been tracked since `4b6a4423` (2026-08-31). `check_doc_claims.py` has no probe for that sentence
  (Grok, F6).
- `examples/README.md` starter table names `examples/corc`, which does not exist
  (`corc_christian`, `corc_gregory` do); "seventeen … each with `.geo`" is wrong for
  `corc_christian` and `pancake`, which ship python generators (Grok, F8).
- `CITATION.cff` `date-released` is 2026-08-14. 89 of 951 `src/` files carry no license header.
  `src/executables/hphirun.cpp` and `hphiTrun.cpp` still ship as sources though not built.
- `bhdata.hdf5` records no source or license in `share/material/README.md`.

## The one ruling

Codex rated `share/material/bscco-2223.hdf5` a blocker: digitised figure data from Turrioni et al.
2008 (AIP Conf. Proc. 986, 451), `meta/reference` present, no `meta/license`, and the `gantry`
example loads it. The embedded `source/` group is digitised CSV point sets, not figure images.
Christian ruled on 2026-08-31 (`todo/cancelled/sch04_history_purge.md`, scope guard) that this
table stays while the Fujikura brochure digitisations were removed. Rights questions are not
settled by vote; relayed, not decided.

## Not checked

The GitLab nightly (pipeline result lives on the server), `make doc`, a configure of a fresh clone
without `nonfree/` (read-only trace by all three voices says the defaults never touch it).

## Applied after Christian's approval (same evening)

- `src/visualizer/visualize.cpp`: no-argument run prints a usage line and exits instead of
  loading a `/home/christian/.../nonfree/` mesh. **Not compiled** — the VTK headers are not on
  this machine and the module is built only under `USE_VTK`; the unchanged includes fail the
  same way before the edit.
- `scripts/scls_env.sh`: `SCLS_PY_PREFIX` is environment-overridable, defaults to
  `$HOME/Applications/python`, and is prepended only when the directory exists.
- `.claude/hooks/ai_exchange_delta.sh`: paths derived from the script's own location.
- `CLAUDE.md` and `doc/coding_philosophy.md`: the "no `.gitlab-ci.yml`" sentence replaced by
  what the file is (scheduled or manual pipeline, runner `belfem-local`, not reproducible from a
  clone). `check_doc_claims.py` 38/38.
- `examples/README.md`: nineteen decks, seventeen with a `.geo`; `corc_christian` and `pancake`
  generate theirs; starter table row `corc` → `corc_gregory`.
- `doc/input_schema.yaml`: `nitsche ghost penalty` consumer is `fn_FEM_ghost_switch.hpp`;
  `has_block_global(` and `reorient_generators()` anchored as functions (the latter in
  `cl_CutFactory.cpp`); `custom` listed as a bare literal so the key inventory knows the refused
  legacy subsection. `belfem-conf drift`: clean, 146/146 anchors. The prose reference needed no
  change — `custom` was already documented as refused at `doc/input_file_reference.md:495`.
- `CITATION.cff`: `date-released` 2026-09-04.
- `examples/corc_christian/input.conf`: the last deck-written `bearing { nodes : 1; }` removed.
  A hand-written bearing replaces the automatic per-component gauge pins rather than adding to
  them, so every unnamed φ component floats; it is an expert function and no shipped deck should
  demonstrate it. `examples/README.md` §"Do not write a `bearing`" already said no deck here
  writes one — that claim is now true. This deck was the only occurrence under `examples/`;
  the remaining hits live in build trees and `tmp/`, which are not shipped.

Left as they were: 89 `src/` files without a license header, `hphirun.cpp`/`hphiTrun.cpp`
sources still in the tree, `bhdata.hdf5` provenance, `bscco-2223.hdf5` (Christian's ruling).
