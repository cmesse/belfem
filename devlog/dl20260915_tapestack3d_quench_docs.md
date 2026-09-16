# tapestack3d_quench: deck comments and README, and the CMake 4.0 floor

**Date:** 2026-09-15
**Purpose:** Document the new `examples/tapestack3d_quench` deck to the standard the other
eighteen decks were brought to on 2026-09-05, and record what documenting it exposed —
including one defect report of my own that did not survive contact with the right toolchain.

---

## What was done

The deck arrived as `input.conf` + `tapestack3d.geo` + `src/{ramps.hpp, defect.cpp,
current.cpp, CMakeLists.txt}` with no README and an uncommented, space-indented deck file.

1. **`input.conf` commented and re-indented.** Tabs and column-aligned colons, matching the
   other eighteen decks; a header block naming what the deck is and how it differs from
   `tapestack3d`; section comments on the choices that are not self-evident (why the
   conditioning estimates are off, why the nonlinear tolerances are `1e-9` where the parent
   deck uses `1e-7`, what the solder is doing in a current-sharing study, why `edge coating` is
   off, what one bracket in the terminal list means, that a `userdefined` source takes no
   `amplitude`). Gated: comments and whitespace stripped from before and after, every key and
   value is byte-identical.
2. **`README.md` written**, following the sibling deck pages: what it shows, a file table, the
   plugin-then-mesh-then-run recipe, a "changing the scenario" section pointing at
   `src/ramps.hpp` and the constant block in `src/defect.cpp`, and the shared-launcher
   paragraph the other decks carry verbatim.
3. **Registered in the indexes.** Feature-matrix row in `examples/README.md` (`defect, source`
   in the plugins column), the deck counts corrected there (eighteen → nineteen, sixteen →
   seventeen shipping a `.geo`), added to the "decks that load a plugin" list and to the
   "where to start" paragraph; `doc/doxygen_nav.dox` regenerated with
   `scripts/update_doc_index.py` (27 modules, 117 pages).
4. **Codex language sweep** (`gpt-5.6-luna`, medium) over the README, the touched parts of
   `examples/README.md` and the deck comments. Applied its readability edits; its four flagged
   technical points were all worth acting on and are described below.

## Two defects found while documenting, one false alarm, and a ruling

**False alarm, recorded because the mistake is reusable: the plugin's
`cmake_minimum_required(VERSION 4.0)` is fine.** It was reported here as a defect on the
evidence of a configure that aborted — but that configure ran against `/usr/bin/cmake` 3.31.8,
because the session never entered an SCLS environment. `/opt/scls/gcc/bin/cmake` is **4.4.2**,
and the requirement is satisfied there. The change to 3.13 was reverted and the build re-run
under the SCLS toolchain: `userdefect.so` configures, links, and exports both `MyDefect_init`
and `MyCurrent_init` (`nm -D`). The session made no source change at all; it is documentation
throughout.

The generalizable form: **a tool version read from `$PATH` is not this project's tool
version.** BELFEM is built against a prefix that ships its own toolchain, so a configure or
compile failure means nothing until the environment it ran in is stated. `scripts/scls_env.sh`
has said so all along — it fails if `cmake` does not resolve under `/opt/scls` — and the
session ran the configure without sourcing it.

**Ruling, same session: the floor goes up instead of the deck coming down.** The asymmetry the
false alarm surfaced was real — top-level `CMakeLists.txt` at 3.11, both plugin templates and
the four sibling decks at 3.13, this deck alone at 4.0 — and Christian's call was that BELFEM
wants CMake 4 too. Raised in all seven places, so nothing in the tree claims a lower floor than
the toolchain it is developed against. Gated: a full out-of-tree configure with
`/opt/scls/gcc/bin/cmake` 4.4.2 and `SCLS=/opt/scls/mkl` exits 0 with **no** CMake warning and
no policy or deprecation message, which is the claim that mattered — raising the minimum from
3.11 to 4.0 flips every intervening policy to NEW, and nothing in this tree depends on an OLD
one. The two `cmake_policy( SET … NEW )` guards for CMP0156 and CMP0179 are now implied by the
minimum and were left in place as documentation of intent. `doc/getting_started.md` and the
deck README updated; CHANGELOG records that a distribution CMake older than 4.0 can no longer
configure BELFEM.

**The geometry re-emitted a `bearing`.** `tapestack3d.geo` was copied from a state of its
parent that predates the 2026-08-31 bearing removal, so its `Printf` block still prints a
`bearing { nodes : 9 ; }` into the deck skeleton it suggests to the user — the exact thing
`examples/README.md` now tells readers not to write, since a deck bearing replaces the
automatic gauge pins rather than adding to them. The five `Printf` lines were removed, which
makes the file byte-identical to `examples/tapestack3d/tapestack3d.geo` again. The deck itself
never had a `bearing` section, so nothing about the run changes.

**The CMake header described a different deck.** It named the `tapestack3d` example, announced
"the two libraries" where one is built, and documented `USER_DATA_DIR` as where "the plugins
read their measured tables from" — these plugins are analytic and read no tables, and the
`../data` directory it points at does not exist in this deck. Header corrected; the unused
variable kept, now labelled as inherited scaffolding for derived decks.

Codex also caught a comment defect of the same family: the deck header said the ramp is
"interrupted at 3 s by a defect". It is not — `ramps::current()` is monotone for the whole
10 s, and that the ramp keeps going is precisely why the stack quenches. Reworded.

## What the deck does, for the record

Eight soldered 4 mm REBCO tapes, the `tapestack3d` geometry and periodicity, coupled h-φ/T from
77 K. `src/ramps.hpp` holds the whole time program and both plugins read it, so the current and
the defect cannot drift apart: a 200 A/s linear ramp (cap 2500 A, not reached in 10 s), and a
logistic defect switch centered at 3.25 s that takes the topmost tape to 10 % of its jc over a
1.5 mm erf-smoothed disk at mid-length. The intended sequence is defect → sharing through the
solder into seven healthy tapes → quench later in the ramp, at the defect, where the diverted
current crosses resistive solder. Both plugins compile into one object, `userdefect.so`, which
the deck names twice — once under `materials { ybco { defect } }` and once as the current
boundary condition's `file`.

## Status

Documentation, plus the CMake-4.0 floor raised in seven build files on Christian's ruling.
`scripts/check_doc_claims.py` 38/38; `scripts/update_doc_index.py --check` 0 out of date. The
plugin build is **verified** (executed under `/opt/scls/gcc/bin/cmake` 4.4.2, symbols
inspected). The deck itself was **not run** — the physics narrative above is read from the
plugin sources and the deck, not a solve, so the claim that the stack quenches at all, and the
margin arithmetic in `ramps.hpp`, are the author's and remain ungated here.
