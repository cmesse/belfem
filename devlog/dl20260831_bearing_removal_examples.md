# Bearing sections removed from the example decks and their generators

**Date:** 2026-08-31
**Purpose:** Record the pre-commit cleanup that took the deprecated `bearing` section out of
every shipped example, out of the two generators that re-emit it, and the doc change that
demotes the key to an expert function

## What changed

Automatic gauge pinning has been the default since 2026-08-28: with no `bearing` in the deck,
`MaxwellFactory` sets one pin per connected φ component and persists it in the `.bfm`. The
example decks still carried hand-written bearings from before that, and a deck bearing
**replaces** the automatic pins rather than adding to them, so every shipped example was
opting out of the feature it should be demonstrating.

Removed the `bearing` section from 11 `input.conf` files:

```
2D_Tapestack  2D_Undulator  RLC_Circuit  circuit  corc_periodic
corc_solder   costheta      gantry       garber   helix
Tape_Quench_obsolete_deleteme
```

Associated comment lines went with the block (gantry's `// node 4 …`, corc_solder's six-line
gauge-node note) along with one separating blank line. CRLF endings in `2D_Tapestack` and
`2D_Undulator` are preserved.

Two generators would have re-emitted the section on the next regeneration, which is why the
deck edit alone was not enough:

- `3D_tapestack/tapestack3d.geo` (and its byte-identical untracked copy `tapestack3d/`) —
  dropped the `Printf("    bearing")` block
- `corc_solder/python/corc/belfem.py` — dropped the emitted block, and reworded the module
  docstring and the comment the generator writes into the deck. Both had justified node 1's
  exclusion from the periodic triple *by* its being the bearing node. Node 1 still stays out,
  but the stated reason is now the real one: three points already fix the affine cap-to-cap
  map, and the cap center lies on the axis. Ran the generator; the deck it produces has the
  periodic block intact and no bearing.

## Why the empty `boundary conditions { }` stays

`circuit` and `RLC_Circuit` had nothing else in that section, so both now carry an empty
block. Deleting the section outright would abort them: `cl_MaxwellFactory.cpp:144` asserts
`section_exists( "boundary conditions" )`. The parser tolerates the empty form —
`Section::num_sections()` returns 0 and the BC loop is skipped.

## Documentation

`bearing` is not removed from the language, so the docs still describe it — with the
recommendation inverted. Three artifacts, kept in lockstep per the input-contract rule:

- `doc/input_file_reference.md` — the old "prefer the automatic pinning unless you have a
  reason" hedge replaced with an explicit *expert function and not recommended*, plus the
  failure mode: a hand-written bearing names one node and silently gives up the pins for
  every φ component the deck does not mention
- `doc/input_schema.yaml` — the same statement in `bearing_note`
- `examples/README.md` — a new `### Do not write a bearing` subsection at the end of §5

No code or key behaviour changed, so nothing was owed on the `src/` side.

## Evidence

`scripts/check_doc_claims.py` passes 37/37. `doc/input_schema.yaml` parses. The corc_solder
generator was executed and its output inspected. `grep -ri bearing examples/` returns only the
two intentional prose mentions in `README.md`.

**Not run:** no deck was executed and nothing was built, so the claim that these examples still
converge on automatic pins is *reviewed*, not verified. The `find_autopins` path was read
(`cl_MaxwellFactory.cpp:924-969`) and confirmed to engage only when no bearing BC is present;
no example ships a stale `.bfm`, so none hits the pre-autopin "potential is unpinned" warning.

**Owed:** the Codex language sweep over the two user-facing docs
(`doc/input_file_reference.md`, `examples/README.md`) has not been run.

## Loose end unrelated to this change

`examples/tapestack3d/` and `examples/tapestack2d/` are untracked directories. The `.geo` edit
in the first was applied on disk so the two copies stay identical, but it will not appear in
the commit — only `3D_tapestack/tapestack3d.geo` is tracked. Whether those directories are
meant to be added or are scratch is undecided.
