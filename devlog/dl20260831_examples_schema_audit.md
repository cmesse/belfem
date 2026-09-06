# Examples audited against the input schema; three checker false positives fixed

**Date:** 2026-08-31
**Purpose:** Record the conformance sweep of `examples/*/input.conf` against
`doc/input_schema.yaml`, the deck defects it found, and the `belfem-conf check`
defects it found — including one false negative that mattered more than any of
the false positives.
**Module:** `python/belfem_conf`, `examples/`

## What was done

All 15 shipped decks were run through `python/belfem-conf check`, and every
reported error was then traced to the C++ consumer before being believed. That
second step is the whole value of the session: of the seven errors the checker
reported, five were the checker's fault.

## Deck defects

**`examples/corc_solder/input.conf` — the `mesh` section header was missing.
FIXED: the block is now byte-identical to `examples/corc_periodic`'s, which is
the deck it shares `corc.msh` with.** With the header restored the deck
validates clean — 58 keys, 17 enum values and 12 unit dimensions, none of which
had ever been reachable, because the parse died before any of it.

Note, not changed: `examples/corc_solder/python/main.py:55` writes the mesh to
`/tmp/corc.msh`, while the deck reads `corc.msh` from the run directory. That
gap is the generator's, not the deck's — `corc_periodic` names the file the
same way.

The original diagnosis, kept because the hazard is general:
The trademark comment block was inserted where `mesh` used to sit, so the deck
opens with a bare `{`. Compare `git show HEAD:examples/corc/input.conf`, the
deck this one was derived from.

This does not fail cleanly. Every line before the `{` is a comment, so after
`InputFile::remove_comments` + `tidy_up` the parse buffer *starts* with `{`
(replicated in Python to be sure). `Section::create_children` then evaluates
`tLabel = mBuffer( k - 1 )` at `k == mStartFlag == 0`
(`src/io/cl_Input_Section.cpp:65`), and `index_t` is `uint64_t`
(`src/core/typedefs.hpp:49`) — unsigned underflow into an out-of-bounds
`Cell<string>` read. Debug asserts; release reads garbage as a `std::string`.

Worth noting as a general hazard, not only as a deck typo: **any** deck whose
first non-comment line is `{` hits this, and it is now the only known way to
reach that read. A guard in `create_children` would be cheap and would turn a
garbage read into a named error; still open.

**`examples/sidecoating/` — MOOT, the deck was deleted during the session.**
The audit had flagged it for loading a `lib/build/libcustom.so` it did not
ship, and for selecting `library : pardiso`, which needs `USE_PARDISO=ON`.
Both findings die with the deck; recorded only so a later reader does not go
looking for a deck that once had them.

**`examples/tape_quench_usermat/input.conf:45` — `target iterations` inside
`nonlinear thermal` is dead.** `Controller::set_params` reads the key only from
the magnetic block (`cl_FEM_Controller.cpp:4231`); the thermal block from
`:4354` reads twelve keys and this is not one. Silently ignored.

Minor, recorded and left alone: a stray `/` above the `homology` header in
sidecoating (survives `remove_comments`, harmless because the header is the
line *immediately* before `{`); `phase : 0` without a unit in 2D_Tapestack
(legal — see below — but the other three phase-setting decks write `0 deg`);
inert `label` keys in RLC_Circuit's topology domains; and `library : pardiso`
in sidecoating, which needs a build with `USE_PARDISO=ON`.

## Checker defects (fixed)

**The false negative, which was the real problem.** `_spec_for` matched deck
sections to schema entries by exact name only. The schema reaches `linear`
through `resolution_chains` and `nonlinear` through `aliases`, so
`linear magnetic`, `linear thermal` and `nonlinear magnetic` matched **no entry
at all** — and a section with no entry is skipped in silence. Every key inside
the magnetic solver configuration of three decks was unchecked. Demonstrated by
planting `another bogus : 7 ;` in `linear magnetic` and watching the deck pass.

Fixed with `_section_aliases` / `_lookup_section` / `_spellings_in_precedence`.
An exact entry beats an alias, so `nonlinear thermal` — which has its own entry
and differs from the magnetic one in three keys — is still validated against
the right list; this is checked in both directions (a thermal-only key in
`nonlinear magnetic` and a magnetic-only key in `nonlinear thermal` are both
reported). `per_field_override` merges `compute conditioning` and
`mumps error analysis` into the field-specific spellings only, because
`cl_FEM_Controller.cpp:4671-4689` reads them from `linear magnetic` /
`linear thermal` and from the `solver` level, never from a bare `linear`.
Coverage on 3D_tapestack went 47 -> 61 keys and 12 -> 17 enum values.

**`builtin : Pb38Sn62` rejected as illegal.** It is a legal solder.
`create_material` falls through to `to_pair` (`core/stringtools.hpp:312`),
whose regex `([A-Za-z]+)(\d+)` matches `Pb`/38 and `Sn`/62 in any order, and
`Alloy` reports `MaterialType::PureMetal`, so the sibling `RRR : 1.5` is
allowed by the DR-146 allow-list too. The schema already said so, but only in
prose (`also_accepts: "alloy formula, e.g. Sn40Pb60"`), which no program can
act on. `also_accepts` is now structured — token pattern, component list,
balance rule — and the checker applies it. `Ybco50Cu50` and `Xx50Cu50` are
refused, `NotAnAlloy` is refused, and `Sn40Pb50` warns rather than errors
because `set_components` overwrites the last fraction with the balance.

**`phase : 0` reported as a runtime abort.** It is not.
`Section::get_value( "phase", "rad" )` ends in `check_unit`, which compares
nothing but the dimension code, and `rad` contributes no exponents
(`stringtools.cpp:955`) — so a bare number matches. The checker was testing the
dimension's *name* against `("dimensionless", "-")` instead of resolving its
signature. Now it resolves the signature; an unresolvable spelling still errors,
and `initial timestep : 5` is still rejected.

## `range` enforcement (added on request, same session)

`range: [0, 8]` with `out_of_range: fatal` was declared on `anderson depth` and
enforced by nobody. It is now enforced, with the severity taken from
`out_of_range:` rather than assumed — a range recorded for documentation must
not silently become a hard error. `constraint:` stays prose and stays unparsed,
deliberately.

Two things fell out of it:

**A schema inconsistency.** The magnetic `anderson depth` carried
`out_of_range: fatal`; the thermal one did not, although both are the same
`BELFEM_ERROR( tDepth >= 0 && tDepth <= 8, … )`
(`cl_FEM_Controller.cpp:4273` and `:4445`). One C++ check was described two
ways. The thermal entry now carries the field too.

**A fourth false positive, in the code path the range check sits in.** The
`int` rule was `re.fullmatch(r"[+-]?\d+")`, so `anderson depth : 3.7` was an
error — but `Section::get_int` is `round( get_real( … ) )`
(`cl_Input_Section.cpp:349`), which accepts it and stores 4. Same mistake as
rejecting `Pb38Sn62`: stricter than the code, on a deck that runs. Now a
warning naming the rounded value, and the range test uses that value. The
rounding is C's `round()` — half away from zero — not Python's `round()`, which
is banker's; `_round_half_away` was checked against `libm` on 8.4, 8.6, -0.4,
-0.6, 0.5, -0.5 and agrees on all six, including the two half-cases where the
Python built-in would have differed.

## Launcher: the concept is current, its wiring is not

Asked whether `Allrun` / `run.conf` went stale when the unified `belfem`
executable landed. They did not — that part was already done. `Allrun`'s
`pick_executable()` sets `EXECUTABLE=belfem` and explains that the deck picks
the physics, `examples/scripts/README.md` (2026-08-18) documents the same, and
`run.conf`'s keys are about gmsh and MPI, orthogonal to which binary runs. Its
one launcher-relevant key, `EXECUTABLE`, is correctly documented as
"empty = decide from `input.conf`".

Three things around it had NOT kept up:

**`examples/README.md` §2 had.** Dated 2026-08-10, it listed `belfem`,
`hphirun` and `hphiTrun` as three co-equal ways to run, never mentioned that
the latter two are retiring, and never mentioned the launcher or `run.conf` at
all. Rewritten: `belfem` first and alone in the code block, the other two named
as retiring halves of what it now decides, `electricalCircuit` explicitly not
an alternative (it reads no input file), and a pointer to
`examples/scripts/`. Its `Allclean` list was also a subset of what the script
removes; now it matches, including the `--mesh` rule.

**`examples/scripts/Allclean:31` never removed the circuit output it exists to
remove.** It named `circuitAnalysis.out`; `electricalCircuit.cpp:42` opens
`circuitAnalysis.txt`, and every other mention in the tree says `.txt`. Found
by the Codex sweep on the README rewrite, verified against the source, fixed.

**`NDIMS` deleted — a knob that could only do harm.** It fed one thing,
`gmsh -$NDIMS`. Measured both directions:

- `-2` versus `-3` on a geometry with no volume is **byte-identical**
  (`2D_tapestack.geo` 6349454 B either way, `costheta.geo` 2016555 B), because
  gmsh meshes dimension by dimension up to the number given and the 3D pass
  then finds nothing: `Done meshing 3D (Wall 3.1202e-05s)`. So `NDIMS=2`, the
  only value any deck ever set, bought nothing.
- `-2` on a 3-D geometry exits 0 and writes a `.msh` with every tetrahedron
  missing: `helix.geo` gives `{0: 130, 1: 1304, 2: 26084}` against
  `{..., 3: 217223}`. `ensure_mesh`'s only guard is "did gmsh produce the
  file?", which that passes — so the wrong value fails silently and reaches the
  solver as a surface-only mesh.

A knob whose intended setting is a no-op and whose wrong setting is silent is
worth deleting. `Allrun` now hardcodes `-3`; the key is gone from the defaults,
the `run.conf` whitelist and the README table. All three shipped `run.conf`
files contained nothing but `NDIMS=2` and were removed with it, so no deck
currently carries one — the mechanism stays, with no users. A leftover local
`run.conf` degrades cleanly: `run.conf: ignoring unknown key 'NDIMS'`, and the
rest of the file still applies (verified with `NTHREADS=4` alongside it).

**`ensure_mesh` no longer assumes gmsh.** `corc_solder` builds its cable in
`python/main.py` and has no `.geo`, so "Neither the mesh nor its geometry is
present" sent the reader hunting for a file that was never meant to exist. It
now looks for a `main.py` in the deck or one level below and names the script
to run; with neither a `.geo` nor a generator it says that instead. Both
branches exit 1 — there is still no mesh — and the gmsh path for decks that do
ship a `.geo` is untouched.

**The per-deck wrappers are 3 of 14.** `examples/scripts/README.md` says "Each
deck holds a two-line wrapper that `cd`s to itself and execs the shared
script", and eleven of them are gone (`git status` shows them deleted); only
`costheta`, `garber` and `helix` still have one, and the three newest decks
never did. So `./Allrun` from a deck directory — the documented usage — works
on three decks. **Not resolved here**, because it is a fork rather than a bug:
either restore the eleven, or drop the wrapper concept and document invoking
`../scripts/Allrun` directly. Left for a ruling.

## Doc indexer: a nested CLAUDE.md aborted the whole nav build

`scripts/update_doc_index.py` collected every `.md` under `examples/` and
requires each to have a `PAGE_GROUPS` entry, so
`examples/corc_solder/python/CLAUDE.md` — the mesh generator's own session
bootstrap, not a documentation page — failed the run outright. The root
`CLAUDE.md` and `AGENTS.md` escape only because they sit outside the scanned
trees, which makes indexing a nested one an accident of location rather than a
decision.

Fixed with a `NOT_DOCUMENTATION` set applied at all three collection sites
(`src/*/doc`, `doc/`, `examples/`). Checked that the exclusion is narrow: a
planted nested `CLAUDE.md` is ignored, while a planted real page in the same
directory still trips the `PAGE_GROUPS` requirement, which is the check's whole
point.

## A concurrent session broke the schema mid-edit

While this work was in progress another session appended to the
`resistivity type` note in `doc/input_schema.yaml` (riva n-source validation,
`set_resistivity_law`). Four of its lines landed with no indentation, which
breaks the block scalar and makes the whole file unparseable — `check`, `drift`
and `roundtrip` all fail on a YAML error, not on anything they are testing.

Repaired here by re-indenting those four lines to the surrounding ten spaces.
**Their wording was not touched.** Flagged rather than merely fixed, because a
read-modify-write on a shared file cannot tell a mangled edit from one it
clobbered itself: if that session had a newer version of the paragraph, this
tree does not have it.

## The header guard, and the audit round that changed it

`Section::create_children` reads the section header as `mBuffer( k - 1 )`. Guard
added, then audited by Codex (gpt-5.6-sol / high) and Grok (grok-4.6 / high)
independently.

**First dispatch produced no audit at all** and is recorded so it is not
mistaken for one: Codex returned "Selected model is at capacity", and Grok was
given `grok-4.2`, which it rejects (`Allowed: grok-4.6 grok-4.5`) -- an error in
the dispatch, not a finding. Redispatched at the depths above.

**The two auditors split on the predicate**, which is the whole value of running
both. Both agreed on everything else: `mKey` is initialised (declaration order
matches the init list, `create_children()` runs in the constructor body),
`BELFEM_ERROR` is the required tier because `Cell::operator()` checks bounds
only under `BELFEM_ASSERT` so an assert would leave release builds unsafe, and
no correctly-parsing deck is rejected.

- Grok: keep `k > 0`. Anything wider is a behaviour change on release day, and
  the message would lie about a nested failure.
- Codex: use `k > mStartFlag`. The invariant is that the header must lie inside
  this section's own scan range; a nested section opening straight onto a brace
  reads the PARENT's brace and names the child `"{"`.

**Codex's reading was adopted**, because its safety claim was checked rather
than taken: a scan of every shipped deck for a brace whose preceding line is
absent, `{` or `}` -- and for stray or unbalanced braces -- returns **zero
hits**, so the stricter predicate cannot reject anything we ship. Grok's
objection to the message was real and was fixed rather than dismissed: the text
now reads "In '%s': an opening brace has no section header on the line above
it", which is true of a file and of a section, and the nested case correctly
reports `In 'solver'` rather than naming the file.

Grok's comment correction was also applied -- the release-build access is
undefined behaviour, not a guaranteed "garbage `std::string`".

Evidence, not inference. Two probes were linked against the debug tree's
`libbelfem.a`, one carrying the `HEAD` version of the parser and one the fixed
version:

| case | before | after |
|---|---|---|
| headerless root deck | `Cell index out of bounds: 4294967295 (expect < 126)`, `cl_Cell.hpp:168` | named error, exit 1 |
| nested `solver { { } }` | parsed clean, built a section named `"{"` | named error, naming `solver` |
| all 14 shipped decks | parse | parse |

`index_t` is `uint32_t` in that build, which is where 4294967295 comes from.
Compiles clean under the module's own flags, `-Wall -Werror -pedantic-errors`
included. `git status` after the Grok round shows no source touched but this
file.

### Left in place, deliberately, by both auditors and me

`--tSectionCount` at `cl_Input_Section.cpp:85` is unconditional and
`tSectionCount` is unsigned, so a stray `}` at depth zero wraps to the maximum
instead of erroring. **Confirmed empirically, and the consequence is silent data
loss:** a deck with one extra `}` after its `mesh` block parses with no error at
all, and `section_exists( "solver" )` then returns false -- every later section
is swallowed. There is likewise no `tSectionCount == 0` check after the loop, so
an unclosed brace is accepted too.

Not fixed today. Unlike the header guard, a `BELFEM_ERROR( tSectionCount > 0 )`
would newly reject decks that currently parse -- badly, but they parse -- and
that is a behaviour change to make deliberately, not on release day.

## Gates

`roundtrip` 15/15 byte-identical; `drift` unchanged at its two pre-existing
findings; `check_doc_claims.py` 37/37. `python/README.md` had a Codex language
sweep (gpt-5.6-luna / medium), which also caught a stale "three decks" in the
schema comment and an overstatement of what the checker independently verifies
— both corrected.

Static work only. **Reviewed, not verified:** no deck was executed, and the
corc_solder out-of-bounds read is reasoned from the source and a Python
replication of the two cleaning passes, not observed under a debugger.

## Left open

- The dead `target iterations` in `tape_quench_usermat`'s `nonlinear thermal`.
- The per-deck `Allrun` wrappers: restore the eleven, or drop the concept.
- The unguarded `mBuffer( k - 1 )` in `Section::create_children`.
- `range: [0, 8]` with `out_of_range: fatal` (`anderson depth`) is declared in
  the schema and not enforced by the checker. `constraint:` is prose and stays
  unenforceable; `range` is mechanical and could be.
- The `custom` material subsection — the pre-0.9.0 spelling of `usermat`,
  refused with a rename message at `cl_MaterialFactory.cpp:252` — is absent
  from the schema's refused-input vocabulary, which is one of `drift`'s two
  standing findings. No deck uses it.
