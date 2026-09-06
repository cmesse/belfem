# DR-146: refusing the input a material shape never reads

**Date:** 2026-08-29
**Purpose:** record the materials-section shape validator, the two defects the executable
control caught that review did not, and what the row still owes
**Module:** `src/physics/materials`, `doc/` input contract

## What was silent

A `materials` subsection selects exactly one shape — a b-h curve, a builtin, or a plugin — and
each shape reads a different set of keys. Everything else was dropped without a word. An earlier
three-AI round inventoried thirteen cases; this session implemented the fix they specified.

The decisive design point, and the reason a "record what we read" sweep was rejected: `RRR` on a
`ybco` **is** read, and then discarded by a constructor that takes no RRR argument. *Was it read*
and *did it do anything* are different predicates, and only the second is the contract.

## The work I did not do

A parallel session had fifty-five uncommitted lines in `cl_MaterialFactory.cpp` — the same file —
renaming `custom` to `usermat` and adding `critical temperature` refusals, which is one of the
thirteen. Their devlog filed DR-146 under "found and NOT fixed" and handed it off. So the row was
mine, but every line number and allow-list had to be read against the working tree including their
edits, not against HEAD. Both auditors were told this explicitly, and Grok returned a table of
five citations in my own plan that were already stale.

## What the audits changed

Neither vendor merely approved. Four corrections changed the outcome:

- **The RRR gate.** I proposed, and Codex confirmed, a list of material names whose constructors
  drop RRR. Grok replaced it with a type test: RRR is legal exactly when
  `tMat->type() == MaterialType::PureMetal`, because `Alloy` forwards it and reports PureMetal
  while HastelloyC276, Magnesia and YBCO all drop it. A name list rots the day someone adds a
  material; a type gate does not. Adopted.
- **The HTS split.** My allow-list held `file` *and* `jc`/`n`/`ec`. Those are an exclusive choice
  in the factory, so a union list would have accepted a deck setting both and silently dropped the
  constants — re-silencing the exact case the row exists to close.
- **The usermat `defect` subsection.** I gated the usermat *keys* on `have(jc)` and forgot the
  subsection, which sits inside the same conditional.
- **Labelled subsections.** `Section` maps only unlabelled children by type, so `defect : mine { }`
  was never consumed by anything. Now refused.

Duplicate material labels became fatal in the same change, which closes both the last-wins
ambiguity and the leaked first `Material*`.

## Two defects the executable control caught, and review would not have

The gate was a probe driving the real `MaterialFactory` over real decks.

**First run: every must-fail case passed.** The probe linked the prebuilt `libbelfem.a`, which
predates the change — the validator was not in the binary at all. This is the same stale-artifact
trap recorded for DR-118 a few hours earlier, walked into from the other direction. Rebuilt with
the modified source compiled in.

**Second run: a legal deck was rejected.** `copper { builtin : copper ; RRR : 50 ; }` failed. The
parser lowercases key names at read time; `key_exists()` lowercases its argument so callers never
notice, but a direct comparison against an allow-list spelling `"RRR"` does not match. Every
copper deck in the tree — 2D_Tapestack, RLC_Circuit, circuit, corc, helix — would have stopped
parsing. CLAUDE.md warns about this exact asymmetry. Both sides are now normalized rather than the
single literal patched.

To be precise about ownership, because the first telling of this was too generous to me: this was
a defect in the allow-list, not in the input layer. Deck-side case-insensitivity already worked
and still does — a deck spelling `BuiltIn`, `JC`, `Resistivity Type` and even the subsection name
`DEFECT` in mixed case behaves identically before and after every change in this session, verified
by running the same deck against both. `get_type`, `get_label` and the key parse all fold to lower
case at read time.

The one genuine asymmetry sits elsewhere and is caller-side only: `section_exists` and both
`section()` overloads did not lowercase their ARGUMENT, while every key accessor does. No in-tree
caller passes a mixed-case section name, so it was latent. Fixed here for consistency
( `cl_Input_Section.hpp`, `cl_Input_Section.cpp` ), with no deck-visible change.

Had the syntax check and the first "ten decks pass" been trusted, this would have shipped as a
validator that did nothing, and then as one that broke half the examples.

## The gate as it stands

Ten of eleven shipped decks pass. `sidecoating` fails on a missing plugin `.so`, raised during
construction before the validator ever runs and reproduced identically on the stale-library
build — environmental, not this change.

Crafted cases all rejected, each naming the material and the shape: `file` beside constants, `RRR`
on ybco and on hastelloy, `resistivity type` on a pure metal, the `builtin iron` misspelling, a
plain `RRRR` typo, a labelled `defect : mine`, and a duplicate label. Curve + bare `builtin ;`,
copper + `RRR`, and HTS-with-constants pass.

## Deck repairs, per Christian's ruling

Fatal was chosen over warning, with both affected decks fixed here rather than left to abort for
whoever ran them next. `examples/2D_Undulator:51` lost its inert `builtin iron ;`;
`cmake-build-debug/tape_quench_usermat` lost the never-opened `file : sst-1.hdf5`, replaced by a
comment saying why — sourcing jc/n from a table while a plugin supplies the rest is not
expressible today, and deleting the line without recording that would erase the intent.

Codex caught that my first edit to the undulator deck had also converted the whole file from CRLF
to LF, turning one deletion into a 667-line diff. Redone at byte level; the diff is one line.

## What the second audit round caught in the LANDING, not the code

Grok returned REQUEST-CHANGES while approving the C++ against the matrix, and both objections
were mine:

- **`doc/input_schema.yaml` no longer parsed.** My `curve_exemption` entry was missing its
  closing quote, so the machine-readable half of the input contract was broken YAML — the one
  failure mode that artifact exists to prevent. Verified broken, fixed, and verified parsing.
- **The contract pair still asserted the behaviour this change made fatal.** Three leftovers:
  the schema preamble saying other shapes' keys "are then never read", the `bh_curve_ferromagnet`
  note saying "other keys in the section are then ignored", and the `RRR` note calling it
  "accepted and silently inert" on ybco/hastelloy/magnesia. All three now state the refusal, and
  the RRR note records the `type() == PureMetal` gate rather than a name list.

Neither was visible in the C++ diff. A doc edit that contradicts the code it documents is the
same defect class as the row itself, one layer up.

## Still owed

- `make check`.
- **The `expk` half of this row is not fixed.** Both vendors agreed it is a design question rather
  than a validator case — `set_sigmoid` stores `log(expk)` while the waveform derives its rate from
  `fuzzyness` alone, which smells like a superseded parameterization — and it touches three
  factories. Documented as inert in all three doc sites; the decision is Christian's.
- Residue Codex flagged outside the settled matrix: duplicate allowed *keys* and duplicate
  unlabelled subsections are still last-wins in the parser itself.
