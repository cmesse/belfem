# BELFEM Python tooling

**Date:** 2026-08-11 (revised 2026-08-13 after the three-AI cross-review)
**Purpose:** Tools for the `input.conf` contract — the drift check, the deck
validator, and the parser round-trip gate; a configurator GUI is planned on
the same schema.
**Module:** `python/belfem_conf`

```bash
python/belfem-conf check <deck>...   # validate a deck against the schema
python/belfem-conf check -v          # ... and list what was NOT checked, and why
python/belfem-conf drift             # check the schema against the C++ parse sites
python/belfem-conf roundtrip         # parse + re-serialise every example deck, byte for byte
```

No installation, no build tree, no `PYTHONPATH`. PyYAML is the only dependency.
That is deliberate: the check has to run in a git hook and on a machine that has
never compiled BELFEM.

---

## Why this exists

The `input.conf` contract is described in two places that must agree:

| Artifact | Audience |
|---|---|
| `doc/input_file_reference.md` | people |
| `doc/input_schema.yaml` | programs |

and both must agree with the C++ that actually parses the deck. The rule is in
`CLAUDE.md` § *The Input Contract: Two Artifacts, One Rule*.

Human discipline keeps the **prose** current — measurably, it has. What it does
not catch is a key added to a factory and to no document, because whoever adds
the key is precisely the person who does not know the document exists. That is
what `drift` is for.

## What `drift` checks

**schema → code.** Every `anchor:` in the schema must still be findable in the
consumer it names. Catches renames and removals. Anchors are searchable tokens,
never line numbers — line numbers in this repository have been observed to be
wrong *in the commit that wrote them*, because a line read from a working copy
is stale before it is committed.

**code → schema.** Every string literal passed to `key_exists( … )`,
`get_real( … )`, `section( … )` and friends, across the input-consuming sources,
must be a key or section the schema knows. This is the direction that catches
the dangerous case. The scan runs over each file as one string, not line by
line — `mesh.unit` is parsed only by a call split across two lines, and a
per-line regex never saw it.

Comparison is whitespace-blind, because real call sites vary between
`key_exists("nodes")` and `key_exists( "nodes" )`; and case-insensitive, because
`Section::get_string` and `key_exists` both call `string_to_lower`, so
`get_reals( "Direction" )` reads the key stored as `direction`.

**enum values.** Several deck values are closed enums that the C++ resolves by
comparing against `to_string( value )` for every enumerator — which makes that
switch the authority and any list in the schema a copy. Rather than transcribe
them, `drift` re-derives each list from `src/` and fails on a difference. A
schema key opts in with `enum_ref: KrylovMethod` next to its `values:`.

This caught a real bug on its first run: `values: [off, blr, automatic]` loads
as `[False, blr, automatic]`, because YAML 1.1 reads a bare `off`/`on`/`yes`/`no`
as a boolean. The schema had silently lost a legal value. Quote those.

**Not checked: doc → code.** A claim can be false while every key and anchor is
correct — the reference once described a shipped feature as unimplemented, and
no grep finds that. It needs a reader.

## Reading the output

The first two lines are coverage, printed before any finding:

```
anchors  124/124 resolve   (declared in file: 124)
keys     111 parsed literals across 18 consumers; schema knows 144
```

`checked` versus `declared in file` is a self-report, and it is there for a
reason. An earlier prototype matched `anchor:` with a line regex, silently
skipped every anchor written inside an inline flow mapping, and reported
"all 53 resolve" for a file that held 115 — passing a check it had never run.
If those two numbers ever disagree, the walker is skipping anchors and the
result means nothing.

A checker that does not say how much it examined is indistinguishable from one
that passes everything.

**Exit codes:** `0` clean, `1` findings, `2` bad usage — including not finding
a BELFEM checkout (the working directory is searched first, then the tool's
own location, so checkout A's copy run from inside checkout B validates B).
Anchors that resolve under a different consumer are reported but do not fail
the run: some rows name the consumer that *uses* a section rather than the one
that looks it up.

## `check` — deck validation

What it checks: required sections **at every nesting level** and required
keys; unknown sections, unknown keys (including on boundary-condition and
topology subsections, whose names are data but whose keys are not); enum
values (including the six re-derived from the C++) and values the parser
accepts but the runtime refuses (`scheme : crank-nicolson`); numeric and
boolean types; **unit dimensions**, with the value shape mirrored from
`create_key`/`get_value` — a bare number on a dimensioned key and a three-word
value are both runtime aborts, and are reported as such; `mesh.unit` as a bare
unit token; topology domain types, including rejecting `cut` and demanding
`material` where the `Domain` constructor does; the `block` xor `blocks` rule;
boundary-condition **type names**, including the refused `background
dirichlet`; the circuit contract (required `topology`, component types, bare
headers, `node +`/`node -`, per-type required keys); `thinshell` ↔ `layers`
pairing; layer materials resolving in `materials` **and layer thicknesses
carrying a length unit**; `periodic` triples being exactly three (ranges like
`1:3` counted expanded, as `get_ids` reads them); duplicate keys and duplicate
sections, which the C++ maps collapse in silence; statements that lost their
`;` and are silently dropped; and the documented landmines — notably a
fractional `background` `direction`, which `std::stoi` truncates in silence.

Unit dimensions are compared the way BELFEM compares them — by **dimension, not
by token** — so `kA` passes where `A` is wanted and `mum` where `m` is. The
table is derived from `unit_to_si` at run time, not transcribed: ~800 lines of
`else if ( tUnit == "V" ) { tMass += tPower; ... }` is a machine-readable table,
and copying it would leave one more thing to rot. The preprocessing is mirrored
too — `µ` → `mu`, `²` → `^2`, `³` → `^3` — and so is the first-slash-only split,
so `µm` and `mm²` pass exactly where the C++ accepts them and `W/m/K` fails
exactly where it aborts. If `unit_to_si` cannot be parsed at all, dimension
checks are skipped and the skip is reported, never silently passed.

**Section names are resolved, not ignored.** Several solver sections have more
than one legal spelling: `nonlinear magnetic` is an alias of `nonlinear`, while
`linear magnetic` and `linear thermal` are the field-specific spellings in the
fallback chains for `linear`. The schema records these two mechanisms as
`aliases:` and `resolution_chains:`, respectively, and `check` now reads them.
It used to read neither. As a result, it called `solver/nonlinear` *missing* on
three shipped decks that spell it `nonlinear magnetic`. Worse, it matched
`linear magnetic` to no schema entry at all, so it checked **nothing inside it**:
a bogus key, an illegal solver name, and a wrong unit all passed in silence.

An exact entry still wins over an alias. Thus, `nonlinear thermal`, which has
its own entry and differs from the magnetic one in three keys, is never
validated against the wrong list. Where the C++ reads a parent key per field —
`compute conditioning` and `mumps error analysis` — the schema's
`per_field_override` merges it only into the field-specific spellings. A bare
`linear` is not a site from which the controller reads those keys.

**An enum is not always the whole contract.** `builtin` falls through to a
formula parser when the name matches no member, so a closed list alone rejected
`Pb38Sn62` — a legal solder — on two shipped decks. The schema's `also_accepts`
now carries the rule in structured form, and `check` applies it. At least one
`<symbol><percentage>` token must match. Each component must be one the schema
lists as a pure metal — the subset whose constructor reports
`MaterialType::PureMetal`, which is what `Alloy::create_component` demands — so
`Ybco50Cu50` is refused. Fractions that do not total 100% produce a warning
rather than an error, because `Alloy::set_components` overwrites the last one
with the balance. That component list is transcribed, not re-derived: a
material's `MaterialType` is fixed by its constructor, which no reader here
parses, so `drift` can only check that each name is also a legal `builtin`.

A prose-only `also_accepts` rule cannot be applied, so `check` passes rather than
fails in that case. When the enum is known to be an incomplete statement of the
contract, reporting a legal deck as broken is the worse error.

**Declared ranges are enforced; prose constraints are not.** A key may carry
`range: [0, 8]`, and `check` now rejects a value outside it. The severity comes
from the schema's own `out_of_range:` rather than being assumed, so a range
recorded for documentation does not silently become a hard error. The sibling
`constraint:` field stays prose — `"> 0, and >= min iterations"` — and is
deliberately left unparsed. Guessing at English would produce exactly the
confident-and-wrong findings this checker exists to avoid.

Integer keys are compared as `Section::get_int` reads them:
`round( get_real( … ) )` — half away from zero, not Python's banker's rounding.
Being fractional is therefore a *warning* that names what the value rounds to,
never an error in itself. The range test then sees the rounded number, and can
still fail on it: `anderson depth : 8.4` draws only the warning, while `8.6`
draws the warning and an error for the 9 it becomes.

**Against the mesh**, when there is one: `check` resolves `mesh { file }`
relative to the deck and verifies that every block and sideset id the deck
names actually occurs in it. Ids come from the **elementary** tag, matching
`GmshReader::create_group_ids` — `$PhysicalNames` is read and discarded by
BELFEM, so checking against it would agree by luck and disagree by design.
The reader mirrors `GmshReader` on purpose: the 2-D/3-D decision is the exact
`min(z) == max(z)` test, only formats 2.2 and 4.1 are parsed (the solver
hard-errors on anything else), binary files are declared unread, and a
truncated mesh reports as unreadable instead of crashing the run.

A missing mesh is **normal**, not an error: the examples ship `.geo` and
generate the `.msh`, so none of the seventeen has a mesh checked in. The validator
says which checks became unavailable and why, rather than passing in silence.
It also notes when a sibling `.bfm` exists, since the factory prefers that when
its stamps match — so the ids that actually run may not be the ones checked.

It also flags **inert keys** — input the parser accepts and no consumer reads.
A key can be tied to a boundary-condition type (`applies_to_bc_type`) or to a
source-function type (`applies_to_function`), and a key on the wrong side of
either is silently ignored at runtime. So a `bearing` carrying `length`,
`direction` and `fuzzyness`, or a `ramp` carrying `phase`, is now reported
rather than accepted.

**It reports what it did not check.** Conditional requiredness (T3) is
enforced; what remains needs data the deck does not contain — mesh-state
gates, bracket-group arity, and the material-shape prose (T4) — and `-v` names
each skip, including a deck-specific one like a missing `materials` section
disabling the layer-material check. A validator that silently skips a class of
check is worse than one that says it skipped it, because the reader cannot
tell "passed" from "never looked".

Severity: `error` fails the run, `warning` and `note` do not, and a deck whose
only diagnostics are notes shows `ok`. A note is for input the C++ accepts and
ignores — a stray `/` above a `homology` header, for instance, has never done
anything.

## `roundtrip` — the parser gate

The parser keeps the **original text as the source of truth**. Every node holds
the byte span it came from, nothing is normalised on the way in, and serialising
an untouched deck returns the original bytes. Edits replace one span, so hand
comments, tabs and layout survive — the example decks are hand-written and
git-tracked, and a tool that reflows them turns a one-key change into a
whole-file diff.

That design makes the naive round-trip check nearly tautological: `dumps()`
returns `self.text` when nothing was edited, so it would pass even if the parser
built no tree at all. So `roundtrip` also **probes every value span**: for each
statement it rewrites that value, requires exactly those bytes to change and no
others, and re-parses to read the new value back. 1367 spans across the
example decks.

Run it after any parser change. If a deck cannot be reproduced byte for byte,
nothing built on the parser can be trusted not to mangle a file on save.

## Layout

```
python/
├── belfem-conf              entry point; no install needed
├── README.md
└── belfem_conf/
    ├── repo.py              locating the checkout, whitespace-blind source search
    ├── schema.py            reads doc/input_schema.yaml as a tree, not as lines
    ├── drift.py             the two drift directions
    ├── enums.py             re-derives C++ enum value lists from to_string()
    ├── units.py             re-derives the unit table from unit_to_si()
    ├── mesh.py              just enough gmsh to check ids against a mesh
    ├── validate.py          the deck checks
    ├── parser.py            lossless span parser for input.conf
    ├── document.py          a deck you can read, and edit without reflowing it
    └── cli.py               argument handling and output
```

`schema.py` walks the parsed YAML rather than the text, which is what fixed the
coverage bug above, and reports findings by **YAML path** rather than line
number — appropriate for a file whose design premise is that line numbers rot.

## When it fires

Add the key to `doc/input_schema.yaml` and `doc/input_file_reference.md`, in the
same session as the code change. If the key changes how the *mesh* is enriched,
also add it to the tag builder in `src/fem/maxwell/fn_mesh_config_tag.hpp` —
see the `mesh_config_tag:` marking in the schema — or a cached `.bfm` will be
reused silently.

## Planned

The wxPython configurator is planned on the same schema; see
`todo/input_conf_configurator_plan.md`. Of the schema's `audit_todo` gaps,
T1–T3 are closed and enforced; T4 (mesh-state gates and the material-shape
prose) still bounds what `check` can see — the skips are reported by `-v`.
