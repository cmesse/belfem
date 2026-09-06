# input.conf Configurator: Schema-First Core with a wxPython Front End

> **DEFERRED 2026-09-03** (todo/ currentness sweep, round 3): the defect-repair half is done; the GUI configurator itself never started and is a sound idea for after the release. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-10
**Purpose:** Author and validate `input.conf` decks from a GUI instead of by hand. The mechanism:
a machine-readable schema of the input contract drives a headless validator, a mesh-introspection
layer that *derives* the `topology` block from the gmsh physical groups, and a wxPython skin
generated from that same schema.
**Module:** new `tools/confedit` (Python; no C++ change proposed)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status (updated 2026-08-11):** IN PROGRESS. The *defect-repair* half is done and reviewed; the
*configurator* half has not started.

| | |
|---|---|
| **Landed + jury-reviewed** | 8 code fixes across 6 files (§4.0), 6 documentation fixes, the two-artifact policy in `CLAUDE.md`, and `doc/input_schema.yaml` covering all 9 top-level sections |
| **Not started** | every step that builds the tool: parser/writer (R1), validator (R4), mesh introspection (R5–R7), wx skin (R9–R10) |
| **Blocked on Christian** | **O6** — validator-first vs GUI-first |
| **Blocked on work** | R4 needs the four `audit_todo` gaps closed in the schema (T1–T4) |
| **Unverified** | nothing has been built or run; the `corc` / `helix` smoke run is the gap between "reviewed by three AIs" and "verified" |

> **The schema is written, audited and self-checking — but it is NOT finished.** Both auditors
> independently concluded that a validator cannot be written from it alone, and four gaps remain
> open (T1–T4, listed in its own `audit_todo`). This is the one status point worth being precise
> about: R4 depends on exactly those gaps, so "the schema is done" would be a costly misreading.

Three audit rounds (Codex + Grok, blind) turned up **six defects in BELFEM itself** — D2, D9, D11,
D12, D15, plus the two the jury found in my own fixes (J1, J2). All are fixed. They were
independent of whether the configurator is ever built, which is the main argument that this
campaign has already paid for itself.

One CRITICAL assumption of the first draft was refuted and corrected: **D1 — deck ids are geometry
(elementary) tags, not physical tags.** That correction is load-bearing for R5–R7.

> **Scope guards:**
> - **OUT of scope:** any change to the C++ parser or to the input contract itself. This tool reads
>   and writes the existing format; it does not extend it.
> - **OUT of scope:** replacing `doc/input_file_reference.md`. The schema is a second, machine-readable
>   view of the same contract; the prose reference remains the human document.
> - **OUT of scope for v1:** `circuit`, `initial conditions`, thermal BCs, 2-D deck specifics.
> - **Compatibility promise:** a deck opened and saved unchanged must be byte-identical.
> - **IN scope:** `mesh`, `solver`, `materials`, `layers`, `homology`, `topology`,
>   `boundary conditions` (maxwell).

---

## 1. Current Behaviour and How It Fails

Decks are hand-written. The format itself is small and regular: `src/io/cl_InputFile.cpp:56-140`
strips comments and splits on `{`, `}`, `;`, while `cl_Input_Section.cpp` builds a nested type/label
tree. The hard part is that **the deck encodes mesh facts the author must remember**, with no checker.

| Failure | Mechanism | Evidence |
|---|---|---|
| Wrong sideset/block id | Ids are bare integers with no names in the deck; the author must remember that `4:9` means the six tapes | `examples/corc/input.conf:93,98` vs `$PhysicalNames` in `corc.msh` |
| `curves` block hand-written | 12 entries, each a sideset intersection `12 @ 4`, all mechanically derivable | `examples/corc/input.conf:101-115`; derivation reproduced 12/12 tonight |
| BC curve lists hand-maintained | `input curves : 1,3,5,7,9,11` must stay consistent with the `curves` numbering above it | `examples/corc/input.conf:128-129` |
| Periodic node triples hand-picked | Three node ids per plane, order-significant | `examples/corc/input.conf:119-120` |
| Deck goes stale after a remesh | Nothing revalidates ids against the mesh; a renumbered sideset silently changes the physics | no existing check |
| Silent unit/type traps | `background`'s `direction` goes through `std::stoi`, so `0, 0.5, 0` becomes `0,0,0` | `cl_Input_Section.cpp:731` (verified: that line is exactly the `std::stoi` push) |
| Bracket groups look like formatting | `1,3,5` is three BCs, `[1,3,5]` is one | `cl_Input_Section.cpp:632-712` |

**Bottom line:** the deck is a set of unchecked text-to-mesh cross-references, each either derivable
from the mesh or checkable against it.

---

## 2. Architecture: Why Schema-First, and Why the Schema Must Not Be Line-Anchored

### 2.1 The measured case against pure md-policy

Christian framed the core discipline as keeping the schema updated, "prescribed as md-policy the same
way we currently update the documentation file." That policy already exists for
`doc/input_file_reference.md`. **Tonight I measured whether it held.** Every documented key in that
file has a `file:line` parse-site citation; I walked all 87
(script: scratchpad `check_citations.py` / `drift_kind.py`, not checked in).

| Outcome | Count | Meaning |
|---|---|---|
| Cited line still correct | 26 | policy held |
| **Line drift** — key parsed in that file, at a different line | **26** | prose true, anchor rotted |
| Apparent divergence | 5 | **4 are false positives** of my heuristic (they cite the *consumer* or *mechanism*, not a quoted key string); 1 (`solver` → `cl_FEM_Controller.cpp:2240`) is probably drift, low confidence |

The drift has a signature: **a constant `+342` line offset across the entire solver block** — 24 of
the 26 rotted citations shift by exactly the same amount. Example: `tolerance` is documented at
`cl_FEM_Controller.cpp:2569-2576`; those lines are `print_line_thermal()`. The real parse is at
`:2911-2913`.

**The decisive detail: the anchors were never correct.** `doc/input_file_reference.md` and
`cl_FEM_Controller.cpp` were modified *in the same commits* (`77011047`, `1d6ef305`, `bc578b5e`) —
the living-document policy was followed to the letter. Yet in commit `bc578b5e` itself:

```
doc/input_file_reference.md:144   | `tolerance` … | `:2569-2576` |
cl_FEM_Controller.cpp:2911        if ( tNonLinear->key_exists( "tolerance" ) )
```

Same commit, 342 lines apart. **The citation was stale before it was committed** — the author read
line numbers from a working copy, then kept editing the file. That commit added 637 lines and
removed 58.

The conclusion is therefore **not** "the policy failed." The prose is accurate and current: the doc
records `absolute tolerance` becoming consumed on 2026-08-09 and the `stall tolerance` retune. Human
discipline held for *semantics*, as designed. The anchor format fails:

- at **authoring time**, not gradually — a fresh citation can already be wrong,
- on **unrelated** edits — nobody touched the input contract; someone added PID timestepping,
- in **bulk** — one commit invalidated 24 citations at once,
- **silently** — nothing checks it,
- **unmaintainably** — no human re-walks 87 citations before every commit.

Corroborating precedent: the 2026-08-09 currentness sweep (`todo/README.md:36`) already had to
re-baseline the line anchors of four todo plans, noting "all had drifted." Hand re-baselining is a
recurring chore that decays within days.

### 2.1b The second, worse tier: semantic divergence (Codex, verified)

Line drift misleads a schema, but not a human reading prose. Codex found a worse tier my line-based
scan could not detect: one that misleads **users**.

> `doc/input_file_reference.md:394-397`: "The current 3-D path still aborts after calling the
> connector builder with `not implemented yet` (`cl_ThinShellFactory.cpp:320`), so treat the key as
> WIP until that stop is removed."

Verified: the string `not implemented yet` **does not exist anywhere in
`cl_ThinShellFactory.cpp`**, and `create_side_connectors` is called unconditionally for 3-D at
`:303-318` and implemented at `:420`. The document tells users a working feature is unusable.

So the drift is two-tier, and the tiers need different defenses:

| Tier | Example | Detected by | Defense |
|---|---|---|---|
| Rotted anchor | `tolerance` cited 342 lines off | line-walk script (26 found) | never store line numbers; anchor by searchable token |
| **Semantic divergence** | `edge coating` documented as WIP, actually shipped | **neither script — needs code→doc key inventory + human review** | code→schema inventory diff (§2.3) + periodic review |

**This argues for the tool independent of the GUI:** the reference doc has at least one user-facing
falsehood today, and nothing would have caught it.

> **Design consequence (the central decision of this plan): the schema anchors keys by
> _searchable token_, never by line number.** The anchor for `tolerance` is the literal parse call
> `key_exists( "tolerance" )` / `get_real( "tolerance" )` in a named file — a string a checker can
> relocate regardless of drift. Line numbers may be *reported* by tooling; they are never *stored*
> as truth.

### 2.2 The spine

1. **`tools/confedit/schema.yaml`** — the input contract: sections, keys, type, unit dimension,
   enum, default, case-sensitivity, required-ness, and a searchable `anchor:` per key.
2. **Core library, no GUI** — schema loader, splice-preserving parser/writer, validator,
   gmsh-backed mesh introspection.
3. **`belfem-conf check <deck>`** — headless CLI. Works over SSH, in a git hook, from vim.
4. **`belfem-conf drift`** — the schema-vs-source checker (§2.3). ~50 lines, and the reason
   md-policy becomes sufficient rather than aspirational.
5. **wxPython GUI** — a thin skin generated from the schema.

**Rejected alternative:** hand-coded wx tabs with the contract embedded in Python. Rejected because
it creates a third copy of the contract under the same discipline that just produced 26 rotted
anchors, with no mechanical check.

### 2.3 What the drift checker actually does

Two directions, and the second is the valuable one:

- **schema → code:** for every schema key, confirm its anchor string still appears in the named file.
  Catches renames and removals.
- **code → schema:** enumerate every `get_real( "…" )` / `get_bool( "…" )` / `key_exists( "…" )`
  literal in the factories and controller, and report any that the schema does not know about.
  **This catches the actual dangerous drift — a key added to the code and to no document.** It is
  also the one failure mode no amount of human discipline catches, because the author of the new key
  is precisely the person who does not know the document exists.

Bonus: because the checker locates the true line for every key, it can mechanically **repair** the
stale citations in `doc/input_file_reference.md`. Done 2026-08-11 — see R8.

> **Both directions are prototyped and both have already earned their keep**, but they live in a
> session scratchpad, not in the repo. Promoting them is the single highest-value remaining step,
> and the one that would make the two-artifact policy hold without a human running a review — see
> §4.0, where that policy failed in my own hands within a day of my writing it.

---

## 3. Gap Table

Per deck section: is it derivable from the mesh, ambiguous, or must it be authored?

| # | Section / state | Needed for | Handled today | Class | Rationale |
|---|---|---|---|---|---|
| 1 | `mesh { file, unit }` | everything | hand | (c) authored | file picker + unit enum; trivial |
| 2 | `solver` (all subsections) | run control | hand | (c) authored | pure schema-driven form; 30+ keys, `cl_FEM_Controller.cpp:2911+` |
| 3 | `materials` | physics | hand | (c) authored | user-named subsections, 4 mutually exclusive shapes (`curve`/`builtin`/`usermat`/fatal), `cl_MaterialFactory.cpp:43-295` |
| 4 | `layers : <tape>` | thin shell | hand | (c) authored, **ordered** | read raw line-by-line, repeats allowed, order = bottom→top, unit mandatory (`cl_MaxwellFactory.cpp:2601-2621` — the doc's `:2515-2534`, which I copied into the first draft, is itself drifted) |
| 5 | `homology { algorithm }` | cuts | hand | (c) authored | single enum |
| 6 | `topology` block/sideset ids | domain tagging | hand | **(a) derivable, corrected** | ids are **geometry/elementary tags**, not physical tags; `$PhysicalNames` is a display overlay only — see O1 (REFUTED) |
| 7 | `topology → curves` | terminals | hand | **(b) → (c) proposal only** | node-set intersection reproduced 12/12 on corc, but `CurveFactory::intersect` is the real algorithm and rejects malformed graphs (`cl_CurveFactory.cpp:825-913`, `:1003-1055`); see O5 |
| 8 | `topology → periodic` | periodic BC | hand | **(b) guarded** | plane pair + rigid transform; helix is not a translation; see O2 |
| 9 | BC `input/output curves` | current path | hand | **(b) proposal only** | follows from curve numbering + an "enters at" choice, but polarity is physics the tool must not silently pick |
| 13 | **Gating: deck sections ignored when the mesh already has them** | correctness of the whole feature | n/a | **(c) explicit** | `curves` is parsed only if `mMesh->curves().size() == 0`, `periodic` only if `!mMesh->has_periodicity()` (`cl_MaxwellFactory.cpp:111-125`). Writing a `curves` block for a `.bfm`-backed run is a **silent no-op** — the GUI must say so |
| 14 | Material section **type** is the material name | materials tab | n/a | **(c) explicit** | `tMatLabel = tMatSection->type()` (`cl_MaterialFactory.cpp:43`) — name is the section *type*, not its label |
| 10 | BC source function keys | excitation | hand | (c) authored | schema-driven; `type` is case-sensitive and effectively required |
| 11 | Bracket grouping in BC lists | BC identity | hand | **(c) explicit** | `1,3,5` ≠ `[1,3,5]`; the GUI must expose grouping, not hide it |
| 12 | Comments / formatting of existing decks | reviewability | n/a | **(c) explicit** | splice-preserving writer; see O3 |

### 3.1 Cross-cutting findings

- **Order is load-bearing in three places** and must never be normalized: id lists are deliberately
  not sorted or uniquified (`cl_Input_Section.cpp:557-628` — sorting broke periodic node pairing);
  `layers` order is physical; `periodic` source/target must stay in *corresponding* order.
  A GUI that round-trips through a `dict` or a `set` breaks all three silently.
- **Derived ≠ authoritative.** BELFEM only consults `curves` / `periodic` when the mesh does not
  already carry them. The tool proposes; the user confirms; the deck stays the record.
- **The tool must refuse to guess** where the mesh is ambiguous (O2), rather than emit a plausible
  wrong answer. A wrong derived `periodic` block is worse than no tool.

---

## 4. Ordered Steps

- [x] **R1** — **LANDED 2026-08-11.** `python/belfem_conf/parser.py` + `document.py`.
      **All six example decks reproduce byte-for-byte.**
      - [x] Spans, not line numbers: every node holds the byte range it came from, nothing is
            normalised on the way in, and an unedited deck serialises to the original bytes.
            Edits replace one span, so O3 (splice vs canonical) is settled in favour of splice —
            the machinery was cheaper than either auditor expected because the parse is already
            span-based.
      - [x] **The gate was strengthened after it first passed, because it was nearly
            tautological.** `dumps()` returns `self.text` when unedited, so the identity check
            would pass even with no tree at all. `roundtrip` now probes **every value span**:
            rewrite that value, require exactly those bytes to change, re-parse and read it back.
            **284 spans across the six decks.**
      - [x] Semantics verified by hand against `corc`: the 7-layer stack keeps its order and its
            duplicate `copper`/`silver` entries; `sidesets : 4:9` splits on the FIRST colon only,
            so the range survives; `1 : 12 @ 4` parses as key `1`; a bare `builtin ;` flag becomes
            the string `"true"`, as the C++ does.
      - [x] Comments are masked before structure scanning, so a `{`, `}` or `;` inside a comment
            is not mistaken for syntax; `newline=""` on read and write, so a CRLF deck is not
            silently rewritten to LF (that would break byte-identity invisibly).
      - [x] Runs on the SCLS Python 3.9 — `Path.read_text(newline=...)` is 3.13+, which the first
            version used and which failed immediately.
      ~~Original wording: parser/writer core with source line spans; round-trip test over
      `examples/*/input.conf`. Superseded by the entry above, which delivers it with byte spans
      rather than line spans and a stronger gate.~~
- [◐] **R2** — `doc/input_schema.yaml`: **written and audited, not finished.**
      - [x] All nine top-level sections, not the three originally planned. Built from a mechanical
            inventory of parse-call literals (110 distinct keys across 10 consumers) rather than by
            transcribing the prose — which is why it caught things the prose had wrong.
      - [x] Audited 2026-08-11 by Codex and Grok blind; ~20 findings applied, every one re-verified
            against the code first. Defaults now come from the members in `cl_FEM_Controller.hpp`,
            not from the prose.
      - [x] `runtime_status` vocabulary adopted (Codex's idea): `live` / `accepted_then_hard_error`
            / `inert` / `silent_noop` / `corrupts` / `rejected`. Encodes that *accepted ≠ usable*,
            which a type-and-enum check alone cannot express.
      - [x] Anchors are bare quoted key literals throughout; 114 of 115 resolve (the one failure is
            prose in the file's own header). All 81 call-form anchors were converted after both
            auditors objected that they only passed because the checker normalises whitespace.
      - [x] Kept in lockstep with the 2026-08-11 code fixes, as its own policy requires.
      - [x] **T1 — CLOSED 2026-08-11, and not by transcribing it.** The C++ resolves these enums by
            comparing the deck value against `to_string( value )` for every enumerator, so the
            switch is machine-readable and hand-copying it would have made a fourth artifact to
            rot. Each key now carries `enum_ref`, and `belfem-conf drift` re-derives the list from
            `src/` and fails on any difference — **six lists checked, not transcribed.**
            Negative-tested all four ways it can fail (value omitted, value invented, unknown enum,
            `enum_ref` with no `values:`).
            **Found by the new check on its first run:** `values: [off, blr, automatic]` was
            loading as `[False, blr, automatic]` — YAML 1.1 reads a bare `off`/`on`/`yes`/`no` as a
            boolean, so the schema had silently lost a legal `compression scheme` value. Now quoted.
            Also recorded: `Preconditioner::GASM` is not deck-selectable, since ASM and GASM both
            stringify to `"asm"` and the lookup takes the first.
      - [x] **T2 — CLOSED 2026-08-11.** `topology.domain_types` enumerates all 24 accepted
            spellings with the DomainType each becomes, whether it reads block or sideset keys, and
            whether `material` is required — the last from the `Domain` constructor, not guessed.
            **Verified exhaustively in both directions: 24 in the code, 24 in the schema, no
            difference either way.** The old `section_type_spellings` block was folded in rather
            than left alongside, so there is one list instead of two that can disagree.
      - [ ] **T3** — structured conditional-requiredness. `required_when` is prose today; the real
            rules needing it are `mesh.unit` (gmsh only), `material` (by domain type), `period`
            (ramp/sigmoid), period XOR frequency (periodic family), HTS `file` XOR (`jc` + `n`).
      - [ ] **T4** — structured gate predicates (`gated_by` is prose: the curves/periodic mesh-state
            gates and the `.bfm`-vs-gmsh split).
      - [ ] `initial conditions` is complete; **`circuit` is still PARTIAL** (component key sets
            sketched, semantics unverified).
- [x] **R3** — `belfem-conf drift` **LANDED 2026-08-11** in `python/belfem_conf`, promoted out of
      session scratch. Runs with no install, no build tree and no PYTHONPATH (PyYAML only), so it
      works in a git hook and on a machine that has never compiled BELFEM.
      - [x] **schema → code**: 119/119 anchors resolve.
      - [x] **code → schema**: 108 parsed literals across 19 consumers, all known to the schema.
      - [x] Coverage self-report (`checked` vs `declared in file`) printed before any finding —
            the direct lesson of the prototype that passed a check it had never run.
      - [x] Failure path exercised deliberately against a gutted schema in a throwaway tree:
            106 unknown keys and 1 unresolved anchor reported, exit 1. **A checker never seen to
            fail is not evidence of anything.**
      - [x] `schema.py` walks the parsed YAML rather than the text — the fix for the inline-flow
            blindness — and reports by YAML path rather than line number, which suits a file whose
            premise is that line numbers rot.
      - [ ] **doc → code** still not built: a claim can be false while every key and anchor is
            correct (the `edge coating` class). Needs a reader, not a grep.
      - [ ] Not yet wired to a hook or `make check`; run manually for now (O7).

      **Found on its first run, and fixed:** the schema described the `boundary conditions`
      `maxwell` and `thermal` subsections only in prose, so a validator reading the schema would
      have rejected a legal deck. Both are now modelled with anchors. Two checker bugs surfaced the
      same way — section names and `resolution_chains` spellings were not being collected — which
      is the tool auditing itself before it audits anything else.
      - [x] **code → doc key inventory** (`key_inventory.py`, `key_sections.py`): 110 keys with
            section attribution from the receiver variable. Found **no missing keys** — the prose
            reference's key coverage is complete. Both auditors independently confirmed this.
      - [x] **schema → code anchor check** (`check_anchors.py`): all 53 anchors resolve. It caught
            **5 broken anchors in the schema on first run**, which is the point. It also produced the
            anchor convention now recorded in the schema header: prefer the bare quoted key literal,
            because `get_value( "simulation time" )` does not match the real
            `get_value( "simulation time", "s")`.
      - [ ] **doc → code**: not built; the `edge coating` class of semantic divergence still needs a
            human or an LLM pass, not a grep.
      - [ ] promote from scratchpad into the tool proper.
- [x] **R4 — LANDED 2026-08-11.** `belfem-conf check <deck>`; all six example decks validate.
      - [x] Structure (required + unknown sections), unknown keys, enum values, numeric/boolean
            types, **unit dimensions**, topology domain types (incl. rejecting `cut` and demanding
            `material` where the `Domain` ctor does), `block` xor `blocks`, `thinshell` ↔ `layers`
            pairing, layer materials resolving, `periodic` triples of exactly 3, and the
            `background direction` `std::stoi` landmine.
      - [x] **Units derived, not transcribed.** `unit_to_si` is ~800 lines of
            `else if ( tUnit == "V" ) { tMass += tPower; ... }` — a machine-readable table. Parsed
            at run time, so the dimension comparison is BELFEM's own and cannot drift. 115 tokens.
      - [x] **Negative-tested.** A deliberately broken deck triggers **16 errors** across every
            check. A validator that has only ever said "ok" is worth nothing.
      - [x] Reports what it did NOT check (`-v`): conditional requiredness (T3), mesh-state gates
            (T4), bracket-group arity. Coverage honesty, same rule as the drift checker.
      - [x] **`applies_to` split by axis, 2026-08-11 — a gap R4 itself exposed.** The field named
            four axes at once, so none was machine-usable: BC type, source-function type, material
            shape, and group kind. Now `applies_to_bc_type`, `applies_to_function` (with
            `function_families` defining `periodic family` = sine/square/triangle/sawtooth,
            confirmed from the factory branch that reads period-or-frequency and phase),
            `applies_to_material_shape` (still prose, still unenforced, and now honestly named),
            and plain notes where it merely repeated `domain_types.groups`.
            **Demonstrated before and after:** a `bearing` carrying `length`, `direction` and
            `fuzzyness`, plus a `ramp` carrying `phase`, validated **clean** beforehand — every one
            of those keys is read by nobody. All four are now reported, and the six real decks
            still pass, so the new check is not firing spuriously.
            Renaming was done by classifying each value against the BC-type and function-type
            vocabularies rather than by line number, so a value could not be misfiled.
      - [◐] **T3 started 2026-08-11 — and the first half needed no new vocabulary at all.**
            `xor_with`, `required_for` and `pairs_with` were already declared in the schema and
            **none of them was ever read**; `validate.py` re-implemented the block/blocks xor by
            hand instead. Wiring them removed that hardcode (and the duplicate report it produced
            once the generic check existed), and turned on `period`/`offset` requiredness for
            ramp and sigmoid. Also added the one genuinely new primitive:
            `required_unless: {sibling: file, ends_with: ".bfm"}`, which covers `mesh.unit` and,
            unexpectedly, the whole `.bfm`-versus-gmsh third of T4.
            Topology domain sections were being key-checked not at all; they are now.
            **Regression caught by the real decks and fixed:** matching an ancestor section by
            name anywhere in the chain made `circuit { topology { … } }` validate its components
            against *mesh* topology keys, reporting 27 false unknown-key warnings. The outermost
            ancestor is what decides.
            Remaining in T3: one-of GROUPS (`period` XOR `frequency`; HTS `file` XOR `jc`+`n`) —
            `xor_with` links pairs, and neither of these is a pair; and section fallback
            (`linear` unless `linear magnetic`), whose predicate is a sibling SECTION rather than
            a key value.
      - [x] **T4 re-scoped, not started.** Working through T3 showed T4 is not one task: its
            `.bfm`-vs-gmsh gate is the same predicate as T3's and is now done, while the
            curves/periodic gates depend on data **outside the deck**. They cannot be expressed as
            a predicate over input, so **T4 is blocked on mesh introspection (R5), not on
            vocabulary design.** Recorded in the schema so nobody designs a predicate language for
            it. `length` required in 2-D moves here from T3 for the same reason.
      - [ ] Not wired to a hook or `make check` (O7, shared with R3).

      **Two parser bugs the round-trip gate could not have caught, found by pointing the validator
      at real decks:**
      1. The section-header span runs from the previous terminator, so the FIRST section in a file
         swallowed the licence banner and parsed as a section named `//---- … mesh`. My earlier
         semantic spot-check missed it because I happened to query `topology` and `layers`, which
         sit after a `}`.
      2. The header must be the LAST line before `{`, as the reference says and the C++'s
         line-split buffer does — not everything since the terminator. `examples/sidecoating` has a
         stray `/` on its own line above `homology`; the C++ ignores it, and now so does this, with
         a `note` rather than in silence.
      Both are cases where byte-identical serialisation was already passing while interpretation
      was wrong — the gate tested the writer, not the reader.

      **Also found, reported not fixed:** `unit_to_si` gives `T`, `G` and the whole flux-density
      family the dimensions of **volts** (`mass·length²·current⁻¹·time⁻³`) rather than tesla
      (`mass·time⁻²·current⁻¹`) — the branch looks copy-pasted from `V`, scale correct, dimensions
      not. Self-consistent for T-vs-T checks, so no deck is affected today, but it makes flux
      density and voltage indistinguishable to `check_unit`. See §8.
- [◐] **R5 / T4 first step — LANDED 2026-08-11.** `python/belfem_conf/mesh.py`, wired into
      `belfem-conf check`.
      - [x] Resolves `mesh { file }` relative to the deck, reads a gmsh 2.x mesh, and **verifies
            every block and sideset id the deck names exists in it** — the check that motivated
            this whole campaign. Ids come from the **elementary** tag, per D1.
      - [x] Dimension from the z extent of the nodes, matching `GmshReader` rather than being
            right in the abstract — a curved 2-D surface in 3-D reads as 3-D there too.
      - [x] Negative-tested: `blocks : 1:5` and `sidesets : 4:99` against the real `corc.msh` are
            both caught, with the id list capped so a mistyped range does not bury the message.
      - [x] Notes when a sibling `.bfm` exists, since the factory prefers it when its stamps
            match — so the ids that actually run may not be the ones checked.
      - [x] **gmsh 4.x parsing — and it was not optional.** Nothing in the examples pins a mesh
            format and gmsh has defaulted to 4.1 since 4.0, so **every freshly generated mesh is
            4.1**: a 2.x-only reader would have checked nothing in practice and worked on `corc`
            only because that `.msh` is an old 2.2 file kept in git. Verified by converting the
            real corc mesh and requiring both formats to give identical blocks, sidesets and
            dimension — they do. Also validated the helix deck against a mesh generated fresh
            from its own `.geo`.
      - [x] **`length` required in 2-D** — `required_when_dimension: 2`, verified both ways on a
            flat mesh generated for the purpose.
      - [ ] `.bfm` not read — HDF5 needs h5py, which is **not installed here** and which this
            tooling avoids depending on so the checks run in a hook and on a machine that never
            built BELFEM. **This is what still blocks the mesh-STATE gates**, since only an
            enriched mesh can carry curves or periodicity.

      **Hazard found while verifying, worth knowing independently of this tooling:** converting a
      mesh from 2.2 to 4.1 through gmsh **drops the type-15 point elements**. The real corc mesh
      has 10; after `gmsh v22.msh -save -format msh41` there are no dim-0 element blocks at all.
      BELFEM builds its vertices from those, and a missing vertex surfaces as a singular matrix or
      a bare "Key not found in map". Generating fresh from `.geo` at 4.1 keeps them (helix: 132),
      so this is a conversion hazard only — do not round-trip a BELFEM mesh through a format
      change.

      **The design constraint, found by looking rather than assuming:** the examples ship `.geo`
      and generate the `.msh`, so **five of the six decks have no mesh at all**. A missing mesh is
      NORMAL and must never be an error — it makes a class of check unavailable, and the validator
      now names which class and why. Had I assumed meshes were present, the check would have
      failed five of six decks on its first run and been switched off.

- [ ] **R5 (remainder)** — Mesh introspection keyed on **geometry/elementary tags**, mirroring
      `GmshReader::create_group_ids` (`cl_Mesh_GmshReader.cpp:771-814`) — not on `$PhysicalNames`.
      Names attach as a display overlay, and only when the physical→geometry map is 1:1; otherwise
      the id is shown bare with a note. Blocks and sidesets live in separate id spaces and may share
      a number (`cl_Mesh.cpp` `block_exists`/`sideset_exists`), so never present a flat "id" list.
      *(after: R1)*
- [ ] **R6** — Cross-reference validation: layer materials resolve, `thinshell` label has a matching
      `layers` block, BC curve ids exist, periodic triples are 3 each, deck ids exist in the mesh.
      *(after: R4, R5)*
- [ ] **R7** — Topology **proposal** (not derivation): surface candidate `curves`, BC in/out lists,
      and `periodic` triples for confirmation, never silent write. **Must fail closed** on: multiple
      disconnected intersection chains, branch degree > 2, intersections with ≤ 1 shared corner node,
      terminal-vs-structural ambiguity, current polarity, and any non-translation periodic map.
      Prefer calling BELFEM's own `CurveFactory::intersect` / `thin_shell_side_curves`
      (`cl_CurveFactory.hpp:51,58,74`) over the Python node-set heuristic — the C++ version already
      builds ordered chains and rejects malformed graphs (`cl_CurveFactory.cpp:825-913`,
      `:1003-1055`). *(after: R5, R6)* — **see O5; this step may reduce to "report candidates".**
- [x] **R8** — Stale citations in `doc/input_file_reference.md` repaired mechanically from the
      located parse sites, 2026-08-11. **50 of 57 walked citations resolve, up from 24.** Two
      passes: a script pass, then a second covering the multi-site rows the jury found it had
      skipped. The remainder are rows whose key is too short to disambiguate (`n`, `t`, `ec`, `jc`)
      or which cite a consumer rather than a parse call; those were left alone rather than guessed
      at, and the header now says so.
      **Method note worth keeping:** the script's first version matched the *first* quoted
      occurrence and would have written new wrong citations (`absolute tolerance` → the thermal
      block; short keys matching comments). It was only safe once it required a parse-call context
      and picked the candidate nearest the existing citation.
- [ ] **R9** — wxPython skin: notebook, schema-generated forms, grid editors for the four
      non-form sections (materials, layers, curves, BC groups). *(after: R4)*
- [ ] **R10** — gmsh pick round-trip via `gmsh.fltk.selectEntities()`. *(after: R9)* Separate
      window, not embedded — FLTK cannot live in a wx panel.
- [ ] **R11** — End-to-end gate: derive `examples/corc` topology from `corc.msh` and diff against
      the committed deck; must reproduce all 12 curves, both BC lists, and the periodic triples.

---

## 5. Open Design Questions (not silently decided)

- **O1 — Are deck ids always gmsh physical tags? → REFUTED 2026-08-10 (Grok, re-verified by Claude).**
  **They are geometry (elementary) tags. `$PhysicalNames` is read and thrown away.**
  `cl_Mesh_GmshReader.cpp:335,338` sets both tags from words 3 and 4, but every grouping site reads
  `geometry_tag()` and nothing ever consumes `physical_tag()`: `create_group_ids` (`:771-814`),
  blocks (`:842-889`), sidesets (`:961-986`). `read_physical_tag` (`:231-241`) only advances the
  buffer cursor past the block. In MSH 4.1 only `set_geometry_tag` is called at all.
  My corc correspondence held **by coincidence**: every nonzero physical group maps 1:1 to its
  elementary tag (`phys N → geom {N}`), which is what gmsh produces when each physical group wraps
  exactly one elementary entity.
  **Design consequence:** the tool must read the *elementary* tag to get deck ids, and may show
  `$PhysicalNames` labels only as a **convenience overlay, valid only when the physical→geometry map
  is 1:1**, saying so when it is not. A mesh built with `Physical Volume("Air") = {1,2,3};` would
  give one name over three deck ids — a GUI offering "Air" as one pickable entity would emit the
  wrong id. Confidence: high (source-verified).
- **O2 — When must periodic derivation refuse?** The corc map is a pure z-translation of 12π, but the
  project's periodic map is a general rigid transform (translate + twist + tilt; helix cases map
  conductor k to k−1). An axis-aligned derivation would be confidently wrong on the meshes that matter most. Options: (a) detect a general rigid transform and refuse if
  residual exceeds tolerance; (b) derive only when the mesh already carries periodicity, else ask;
  (c) always require confirmation. Leaning (a)+(c).
- **O3 — RESOLVED 2026-08-11 by building it: splice, and it was cheap.** Codex was right that
  editing git-tracked, hand-commented decks needs splice; Grok was right that a bolted-on splice
  layer would be a research project. Both were reasoning about splice as something added *after* a
  normalising parse. Making the parse span-based from the start means splice is not a layer at all
  — the original text IS the document, and an edit replaces one byte range. The whole writer is
  ~30 lines (`Document.dumps`), and byte-identity is provable rather than tested for.
  *Superseded discussion follows.*
- **O3 (superseded) — Splice-preserving writer, or canonical formatter? THE AUDITORS SPLIT.**
  - **Codex: splice is worth it.** The C++ parser destroys comments and trivia on load
    (`cl_InputFile.cpp:55-68,73-137`), layer order and repeated names are meaningful
    (`cl_MaxwellFactory.cpp:2601-2611`), and reflowing examples on every GUI save produces noisy
    diffs. Canonical formatting is fine for *newly generated* decks or an explicit `format` command.
  - **Grok: splice is "a research project dressed as a feature."** Costs cited: multi-statement
    lines (`tidy_up` splits on `;`, `:104-120`), layers as a raw buffer, keys containing spaces.
    Recommends validate-only first, and if writing is needed, deterministic rewrite of whole
    sections rather than line-span splicing.
  - **My read:** they are answering different questions. Codex is right that *editing existing
    git-tracked decks* needs splice; Grok is right that it is not needed if v1 never writes. The
    decision therefore collapses into O6 — if v1 is validate-only, O3 does not arise yet.
- **O4 — Where does the schema live?** `tools/confedit/schema.yaml` (with the tool) or `doc/`
  (with the contract it mirrors)? The md-policy line must name one path.
- **O5 — Is the intersection heuristic safe in general?** It worked on one mesh. Tangent surfaces,
  non-manifold junctions, and terminals that are a strict subset of an intersection are untested.
  *Dispatched as Q3.*
- **O6 — RESOLVED 2026-08-11, Christian: validator first, GUI after.** Module home is `./python`
  (not `tools/confedit` as originally drafted, and not `scripts/`, which is dev/build tooling wired
  into CMake targets). Suggested package path `python/belfem_conf/`.
  *Superseded discussion follows.*
- **O6 (superseded) — the one thing gating further work.** Update 2026-08-11: the argument
  has shifted since it was written. The defect-repair half is finished, so the release-risk concern
  that motivated "validator first" is largely spent — those fixes are in and reviewed. What remains
  is genuinely new work in either order, and R1/R2/R4 are prerequisites for the GUI regardless.
  The practical question is now narrower: **close T1–T4 and build the validator, or start the wx
  skin against the schema as it stands and accept that its enum coverage is incomplete?**
  Original framing follows.
- **O6 (original) — Does this compete with the end-of-August release? → BOTH AUDITORS SAY YES, INDEPENDENTLY.
  This decision needs Christian; it is a scope question, not a go/no-go.**
  Neither says "don't build it"; both say **build the validator first and the wx skin after**.
  - Codex: scope to (1) `belfem-conf check`, (2) schema inventory + drift test, (3) mesh id /
    topology diagnostics, (4) "wx GUI later, generated from the proven schema."
  - Grok: "Do not build the wxPython multi-tab author before the August release… GUI after release,
    only if `check` earns trust." Reasons: schema truth is unsolved (§2.1), the mesh-id model was
    mis-stated (O1), auto-topology is unsafe without refuse-to-guess discipline (O5), and — the
    sharpest point — **release work is actively churning the very files the schema tracks**
    (`Controller::set_params` grew 637 lines today), so every release fix re-breaks line citations.
  - **Counterpoint for the record:** Christian is building wxPython tonight and wants to write the
    GUI tomorrow. The auditors' ordering does not forbid that — R1/R2/R4 (parser, schema, validator)
    are prerequisites for the GUI *anyway*, since the GUI is generated from the schema. Building
    those first is not a detour; it is the first half of the same work. The only thing genuinely
    deferred is R9/R10 (wx skin, gmsh picking).

- **O7 — Should the checker also gate CI?** Both auditors propose failing a build on unknown parse
  keys. The repo has no CI (per the build-flags note) and `USE_TEST` is off by default in the shared
  tree, so "fail CI" has no home today. Options: a `make check` test, a post-commit hook alongside
  the existing autoreview hook, or a manual `belfem-conf drift` run. Leaning: manual first, hook
  second — a gate nobody runs is not a gate.

---

## 6. Schema Sketch

```yaml
# tools/confedit/schema.yaml  (shape only — not the final key set)
solver:
  required: true
  consumer: cl_FEM_Controller.cpp
  sections:
    nonlinear:
      aliases: [nonlinear magnetic]        # alias wins if present
      required: true
      keys:
        tolerance:
          type: real
          default: 1.0e-6
          aliases: [relative tolerance]    # consulted only if `tolerance` absent
          anchor: 'key_exists( "tolerance" )'   # searchable — NEVER a line number
        anderson depth:
          type: int
          range: [0, 8]                    # else fatal
          default: 0
          note: absent → timestep{anderson stabilization:true} fills 3
          anchor: '"anderson depth"'
topology:
  sections:
    periodic:
      keys:
        source: {type: id_list, exact_count: 3, ordered: true, pairs_with: target}
        target: {type: id_list, exact_count: 3, ordered: true}
```

Three properties carry the design: `anchor` is a string not a line; `ordered: true` forbids the
normalization that would break periodic pairing; `pairs_with` encodes a cross-key constraint the
validator can enforce.

---

## 7. Definition of Done

- [ ] Every gap-table row maps to a step or an open question
- [ ] Each claimed gap carries a citation or is marked an assumption (O1 is marked)
- [ ] Ordered steps with dependencies
- [ ] Open questions logged, not decided
- [ ] R1 byte-identical round-trip passes on all six example decks
- [ ] R11 derives corc topology matching the committed deck exactly
- [ ] `belfem-conf drift` reports zero unknown keys, or each is filed

---

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/input_conf_gui_plan.md` (brief), with per-voice replies in
  `input_conf_gui_plan_codex.md` and `input_conf_gui_plan_grok.md`.
- Brief dispatched 2026-08-10 with six questions; Q1 (is md-policy sufficient?) is Christian's
  stated core concern and was additionally measured directly — see §2.1.
- Both auditors ran blind (separate slugs, no visibility into each other's replies).
  **Every finding below was re-verified by Claude against the cited code before inclusion**, per the
  third-voice rule; the verification commands are in this session's transcript.

### 4.0 Fixes landed 2026-08-11 (approved by Christian, jury review pending)

Seven code fixes and six documentation fixes. **No behaviour change for any deck in the repository**
— verified per item. The schema was updated in the same session, as its own policy requires.

- [x] **D9 → fixed.** `cl_MaxwellFactory.cpp` now dispatches on
      `domain_type( tSection->type() ) == DomainType::ThinShell` at both the counting and the
      building loop, so `tape` and `shell` are true aliases. Verified no raw `"thinshell"` compare
      survives anywhere in `src/`.
- [x] **D9b → fixed.** Singular `curve { }` now falls back correctly; plural still wins if both exist.
- [x] **D9c → fixed.** `"cut"` removed from `domain_type()`; it now fails with
      `Unknown Domain Type: cut`. `DomainType::Cut` stays live for internally generated cuts and the
      enum→string direction is untouched. Brace balance of the else-if chain verified from the diff.
- [x] **D12 → fixed at the root.** The BC creation switch's `default:` is now a hard error naming
      the type. **Verified safe:** `boundary_condition_type()` maps only `neumann`, `dirichlet`,
      `bearing`, `gauge`, `current`, `voltage`, `background` — all of which have push cases — plus
      `background dirichlet` and `UNDEFINED`, both of which should error. `CircuitCurrent` and
      `CircuitVoltage` are set programmatically by the circuit factory and are **not reachable from
      a deck string**, so no legitimate deck can hit the new error.
- [x] **D11 → fixed.** The 2-D branch now also tests `CircuitVoltage`. Behaviour change confined to
      a 2-D terminal pair omitting the output list, which now requires `length`; `examples/circuit`
      is 3-D and supplies both lists.
- [x] **D15 → fixed.** `Gauge` joins `Dirichlet` in the thermal units switch (`K`), matching what
      `cl_FEM_PhysicalBoundaryCondition` already does — a gauge imposes a Dirichlet temperature.
- [x] **D7 → fixed.** Both group checks are `BELFEM_ERROR`. Also repaired a latent format-string
      bug found while editing: the xor message had two `%s` and three arguments, printing `mLabel`
      where `tGroup` belonged.
- [x] **D2, D10, D13, D14 → documented.** Reference corrected: `edge coating` is live; source
      functions exclude `bearing`; the third `builtin` path; `period` wins over `frequency`.
- [x] **B4 → documented.** The domain-type table now states that `tape`/`shell`/`curve` are true
      aliases as of 2026-08-11 and that `cut` was removed, with the reason.
- [x] **B6 → repaired mechanically.** 24 citations re-derived from the parse sites; **42 of the
      walked citations are now correct, up from 24.** The repair script required a parse-call
      context and picked the candidate nearest the original citation — its first version, matching
      the first quoted occurrence, would have written *new* wrong citations (`absolute tolerance`
      → the thermal block, and short keys like `n`/`t`/`ec` matching noise). Those rows and all
      multi-site rows were deliberately left alone rather than guessed at, and the header now says so.
- [ ] **Smoke run on `corc` and `helix`** — the only fix touching thin-shell setup is D9; a run is
      wanted before this is called done. Not run (builds are Christian's).

#### Jury review of the fixes (2026-08-11, Codex + Grok, blind, pre-registered)

**Verdict on the code: all seven fixes sound.** Grok verified F1–F7 individually at 90–95 %
confidence with citations, independently reproducing my F4 enum-coverage argument (only
`background dirichlet` and `UNDEFINED` reach `default`; `CircuitCurrent`/`CircuitVoltage` are not
returned by the string parser) and confirming no third `"thinshell"` compare exists. Both auditors
agreed the **real defect was my documentation pass, not the C++**.

- [x] **J1 — CRITICAL, a bug my own fix created (Codex).** Activating the 2-D `length` branch
      exposed that it divides by `length` with no positivity check: zero gives an infinite scale,
      negative silently flips the sign of the imposed voltage. Guarded in the circuit path **and**
      in the Maxwell path, which had the identical unguarded division all along.
- [x] **J2 — my `Voltage ||` arm was dead code (Grok, R4).** `tType` is initialised to
      `CircuitVoltage` and never reassigned, so the disjunct could never fire. Reduced to the single
      live test with a comment.
- [x] **J3 — the prose pass was incomplete and self-contradictory (both).** §13 still called the
      block/sideset rule "assert-level" *after* F7 made it an error; §9 and §13 still described
      `background dirichlet` as harmlessly inert *after* F4 made it fatal; the §9 amplitude table
      still said "others dimensionless" *after* F6 gave thermal gauge K; and
      `circuit_usage_guide.md` still called the `length` branch unreachable *after* F5 reached it.
      All corrected. **This is the two-artifact policy failing in my own hands within a day of my
      writing it** — I changed code and updated the schema, but propagated to only part of the prose.
- [x] **J4 — my "citations re-derived" banner was overstated (Grok, R1).** Nine further rows were
      still wrong, all of them the multi-site rows my script deliberately skipped. Each corrected
      value was verified against the source before applying. **50 of 57 walked citations now
      resolve, up from 24**; the banner now states the real number and why the rest were left.
- [x] **J5 — `CLAUDE.md` contradicted its own convention (Grok, R6).** The anchor example used the
      call form `'key_exists( "tolerance" )'` while the rule mandates the bare literal. Fixed.
- **J6 — noted, not acted on (Grok, R5).** F4 refuses `background dirichlet` rather than
      implementing it; `impose_bc` still carries a BackgroundDirichlet branch. The feature is
      unfinished, not merely disabled — recorded in §13 so the distinction is not lost.
- **J7 — noted, out of scope.** The thermal creation switch has the same `default: break` shape as
      the Maxwell one, but errors later before `set_function`, so it cannot corrupt a neighbor.
      Grok flags a possible `SourceFunction` leak on that abort path. Not investigated.

### Defect tracker (audit round 2026-08-10; Codex + Grok, blind)

- **D1 — CRITICAL, my error.** "Deck ids are gmsh physical tags" is false. Both auditors refuted it
  independently with different citations; verified: grouping reads `geometry_tag()` everywhere
  (`cl_Mesh_GmshReader.cpp:771-814,842-889,961-986`), `read_physical_tag` (`:231-241`) discards the
  names, MSH 4.1 never calls `set_physical_tag`. corc holds only because its physical→geometry map
  is 1:1 (measured). **Fixed in plan** — O1, gap row 6, R5 rewritten. Found: Grok + Codex.
- **D2 — HIGH.** `doc/input_file_reference.md:394-397` states 3-D `edge coating` aborts with
  `not implemented yet`; that string is absent from `cl_ThinShellFactory.cpp` and
  `create_side_connectors` is live (`:303-318`, impl `:420`). A user-facing falsehood in the
  reference today. Found: Codex. Verified: Claude. **→ needs a doc fix independent of this tool.**
- **D3 — MEDIUM, my error.** I cited the layers loop at `cl_MaxwellFactory.cpp:2515-2534` (copied
  from the reference doc); the real loop is `:2601-2621`. Illustrates the drift being discussed.
  Found: Grok. Verified: Claude. **Fixed in plan** — gap row 4.
- **D4 — MEDIUM.** Deck `curves`/`periodic` are parsed **only** when the mesh lacks them
  (`cl_MaxwellFactory.cpp:111-125`). A tool that writes a `curves` block for a `.bfm`-backed run
  produces a silent no-op. Found: Grok. Verified: Claude. **Added** — gap row 13.
- **D5 — LOW.** Material name is the section **type**, not its label
  (`cl_MaterialFactory.cpp:43`). Found: Grok. Verified: Claude. **Added** — gap row 14.
- **D6 — MEDIUM, scope.** The brief's tab list omits `initial conditions` and `circuit`, both real
  top-level sections. Found: Grok. **Acknowledged** — already declared out of scope for v1 in the
  scope guards, but the schema must not imply they do not exist.
- **D7 — MEDIUM. VERIFIED 2026-08-10.** Both the "group defined at all" check *and* the
  `block`/`sideset` xor check are `BELFEM_ASSERT`, so **a release build silently accepts a domain
  with neither key (empty group list) or with both**. Found: Grok. Verified: Claude. The validator
  must error on both regardless of what release builds tolerate.
- **D8 — LOW. VERIFIED 2026-08-10.** `//unique( tGroupIDs );` is commented out in the id-list
  reader — lists are genuinely never uniquified, confirming the ordering rule the schema encodes.
  Found: Grok. Verified: Claude.

### Audit round 2 — reference-document currentness (2026-08-10, Codex + Grok blind, split coverage)

Dispatched with complementary priorities (Codex §1-§8 first, Grok §9-§14 first) so the union covers
the document. **Both independently found the same top defect.** All findings below re-verified by
Claude against the cited code.

- **D9 — CRITICAL. `tape` / `shell` / `curve` / `cut` are inert section types.**
  `domain_type()` (`src/mesh/en_DomainType.cpp`) accepts spellings that no consumer path handles:
  `tape` and `shell` map to `DomainType::ThinShell` but `read_thin_shell_data()` gates on
  `type() == "thinshell"` exactly, so **no Protoshell is built — the deck looks right and there is
  no thin shell**. Singular `curve { }` never reaches `create_curves()` (only `curves` is looked
  up). `cut` has no read path at all. Found: Grok (tape/shell) and Codex (all four, rated CRITICAL).
  Verified: Claude. **In the schema as `topology.inert_section_types`; this is the single
  highest-value validator check.**
- **D10 — HIGH. Source-function keys do not apply to `bearing`.** The doc says they apply to all BC
  types; `bearing`'s second switch only re-reads `nodes`, in both the Maxwell and thermal factories.
  Found: Grok. **In the schema** as `does_not_apply_to`.
- **D11 — HIGH. Circuit `terminal pair` 2-D `length` is dead code.** The pair sets
  `BoundaryConditionType::CircuitVoltage` but the 2-D fallback tests plain `Voltage`, so `length` is
  never required or read and the scale stays 1.0 — while the doc promises §9 rules including
  1/length scaling. Found: **Grok and Codex independently.** In the schema as `circuit.known_defect`.
- **D12 — HIGH, live defect (not just a doc bug). `background dirichlet` corrupts a neighboring
  BC.** The creation switch falls through `default: break` without pushing, but `BackgroundDirichlet`
  is in the source-function switch, so the loop calls `set_function` on
  `mPhysicalBoundaryConditions( tCount - i )` — a *previously created* BC, or an unsigned underflow
  if it is the first. The doc's "parses but creates no BC object" understates it. Found: Grok.
  Verified: Claude. **→ this is a code fix, not a doc fix.**
- **D13 — MEDIUM. Undocumented third `builtin` selection path.** Besides `builtin : <type>` and a
  section labeled `builtin`, `cl_MaterialFactory.cpp` also accepts a key at the `materials` level
  named like the material subsection with value `builtin`. Found: Codex. Verified: Claude.
  **In the schema.**
- **D14 — MEDIUM. `period` wins over `frequency`** when both are present (`if period … else if
  frequency`); the doc says "either". Found: Grok. Verified: Claude. **In the schema.**
- **D15 — MEDIUM. VERIFIED 2026-08-10.** The thermal units switch sets `tUnits` for `Dirichlet`
  (`"K"`) and `Neumann` (`"W/m^2"`) and falls to `default: break` otherwise — but `Gauge` **is** in
  the following switch that reads `type` and `amplitude`. So a thermal gauge with a source function
  validates its amplitude against an **empty** required unit. Found: Grok. Verified: Claude.
  The same read confirms D10 on the thermal side: `Bearing` only re-reads `nodes`.
- **D16 — corrected defaults, all verified against `cl_FEM_Controller.hpp`** and now in the schema:
  thermal `tolerance switch` 1e-3 (not the magnetic 1e-4), thermal `min relaxation` 0.1 (not 0.001),
  `update gate` `BELFEM_REAL_MAX` = disabled, `coupling` fully-coupled, `floor retries` 20,
  `restart` **true**, `scheme` BDF1, `minimum timestep` 1e-10. Found: Grok.

**Missing keys: none.** My mechanical inventory and both auditors independently agree that every
string key the parser consumes appears somewhere in the reference document. The document's *key
coverage* is complete; its *behavioural claims* and *line anchors* are what drift.

### Audit round 3 — the schema itself (2026-08-10, Grok; Codex still running)

Grok audited `doc/input_schema.yaml` v1 and found it "a strong inventory scaffold, not a finished
contract." All findings below verified by Claude; **all are now fixed in the schema.**

- **D17 — CRITICAL. `mesh { unit }` is effectively REQUIRED**, not optional: on the gmsh path
  `get_string( "unit" )` is called with no `key_exists` guard, so a missing `unit` is fatal. Only
  the `.bfm` path skips it. My schema said `required: false`.
- **D18 — CRITICAL. The `scheme` enum was incomplete.** Also accepted: `explicit`, `crc`,
  `crank-nicolson`, `galerkin` — not just `bdf1`…`bdf5`.
- **D19 — CRITICAL. Singular `block` was missing.** Volume groups accept `block` OR `blocks`, just
  as sidesets accept `sideset` OR `sidesets`. A valid deck would have failed the validator.
- **D20 — HIGH. `linear` / `linear magnetic` / `linear thermal` are three distinct sections with a
  fallback chain, not aliases of one.** And `linear` is required unless the field-specific section
  exists, because the fallback calls `section( "linear" )` unguarded.
- **D21 — HIGH. Five linear-solver keys were typed as free `string`; they are closed enums**
  resolved through `en_SolverEnums` converters, unknown value fatal. A GUI would have offered free
  text and produced decks the C++ rejects. Value lists still to transcribe (schema `audit_todo.T1`).
- **D22 — HIGH. `coupling` is a closed, case-sensitive enum**: `fully coupled` | `segregated`.
- **D23 — HIGH. `period` is REQUIRED for `ramp` and `sigmoid`** — no frequency fallback there.
- **D24 — MEDIUM. `anderson stabilization` fills depth 3 for magnetic but 1 for thermal**, not 3
  for both. Plus missing defaults: `adapt timestep` true, `maximum timestep` `REAL_MAX`,
  `coupling factor` 1.
- **D25 — MEDIUM. `initial conditions` fully specified now:** only `t` / `temp` / `temperature`
  exist, in K, setting `gTbulk`; any other key is silently ignored; absent section leaves NaN.
- **D26 — MEDIUM. `circuit` component headers must be bare types** (`resistor { }`), because the
  type comes from the section key — `resistor : r1 { }` does not resolve.
- **D27 — FALSE POSITIVE (retracted 2026-08-10).** Grok rated "`#` is not a comment token" as
  CRITICAL. **Refuted:** `clean_string` does `else if ( tChars[c] == tHash ) { break; }` and
  `remove_comments` calls it on every line. Grok read only `cl_InputFile.cpp`. Recorded in the
  schema's `known_false_positive` block so it is not re-raised.

**Codex's independent schema audit (same round, arrived later).** Confirmed Grok on `mesh.unit`,
the `linear` fallback model, the `scheme` enum, singular `block`, the homology and `coupling` enums,
`coupling factor` default 1, `initial conditions`, and `circuit` requiredness — **eight independent
confirmations**, all already applied. New findings:

- **D29 — HIGH. `a:b` ranges always expand ASCENDING.** `9:4` gives `4,5,6,7,8,9`, not a descending
  list — the reader branches on `tA < tB` and iterates upward either way. My schema's condensed
  phrasing ("inclusive range in either direction") could be read as preserving descent, which
  matters precisely because order is load-bearing everywhere else in this contract. Wording fixed.
- **D30 — MEDIUM. A bare number IS a legal `value`** (treated as already-SI and dimensionless); the
  unit is not globally mandatory. It *is* mandatory per-key in `layers`. My schema stated the strict
  rule globally. Fixed with `unit_mandatory:` marked per key.
- **D31 — MEDIUM. `materials.builtin` is `enum_or_pattern`**, not enum-plus-prose: alloy formulas
  go through a structured formula parser, so a validator needs both the list and the pattern.
- **D32 — HIGH, and it invalidated a fix I had just made.** I had added `crank-nicolson` and
  `galerkin` to the `scheme` enum as plain valid values. They parse into `EulerMethod` and then hit
  `BELFEM_ERROR( false, … )` in `cl_IWG_Timestep` — the Newton correction is exact for the BDF
  family only, and they are deliberately not aliased onto bdf1. Accepted ≠ usable.
- **D33 — HIGH. Amplitude dimension is context-dependent**, not a flat map: thermal `dirichlet` is
  K, while Maxwell `dirichlet` is dimensionless. Now keyed by consumer.
- **D34 — the best structural idea of the night: `runtime_status`.** The parser accepts strings the
  rest of the code then refuses, ignores, or mishandles, and a validator checking only "is this a
  legal value" passes decks that cannot run. Adopted as schema vocabulary:
  `live` / `accepted_then_hard_error` / `inert` / `silent_noop` / `corrupts`, with the policy
  *live passes, silent_noop warns, everything else errors*. It unifies findings that were scattered:
  CN+Galerkin (`accepted_then_hard_error`), `tape`/`shell`/`curve`/`cut` (`inert`), `curves` on a
  mesh that already has them (`silent_noop`), and `background dirichlet` (`corrupts`).
- **D35 — FALSE POSITIVE (retracted 2026-08-10).** Codex rated several anchors CRITICAL-unresolvable
  (`section_exists( "circuit" )`, `key_exists( "number of nodes" )`, …). They resolve: the source
  writes `section_exists("circuit")` without inner spaces and the checker normalises whitespace, as
  the schema header declares. Codex tested a stricter literal-grep rule than the schema states.
  **But the underlying style criticism was right and both auditors raised it** — the anchors should
  not depend on normalisation at all. All 81 call-form anchors were mechanically converted to bare
  quoted key literals; the 5 remaining expression anchors are genuine gates/functions.

**D28 — my own tooling was wrong, and this is the lesson of the round.** `check_anchors.py`
reported "all 53 anchors resolve." The file actually contains **115**; the script's regex matched
only block-style `anchor:` at line start and silently skipped every anchor written inline inside a
flow mapping. It passed a check it had never run. Grok independently predicted the one broken
anchor (`key_exists( "blocks" )` — the code reads that key through the variables `tGroup`/`tGroups`,
so the literal never appears). After fixing the parser: **116 of 117 resolve**, the one failure
being prose in the file's own header. *A checker that does not report its own coverage is
indistinguishable from a checker that passes everything.*

### Where the auditors disagreed

Only one substantive split: **O3 (splice vs canonical writer)** — Codex for splice, Grok against.
Recorded in O3 rather than resolved. Everything else was independent agreement, most importantly
D1 and the Q1 verdict that md-policy alone is insufficient.
