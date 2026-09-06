# BFM Stale-Cache Detection: Stamp the Processing Options, Not Just the Geometry

**Date:** 2026-08-10
**Purpose:** A `.bfm` file caches an *enriched* mesh — cuts, thin-shell layers, edge-coating
walls, periodicity, hanging entities. The reuse test compares only a checksum of the **raw**
geometry, so changing an input option that alters the enrichment (a layer thickness, a coating
width, a sideset's domain type) reloads the stale cache **silently** and the new options are
ignored. Fix: fold a small tuple of the mesh-defining numbers ( thin-shell count, layers per
shell, thicknesses, coating flag, cut algorithm ) into a second tag stored beside the existing
checksum, and require both to match — warn and rebuild when they do not.
**Module:** `src/mesh` (BfmFile), `src/fem/maxwell` (MaxwellFactory), `src/io` (Section)
**AIs involved:** Claude (exploration + plan), Codex + Grok (two blind jury rounds), Christian (design ruling)
**Status:** **CLOSED 2026-08-11** (currentness sweep) — R1–R8 all done and the verification gate
has run. What verified it: after Christian's build, `hphirun` exercised the write path end to end
on a real deck — the first run stored `meta/config` + `meta/config_text` with a tag identical to
the standalone builder's value for the same deck (round-trips through HDF5 unchanged), a
1.6 → 1.7 µm ybco thickness edit rebuilt and *named the changed value* before recomputing cuts,
and restoring the deck reused the cache silently (commit `2c0a7382`). The old-file branch — a
`.bfm` carrying `meta/checksum` and no `meta/config` — was then closed too against a genuinely
pre-feature file: one rebuild, one message, and a stamped file afterwards (§6 DoD).
R8's schema half, outstanding when the plan was last
touched, landed in a parallel session and is verified present in this sweep (§4, R8).
Register row `DR-21` is struck.
*(Superseded status, kept for the record: IMPLEMENTED 2026-08-11 on Christian's go-ahead — R1-R7
done, R8 partial.)* Compiles clean under the project's own `-Wall -Werror -pedantic-errors`; behaviour
verified against three real decks (§7). Audited by two blind jury rounds; design ruled by Christian (tuple of numbers, not
a recipe text), O4 resolved (defined byte algorithm over fixed-precision text), O5 closed as
out of scope, R7 reframed (warn, never block). Threads: `tmp/ai_exchange/review_bfm_stale_cache.md`
(design audit), `tmp/ai_exchange/review_bfm_cache_method.md` (method choice).
**Register row:** `debt_register.md` DR-21 (P1, blocking 1.0)

> **Scope guards**
> - **Do not widen `Mesh::compute_checksum`.** It is also the memdump guard
>   (`cl_Mesh.cpp:3014`, "Memdump checksum mismatch"); changing what it covers changes restart
>   semantics. The new stamp lives *beside* it.
> - **Do not invalidate on unrelated input changes.** A solver tolerance or a timestep setting
>   must not force a cohomology recomputation — that is the expensive thing this cache exists
>   to avoid.
> - No change to the `.bfm` payload format beyond adding one dataset to the existing `meta`
>   group. Old files must still load.
> - This adds **no new `input.conf` key**, so the Input Contract (CLAUDE.md) is not triggered.

---

## 1. Current behaviour and how it fails

`MaxwellFactory` decides to reuse a cached mesh here (`cl_MaxwellFactory.cpp:279-306`):

```cpp
size_t tChecksum = aMesh->checksum();          // raw .msh only
if ( std::filesystem::exists( tBfmFilePath ) )
{
    mesh::BfmFile tBfmFile( tBfmFilePath );
    if ( tChecksum == tBfmFile.checksum() )    // <-- the whole test
    {
        delete aMesh ; tBfmFile.load(); mComputeCohomologies = false ;
        return tBfmFile.get();
    }
}
```

and `Mesh::compute_checksum` (`cl_Mesh.cpp:2379-2417`) hashes exactly: number of dimensions,
number of nodes, every node coordinate, number of elements, and every element's node ids.

**Nothing about processing.** But between that test and the `mMesh->save( mMeshPath )` that
writes the file (`cl_MaxwellFactory.cpp:513`), the mesh is transformed by `create_cuts()`,
`create_edges_and_faces_on_mesh()`, `create_thinshells()`, periodicity update,
`create_hanging_edges_and_facets()` and an RCM renumbering — all driven by `input.conf`.

**Reproducer (the register's):** run a deck, edit a layer thickness or an `edge coating width`,
rerun. The `.msh` is untouched, so the checksum matches, so the old `.bfm` loads and the edited
value never reaches the mesh. No warning is printed. The run appears to succeed.

**Partial mitigation already in tree, and why it is not enough.** `save_meta_data` writes a git
provenance stamp with the comment *"the checksum only covers the base mesh, so this stamp is
what identifies a file written by an older construction state"* (`cl_Mesh_BfmFile.cpp:177-184`).
That catches *code* drift between writer and reader. It cannot catch an **input** change made
with the same binary, which is the common case.

## 2. Architecture: two orthogonal stamps, and a tuple of numbers

**Why the existing checksum stays exactly as it is** (Christian, 2026-08-10). The two stamps
answer independent questions and neither implies the other:

| Stamp | Question | Blind to |
|---|---|---|
| `checksum` (unchanged, **is** the `.msh` checksum) | is this `.bfm` built from the same **base mesh**? | every processing option |
| config tag (new) | is it built from the same **derivation**? | any change to the base mesh |

Cache validity is the conjunction. A config tag alone would reload a mesh derived from a
*different geometry* whenever `input.conf` happened to be unchanged; the checksum alone is the
present bug. The checksum must also keep equalling the `.msh` value, because that identity is
what makes the positive reuse decision legitimate.

**The ordering constraint that fixes the shape of the fix** *(corrected after audit — the
first version of this paragraph overstated it)*. `read_mesh()` — which contains the reuse test
— runs at `cl_MaxwellFactory.cpp:69`. `mInputFile` is constructed at `:60`, so the **`Section`
tree is already parsed** at decision time; what is *not* available is the mesh **binding**,
because `read_domain_types()` (`:110`) dereferences `mMesh->block( tID )` (`:338-352`). The
correct constraint is therefore: **the tag must be derivable from the input tree without
mesh binding.** A structured extractor is fine; hashing the parsed `Protoshell` / domain-type
objects is not, since those exist only after binding.

**One consequence for the extractor:** it cannot be built on the key API alone.
`read_thin_shell_data` parses the layer stack from raw buffer lines with word splitting
(`cl_MaxwellFactory.cpp:2601-2619`), not via `num_keys()` / `key( i )`, so a key-only walker
would silently omit exactly the thicknesses and materials this plan exists to catch.

**Two properties any version of this tag must hold**, carried over from the superseded drafts:

1. **Semantic** — values parsed and unit-converted to SI, so `100 um` and `0.1 mm` compare
   equal and formatting costs nothing.
2. **Deterministic** — folded in a fixed order, never `Map` iteration order. A tag that varies
   run to run rebuilds on every reload; that is the one defect that would make the feature
   worse than useless.

*(Two earlier designs are recorded in §8: a whitelist of section **text**, and jury round 2's
blacklist-filtered recipe with the text stored for per-setting diffs. Both are superseded by
the decision below, though the stored canonical string in O4 keeps a readable mismatch
message.)*

**DESIGN DECIDED BY CHRISTIAN, 2026-08-10: hash a small tuple of numbers, not a recipe text.**

> *"Maybe just check the numbers? Same amount of thin shell configurations (there can be more
> than one!) same amount of layers per shell, same thicknesses, side coatings on or off.
> That's just a handful of numbers to be checked from the input file."*

This supersedes both the first draft (whitelist of section **text**) and jury round 2's
recommendation (blacklist-filtered recipe with stored text). What it buys:

- **No unit normalization problem.** Values are parsed to SI and compared as numbers, so
  `100 um` and `0.1 mm` are equal by construction — the text designs had to work for this.
- **No stored deck.** A tag beside the existing checksum, plus a couple of hundred bytes of
  canonical string for the mismatch message (O4).
- **Transparent.** The whole check is one function anyone can read in a minute.
- **Computable at decision time**, which is the binding constraint (see the ordering note
  below): every value comes from the input tree with no mesh binding.

**The tuple** — Christian's list, plus two entries the jury verified do change the cached mesh:

| Contributor | Why it is in |
|---|---|
| number of `thinshell` sections | Christian: there can be more than one |
| per shell: layer count | sets how many layer blocks are built |
| per shell: thicknesses, in SI | layer geometry |
| per shell: material labels | **not cosmetic** — `hasDuplicates = ( tA != tB ) && tA->have( rho ) && tB->have( rho )` (`cl_ThinShellFactory.cpp:132-143`); a material swap flips inter-layer node duplication with `layers` otherwise unchanged |
| per shell: sideset ids | which surfaces become shells |
| side coating on/off, and its width | Christian |
| `homology { algorithm }` | read only on the recompute path (`cl_MaxwellFactory.cpp:789-805`), passed to `CutFactory` (`:837-838`), which switches behaviour on it (`cl_CutFactory.cpp:271-344`) — cuts are baked into the `.bfm` |
| terminal domain ids | feed the suggested homology (`cl_MaxwellFactory.cpp:1027-1108`, `:852-855`; `cl_CutFactory.cpp:221-224`, `:245-248`) |
| periodic source / target ids | node pairing and duplication |

**Accepted trade-off, recorded so it is a choice and not an accident.** This is a whitelist,
and whitelists fail *unsafe*: a mesh-affecting key added later is silently omitted until
someone extends the list, and the symptom is a stale mesh rather than a wasted rebuild. Jury
round 2 argued for a blacklist on exactly that asymmetry, and this plan's own first draft
leaked two entries within an hour of being written. **Christian's ruling stands** — the list is
small, lives in one function, and covers everything found to change the mesh today. Mitigation:
keep the tuple builder in one place and reference it from `doc/input_schema.yaml`, so the
code→schema check has something to compare against.

**Also in the tuple, from tracing the enrichment** (still "just numbers" — ids, flags, counts):

| Contributor | Why it is in |
|---|---|
| domain type per block / sideset id | decides Conductor / Air / Ferro / ThinShell, which drives both cuts and enrichment |
| the side-connector flag | whether connector blocks exist at all |
| `topology`: `curves` ids | curve entities are part of the stored topology |

**Explicitly NOT in the tuple** — the user must stay free to change these against a cached
mesh: solver library and settings, tolerances, timestepping, boundary-condition **amplitudes**
(but *not* terminal domains, which feed cut construction), and output options.

**Note the one that looks safe to exclude and is not:** material *property* values. An earlier
draft excluded them. `ThinShellFactory` decides topology from material identity
(`cl_ThinShellFactory.cpp:132-143`), so the layer material **labels** are in the tuple. Whether
a change *inside* a material definition — one that adds or removes `rho` — must also invalidate
is O5 below.

**Deliberately out of scope:** `mesh { unit }` is already covered, because `scale_mesh` runs
*before* `aMesh->checksum()` (`cl_MaxwellFactory.cpp:277-280`) — a unit change moves every node
coordinate and the geometric checksum already catches it. Worth stating so nobody adds it twice.

**Open question O1** covers element order and enrichment, which are borderline.

## 3. Gap table

| ID | Gap | Evidence |
|---|---|---|
| D1 | Reuse test ignores every processing option | `cl_MaxwellFactory.cpp:296` |
| D2 | No option stamp exists in the file to compare against | `cl_Mesh_BfmFile.cpp:165-187` — `meta` holds dimensions, entities, groups, checksum, belfem, git, branch |
| D3 | A stale hit is silent — no message on the reuse path | `cl_MaxwellFactory.cpp:296-305` returns without printing |
| D4 | Old `.bfm` files carry no stamp, so "absent" must mean "rebuild", not "error" | pattern already used for the provenance stamps, `cl_Mesh_BfmFile.cpp:200-208` |
| D5 | **The direct `.bfm` path performs NO validation at all** — if `mesh { file }` names a `.bfm`, `read_mesh()` loads and returns before any checksum exists to compare | `cl_MaxwellFactory.cpp:242-253` (found by Codex, jury round 2; **pre-existing, not introduced here**) |
| D6 | `Mesh::save()` builds its own `BfmFile` internally, so a stamp cannot be handed to the writer by a setter on `BfmFile` | `cl_Mesh.cpp:353-356` (jury round 1) |
| D7 | `load_meta_data()` runs *inside* `load()`, but the decision must precede `load()` | `cl_Mesh_BfmFile.cpp:101-121` vs the gate at `cl_MaxwellFactory.cpp:291-296` (jury round 1) |

## 4. Ordered steps

- [x] **R1 (**DONE** — `src/fem/maxwell/fn_mesh_config_tag.hpp`; behaviour verified against three real decks, see §7) — the tuple builder.** One free function (proposed
  `src/fem/maxwell/fn_mesh_config_tag.hpp`) taking the `InputFile` and returning the canonical
  **string**: the §2 tuple in fixed order, each value formatted `%.12g` (O4). A second tiny
  function FNV-1a's those bytes into the tag. **Not `Hash`/`std::hash`** — see O4. Reads the input tree only — **no mesh binding**, which is
  what makes it callable at the reuse decision. Note it must read the layer stack from the raw
  section buffer, not `key( i )`: `read_thin_shell_data` parses layers by word-splitting
  buffer lines (`cl_MaxwellFactory.cpp:2601-2619`), so a key-API walker would silently omit the
  thicknesses. Iterate shells in a deterministic order (sort by label) — never `Map` order, or
  the tag varies run to run and every reload rebuilds.
- [x] **R2 (**DONE** — `meta/config` + `meta/config_text` in `BfmFile::save_meta_data`. D6 solved without restructuring the save path: the tag rides on the `Mesh` (`set_config_tag`), which `Mesh::save()` already carries into its own `BfmFile`) — write both.** In `BfmFile::save_meta_data`, add `meta/config` (the tag) and
  `meta/config_text` (the canonical string, a couple of hundred bytes) beside the existing
  `checksum`. **D6 blocks the obvious route:** `Mesh::save()` constructs its own `BfmFile`
  (`cl_Mesh.cpp:353-356`), so a setter on `BfmFile` is unreachable from
  `cl_MaxwellFactory.cpp:513`. Either save through a `BfmFile` the factory owns, or thread the
  value through `Mesh::save`. This is the one structural change; decide it before coding.
- [x] **R3 (**DONE** — `BfmFile::config_tag()` / `config_text()`, standalone opens mirroring `checksum()`, `dataset_exists`-guarded) — read it, *before* `load()`.** D7: `load_meta_data()` runs inside `load()`, too
  late for the decision. Add a lightweight probe mirroring `checksum()`
  (`cl_Mesh_BfmFile.cpp:44-56`) that opens `meta` and returns the value without loading the
  payload, guarded by `hdf5::dataset_exists` (D4).
- [x] **R4 (**DONE** — both stamps required at `cl_MaxwellFactory.cpp`; absent tag = mismatch with its own message; `report_config_difference` prints the changed lines, capped at 8) — extend the reuse test.** At `cl_MaxwellFactory.cpp:296`, require **both** the
  checksum and the config tag to match. Absent tag (every existing `.bfm`) counts as a
  mismatch. On mismatch, print a Default-level line naming the file and saying the mesh
  configuration changed, then fall through to the existing recompute path, which already does
  the right thing. Because the canonical string is stored (O4), the message can print the
  stored and current strings so the changed value is visible, rather than only saying the tag
  differs.
- [x] **R5 (**DONE** — Verbose line on the reuse hit) — say something on the *hit* too (D3).** One Verbose-level line confirming the cache
  was reused, so "did it use the cache?" stops being a guess.
- [x] **R6 (**DONE by construction** — `read_mesh()` already runs on rank 0 only and the decision travels in the existing `mComputeCohomologies` broadcast; no new communication added) — MPI placement.** Compute and compare on rank 0, then broadcast the resulting cache
  decision exactly as `mComputeCohomologies` is already broadcast
  (`cl_MaxwellFactory.cpp:65-104`). Never let ranks decide independently.
- [x] **R7 (**DONE** — direct `.bfm` path warns and continues) — the direct `.bfm` path (D5): warn, never block.** *Reframed 2026-08-11 — this is
  not the hole it looked like.* When the deck names a `.bfm` directly, `read_mesh()` loads and
  returns at `cl_MaxwellFactory.cpp:242-253` with no validation, and the reason is structural:
  there is no `.msh` to validate against. **That is the "someone deliberately named a prepared mesh" path.**
  There is no `.msh` *in that code path* by construction — the deck names the `.bfm` — so the
  geometric checksum cannot be formed there regardless of what happens to sit on disk. A user
  who asked for that file explicitly should get it. So: compare the config tag if the file
  carries one, emit a **Default-level warning** naming the mismatch, and **continue**. Never
  abort, never rebuild on this path. *(Note that in Christian's transfer scenario the recipient
  also receives the `.msh` and deck, so they would normally point at the `.msh` and take the
  ordinary sidecar path, where the full check applies.)*
- [x] **R8 (**DONE** — `doc/input_file_reference.md` §3 documents which settings invalidate a cache;
  the schema half landed in a parallel session on 2026-08-11 and is verified present in this sweep:
  `doc/input_schema.yaml` carries the `mesh_config_tag` vocabulary — `wholesale` / `included` /
  `filtered` / `excluded` — the load-bearing rule that a new enrichment key outside a `wholesale`
  section must also be added to `src/fem/maxwell/fn_mesh_config_tag.hpp`, and the recorded
  accepted gap on material *bodies*. The file is on disk and still untracked, so it is not yet in
  history) — docs.** `doc/input_file_reference.md` gains a note naming which settings invalidate
  a `.bfm`; `doc/input_schema.yaml` records which keys feed the tag (the single source of truth the R1
  tuple mirrors); the mesh module doc gains the stamp's meaning. Codex prose pass on the touched
  sections.

### 4.0 Implementation progress

Not started — this file is the plan.

## 5. Open design questions

- [ ] **O1 — element order and enrichment in the stamp?** `mElementOrder` is taken from the
  mesh itself (`mMesh->max_element_order()`, `cl_MaxwellFactory.cpp:86`), not from input, so it
  is implied by the geometry. `mUseEnrichment` already hard-errors rather than saving
  (`cl_MaxwellFactory.cpp:510`). **Claude's reading: neither needs to be in the stamp.** Worth
  one auditor's eye, because if the order can ever be raised from input the stamp must cover it.
- [x] ~~**O2 — normalization strength.** Should a changed *comment* invalidate the cache?~~
  **CLOSED by jury round 1: the question is moot.** `InputFile` calls `remove_comments()` and
  `tidy_up()` before the section tree is built (`cl_InputFile.cpp:18-28`, `:55-68`, `:73-140`),
  so comments and blank lines never reach anything the recipe can see.
- [x] ~~**O4 — which hash for a stamp read by a later process?**~~ **RESOLVED 2026-08-11: use a
  defined byte algorithm (FNV-1a 64) over fixed-precision text — NOT `std::hash`.**

  *Rationale corrected the same day, after Christian clarified the transfer case.* The first
  version of this entry argued that a false mismatch on another machine is a dead end because
  the recipient has nothing to rebuild from. **Christian: "I assume that the user will also get
  an input file and a gmsh file."** With the `.msh` and deck in hand a false mismatch costs a
  rebuild, not a failure, so that argument is void and is withdrawn.

  The decision stands on weaker but sufficient grounds: the defined algorithm costs ~20
  auditable lines, guarantees the **new** check adds no portability risk of its own, and avoids
  a check that can spuriously fire — which teaches users to ignore warnings, costing more than
  it saves. Priority order given was stability, then user friendliness, then maintainability.

  **The serialization matters more than the hash function.** Do not hash raw IEEE bytes: that
  inherits endianness and last-ulp noise, so `100 um` and `0.1 mm` could tag differently for
  the same physical mesh through different unit-conversion paths. Format every value with
  `%.12g`, concatenate in fixed order, FNV-1a the bytes. Twelve significant digits absorbs
  conversion noise (1 ulp at 1e-4 is ~1e-16 relative) while still separating any thickness a
  user would type, and text has no endianness.

  **Store the short canonical string too.** Rejected earlier when it meant storing a deck; with
  a tuple it is a couple of hundred bytes and buys the user-friendliness ranked second — print
  both strings on mismatch and the difference is visible at a glance.

  > **Pre-existing property of the transfer workflow, found while resolving this and NOT in
  > scope to change.** `Mesh::compute_checksum` hashes every node coordinate as a `double`
  > through `Hash`, and libstdc++ implements `std::hash<double>` as `_Hash_impl::hash( __val )`
  > — MurmurHash over the raw bytes (`/usr/include/c++/11/bits/functional_hash.h:244-252`);
  > libc++ uses a different scalar hash. **So the existing base checksum is already
  > library-dependent:** a `.bfm` written by a GCC/libstdc++ build and read by a Clang/libc++
  > build fails its base checksum today, independently of anything in this plan. It is
  > invisible in practice — both ends run the same SCLS GCC stack, and a recipient holding the
  > `.msh` merely rebuilds and notices only that it was slow. Do **not** re-base it: it is the
  > memdump guard, and changing it would invalidate every existing `.bfm`. Recorded so the next
  > person to hit a "why did it rebuild on my colleague's machine" question finds the answer.

- [x] ~~**O5 — does a change *inside* a material definition need to invalidate?**~~ **CLOSED as
  out of scope (Christian, 2026-08-11):** *"changing a user defined material or boundary
  condition on the code side feels also out of scope of what we want to check."* The tuple keeps
  the layer material **labels**, so `copper` → `silver` is still caught; editing a custom
  material so it gains or loses `rho` is knowingly **not** caught, and the plan documents that
  rather than chasing it. *Claude's reading of the ruling, flagged for correction: "boundary
  condition on the code side" is taken to mean BC **implementations**, not the terminal ids in
  `input.conf`, which stay in the tuple because they feed cut construction.*
- [x] ~~**O3 — is the hash worth making portable?**~~ **RESOLVED before the audit (2026-08-10):
  no new assumption is introduced, so reuse `Hash`.** `Hash::operator+=` folds
  `std::hash<T>{}( aValue )` (`cl_Hash.hpp:64-70`), which libstdc++ implements without a
  per-process seed — stable run to run, not guaranteed across standard-library
  implementations. That is acceptable here because **`Mesh::compute_checksum` already hashes
  `double` node coordinates with the same utility and already compares the result across runs**
  (`cl_Mesh.cpp:2379-2417` vs the reuse test). The new stamp is therefore exactly as portable
  as the test it extends. Worst case after a toolchain change: the stamp mismatches, the user
  gets the warning, and the mesh is rebuilt once — the same safe outcome as an old file with no
  stamp.

## 6. Definition of done

- [x] Editing a layer thickness and rerunning **rebuilds** the mesh and says why. *(2026-08-11,
  corc: message named `layers:tape.layer[2].thickness = 1.7e-06` was `1.6e-06`, then
  "Creating cuts ...".)*
- [x] Editing a solver setting and rerunning **reuses** the cache. *(Tag-level: `mumps` →
  `strumpack` leaves the tag identical, and an identical tag is what the gate compares —
  the same condition under which the restored-deck run was observed to reuse.)*
- [x] An existing `.bfm` with no stamp triggers exactly one rebuild, with a message.
  *(2026-08-11, using a genuinely pre-feature `.bfm` — `h5dump` confirms it has `meta/checksum`
  and no `meta/config`: "corc.bfm predates the mesh-configuration check and is being rebuilt",
  followed by a real rebuild. It then writes a stamped file, so the next run reuses.)*
- [x] Reindenting a `topology` section does **not** rebuild. *(2026-08-11: 93 changed lines of
  whitespace plus an added comment — tag unchanged. `InputFile` strips comments and blank lines
  before the section tree, and values are canonicalized, so layout cannot leak in.)*
- [x] `debt_register.md` DR-21 updated; devlogs written
  (`dl20260811_bfm_stale_cache_detection.md`, §6 end-to-end and §7 audit).

**All five met. The only unticked step in this plan is R8's schema half** — see §4.

## 7. Verification (2026-08-11)

**Compile:** `cl_MaxwellFactory.cpp`, `cl_Mesh.cpp` and `cl_Mesh_BfmFile.cpp` all pass
`-fsyntax-only` with the module's real flags, including `-Wall -Werror -pedantic-errors`.

**Behaviour:** the tag builder was linked into a standalone driver and run against
`examples/corc`, `examples/helix` and `examples/sidecoating` — three structurally different
decks (layer stack + curves; conductor blocks + terminals; edge coating). Discrimination
tests on corc, each a separate process, so cross-run determinism is covered:

| edit | tag | expected |
|---|---|---|
| none (rerun) | unchanged | unchanged — deterministic |
| `ybco 1.6 -> 1.7 mum` | **changes** | changes — the headline case |
| extra layer added | **changes** | changes |
| terminal id `12 -> 14` | **changes** | changes — cuts depend on it |
| `hastelloy 50 mum -> 0.05 mm` | unchanged | **unchanged** — same mesh, different unit |
| `solver library mumps -> strumpack` | unchanged | **unchanged** — must not force a rebuild |

**Two real defects were found by running it, neither of which compiling could have caught:**

1. `input curves : 1,3,5,7,9,11` collapsed to `1`. `looks_numeric` tested only the first
   character, so an id list went through `to_real`, which stops at the first separator —
   `...,11` and `...,13` would have produced the SAME tag. Now the whole token must parse as a
   number, and anything with a separator stays a string.
2. The terminal filter matched only `terminal`, but corc writes `input curves` /
   `output curves` for the same thing. corc's cut-defining ids were absent from the tag
   entirely. Both spellings are now matched.

**END-TO-END, 2026-08-11 (after Christian confirmed `hphirun` builds).** Run in a scratch dir
holding only `corc.msh` + `input.conf`:

1. First run wrote `corc.bfm` carrying both `meta/config` and `meta/config_text` (`h5dump`),
   and the stored tag `11665186923197119937` is **identical** to the value the standalone
   builder computes for that deck — so the tag round-trips through HDF5 unchanged.
2. `ybco 1.6 -> 1.7 mum`, rerun:
   ```
   corc.bfm was built with a different mesh configuration and is being rebuilt
       now : layers:tape.layer[2].thickness = 1.7e-06
       was : layers:tape.layer[2].thickness = 1.6e-06
   Creating cuts ...
   ```
   The message names the exact changed value and the mesh is genuinely rebuilt.
3. Deck restored, rerun: **no** rebuild message, **no** "Creating cuts", straight through to
   element connectivities — the cache is reused.

That exercises the write path, the pre-`load()` probe, the comparison, the diff and both
outcomes. **Not yet exercised:** a file whose `config` dataset is absent (the old-file branch);
its message differs but the code path is the same `== 0` test.

## 8. Audit trail

- 2026-08-10 — Claude: traced the reuse path, `compute_checksum`, the enrichment sequence and
  the `layers`/`topology` parsing; wrote the first draft (whitelist of mesh-affecting sections,
  hashed as text).
- 2026-08-10 — **jury round 1** (Codex + Grok, blind, on the draft plan;
  `review_bfm_stale_cache.md`). Verdict: do not approve as written. Both found
  `homology { algorithm }` missing from scope; Codex additionally found the terminal domains
  under `boundary conditions`, plus D6 (`Mesh::save` owns its `BfmFile`) and D7
  (`load_meta_data` runs too late). O2 closed as moot.
- 2026-08-10 — Christian: keep the geometric checksum unchanged (it is the `.msh` identity and
  the memdump guard); add a second tag for the layer/cut configuration; solver settings must
  stay freely changeable. Asked the jury to choose between an enriched-mesh checksum, a
  line-by-line deck backup, and the recipe.
- 2026-08-10 — **jury round 2** (Codex + Grok, blind, on the method question;
  `review_bfm_cache_method.md`). Both independently: **Method C persisted, built with a
  blacklist, recipe text stored**; Method A cannot decide validity but is a legitimate
  post-load integrity guard. Codex additionally found **F6** (`materials` must not be
  blacklisted — it changes thin-shell topology) and **F12/D5** (the direct `.bfm` path has no
  validation at all), and corrected this plan's ordering claim (the input tree *is* parsed; only
  mesh binding is unavailable). All verified at source by Claude before folding in.
- **Reversal recorded honestly:** the whitelist was Claude's own recommendation and is now
  withdrawn. Its first draft leaked two mesh-affecting settings within the hour, which is the
  clearest available evidence about how whitelists age.
