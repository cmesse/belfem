# Devlog 2026-08-11 — BFM Stale-Cache Detection (DR-21)

**Date:** 2026-08-11
**Topic:** A `.bfm` cached the enriched mesh but was reused on a checksum of the raw geometry
alone, so editing a layer thickness, a cut algorithm or the coating flag silently reloaded the
old mesh. Adds a second, orthogonal stamp.
**Module:** `src/fem/maxwell`, `src/mesh`
**AIs involved:** Claude (trace, plan, implementation), Codex + Grok (two blind jury rounds),
Christian (design ruling, scope, and the use case that settled the hash choice)
**Claude Confidence:** high on the mechanism (behaviour tested against three real decks); the
untested surface is the write path — no `.bfm` was actually written and reloaded in this session
**Verification:** compile — all three touched TUs pass `-fsyntax-only` under the project's own
`-Wall -Werror -pedantic-errors`; behaviour — standalone driver over `examples/{corc,helix,sidecoating}`
with six discrimination tests (§3). **Not built into the tree, not run end to end.**
**Plan:** `todo/bfm_stale_cache_detection.md` · **Register:** DR-21

## Summary

The reuse test was one equality:

```cpp
if ( tChecksum == tBfmFile.checksum() )   // raw .msh identity only
```

against a checksum that hashes dimensions, node coordinates and element node ids — while
everything the file actually caches (cuts, thin-shell layers, coating walls, periodicity,
hanging entities, RCM renumbering) is produced *afterwards* from `input.conf`. Edit a
thickness, rerun, and the edit never reaches the mesh. Nothing is printed.

Now a second stamp travels in the file's `meta` group and both must match.

## 1. Why the old checksum stays untouched

Christian's framing, which shaped the design: the two stamps answer independent questions.
The checksum **is** the `.msh` identity, so it catches a changed base mesh under an unchanged
deck — something no options tag can see. The new tag catches a changed deck under an unchanged
base mesh — something the checksum cannot see. Validity is the conjunction. It is also the
memdump guard, so widening it would have moved restart semantics.

## 2. What the tag is, and the two rejected alternatives

**Rejected — a checksum of the completed mesh.** Christian put it best: the point of the file
is to *avoid* rerunning cohomology and thin-shell generation, so a stamp computable only after
rerunning them has nothing left to save. It survives only as a possible future *integrity*
check, which is a different feature.

**Rejected — hashing the deck text.** Two earlier drafts (a whitelist of section text, then a
blacklist-filtered recipe) died on the same rock: `100 um` and `0.1 mm` are the same mesh but
different text, and a reindent is the same mesh too.

**Built — a tuple of numbers**, Christian's proposal: thin-shell count (there can be several),
layers per shell, thicknesses, coating on/off, plus the entries the jury proved also change the
cached mesh — `homology { algorithm }`, the terminal ids, domain types, periodic ids, curves.
Parsed, converted to SI, formatted `%.12g`, sorted, and hashed with **FNV-1a** rather than
`std::hash`, because a `.bfm` may be mailed to another machine and `std::hash` is unspecified
across standard libraries. The canonical text is stored beside the tag, so a mismatch names the
line that changed instead of asserting that two numbers differ.

`Mesh::save()` constructs its own `BfmFile`, so the writer could not be handed the tag through
a setter (the plan's D6). Resolved without restructuring the save path: the tag rides on the
`Mesh`, which `Mesh::save()` already carries into the file. The reuse decision happens *before*
`load()`, so `BfmFile` gained `config_tag()` / `config_text()` probes that open the `meta` group
standalone, mirroring the existing `checksum()`.

## 3. Running it found two defects that compiling could not

This is the part worth keeping. Both TUs compiled clean, and both of these were live:

1. **`input curves : 1,3,5,7,9,11` collapsed to `1`.** `looks_numeric` tested only the first
   character, so an id list went through `to_real`, which stops at the first separator.
   `...,11` and `...,13` would have produced the **same tag** — a silent miss in exactly the
   feature meant to prevent silent misses. Now the whole token must parse as a number, and
   anything carrying a separator stays a string compared verbatim.
2. **The terminal filter matched only `terminal`.** helix writes `input terminals`, but corc
   writes `input curves` for the same cut-defining ids, so corc's terminals were absent from
   the tag entirely. Both spellings now match.

Discrimination results on corc, each a separate process (so cross-run determinism is covered):

| edit | tag | wanted |
|---|---|---|
| rerun, no edit | unchanged | unchanged |
| `ybco 1.6 → 1.7 mum` | **changes** | changes |
| extra layer | **changes** | changes |
| terminal id `12 → 14` | **changes** | changes |
| `hastelloy 50 mum → 0.05 mm` | unchanged | **unchanged** (same mesh) |
| solver `mumps → strumpack` | unchanged | **unchanged** (must not rebuild) |

One test failure was my own: a first attempt reported "thickness change does not move the tag",
which was a `sed` that never matched — the deck writes `mum`, not `um`. Worth recording because
the alarming result was the test's fault, not the code's, and I nearly acted on it.

## 4. Scope, decided by Christian

In: everything that yields a different mesh configuration. Out: solver, timestepping,
tolerances, output, BC amplitudes — a user must stay free to retune a run against a cached
mesh. Also out by ruling: changes *inside* a user-defined material or boundary condition. Layer
material labels are in, so `copper` → `silver` rebuilds; a custom material gaining or losing
`rho` does not, and that gap is documented rather than chased.

The direct `.bfm` path (`mesh { file }` naming a `.bfm`) **warns and continues**. It has no
`.msh` by construction, and it is the "someone handed me a prepared mesh" path — a hard check
there would break the workflow it exists for.

## 5. Handoff

Not built into the tree; `make` and the end-to-end reproducer are Christian's: edit a thickness
in a directory that already holds a `.bfm`, rerun, confirm the rebuild message names the changed
line, then confirm a solver edit still hits the cache.

`doc/input_schema.yaml` was deliberately not touched — it is uncommitted work in progress — so
recording the tag's key set there is still owed (plan R8).

## 6. End-to-end verification (2026-08-11, after Christian's build)

`hphirun` builds in the real tree, so the write path — the one surface the standalone driver
could not reach — was finally testable.

- The first run wrote `corc.bfm` with both `meta/config` and `meta/config_text` present
  (`h5dump`), and the stored tag `11665186923197119937` is **identical** to the value the
  standalone builder computes for the same deck. The tag round-trips through HDF5 unchanged.
- Editing `ybco 1.6 → 1.7 mum` and rerunning produced exactly the intended behaviour:

  ```
  corc.bfm was built with a different mesh configuration and is being rebuilt
      now : layers:tape.layer[2].thickness = 1.7e-06
      was : layers:tape.layer[2].thickness = 1.6e-06
  Creating cuts ...
  ```

- Restoring the deck and rerunning reused the cache: no rebuild message, no cut recomputation,
  straight through to element connectivities.

Both outcomes, the probe, the comparison and the diff message are therefore confirmed against a
real deck. The remaining untested branch is a `.bfm` whose `config` dataset is absent — the
old-file case — which shares the code path and differs only in wording.

## 7. Implementation audit (2026-08-11) — one regression and one real hole

Christian asked for a jury round on the committed change. It was worth it: two findings were
things my own end-to-end testing had passed straight over.

**A regression I introduced and did not see.** `BfmFile::load()` restores the checksum but
never restored the new stamp, so a mesh loaded from a `.bfm` and saved again wrote an
**untagged** file — which the reader then treats as a permanent mismatch, rebuilding on every
run and never becoming valid. I had pre-registered the write path as my top risk but guessed
the wrong mechanism (un-enumerated `Mesh::save()` callers); the hole was the *load* path
leaving the mesh in that state. Fixed.

**A stale-cache hole of exactly the species I predicted, and still missed.**
`examples/circuit/input.conf` declares its cut-defining ids as `input curves : [1,2]` under
`circuit/topology/terminal pair` — **not** under `boundary conditions`, which is the only tree
my walker visited. Those pairs become CircuitVoltage BC domains and reach
`CutFactory::set_terminals`, so editing a circuit curve id changed the cuts while the tag stood
still. Confirmed by execution both ways: before the fix the edit left the tag unchanged, after
it the tag moves. One recursive collector now serves the BC tree and the circuit tree, which
also closes the nested `boundary conditions/maxwell/...` form the auditors flagged separately.

**Found by me while fixing those:** two unlabelled `terminal pair` sections shared a path, so
swapping ids *between* pairs left the line set unchanged. A section index now separates
same-type siblings — order sensitivity being the safe trade, since reordering costs a rebuild
while a collision costs a stale mesh.

Also fixed: boolean spellings (`on`/`true`/`yes`) produced different tags for one flag; `a`
prefixes on locals, which belong to arguments; `Mesh::memory()` did not count the new string;
and `bfm_file_format.md` documented neither the new datasets nor the provenance stamps it had
already been missing.

Accepted rather than fixed, with the auditors' blast-radius analysis on record: material
definitions stay untracked (Christian's ruling), `1:3` versus `1,2,3` tags differently (false
rebuild only), and the probes still abort on a corrupt `meta` group exactly as the pre-existing
checksum probe does. MPI was explicitly cleared by both auditors.

**Note for the first run after this:** the tag values changed (the section index altered the
canonical text), so any `.bfm` written by the previous commit rebuilds once.
