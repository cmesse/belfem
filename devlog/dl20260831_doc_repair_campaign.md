# Devlog 2026-08-31 — Documentation Repair Campaign (R1–R7)

**Date:** 2026-08-31
**Purpose:** Record of the seven-batch repair of the defects found by the documentation sweep,
run as plan+audit → correct → audit with Codex and Grok on every batch
**Module:** cross-cutting

---

## Summary

The sweep (`dl20260831_doc_currentness_sweep.md`) found 19 defect groups across the 96 markdown
files in `doc/` and `src/*/doc/`. This is the repair. Christian's instruction was to treat it like
code — plan, audit the plan, correct, audit the correction — and to batch it.

**89 files, 1295 insertions, 662 deletions. Seven batches, 16 auditor rounds** (a plan round plus
both auditors on each batch). `todo/doc_currentness_fixes.md` grew to 1411 lines as the working
record.

**The headline is not the repair. It is that every single batch shipped at least one new false
statement, and the audit caught all of them.** Two were only caught two batches later. If this had
been run as a straight edit pass — even a careful one — a substantial fraction of the corrections
would have been wrong in the same way the originals were.

## The tally, batch by batch

| Batch | What I got wrong |
|---|---|
| R1 | A false universal: `belfem` names its Exodus output after the mesh — true except on the segregated coupled path, which hardcodes `hphi_results.e-s` (`belfem.cpp:303`). I removed one false claim about the output filename and introduced another |
| R2 | A **mathematical error** — `dmudH` is dμ/dH, not ∂B/∂H; the differential permeability is μ + H·dμ/dH and is returned by nothing. Also taught an unsafe MPI pattern and missed a third valid pairing |
| R3 | **A regression.** Replacing the non-existent `create_local_mesh()` with `run()`+`partial_mesh()` on every rank turned a compile error into a **runtime abort on rank 0** — `partial_mesh()` is worker-only. The original at least failed loudly at build time |
| R4 | A flag that manufactured an open question the tree had **already answered** — `nedelec_derivation.md` settles the TRI6/TET10 circulation split explicitly, and my flag repeated the very conflation it was flagging |
| R5 | Four defects created by my own corrections, including a rename that left `propane.v(...)` calling a variable I had renamed to `methane` |
| R6 | A false repository-wide negative, caused by tooling — see below |
| R7 | A non-compiling example (`Objective` has a required ctor and a non-virtual `dimension()`), a descending line range, and an unchecked replacement claim |

**Two were found late.** R3's blocker surfaced in R6's audit, three batches after it shipped and
after R3's own audit had passed it. R6's tooling defect invalidated a class of claim I had been
making since R2.

## The tooling defect, because it outlives this campaign

`grep` in this environment is a **shell function wrapping `ugrep --ignore-files`**, which honours
`.gitignore`. This repository ignores `archive/`, `nonfree/`, `literature/` and `tmp/`. Every
`grep -r … .` from the repository root therefore returned a **false negative** for those trees.

Codex found it by locating `archive/graph/fn_Graph_tarjan.{hpp,cpp}` — a real Tarjan
articulation-point pocket detector — after I had asserted in a module document *and in a handoff to
the module owner* that no such implementation existed anywhere. A stored memory said it did. I
trusted the empty search over the note and did not investigate the discrepancy.

I re-ran every load-bearing negative from R2–R6 with `command grep`. **All of them held**; only the
Tarjan claim was wrong. The mechanism is now a memory
(`feedback_grep_ignores_gitignored`): use `command grep` for negatives, or search excluded trees by
explicit path.

## Two code defects, both approved mid-campaign

**`cl_Graph_Vertex.hpp:42` — `mOwner` initialised to the wrong sentinel.** Found while documenting
the graph vertex table; Christian ruled it a defect rather than a convention. `gNoID` is
`numeric_limits<id_t>::max()` over `unsigned int`, `gNoOwner` is `numeric_limits<proc_t>::max()`
over `int`; the narrowing made a fresh vertex's owner **−1**, so `owner() == gNoOwner` was false.
Confirmed by compiling and running the conversion.

A dedicated agent applied the one-line fix and traced the consequences, correcting two things I had
said:

- **My "−1 underrun" framing was wrong.** `Vector<T>::operator()` takes `size_t`, so the old value
  converted to 2⁶⁴−1 — always faults. The new sentinel gives 2³¹−1, *marginally less certain* to
  fault. The change is right, but it is a small regression in failure **loudness**.
- **My severity was overstated.** `mesh::Basis::Basis()` sets owner to 0
  (`cl_Mesh_Basis.cpp:25`), so **no mesh entity ever observes the default**. The blast radius is a
  handful of direct `graph::Vertex` users.

The strongest evidence that `gNoOwner` was always intended is not the naming but the algorithms:
the `std::min` ownership sweeps at `cl_Mesh_Partitioner.cpp:253-260` and `cl_FEM_Kernel.cpp:427`
only work if the sentinel behaves as **+∞**. Under the old value it was an absorbing element and
every facet would have ended unassigned. Several sites are **repaired** by the change.

**`tests/math/test_GraphVertex.cpp:27` pinned the defect** — `EXPECT_EQ( (id_t) tV.owner(), gNoID )`,
where the cast to `id_t` was the tell that the test matched the buggy initialiser rather than the
contract. Corrected, and the now-dead `fn_Graph_clear.hpp` include removed after checking the one
thing that could break it: the header also supplies `cl_Cell.hpp`, but `cl_Graph_Vertex.hpp:19`
supplies it too.

## The generator, and O5

R7e's new `src/numerics/opt/doc/README.md` **broke `update_doc_index.py --check`** — every collected
page needs a `MODULE_GROUPS` entry. Fixing it on approval surfaced a **second, pre-existing failure
the first was masking**: two `examples/*/data` pages from the peer session's usermat migration.

Fixing the generator also **resolved O5**. The stale `hphirun`/`hphiTrun` names lived in
`update_doc_index.py`, which *generates* `doxygen_nav.dox` and `groups.dox` — so editing the `.dox`
files (as R1 could not) would have been reverted on the next run. Fixed at the source and
regenerated. `--check` now exits 0, 0 files out of date, 27 modules, 101 pages.

## Scope

Documentation and comments only, and it is checkable rather than asserted: of the 15 non-markdown
files in the diff, each was filtered to its non-comment lines and **eight come back empty** — the two
example `matlib.cpp` files, `examples/scripts/Allrun`, `src/executables/CMakeLists.txt`,
`fn_Graph_METIS.hpp`, `UserMaterialTemplate.cmake`, `cl_Material.hpp` and
`cl_Material_UserDefined.hpp`. Those are comment-only by measurement.

Seven carry real changes. Six were explicitly approved: the `gNoOwner` fix, its test,
`update_doc_index.py`, the two `.dox` files that script generates, and `doc/input_schema.yaml` as the
input contract's machine-readable twin, updated in lockstep with its prose half. The seventh,
`.claude/ai_exchange_pos.txt`, is the exchange cursor the tooling maintains on its own.

**No file under the §7.1 cohomology-core ban was touched**, checked after every batch and again at
close.

## Non-documentation findings — handed off, not fixed

Four things the campaign found were outside its scope, and one is an error of mine. They are written
up with their evidence in **`todo/handoff_20260831_open_defects.md`**:

| ID | Finding |
|---|---|
| F1 | `graphtools.hpp:142-146` guards only the `comm_size()` marker, so a vertex holding `gNoOwner` indexes a `comm_size()`-length `Vector` at 2147483647 — a debug assert, an **out-of-bounds increment in release** |
| F2 | `parmetis_nd`, `ptscotch_nd` and `build_pargraph_adjacency` are **defined and called from nowhere** — verified across `src/`, `tests/`, `archive/`, `nonfree/`. Needs a wire/archive/delete ruling; F1 only matters if it is wired |
| F3 | The `mOwner` fix is right but made failure **quieter**, and two claims made during the campaign are corrected downward there |
| F4 | **My error.** I ticked D6 while `homology_usage_guide.md:1329-1350` and `:335` still taught `manifold_filter_3d()` and `check_surface()` — both **zero hits tree-wide** — as the live cleanup pipeline. Box unticked at close |
| F5 | The `grep` false-negative exposure in shipped docs, **measured**: five absolute-negative claims, the two load-bearing ones re-verified and holding. One wrong claim total, fixed |
| Q1 | Whether the absent AC-loss postprocessing (`loss`: **zero hits** under `src/fem/maxwell/`) is deliberate — Christian's ruling |

## Open Questions

- **O1** — whether `maxwell_usage_guide.md`'s static `∇·(μ∇φ) = 0` air form is stale or a
  deliberate modelling choice. Formulation; Christian's.
- **O2** — the static-condensation attribution. `messe2023.txt:509-517` says the paper implemented
  Lagrange multipliers and credits condensation to Alves et al. Christian's paper, Christian's call.
- **`src/homology/doc/handoff_for_gregory_20260831.md`** — twelve API-drift rows, three questions,
  for Gregory Giard. Deliberately untracked.
- `examples/Tape_Quench_obsolete_deleteme` still exists, named in no README.

## Files Updated

89 files. Full per-defect record in `todo/doc_currentness_fixes.md` (19 defect groups `D1`–`D19`,
one code defect `C1`, five open questions, seven executed batches).

## Status

**Reviewed, not verified**, with two exceptions. Nothing was built or run except two standalone
probes (the sentinel narrowing) and one real compile gate: `-fsyntax-only -std=gnu++17 -Wall
-Werror` on `test_GraphVertex.cpp` after the include removal, which passed clean.

**Owed:** a build, and `make check`. `check_doc_claims.py` held 37/37 throughout; citations went
from 564 resolving to 665.
