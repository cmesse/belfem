# db2exo: Database → Exodus Converter Audited and Repaired

**Date:** 2026-08-28
**Purpose:** Audit the new `db2exo` tool, fix what it got wrong, then run a jury round on the repair and fix what that found.
**Module:** physics/database ( new executable ), io ( `hdf5_tools.hpp` ), config ( install list )

## What the tool is

`db2exo <hdf5file>` opens a BELFEM material database, walks its top-level HDF5
groups, and for every group carrying the full record written by
`Database::save` — `order`, `origin`, `points`, `step`, `values`
( `cl_Database.cpp:136-140` ) — reconstructs the tensor mesh and attaches the
table as a node field, then writes one `.exo` for ParaView.

## The defect that mattered

The draft's grid-compatibility test could never pass:

```cpp
tPoints -= tMesh->tensorconf()->num_nodes() ;   // total node count, not the per-dim vector
```

`num_nodes()` takes a defaulted `aDimension = BELFEM_UINT_MAX` and falls through
to `return mNumNodes` ( `cl_TensorMeshConfig.hpp:101-116` ); the per-dimension
vector is `num_nodes_vector()` ( `:60-63` ). So a scalar of order 1e5 was
subtracted from each entry of a `Vector< index_t >`. `index_t` is unsigned, so
on `sp-ap_new.hdf5` ( `points = {23,31,181}` ) every entry wrapped to ~4.29e9,
and the following `for ( real x : tPoints )` widened the wrapped value into a
`real` where `std::abs( x ) > BELFEM_EPSILON` was trivially true.

Every table after the first was therefore dropped, silently. `sp-ap` and `sst-1`
both hold `jc` **and** `n` on an identical grid, so they converted with `jc`
only. Christian hit exactly that symptom — `n` missing from the output — which
is what puts this finding on the reproducer rung rather than source trace.

Five more in the same file: length-mismatched `operator-=` aborting instead of
skipping on a mixed-dimension file; no guard that `values` length equals the node
count, where the only downstream check is a `BELFEM_ASSERT` that compiles out in
release ( `cl_Mesh_ExodusWriter.cpp:572-574` ); `data( 1 )` taken blindly even
though `Arguments` leaves `-v` in the list; the extension stripped from the whole
path rather than the basename; and a silent zero-exit when no table was found,
which is what kept the missing-`n` bug invisible.

## The jury round

Pre-registered in `tmp/ai_exchange/review_db2exo.md` before dispatch, including
a Part B of things checked and found correct ( node ordering, finalize, CMake
wiring, the `get_groups` group filter ) so the auditors would not re-litigate
them, and a Part C self-review of the rewrite.

**Codex** raised six, all confirmed. Two found real bugs in the repair, not the
draft: an input already named `.exo` is destroyed by its own conversion
( `ex_create( ..., EX_CLOBBER, ... )`, `cl_Mesh_ExodusWriter.cpp:158` ) and — the
half Part C missed — a hidden file `.database` collapsed to `.exo`. The sharpest
was `order == 0`: `TensorMeshConfig` evaluates `( points - 1 ) % order` while
initializing `mNumElementsPerDim`, which is **declared before** `mElementType`
( `cl_TensorMeshConfig.hpp:29` vs `:39` ), so the clean "Invalid element order"
error at `:87` is unreachable — you get SIGFPE first.

**Grok** was blocked for thirteen minutes by a host fault, not an API one:
`sandbox profile resolve failed: could not resolve runtime-socket deny path
/run/podman/podman.sock: Permission denied`. It recovered on retry, read the
**live** tree rather than the dispatch snapshot, and so audited the fixes
themselves. That produced the round's best finding, which neither Codex nor the
pre-registration had:

> `collect_names` does not compile with `USE_HDF5=OFF`.

The body was `#ifdef BELFEM_HDF5`'d; the **signature** was not, and it names
`H5L_info2_t`, which `hdf5_types.hpp` does not stub. Since `cl_HDF5.hpp` includes
the header and `src/io` always compiles `cl_HDF5.cpp`, that was a tree-wide break
of a documented build option. Reproduced as a compile gate:

```
src/io/hdf5_tools.hpp:80:15: error: 'H5L_info2_t' does not name a type
```

Grok also corrected a comment I had written asserting that comparing `tName` to
`tPath` was sufficient to catch self-overwrite. It is not: `file.EXO` maps to
`file.exo`, the strings differ, and on the case-insensitive volumes this project
also ships on that is the same file. The guard is now on the extension, using the
`string_to_lower( filetype( ... ) )` idiom `Mesh::save` already uses.

## Fixed

Grid comparison rewritten to compare per-dimension and tolerate a length
mismatch; always-active validation of order, point counts and finite positive
steps before `new Mesh`; `values` length checked against the node count;
first-non-flag-argument selection; basename-only extension replacement with
dotfile handling; self-overwrite refused on the extension; non-2D/3D tables
skipped with a message instead of aborting the run; skipped and converted tables
both reported. In `hdf5_tools.hpp`, `collect_names` fully guarded,
`H5_ITER_NATIVE` → `H5_ITER_INC` so "the first table" is deterministic, and an
HDF5 1.12 `#error` stating the minimum that `belfem_find_package` cannot express.
Target gated on `USE_HDF5 AND USE_EXODUS`; `db2exo` added to
`BELFEM_INSTALL_EXECUTABLES`.

## Confirmed but deliberately not fixed

Soft-linked tables vanish without a message — emitting from inside an HDF5
callback would mean threading a logger through `void * aData`. The
`( points - 1 ) % order` check is not hoisted; the callee still errors on the
always-active tier, just without naming the table. `collect_names` keeps C++
linkage as a C callback; `extern "C"` on an `inline` in a namespace would export
an unqualified symbol for a risk Grok itself rates low.

## Adjudicated: order 3 degenerates ( ruling 2026-08-29 )

Order-3 tables reach Exodus with corner connectivity only — `fix_num_nodes` maps
QUAD16 → 4 and HEX64 → 8 ( `cl_Mesh_ExodusWriter.cpp:1055`, `:1076` ) — so
interior samples are written but connected to nothing. Christian accepted the
degradation, and a warning was added rather than silence. Codex's alternative is
to build the visualization mesh at order 1 always: node count and ordering are
unchanged so the field still aligns, every sample becomes a connected corner, at
the cost of showing trilinear where the solver interpolates quadratically. One
line either way. Every database in the tree is order 2, which renders fully as
HEX27, so nothing is currently affected.

**Ruled: degenerate.** Christian confirmed on 2026-08-29, after Codex's
reframing of the impact ( 8 of 64 nodes connected per HEX64 element, not merely
"corner connectivity" ). The order-1 remesh is declined rather than deferred:
the converter renders what the database is, and substituting a trilinear mesh
everywhere would misrepresent the faithful order-2 case to rescue a case that
does not occur in the tree. The warning stands as the implementation of that
ruling; no code change followed.

## Evidence

`exodus_path` probed 10/10 including both dotfile forms. `db2exo.cpp`,
`cl_Database.cpp` and `cl_HDF5.cpp` syntax-clean under the tree's own
`-Wall -Werror -pedantic-errors`; `cl_HDF5.cpp` and `cl_Mesh.cpp` also clean with
`BELFEM_HDF5` and `BELFEM_EXODUS` undefined, which is the regression gate for
Grok's P1. `check_doc_claims.py` 34/34.

**Not run:** `db2exo` has never been built or executed against a real database.
The gate still owed is `make check-fast` plus one conversion of `sp-ap_new.hdf5`
confirming both `jc` and `n` land as node fields.
