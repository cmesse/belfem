# DR-177 finished: src/homology compiles under BELFEM_INT64, with a sweep of every key-holding variable

**Date:** 2026-10-09
**Purpose:** The blind audit of f3841bc (DR-177, the `Chain` key type) cleared the code for the
default build but refuted its INT64 claim: `cl_CutFactory.cpp:526` was ill-formed and four core
sites still narrowed simplex keys into `uint`. Christian ruled "Please fix the 64 bit issue".
This entry records the compile gate, the sweep, the fixes, and the guide corrections the audit
asked for.
**Module:** `src/homology` (`cl_SimplicialComplex.{hpp,cpp}` by Christian; `cl_CutFactory.cpp`,
`doc/homology_usage_guide.md`), `todo/debt_register_closed.md`
**AIs:** Claude Opus 5.5 ("tests" session): sweep, non-core fixes, guide. Audit of f3841bc:
Codex `gpt-6-astra`/xhigh and Grok `grok-4.7`/xhigh, blind (`tmp/ai_exchange/f3841bc_audit_*.md`).
Codex `gpt-5.6-luna`/medium language sweep of the guide corrections.
**Status:** COMPLETE within the scope Christian ruled on 2026-10-10: simplex keys and mesh-scale counts are 64-bit-safe; generator, cut, sideset and Smith-matrix cardinalities stay 32-bit (Open). Landed in 7011f44, fc41443 and 021dbc8
**Artifacts:** scratchpad only: INT64 compile tree, `-Wconversion` logs of both builds

## Gate: compile under BELFEM_INT64

Configure with `-DUSE_MKL_64BIT_API=ON` (defines `BELFEM_INT64`), release, mkl flavor. There is
no 64-bit SCLS stack, so this gate is compile-only: nothing was linked or run under INT64.

- f3841bc: `belfem_homology` fails in one file, `cl_CutFactory.cpp`: `:526` (an
  `OrderedMap< id_t, int > &` bound to the cochain's `OrderedMap< index_t, int >`) and
  `:744`/`:745` (`-Werror=overflow`: the `gNoIndex` sentinel stored in a `Vector< uint >`). The
  eight `tests/homology` objects compile.
- After the fixes: `belfem_homology` and the eight test objects compile with 0 errors. The
  default release `belfem` builds with 0 warnings.
- Other libraries, report only: compiling every other object library under INT64 gives 63
  errors in 14 files outside `src/homology`: `src/mesh` (`cl_TensorMeshFactory`,
  `cl_ProtoMesh`, `cl_Mesh_Distributor`), `src/sparse` (STRUMPACK, SuperLU, PETSc, the dist
  matrix, `petsctools`), `src/fem/kernel` (three DofMgr units), `src/fem/maxwell`
  (`cl_Maxwell_FieldList`), `src/fem/iwg` (`cl_IWG`), `src/physics/database`. Not widened here.

## Sweep: every variable that holds a simplex key

Two passes. (1) `-fsyntax-only -Wconversion -Wsign-conversion` on every `src/homology` source,
once with and once without `BELFEM_INT64`. Warnings that exist only under INT64 (compared by
line and message, so a line that warns in both builds for different types still shows) are
the 64→32 narrowings. (2) Every call that takes a key (`add_simplex_to_*`, `set_coefficient`,
`get_coefficient`, `key_exists`, `erase_key`, `reduce_pair`, `coreduce_pair`, the chain and
cochain map lookups): the declared type of its key argument, found by walking back to the
declaration. Pass 2 catches a key carried in a `uint` that never narrows, which the compiler
cannot see. 81 call sites: 53 `index_t`, 20 structured bindings over `index_t` maps, and the
rest listed below. Maps keyed by `id_t` are keyed by node, block and sideset ids; maps keyed by
`int` in `Cohomology` and `Homology` are keyed by the dimension.

| site | holds | before | verdict |
|---|---|---|---|
| `cl_CutFactory.cpp:526` | cochain key map | `OrderedMap< id_t, int > &` | key: **fixed**, `index_t` |
| `cl_CutFactory.cpp:744-745` | local edge number of an element, sentinel `gNoIndex` | `Vector< uint >` | small count: sentinel **fixed** to `BELFEM_UINT_MAX`; never compared, same value on the default build |
| `cl_SimplicialComplex.cpp:640`, `.hpp:95` | `reduce_pair` parameter `a` | `const uint` | key: **fixed by Christian** |
| `cl_SimplicialComplex.cpp:668`, `.hpp:100` | `coreduce_pair` parameter `a` | `const uint` | key: **fixed by Christian** |
| `cl_SimplicialComplex.cpp:948` | `reduceOmit` `tKey` | `uint` | key: **fixed by Christian** |
| `cl_SimplicialComplex.cpp:1244` | `coreduceOmit` `tKey` | `uint` | key: **fixed by Christian** |
| `cl_BeltedTree.cpp:332` | cochain key passed to `Mesh::element( id_t )`, an id lookup | — | key used as an id: latent defect, open (see below) |
| `cl_CutData.cpp:139` | local side number 0..2 passed to `set_master` | `index_t` | small count, fine |
| `cl_CutData.cpp:32` | cut id = generator index + 1 | `id_t` | id, fine |
| `cl_CutFactory.cpp:1507` | facet ids | `id_t` | id, fine |
| `cl_CutFactory.cpp:443` | the original-sideset count passed to the `uint` parameter of `InterfaceProcessor` (`max_node_id()` is already `id_t`) | `uint` parameter | count, fine (corrected after the delta audit; first described as the node id) |
| `cl_CutFactory.cpp:368-369`, `:372`, `:394` | global node and edge counts and the loops over them, in the caller-less debug writer | `id_t`, `uint` | count and global index: **widened to `index_t`** after the delta audit |
| `cl_CutFactory.cpp:2615` | neighbor count into `allocate_element_container` | `uint` | count, fine |
| `cl_CutSet.cpp:46`, `:56` | number of cuts a set touches, and the loop over the abstract nodes | `uint` | count: **widened to `index_t`** after the delta audit (unbounded in principle) |
| `cl_CutProcessor.cpp:1039` | hash key, `size_t` into a `Map< size_t, … >` | `size_t` | not a narrowing at all; listed by the census, not by the compiler |
| `cl_InterfaceProcessor.cpp` (`duplicate_nodes`, line moved by the revert in 99f8a8f) | node id | `id_t` | id, fine |

After 7011f44 the INT64-only warning list held the id and count rows above and
`cl_BeltedTree.cpp:332`. The delta audit found what a warning comparison cannot see: values
that are `size_t` in both builds and narrow into `uint` in both (the core count accessors), and
global index loops that never narrow because they start as `uint`.

## Changes

- 0579d3a (Christian, core): `reduce_pair`/`coreduce_pair` take `const index_t a`;
  both `tKey` are `index_t`; the two declarations follow. Applied from his working tree verbatim.
- `cl_CutFactory.cpp`: the key map at `save_edges` is `OrderedMap< index_t, int > &`; the
  master/slave local-edge vectors are filled with `BELFEM_UINT_MAX`.
- Guide (`homology_usage_guide.md`), the audit's findings: `Face::index()` for 3D triangles,
  not `Facet::index()`; the constructor flags describe ownership correctly (two levels, with
  the dimension-0 and dimension-3 exceptions); only `operator+`, `operator-` and `operator*`
  modify `this`, `operator()` returns the pairing; `add_chain_to_chain` adds the stored
  boundary and coboundary, so the source must own the members; `Cochain` returns early on a
  zero coefficient, `Chain` does not; scaling a cochain means adding it into an empty one; the
  `BeltedTree` writer signatures match the header.
- `todo/debt_register_closed.md`: DR-177's resolution states what the INT64 gate is.

## Open

- `cl_BeltedTree.cpp:332`: `create_cohomology_field` looks up `aEdgeMesh->element( tID )`, an
  id lookup, with a cochain index. Latent; the writer has no callers. Not fixed, on Christian's
  word (0.9.3 without new scope).
- INT64 is compile-only for `src/homology`; nothing was linked or run under INT64, and 63
  errors remain in other libraries.

## Delta audit of 202e62a..7011f44 and the follow-up

Blind, same brief, separate worktrees: Codex `gpt-6-astra`/xhigh and Grok `grok-4.7`/xhigh
(`tmp/ai_exchange/delta_7011f44_*.md` in the shared checkout's exchange). Both verify the edited hunks (core retype,
`CutFactory:526`, the sentinel), the guide, and default-build identity. Grok: CLEARS (~85 %),
every stored key is `index_t`, the only 32-bit pass-through is `BeltedTree:332`. Codex: does
not clear the completeness claim, citing cardinalities rather than keys: the caller-less debug
writer's global counts and loops (`CutFactory:368-394`), the cut count in `CutSet:46`, and the
core's five count accessors returning `uint` from a map size (`cl_SimplicialComplex.hpp:174-195,
262-312`), which would make the omit sweeps at `.cpp:949` and `:1245` see zero simplices at
2³² or more. Both found three table rows misdescribed (corrected above). Every citation was
checked against the code.

Follow-up: outside the core, the debug writer and `CutSet` counters are `index_t` (021dbc8).
Inside the core, the five count accessors return `index_t` (fc41443): proposed as a patch,
reviewed and applied by Christian.
The log casts to `unsigned int` and the `int` sizes of the dense Smith matrices stay; they are
display and dense-matrix bounds far below 2³².

Gates: `make check` on 7011f44 (default release, mkl flavor, Armadillo) 17/17. The five decks
(tapestack2d_layered, tapestack3d, helix, gantry, corc_periodic_bc) set up on 7011f44 with logs
and `.bfm` identical to 337a318 apart from `/meta/git`.

## Second delta audit (f6dca70..021dbc8) and the scope of the INT64 claim

Blind, same seats (`tmp/ai_exchange/delta_021dbc8_*.md`). Both verify the change: the accessor
declarations and definitions agree, every caller behaves as before on the default build, the
widened debug-writer and `CutSet` values reach only `index_t`/`size_t` sinks, and the default
build is type-identical. Both refute a module-wide completeness claim with a census of
pre-existing 32-bit counters: the `uint` ordinals in `createMatrixFromBoundaryMap` and
`createMatrixFromCoboundaryMap` and the `int` captures in the latter, `generatorsOfHomology` /
`generatorsOfCohomology` counters and the `Map< int, uint >` Smith ranks, the `uint` dimensions
throughout `fn_Smith`, the cut counts (`CutProcessor::mNumberOfCuts`, `CutData::mNumberOfCuts`),
the sideset counts into `InterfaceProcessor` (capped by `Mesh::number_of_sidesets()` returning
`uint`), per-node duplicate multiplicities into `Node`'s `int` counter, BeltedTree traversal
counters, a vector length serialized in `Topology`'s `Vector< id_t >` buffer, and the Exodus
writer's `int` connectivity and `uint` counts. Grok rates all of them P2 (no blocker, they need a
reduced complex or a generator, cut or sideset count of 2³¹ to 2³², which no dense
`Matrix< int >` Smith factorization reaches); Codex rates many P1 "unless explicitly bounded".

Christian ruled 2026-10-10: push fc41443 and 021dbc8, and scope the claim instead of another
code round. **The INT64 claim for `src/homology` is: it compiles, and simplex keys and
mesh-scale entity counts are carried in `index_t`. Generator, cut, sideset and
reduced-complex/Smith cardinalities stay 32-bit and are not supported above 2³¹.** The census
above is the Open list for that work; a size guard before the dense Smith matrices are sized
would be its first step.

Gates on 021dbc8: INT64 `belfem_homology` and the eight `tests/homology` objects 0 errors;
default release `belfem` 0 warnings; `make check` 17/17 (release, mkl flavor, Armadillo).
