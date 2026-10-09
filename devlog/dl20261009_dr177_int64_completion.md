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
**Status:** fixes landed; delta audit, `make check` and the deck comparison on the pushed state follow
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
| `cl_CutFactory.cpp:443` | max node id into `InterfaceProcessor` | `id_t` parameter | id, fine |
| `cl_CutFactory.cpp:368-369` | node and edge counts in the debug writer | `id_t` | count, fine |
| `cl_CutFactory.cpp:2615` | neighbor count into `allocate_element_container` | `uint` | count, fine |
| `cl_CutSet.cpp:46` | bitset count | `uint` | count, fine |
| `cl_CutProcessor.cpp:1039` | hash key | `size_t` | not a simplex key, fine |
| `cl_InterfaceProcessor.cpp:216` | node id | `id_t` | id, fine |

After the fixes, the INT64-only warning list holds exactly the nine id/count rows above and
`cl_BeltedTree.cpp:332`.

## Changes

- `0579d3a`-equivalent (Christian, core): `reduce_pair`/`coreduce_pair` take `const index_t a`;
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
