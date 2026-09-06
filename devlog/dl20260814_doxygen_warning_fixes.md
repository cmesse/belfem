# Doxygen Warning Elimination: 20 → 0

**Date:** 2026-08-14
**Purpose:** Eliminate all 20 warnings in the `make doc` log at the source and lock the improvement in with a zero baseline.
**Scope:** 14 source/doc files + `cmake/report_doxygen_warnings.cmake`. No runtime behavior change.

## What was done

`cmake-build-debug/doc/doxygen_warnings.log` carried 20 warnings, previously
accepted as a 16-warning "parser limitation" baseline plus 4 markdown warnings.
All 20 turned out to be fixable at the source without changing behavior:

| root cause | count | fix |
|---|---|---|
| Anonymous comparator structs (`belfem::@N` unmatched `operator()`) | 9 | name the struct types (`OpVertexDegree`, `OpVertexID`, `OpVertexIndex`, `OpVertexLevel`, `OpVertexOwner`, `OpNodeIndex`, `OpSegmentIndex`, `OpHeatPoly`, `OpTransportPoly`) |
| C++11 extended friend declarations misparsed | 4 | `friend mesh::GmshReader;` → `friend class mesh::GmshReader;` (cl_Mesh.hpp ×3, cl_GT_GasData.hpp ×1) |
| `**Author, *Title***` nested-emphasis cascade in `doc/coding_philosophy.md` | 3 | inner italics switched to `_…_` on the Ousterhout and Lakos lines (798-799); the warnings surfaced 30 lines downstream (828/830/841) |
| `` `Cell<T>::operator()` `` parsed as unresolvable explicit link | 1 | reworded to `` `Cell<T>`'s `operator()` `` (coding_philosophy.md:587) |
| Duplicate `@param normJ`/`@return` merged across declaration + definition | 2 | dropped the header duplicates in `cl_Material.hpp`; the richer `powerlaws.hpp` blocks (formula, Rhyner 1993 / Messe 2023 citations) survive |
| Pointer-to-member-function member unparseable by doxygen 1.9.1 | 1+1 | `using SwapFunction = void ( MeshChecker::* )( mesh::Element * );` alias — applied to **both** `cl_MeshChecker.hpp` copies (see below) |

`cmake/report_doxygen_warnings.cmake`: `BELFEM_DOXYGEN_WARNING_BASELINE` 16 → 0;
the stale artifact list in its header comment replaced. Every future doxygen
warning is now a regression by definition.

## The one semantic delta (declared, verified harmless)

Naming a previously anonymous struct promotes its header-defined variable from
internal to external linkage. A 2-TU reproducer confirmed that naming without
`inline` breaks the link (`multiple definition`), so the three comparators that
lacked it (`op_Graph_Vertex_Degree.hpp`, `op_Node_Index.hpp`,
`op_Segment_Index.hpp`) gained C++17 `inline`, matching the four graph headers
that already carried it. The `static` gastables pair keeps internal linkage
unchanged. All comparators are empty and stateless and nothing in the tree
observes their address, so runtime behavior is identical; only mangled names
change (`belfem::@N` → named), which a normal incremental rebuild absorbs.

## Cross-review (jury round, exchange `review_doxygen_warning_fixes`)

Codex: safe to apply as-is (high). Grok: every hunk behavior-safe, but refused
"apply as-is" on two confirmed completeness gaps, both incorporated:

1. **The patched `src/fem/postproc/cl_MeshChecker.hpp` is the orphaned copy** —
   only `src/fem/kernel/cl_MeshChecker.cpp` is in a library
   (`src/fem/kernel/CMakeLists.txt`), and `src/fem/kernel/cl_MeshChecker.hpp`
   carried the identical unparseable member. Both copies now have the alias.
2. **The warning baseline would have silently stayed at 16.** Now 0.

Both auditors also caught the pre-registration miscounting the anonymous-struct
bucket (9, not 11 — two friend warnings were misallocated), and Grok correctly
downgraded an over-read `nm` evidence bullet; the 2-TU reproducer is the
load-bearing linkage evidence.

**Pre-existing debt noted, deliberately not touched:** two `MeshChecker`
headers/implementations share one include guard (`CL_MESHCHECKER_HPP`) and have
already diverged (kernel ctor carries parallel/edges/faces guards the postproc
copy lacks); the postproc `cl_MeshChecker.cpp` is compiled by no target.
Merging or deleting is a separate decision.

## Verification

- doxygen 1.9.1 over the applied tree (same INPUT as `make doc`): **0 warnings**.
- All ~20 affected/consumer TUs compile clean with the production
  `compile_commands.json` flags (`-Wall -Werror -pedantic-errors -std=gnu++17`,
  Blaze backend).
- 8 comparator-consumer TUs object-compiled; `ld -r` partial link clean.
- NOT run: full `make` and the real `make doc` (user runs builds) — the next
  `make doc` is the closing gate and should report 0.
