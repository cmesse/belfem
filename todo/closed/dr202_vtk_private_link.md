# DR-202: Link VTK Privately, So Consumers Stop Inheriting VTK's Bundled HDF5

**Date:** 2026-11-07
**Purpose:** With `USE_VTK=ON`, every target linked against `libbelfem` compiles `<hdf5.h>` from
VTK's renamed, bundled HDF5 rather than the SCLS one, and every HDF5 test fails. Link VTK to the
project library `PRIVATE`, and link the executables that really include VTK headers against the VTK
libraries without VTK's usage requirements (they reach the VTK headers through the same global
`include/vtk-9.6` path the visualizer objects use).
**Module:** build system (`config/scripts`, `src/visualizer`)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** ✅ COMPLETE (2026-11-08). R1–R3 landed: `Add_BelfemLibrary.cmake` links VTK `PRIVATE`; `visualize` and the nonfree VTK executables link VTK `$<LINK_ONLY:…>`. Plan jury and code jury (Codex `gpt-5.6-terra`/high + Grok `grok-4.7`/high) no P0/P1; the code round's one P2 (R1 comment scope) applied on approval. Gates G1–G6 confirmed by Christian 2026-11-08; G1 also green at configure in a scratch tree. Residue, not owed by DR-202: O2 (libc++/libstdc++ mix) and O3 (VTK include path global) stay open here

> **Scope guards:**
> - In scope: how VTK's usage requirements reach targets. Nothing in `src/` C++ changes.
> - Out of scope: VTK module autoinit (`vtk_module_autoinit` is not called anywhere today, so this
>   plan neither adds nor removes it); the libc++/libstdc++ mix between Apple-Clang-built VTK and
>   GNU-built BELFEM (O2); DR-203 (mesh, unrelated).
> - The nonfree tree's VTK consumers get the same change under the nonfree rule (no exchange, no
>   jury), in step R3. They are not named or reviewed here.

---

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| Every HDF5 read in a test fails, `H5Tget_class: not a datatype`, class −1 | the test TU's `<hdf5.h>` comes from VTK's renamed copy (`vtk_hdf5_mangle.h` renames `H5Tget_class` → `vtkhdf5_H5Tget_class`); `libbelfem.a` was compiled against SCLS HDF5; the two libraries exchange type ids in one process | `make check` 2026-11-07 20:12, VTK ON: 51 cases in io, circuit, physics and `BfmFile` (in `test_mesh`, `tests/mesh/CMakeLists.txt:18`), log line `vtkhdf5_H5Tget_class(): not a datatype` (log since overwritten). **Reproducer:** the 20:26 run with `USE_VTK=OFF` passes all of them and fails only DR-203's three cases |
| VTK's HDF5 include directories come before SCLS's | `target_link_libraries( ${LIBTARGET} PUBLIC ${VTK_LIBRARIES} )` (`config/scripts/Add_BelfemLibrary.cmake:40-41`) passes the interface include directories of every VTK target on to every consumer. `VTK::vtkhdf5_src` carries `include/vtk-9.6/vtkhdf5;include/vtk-9.6/vtkhdf5/src` (`/opt/scls/lib/cmake/vtk-9.6/VTK-targets.cmake:3100`) | VTK-ON snapshot of `tests/io/CMakeFiles/test_io.dir/flags.make` (one token per line): `vtkhdf5`, `vtkhdf5/src`, `vtkhdf5/hl/src` at 98/100/102, `/opt/scls/include` at 161; overwritten by the VTK-OFF reconfigure. `src/visualizer/CMakeFiles/visualize.dir/flags.make` still holds the VTK-ON flags |
| Library objects are clean | module OBJECT libraries do not link VTK; they see only `include/vtk-9.6` through `BELFEM_INCLUDES` (`config/system/find_vtk.cmake`, the `list( APPEND BELFEM_INCLUDES …)` line) | `src/visualizer/CMakeFiles/belfem_visualize.dir/flags.make`: `/opt/scls/include/vtk-9.6` only, no `vtkhdf5` |
| CI never sees it | `USE_VTK` defaults OFF (`CMakeLists.txt:123`); `.gitlab-ci.yml` does not set it | — |

Who actually includes VTK headers (`grep '#include <vtk\|#include "vtk[A-Z]'`): `src/visualizer/`
(`cl_VTK_Curve.cpp`, `visualize.cpp`, `vtktypes.hpp`) and one nonfree module. `src/mesh/vtktools.*`
and the other "vtk" names write the VTK file format by hand and use no VTK library code.
The solver is on the same path: `src/executables/belfem.cpp:21` includes `cl_Mesh.hpp`, which
reaches `<hdf5.h>` through `hdf5_types.hpp:16`; R1 fixes it with no extra step.

Downstream consumers of an installed BELFEM are not affected either way: the headers are installed
(`CMakeLists.txt:370-376`), but no `install( EXPORT … )` or package config exists, so no usage
requirement reaches them under `PUBLIC` or `PRIVATE`.

**Bottom line:** a library that needs VTK only for its own visualizer objects advertises VTK's
whole interface, bundled HDF5 included, to every consumer. (high confidence)

## 2. Mechanism

`PRIVATE` on a **static** library: CMake still exports the VTK libraries on the consumer's link
line (as `$<LINK_ONLY:…>` in `INTERFACE_LINK_LIBRARIES`), so the visualizer objects inside
`libbelfem.a` keep resolving, but no include directory or compile definition of VTK propagates.
`PRIVATE` on a **shared** library (`USE_SHARED_LIBS=ON`): `libbelfem.dylib` links VTK itself;
consumers do not link VTK at all, which is correct for everyone but the VTK-header executables,
and those link it explicitly (R2).

VTK's bundled HDF5 is symbol-renamed, so linking both HDF5 libraries into one test binary is
harmless once the headers agree; both are already on the VTK-ON link line of `visualize`
(`visualize.dir/link.txt`). (high confidence)

Rejected: forcing `/opt/scls/include` ahead of VTK's directories. Both `hdf5.h` copies stay on
the path and correctness hangs on search order, which any later include change can flip again.

## 3. Steps

- [x] **R1** `config/scripts/Add_BelfemLibrary.cmake:40-41`: `PUBLIC ${VTK_LIBRARIES}` →
  `PRIVATE ${VTK_LIBRARIES}`, comment states why (only the VTK-based modules need VTK, and VTK's
  bundled HDF5 must not reach the consumers). Done 2026-11-07; comment reworded from "only the
  visualizer objects" after the code jury (all three voices), approved 2026-11-08.
- [x] **R2** `src/visualizer/CMakeLists.txt`, after the `Add_Executable.cmake` include: link
  `visualize` against VTK **without** VTK's usage requirements:
  `foreach( tVtkLib ${VTK_LIBRARIES} ) target_link_libraries( visualize $<LINK_ONLY:${tVtkLib}> ) endforeach()`.
  Two constraints shape this line (amended before dispatch, Claude):
  - **Link only, not the full interface.** `visualize.cpp` includes `cl_Mesh.hpp`
    → `hdf5_types.hpp:16` → `<hdf5.h>`. A plain `target_link_libraries( visualize ${VTK_LIBRARIES} )`
    would hand this TU the `vtkhdf5` include paths again and recreate the DR-202 mismatch in
    `visualize` itself (it has it today: `visualize.dir/flags.make` contains `vtkhdf5`). The VTK
    headers it does include are reached through `include/vtk-9.6` in `BELFEM_INCLUDES`, exactly as
    for the visualizer objects. The interface compile definitions VTK's targets carry (`FMT_SHARED`,
    `SCN_*`, `kiss_fft_scalar`, `DIY_*`, `VTK-targets.cmake:131,212,371,1758,1872`) belong to
    third-party modules `visualize` does not include. (medium confidence; G4 decides)
  - **Plain signature, not `PRIVATE`:** `Add_Executable.cmake` calls
    `target_link_libraries( ${EXECNAME} ${LIBTARGET} )` without a keyword, and CMake refuses to mix
    the plain and keyword forms on one target.
  For a static `libbelfem` R2 is redundant (R1 already exports the libraries link-only); it is
  needed for `USE_SHARED_LIBS=ON`, and kept unconditional so both builds read the same.
  (after: R1)
- [x] **R3** (nonfree, outside this round) the same one-line addition for the nonfree
  executables that include VTK headers. (after: R2)
- [x] **R4** gate, run by Christian (§4). (after: R1–R3)

## 4. Gate

| ID | Check | Pass |
|---|---|---|
| G1 | reconfigure `cmake-build-debug` with `-DUSE_VTK=ON` | `grep -l vtkhdf5 tests/*/CMakeFiles/*/flags.make src/executables/CMakeFiles/*/flags.make src/visualizer/CMakeFiles/*/flags.make` prints nothing |
| G2 | after `make`: `for f in test/test_io test/test_mesh test/test_circuit test/test_physics bin/belfem; do nm -u $f \| grep -c vtkhdf5_H5; nm -u $f \| grep -c ' _H5Tget_class$'; done` | first count 0 for every binary; second count ≥ 1 for `test_io` (the positive control: the binary still calls SCLS HDF5) |
| G3 | full `make` with `USE_VTK=ON`, then `make check` | everything builds (closes O1); io, circuit, physics, `BfmFile` green; mesh fails only in `TensorMeshBSpline.*` until DR-203 lands |
| G4 | `make visualize` (and the nonfree VTK executables when nonfree is on), then one `visualize` run that writes a PNG | compiles, links, and the run produces the image (catches a runtime factory miss the link would not) |
| G5 | a `USE_VTK=OFF` tree | configure and `flags.make` unchanged: R1 sits in `if( USE_VTK )`, and R2's directory is added only under VTK (`src/CMakeLists.txt:20-21`) |
| G6 | `USE_SHARED_LIBS=ON`, `USE_VTK=ON` | `visualize` links; `make check` as G3 |

## 5. Open Questions

- **O1** Does anything other than `visualize` (and the nonfree executables) compile against VTK
  headers through the inherited interface today? The grep in §1 says no; G4 plus a full `make`
  with `USE_VTK=ON` is the check.
- **O2** The libc++ (Apple Clang VTK) versus libstdc++ (GNU BELFEM) mix. Any VTK call that passes
  a `std::` type across the boundary cannot work. Usually that is an unresolved symbol (the two
  libraries mangle `std::` differently), but an inline header function or a C signature that
  carries a `std::` object can break layout silently. Independent of this fix; logged, not decided.
- **O3** R2 keeps VTK's include path global (`find_vtk.cmake:56` → `BELFEM_INCLUDES`) rather than
  attaching it to the visualizer targets. That is the existing mechanism, kept on purpose to keep
  DR-202 small; moving it is a separate cleanup (Codex, plan round 1).

## 6. Audit Trail

- `tmp/ai_exchange/review_dr202_vtk_private_link.md` — code round (same slug, diff `review_dr202_vtk_private_link.diff`, own files only): no functional defect, one P2 comment rewording (Claude, Codex, Grok); Grok's open risk on `visualize` include suppression resolved by the scratch configure. First code dispatch failed on the `cross_review.sh:95-96` literal-`mktemp` collision, re-run with a private `TMPDIR`.
- same file — plan jury round 1 (Codex gpt-5.6-terra/high, Grok grok-4.7/high); verification and reconciliation appended. v2 folds in: purpose wording (both), named G2 with positive control (both), G1 over all test, solver and visualizer targets (Grok), full VTK-ON build gate (Codex), PNG smoke in G4 (Grok), O2 wording (Grok), O3 (Codex). Codex's installed-header finding refuted: no `install( EXPORT )` exists.
