# VTK linked privately, so tests and the solver stop compiling against VTK's bundled HDF5; DR-201 filed

**Date:** 2026-11-08 (session opened 2026-11-07)
**Purpose:** Christian asked for "DR-201", which did not exist yet. It turned out to be the nonfree
tests that `make check` runs without building them (filed, not fixed). On the way, his VTK-ON
`make check` showed 54 failures; 51 of them traced to DR-202 (filed the same evening by the linalg
session), which this session planned, fixed and closed. Christian also asked whether BELFEM can
run with VTK on this Mac, given that VTK is built with Apple Clang and BELFEM with GNU.
**Module:** build system (`config/scripts`, `src/visualizer`; one nonfree CMake file)
**AIs involved:** Claude Opus 5.5 (primary); Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high
(plan jury); Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high (code jury). Peer session
belfem-de (linalg, DR-200/DR-203) coordinated by message throughout.
**Status:** COMPLETE for DR-202; DR-201 filed and open; `cross_review.sh` mktemp defect reported
**Verification:** reproducer + Christian's gates. Diagnosis: `make check` 2026-11-07 20:26 with
`USE_VTK=OFF` in `cmake-build-debug` — io, circuit, physics, `BfmFile` green, only DR-203's three
cases red. Fix: configure rung run by Claude (scratch tree, `-DUSE_VTK=ON -DUSE_NONFREE=ON`: no
`vtkhdf5` in any compiling target's `flags.make`); gates G1–G6 of the plan run and confirmed by
Christian 2026-11-08 ("gates confirmed")
**Exchange:** `tmp/ai_exchange/review_dr202_vtk_private_link.md` (+ `.diff`)

## Summary

`target_link_libraries( belfem_lib PUBLIC ${VTK_LIBRARIES} )` handed every consumer of the library
the interface include directories of every VTK component. `VTK::vtkhdf5_src` carries
`include/vtk-9.6/vtkhdf5;…/vtkhdf5/src` (`/opt/scls/lib/cmake/vtk-9.6/VTK-targets.cmake:3100`;
`vtkhdf5_hl_src` adds `hl/src`). In the VTK-ON `test_io` flags those came at tokens 98/100/102,
ahead of `/opt/scls/include` at 161. Both `hdf5.h` copies use the guard `HDF5_H` (Grok), so the
first include wins for the whole translation unit. The test files therefore compiled against VTK's
renamed HDF5 (`vtk_hdf5_mangle.h:3202`, `H5Tget_class` → `vtkhdf5_H5Tget_class`), while
`libbelfem.a` used SCLS's. Two HDF5 libraries then exchanged type ids in one process, and every
read failed with `vtkhdf5_H5Tget_class(): not a datatype`.

Library objects were clean: module OBJECT libraries do not link VTK and see only
`include/vtk-9.6` through `BELFEM_INCLUDES` (`find_vtk.cmake:56`). (high)

## What landed

- **`config/scripts/Add_BelfemLibrary.cmake`.** VTK now links `PRIVATE`, with a comment stating why.
  A static library still exports a private dependency as `$<LINK_ONLY:…>`, so every consumer keeps
  the 160 VTK libraries on its link line and gets none of VTK's compile usage requirements.
- **`src/visualizer/CMakeLists.txt`.** `visualize` links each VTK library as `$<LINK_ONLY:…>`, one
  per call; an unquoted `$<LINK_ONLY:${VTK_LIBRARIES}>` would wrap only the first element (Grok).
  The plain signature is required, because `Add_Executable.cmake:45` links without a keyword and
  CMake refuses to mix the two forms. The loop is redundant for the static library and needed for
  `USE_SHARED_LIBS=ON`.
- **Nonfree VTK executables.** The same loop. Under the nonfree rule it went through no exchange
  and no jury, so it is not described further here.
- **Plan.** `todo/closed/dr202_vtk_private_link.md` v2. The register row moved to
  `debt_register_closed.md`.

**My own first R2 was wrong, and I caught it before dispatch.** Linking `visualize` against
`${VTK_LIBRARIES}` with the full interface would have recreated the bug in that executable:
`visualize.cpp:18` includes `cl_Mesh.hpp`, which reaches `<hdf5.h>` through
`cl_Mesh.hpp:34` → `hdf5_types.hpp:16`. Before the fix, `visualize` was in fact already
compiling against the bundled headers. The solver is on the same path (`belfem.cpp:21`), and
R1 alone fixes it.

## What the juries changed

**Plan round.** No P0.
- Both auditors: the purpose sentence still described the rejected full-interface R2, and G2 was
  not executable as written. G2 now names its binaries and carries a positive control (`_H5Tget_class`
  ≥ 1 in `test_io`).
- Grok: G1 skipped `tests/mesh`, where `BfmFile` lives (`tests/mesh/CMakeLists.txt:18`), and the
  solver. The flags files and the 20:12 log I had cited were overwritten by Christian's VTK-OFF
  reconfigure; both are now cited as snapshots. The O2 wording also overclaimed.
- Codex: O1 needed a full VTK-ON `make` gate.
- **Refuted:** Codex's "installed headers lose VTK's usage requirements". The tree has no
  `install( EXPORT … )` and no package config, so no usage requirement reaches a downstream consumer
  under either `PUBLIC` or `PRIVATE`.

**Code round.** No functional defect, no P0/P1.
- One P2, raised by all three voices: the comment's "only the visualizer objects" was too narrow. It
  now reads "only the VTK-based modules", applied on Christian's approval.
- Grok's one open risk was that `$<LINK_ONLY>` might still apply VTK's include directories to
  `visualize`. The scratch configure had already settled it: `visualize.dir/flags.make` and the
  nonfree executables' flags hold `include/vtk-9.6` only.
- Off-topic: Grok spent a paragraph on rpm specs and `libscotchmetis`. That was prompted by a global
  instruction, not by the diff, and was ignored.

## Answer to "can VTK run on this Mac?"

Yes, and the failure was not a compiler-mix problem: Apple Clang throughout would have failed
the same way. What the mix does risk is plan O2, left open. VTK is built against libc++ and
BELFEM against libstdc++, so any `std::` object crossing the boundary is broken. Usually that is
an unresolved symbol, but an inline header function or a C signature that carries such an object
can fail silently. The visualizer sources pass no `std::` types (`SetFileName( … c_str() )`). (medium)

## DR-201, filed, not fixed

With `USE_NONFREE=ON`, ctest runs `nonfree-suite-a`, `nonfree-suite-b` and `the-nonfree-python-extension`, but the `check` target's
dependency loop lists only the open-source suites. The 2026-11-07 19:11 run therefore showed
"Unable to find executable" for the first two and "No module named 'the-nonfree-python-extension._the-nonfree-python-extension'" for the
third. With nonfree off, none of them is registered, which is correct.

Proposed fix, awaiting approval: a global property filled by `Add_Test.cmake` for test sources
under `nonfree/`, plus by `nonfree/python/CMakeLists.txt`, that the root `check` block depends on.
No nonfree name would appear in the open-source `CMakeLists.txt`.

## Tooling defect found

`scripts/cross_review.sh:95-96` passes `crossrev_subject_XXXXXX.md` to `mktemp`, and BSD mktemp
replaces only trailing X's, so the name is literal. My first code-round dispatch failed with
"mkstemp failed … File exists" because the peer's DR-203 round held the file. Workaround: a private
`TMPDIR`. The EXIT trap on `:97` removes the file after a normal finish, so only concurrent rounds
collide. Not fixed here: a tooling script owes its own jury round (ruling 2026-09-22).

## Open Questions

- DR-201: which fix, and approval to edit.
- Plan O2 (libc++/libstdc++) and O3 (VTK's include path stays global in `BELFEM_INCLUDES`) remain
  open in the closed plan; neither is owed by DR-202.

## Files Updated

`config/scripts/Add_BelfemLibrary.cmake`, `src/visualizer/CMakeLists.txt`, nonfree VTK executables'
`CMakeLists.txt` (nonfree repository), `todo/closed/dr202_vtk_private_link.md` (new),
`todo/README.md`, `todo/debt_register.md` (DR-202 out, DR-201 in), `todo/debt_register_closed.md`,
`devlog/README.md`.
