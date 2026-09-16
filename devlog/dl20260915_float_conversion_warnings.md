# Devlog 2026-09-15 — Clearing the `-Wfloat-conversion` warnings; third-party headers as system headers

**Date:** 2026-09-15
**Topic:** The 60 warning lines a post-0.9.1 `make` printed (`-Wfloat-conversion`, added in `8465284`): 48 from STRUMPACK headers under `/opt/scls/mkl/include`, 12 from BELFEM sources, 2 pre-existing gfortran rank-mismatch lines
**AIs involved:** Claude (Fable 5.1), Codex (gpt-5.6-terra, high), Grok (grok-4.6, high) — plan jury and code jury, both blind
**Claude Confidence:** high on the conversion sites and the mechanism; medium on the `-isystem` builtin-directory guard on hosts not available here (Ubuntu 24.04 / GCC 13, Apple Clang)
**Codex Audit Confidence:** high (plan: "do not approve unchanged", 6 findings; code: "request changes", 3 findings — all adjudicated below)
**Grok Audit Confidence:** high (plan: 7 challenges; code: "ship the code, fix the prose", 4 refutations)
**Literature References:** none — build system and mechanical conversions
**Verification:** probe level — all nine affected translation units (plus `test_StringTools.cpp` after the guard round) compile with `-fsyntax-only` under their module's real `flags.make` (`-Wall -Werror -Wfloat-conversion`) with zero warnings; the new CMake include block exercised with `cmake -P` on a synthetic list (duplicates, trailing slash, a builtin dir, empty list). `check_doc_claims.py` 38/38. **Then verified — focused regression:** Christian reconfigured and rebuilt `cmake-build-debug` (GCC 11, `USE_DEBUG=ON`, mkl flavor) at `14bb703` and ran `make check`: **17/17 passed** (fast 8, mpi 4). `flags.make` of the rebuilt tree carries `-isystem` with the SCLS dir once and no per-test `BEFORE` include, so the new include block was in effect; `test_core` was rebuilt after the commit and `UnitToSiExponentGuard` appears in `LastTest.log`, so the new guard test ran. Not reported: the rebuild's warning count, which is the direct falsifier for the `-isystem` claim. Not run: any other platform.

## Summary

Two classes of warning. The 48 foreign ones came through `config/compiler/finalize_compiler.cmake`, which spliced every third-party include directory into `CMAKE_CXX_FLAGS` as `-I`; they are now `-isystem`, which GCC and Clang exclude from diagnostics, with a guard that keeps `-I` for any directory already in `CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES` (so `/usr/include` can never be moved ahead of the standard library headers). No foreign file was touched. The 12 BELFEM sites were all intended narrowings of `ceil`/`round`/scale results and are now explicit casts of the *same* expressions; two more of the same class that only the test build compiles (`cl_Genome.hpp`, Grok's find) are included. The one real API change: `spline::create_helpmatrix` took its point count as `const real &`; it is now `const index_t aSize, const real aDeltaX`, both by value (Christian's decision; house convention is 399 by-value to 30 by-reference scalar parameters in `src/**/*.hpp`).

## Key Findings

- **The flag is doing its job, and the jury round was worth its cost.** Of my original plan, Codex refuted three rows as *not* behavior-preserving: `std::stoi` for a unit exponent (`stod("1e2")` = 100, `stoi` = 1), and `static_cast<size_t>(ceil(x)) + 1` versus the original `size_t(ceil(x) + 1)` (differs, with UB, for `ceil(x) < 0`). Both would have shipped as "no behavior change".
- **My caller census was truncated** — the tests grep had been piped through `head` while the plan claimed an enumeration (Iron rule 4). Full census: 40 calls in 9 files, every first argument integer-typed (`uint`, `size_t`, `int`, `index_t`, literals).
- **`config/scripts/Add_Test.cmake:208`** gave every test `BEFORE PRIVATE $ENV{SCLS}/include`. Under `-isystem`, GCC drops a duplicate `-I`, so tests would have searched SCLS after BELFEM's headers — the order the library always used. Probe: zero basename collisions between BELFEM's 170 include-dir headers and SCLS's 611 top-level entries; gtest exists only under SCLS on this host. The line was redundant in both `$SCLS` branches (`find_scls.cmake:18` already adds it; unset gives `-I/include`) and is removed.
- **`make check-fast` does not run `test_physics`** — it is `ctest -L fast` and `tests/physics/CMakeLists.txt` sets no label. The gate for the material rows is `make check`.
- **The two gfortran "Rank mismatch" warnings at `pardisotools.f90:369,424` are by design**: `-fallow-argument-mismatch` (`config_gcc.cmake`, both configurations) demotes gfortran ≥ 10's error on PARDISO's F77 scalar-dummy idiom to a warning. Christian asked for this to be documented; a rationale comment now sits above the `USE_DEBUG` block. In this tree only the PARDISO wrapper trips it (MUMPS goes through `dmumps_struc.h`, ARPACK's calls are rank-consistent).
- **My first `-isystem` comment and CHANGELOG bullet told a false story** (Grok): they attributed 0.9.0's Ubuntu build failure to a vendor header, but 0.9.1's record says it was `cl_Logger.hpp` and `fn_sprint.hpp` — BELFEM's own. Both rewritten as preventative. L-08: a correct change with a false reason poisons the next design.
- **`cmake -P` probe caught my first include block** emitting `/opt/scls/mkl/include` twice: it deduplicated before stripping trailing slashes. Fixed to normalize → dedup → emit, and the deduplicated list is written back to `BELFEM_INCLUDES` so the Fortran loop shares it.
- **Clang qualifier** (Grok): `config_gcc.cmake:76` has `-Werror=uninitialized` only; the "`-Werror` in both configurations" statement is GCC-only, and both are gated by `USE_WARNINGS` (default ON), not `USE_DEBUG`.
- **Follow-up guards (Christian's ruling: input parsing is setup-tier, `BELFEM_ERROR` is right).** `Section::get_int` and the unit-exponent parse in `unit_to_si` cast an unchecked `double` to `int` (UB for non-finite or huge values). The first draft of the guards went through its own blind jury (Codex terra/medium, Grok 4.6/high) and lost twice: `std::abs( x ) <= INT_MAX` rejects `INT_MIN` (both auditors) — replaced by inclusive `min()..max()`; and an exponent bound of `INT_MAX` still lets `tPower * 2` / `* 3` / `* 4` overflow (Codex, Grok; the seven dimension accumulators are `real`, only those products are `int`) — capped at `|e| ≤ 32`. Grok added: `std::stod` throws on `m^abc` before any guard runs — replaced by the file's own `to_real()` (NaN on no-conversion); a fractional exponent was silently truncated (`m^1.5` → `m^1`) — now rejected as an input defect; the error text "not a valid integer" overclaimed since `get_int` still rounds — reworded. `StringTools.UnitToSiExponentGuard` locks the seven rejected spellings and three boundary cases.

## Changes Made

- `config/compiler/finalize_compiler.cmake` — third-party dirs normalized, deduplicated, `-isystem` unless in `CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES`; Fortran loop unchanged (`-I`, deduplicated list)
- `config/scripts/Add_Test.cmake` — per-test `BEFORE` SCLS include removed
- `config/compiler/config_gcc.cmake` — `-fallow-argument-mismatch` rationale comment (Christian's request)
- `src/numerics/spline/cl_Spline.{hpp,cpp}` — `create_helpmatrix( const index_t aSize, const real aDeltaX, … )`
- `src/physics/materials/cl_Material_Alloy.cpp`, `cl_Material_SplineLookupTable.cpp`, `cl_Material_Metal.cpp` — `static_cast< size_t >( std::ceil( … ) + 1 )` / `static_cast< int >( std::ceil( … ) )`; explicit `static_cast< index_t >` at the four non-`index_t` `create_helpmatrix` callers
- `src/numerics/integration/fn_intpoints.cpp` (4 sites), `src/io/cl_Input_Section.cpp` (1), `src/containers/cl_Genome.hpp` (2) — `static_cast` + `std::` qualification, `#include <cmath>` added to all three
- `src/core/stringtools.cpp` — `reserve( length + length / 2 )`; `1 : -1`; unit exponent via `to_real()`, `BELFEM_ERROR` unless finite, integral and `|e| ≤ 32`, then `static_cast< int >`
- `src/io/cl_Input_Section.{hpp,cpp}` — `get_int` rounds, then `BELFEM_ERROR` unless finite and within `int`; header doc updated; `src/io/doc/io_usage_guide.md:588` one-line note
- `tests/core/test_StringTools.cpp` — `UnitToSiExponentGuard`
- `CHANGELOG.md` — `[Unreleased]`: Changed ×2 (system headers; `create_helpmatrix` signature), Fixed ×1
- `tests/doc/tests_09_spline.md:58` — gotcha line updated to the new signature

## Open Questions

- **Gate passed:** `make check` 17/17 at `14bb703` (see Verification). Still to confirm from the same build: zero `warning:` lines from `src/` and `/opt/scls/` (the two gfortran rank-mismatch lines remain by design).
- **SemVer:** the `create_helpmatrix` signature change is a public API break in an installed header; by SemVer the next release carrying it would be 0.10.0 (Codex). **Christian's ruling (same day): the next release stays 0.9.2** — pre-1.0, static library by default, the CHANGELOG line that consumers must recompile is sufficient.
- **Two product choices in the exponent guard for Christian to veto:** the cap `|e| ≤ 32` (any bound ≪ INT_MAX/4 closes the overflow; 32 is far above any physical unit), and rejecting fractional exponents instead of truncating them (no deck or test in the tree uses one; supporting them would mean making `tPower` a `real`, a feature, not this fix).
- **Still unchecked, out of scope:** `Section::create_key` admits `inf` into the real keys (`!isnan`, not `isfinite`, `cl_Input_Section.cpp:161-171`), so `get_real` can still return inf; `get_int` now refuses it. No `Section` test fixture exists, so the `get_int` abort path has no test.
- **Unverified platforms:** the builtin-directory guard is exact string match against CMake's implicit list; a TPL path spelled differently (symlink, `..`) for the same directory would get `-isystem`. Ubuntu 24.04 / GCC 13 and Apple Clang not run here.
- **Coverage:** `make material` compiled all of `libbelfem` but no tests or other executables; `cl_Genome.hpp` was found by Grok reading the tree. Christian's `make check` build is the first full-tree pass with the flag.
- **Shared checkout:** a peer session's WIP (tapestack3d_quench docs, copper-Poisson devlog/todo) sits in the same tree, including `devlog/README.md`; stage explicit paths.

## Files Updated

- CHANGELOG.md
- config/compiler/config_gcc.cmake
- config/compiler/finalize_compiler.cmake
- config/scripts/Add_Test.cmake
- src/containers/cl_Genome.hpp
- src/core/stringtools.cpp
- src/io/cl_Input_Section.cpp
- src/io/cl_Input_Section.hpp
- src/io/doc/io_usage_guide.md
- tests/core/test_StringTools.cpp
- src/numerics/integration/fn_intpoints.cpp
- src/numerics/spline/cl_Spline.cpp
- src/numerics/spline/cl_Spline.hpp
- src/physics/materials/cl_Material_Alloy.cpp
- src/physics/materials/cl_Material_Metal.cpp
- src/physics/materials/cl_Material_SplineLookupTable.cpp
- tests/doc/tests_09_spline.md
- devlog/dl20260915_float_conversion_warnings.md, devlog/README.md
