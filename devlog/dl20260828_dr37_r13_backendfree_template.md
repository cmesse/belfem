# DR-37 closed, then the templates audited and repaired

**Date:** 2026-08-28
**Purpose:** close DR-37 entirely. (b)/R13 by deleting the dead backend-detection block in
`UserMaterialTemplate.cmake` and running its named standalone-plugin gate; (a)/R12 by executing
the owed in-tree compile gate in an independently configured tree, on both linear-algebra
backends, with the guards proven non-vacuous. Also retracts a false positive the same dead code
had produced in `shared_library_and_install_plan` O9
**Module:** `src/physics/materials`, `todo/`

## What the row asked for, and why the answer was deletion

DR-37 (b) recorded `UserMaterialTemplate.cmake` backend detection as dead code
(`BELFEM_CACHE_LOCATIONS` set at `:69`, never read; no `target_compile_definitions`, no TPL
includes) and framed the fix as *implementing* the missing "template-local compile-definition
import".

That framing is wrong, and R12 is what makes it wrong. The user-material API is backend-free by
construction:

- `cl_Material.hpp:24-25` carries the standing instruction: *"NEVER include cl_Vector.hpp or
  cl_Matrix.hpp here, or any class that uses it. Doing so would break the API for the user
  defined materials."*
- No `Vector<` or `Matrix<` appears anywhere in the header. Both
  `set_user_defined_polynomial` overloads take `Cell<real>` / `std::vector<real>`
  (`:1271`, `:1276`). `material::BhCurve` (`:90`, `:508`, `:1122`) appears only as a
  forward-declared pointer.
- R12's own gate, `tests/physics/backendfree/`, compiles that closure with the backend macro
  stripped and `#error`-guards both the backend macros and the linalg include guards.

So there is no backend for a plugin to detect, and nothing for a compile definition to import.
Implementing the block would have added machinery for a requirement the tree actively forbids.
It is removed and replaced by a comment that states the contract and names the gate that pins
it, so the block is not reintroduced by the next reader who notices the absence.

## R13's gate was executed, not reviewed

R13's reproducer column named "standalone plugin build". It had never been run. It has now,
against both the pre-edit and post-edit template, in a scratch tree (no shared build tree
touched):

```
cp src/physics/materials/UserMaterialTemplate.cmake   <scratch>/CMakeLists.txt
cp src/physics/materials/example_user_material.cpp    <scratch>/myalloy.cpp
cmake -DBELFEM_DIR=/home/christian/codes/belfem ..    # exit 0, no warnings
make                                                  # links libmyalloy.so
```

Result, identical before and after the edit:

| check | result |
|---|---|
| configure | exit 0, no errors or warnings |
| compile + link | `libmyalloy.so`, 19152 bytes |
| exported init symbol | `T MyAlloy_init` |
| Blaze / Armadillo symbols | **0** |
| host symbols required | exactly one: `belfem::Material::set_constant` |

The single undefined BELFEM symbol is the point. Everything else on the plugin's surface
(`set_user_defined_function`, `rho`) is header-inline, and none of it reaches a backend type.

This moves R13 from *reviewed* to *verified* on the register's evidence ladder: a named gate
that can observe the defect actually ran.

## R12: the in-tree gate, executed

R12 landed 2026-08-27 as `tests/physics/backendfree/` but its status read *"reviewed, not
verified, until Christian's configure + `make` compiles it in the shared tree."* That gate has
now run, in a tree configured from scratch under the scratchpad rather than the shared one, with
options matched to `cmake-build-debug` (Debug, `USE_TEST=ON`, `USE_SHARED_LIBS=ON`,
`USE_MKL=ON`, `SCLS=/opt/scls/mkl`, system GCC 11.5 per the Linux compiler policy).

### 1. The strip is exact

The whole design rests on `remove_definitions()` taking away the backend macro and nothing else.
Checked in `flags.make`, not `compile_commands.json`. Diffing the gate target's `CXX_DEFINES`
against a normal sibling target (`test_physics`) gives **exactly one** entry of difference:

| tree | difference |
|---|---|
| Armadillo | `-DBELFEM_ARMADILLO` only |
| Blaze | `-DBELFEM_BLAZE` only |

`ARMA_DONT_USE_SUPERLU`, `BELFEM_{ARPACK,EXODUS,GCC,HDF5,METIS,MKL,MPI,MUMPS,NLOPT,PARDISO,
PARMETIS,PARPACK,PETSC,PTSCOTCH,SCOTCH,STRUMPACK,SUPERLU}`, `DEBUG` and `OMP` all survive, which
is exactly what the directory's CMakeLists claims. Note `-DARMA_ALLOW_FAKE_GCC/CLANG` never
appears in either tree: it only accompanies Intel builds, so that clause of the
`remove_definitions()` call is a no-op here and remains untested.

### 2. It compiles

`make test_material_backendfree` exits 0 in both trees and emits both objects
(`test_MaterialBackendFree.cpp.o`, `example_user_material.cpp.o`).

### 3. It is not vacuous

This is the part that mattered, and the part a passing build cannot tell you. Both `#error`
guards were driven, using the real in-tree compile command lifted from `flags.make`:

| probe | change | result |
|---|---|---|
| control | none | exit 0 |
| guard 1, Armadillo tree | append `-DBELFEM_ARMADILLO` | `#error "backend define leaked into the backend-free material gate"` |
| guard 1, Blaze tree | append `-DBELFEM_BLAZE` | same `#error` |
| guard 2 | inject `#include "cl_Vector.hpp"` after the `cl_Material.hpp` include | `#error "cl_Material.hpp pulled a linalg/spline header into the backend-free material gate"` |

Guard 1's probe also **confirms the row's own accepted residual**: a `-D` arriving through
`CMAKE_CXX_FLAGS` does bypass `remove_definitions()`, and guard 1 catches it anyway. That was
written as a claim; it is now a measurement.

The other accepted residual stands and is inherent: a direct `#include <armadillo>` or
`<blaze/...>` in a plugin bypasses both guards, because no guard covers it.

### 4. The three wiring claims hold

Read out of the generated `CMakeFiles/Makefile2` rather than inferred:

```
all -> tests/all -> tests/physics/all -> tests/physics/backendfree/all
                                      -> test_material_backendfree.dir/all
```

plus explicit edges `CMakeFiles/check-fast.dir/all -> …test_material_backendfree.dir/all` and
`tests/physics/CMakeFiles/test_physics.dir/all -> …test_material_backendfree.dir/all`. So ALL,
`make check` and `make check-fast` each pull the gate, as claimed.

### What was not run

The full suite. Only the gate target was built, on a one-core box; nothing else in the tree was
compiled, and `make check` was not executed. The gate is compile-only by design (OBJECT library,
never linked, never run), so building the target *is* the whole test.

## The false positive this dead code had already caused

`shared_library_and_install_plan` O9 (found 2026-08-27) inferred from the same block that a
plugin "is compiled without `BELFEM_ARMADILLO`/`BELFEM_BLAZE` and silently disagrees with the
framework about the matrix ABI."

Refuted. There is no matrix ABI on this surface to disagree about, and the probe shows zero
backend symbols in the produced `.so`. The claim is struck in place in O9 with the retraction
written beside it, per the register's standing trap 1 ("Retract a false positive in the row; do
not quietly delete it"), rather than removed.

**O9's other half survives and is untouched:** both plugin templates put `${BELFEM_DIR}/include`
alone on the include path while R6 installs headers as `include/belfem/<module>/…`, so no plugin
can compile against an installed tree. That is a genuine defect, it is unrelated to the backend
question, and it stays with R11. O9 therefore stays unticked.

## State after this session

- **DR-37 struck.** Both halves are finished with executed gates, so the ID and description are
  struck and the status column keeps the evidence, per the register's convention. The `[P]`
  count in the preamble drops 24 → 23. The row is not yet moved to `debt_register_closed.md`;
  archiving has been running a day behind striking and is left to the next sweep.
- One deviation worth Christian's eye: R12's reproducer said *Christian's `make`*, and this ran
  in an independently configured tree instead. The options were matched to `cmake-build-debug`
  and the gate was additionally exercised on Blaze, which the shared tree does not cover, so the
  evidence is broader than the row asked for but it is not literally the shared tree.
- `todo/shared_library_and_install_plan.md` O9: dead-code half closed and its ABI inference
  retracted; include-path half open, still with R11.

## Not done

- **No three-vendor audit round.** The standing rule is plan+audit → code+audit with both Codex
  and Grok for code work. The only source change this session is a deletion plus a comment in a
  template that is not part of any build, backed by an executed gate, so the round was skipped
  pending Christian's call. Raised in chat.
- **`make check` was not run**, and no full build was attempted.
- The `-DARMA_ALLOW_FAKE_GCC` / `-DARMA_ALLOW_FAKE_CLANG` clause of the gate's
  `remove_definitions()` is untested: those defines only exist on Intel builds.
- The template's include-directory list (`src`, `core`, `containers`, `math/graph`, `io`,
  `physics/materials`) is wider than the three directories R12's gate proves are needed
  (`core`, `containers`, `physics/materials`). Harmless, out of scope for this row, not changed.

## Three-AI round on the template, and the repairs that followed

With DR-37 closed, the edited template went through a `--jury` cross-review
(`tmp/ai_exchange/review_user_material_template.md`, pre-registration frozen before dispatch).
Both auditors endorsed the R13 deletion — Grok unprompted: *"Not adding
`BELFEM_ARMADILLO`/`BELFEM_BLAZE`. That comment is the load-bearing part of DR-37 R13. Re-adding
detection would be the regression."* Then they took the rest of the file apart.

Two findings were promoted from reasoning to measurement during verification, and both were
Grok's, single-raiser:

**The template's own documented workflow never worked.** `BELFEM_DIR` was
`set(... CACHE PATH ...)` without `FORCE`, so the placeholder written on the first configure
shadowed every later edit. Edit the file as the USAGE block instructs, reconfigure, and you get
the same error still quoting `/path/to/belfem`. This also indicts the R13 gate run earlier the
same day: it passed only because it used `-DBELFEM_DIR=`, the one workaround the file never
mentions. Grok's suggested fix — empty `CACHE PATH` plus a not-set check — was applied, tested,
and **failed**: an empty cached value is still cached. `BELFEM_DIR` is now not a cache variable
at all.

**Including `cl_Material.hpp` silently disabled warnings for the rest of any plugin TU.** The
header pushes diagnostics on `__clang__`/`__GNUC__` and pops on `BELFEM_CLANG`/`BELFEM_GCC`. A
plugin defines neither, so it pushed and never popped. Measured under the template's own
`-Wall -Wextra`: plugin 0 unused-parameter warnings, in-tree 2, control 2. Fixed in the header
rather than by defining the identity macro in the template, because the header is the root cause.
In-tree behaviour is provably unchanged (with `BELFEM_GCC` defined the old pop already fired).
The Intel branch also had `#pragma warning pop` against `#pragma warning(push)`; corrected,
untested.

### The sibling template had the mirror-image defect

`UserLibraryTemplate.cmake` was outside the review's scope but shares every structural defect.
Fixing it surfaced what the round did not see: its include set reaches `src/linalg`, and
`cl_Vector.hpp:16-20` is `#ifdef BELFEM_ARMADILLO / #elif BELFEM_BLAZE`. Without a backend macro
those headers collapse to forward declarations, so that template could never have compiled
anything using `Vector`/`Matrix`. It was missing precisely the `target_compile_definitions` and
TPL includes DR-37 R13 originally asked for — **the dead `BELFEM_CACHE_LOCATIONS` block had been
in the wrong file all along.** It now takes a validated `BELFEM_BACKEND`, sets the define, adds
the backend wrapper directory, and exposes a `BELFEM_TPL_INCLUDE_DIRS` hook. A source using
`Vector<real>` builds under both backends.

### Also fixed

Installed-tree include path (both templates, header-probe rather than directory-probe, one shared
`BELFEM_HEADER_ROOT`); `SHARED` → `MODULE` for a dlopen-only artifact; Windows `.dll` claim
removed since both loaders are POSIX `dlopen`; `GNUInstallDirs` for the install destination; four
false comments; `CMAKE_BUILD_TYPE` defaulted to Release; the ABI-relevant defines (`NDEBUG`/`DEBUG`,
`BELFEM_INT64`) documented.

**Not applied:** Codex's `CMAKE_CXX_EXTENSIONS OFF` — the host leaves it unset too, so it would
create the divergence it appears to prevent. Linking `libbelfem` stays R11's design question.

**Codex fabricated a citation.** Its ABI finding cited `doc/lessons_learned_evidence.md:408` for a
plugin landing in wrong vtable slots "after eleven insertions". Line 408 is INC-316 about `ybco`
Buffer classification; neither phrase occurs in that 922-line file. The finding stands on its own
trace, the precedent does not exist, and DR-132 records this so it is not repeated.

### What is gated and what is not

Run: both `BELFEM_DIR` routes including a dirty-build-dir reconfigure, the installed layout
against 626 staged headers, MODULE symbols, the library template on both backends, the warning
leak, the in-tree R12 gate in both trees with both guards re-proven, CMake 3.31.8 and 4.4.2.

Not run: **macOS** — the `-Wl,-undefined,dynamic_lookup` fix and the `MODULE` rename are reviewed
only, and the rename is the risky half because it changes the artifact name on Darwin. That is
**DR-131**. Intel's `#pragma warning(pop)` is likewise untested. `make check` was not run.

O9 is closed in `shared_library_and_install_plan`; DR-131 and DR-132 filed.
