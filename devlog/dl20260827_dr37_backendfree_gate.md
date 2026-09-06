# Devlog 2026-08-27 — DR-37/R12: Checked-in Backend-Free Compile Gate

**Date:** 2026-08-27
**Topic:** Promote the ad-hoc backend-free compile check from the 2026-07-02 SplineLookupTable
refactor into a checked-in build gate (`tests/physics/backendfree/`)
**AIs involved:** Claude (plan, implementation, probes); Codex + Grok (plan audit and code audit,
`tmp/ai_exchange/dr37_backendfree_gate.md`)
**Claude Confidence:** high on the gate mechanics (four executable g++/CMake probes); the in-tree
build gate is still owed
**Literature References:** N/A (build-system gate)

## Summary

DR-37 bundles two residuals from `dl20260702_spline_lookup_table_refactor.md`: R12 (check in the
backend-free compile check) and R13 (standalone user-material CMake import). This session lands
R12 only; R13 stays open, split between `todo/shared_library_and_install_plan.md` G10/O4 and a
template-local half those steps do not cover. Christian approved R12 explicitly; the 2×2-element
runtime-test idea floated mid-session was withdrawn.

Along the way: `src/physics/materials/example_user_material.cpp` — the shipped plugin example —
was in no CMake sources list at all. It is now compiled by the gate, so it can no longer rot
silently.

## What landed

- **`tests/physics/backendfree/test_MaterialBackendFree.cpp`** (new): compile-only TU.
  Guard 1 `#error`s if `BELFEM_ARMADILLO`/`BELFEM_BLAZE` reaches the TU (the gate must fail
  loudly when its premise breaks, not pass vacuously). Guard 2, after
  `#include "cl_Material.hpp"`, `#error`s if the include-guard macro of `cl_Vector.hpp`,
  `cl_Matrix.hpp`, `cl_Spline.hpp`, `cl_SpMatrix.hpp`, or the native wrappers
  (`cl_BZ_Vector/Matrix.hpp`, `cl_AR_Vector/Matrix.hpp`) is defined. The body exercises the
  user-material surface: `MatFunc1` property functions, `set_constant`,
  `set_user_defined_function`, and the backend-neutral
  `set_user_defined_polynomial( property, std::vector<real> )`. No `main`, external linkage
  (avoids `-Wunused-function` under `-Wall -Werror`).
- **`tests/physics/backendfree/CMakeLists.txt`** (new): `remove_definitions` of the backend
  selectors (+ `ARMA_ALLOW_FAKE_*`, which only ever accompany them) in its own directory scope —
  doing it in `tests/physics/` would strip `test_physics`, whose YBCO/Jc tests need the backend.
  Explicit `include_directories` for `core`, `containers`, `physics/materials`; OBJECT library
  over the gate TU and the example; `add_dependencies( test_physics ... )`.
- **`tests/physics/CMakeLists.txt`**: `add_subdirectory( backendfree )` after the
  `USE_GASMODELS` block, so the gate exists with gas models off.
- **`CMakeLists.txt`**: `TARGET`-guarded `add_dependencies( check-fast
  test_material_backendfree )` after the check-fast foreach — that foreach omits `physics`, so
  without the hook the documented fast loop would never build the gate.

Coverage: plain `make` (OBJECT libraries are in ALL), `make check`, `make check-fast`. A named
target like `make hphirun` does not build it; a materials-header recoupling therefore fails the
next full build, which is the intended contract-gate behaviour.

## Why the guards are shaped this way

A bare `#include "cl_Vector.hpp"` without a backend define **compiles clean**: the header defines
its include guard, skips the backend wrapper, and the unconditional operator headers only carry
uninstantiated templates over the forward-declared `Vector`. So macro-absence alone is not a
detector — the include-guard `#error`s are the part that catches a vacuous recoupling. (Grok
initially claimed the opposite in the plan round and retracted it in the code round after
re-reading `op_VectorPlus.hpp`.) The native wrappers are the converse hole: they include
Blaze/Armadillo with **no** macro gate and sit on the inherited include path, so their guards are
checked too.

Accepted residuals, deliberately not claimed: a backend `-D` injected through
`CMAKE_CXX_FLAGS` bypasses the strip (guard 1 still fires); a direct
`#include <blaze/Math.h>`-style include bypasses both guards.

## Audit trail

Three-round exchange in `tmp/ai_exchange/dr37_backendfree_gate.md`: pre-registered plan, Codex +
Grok plan audits (both "revise"), reconciliation, implementation, Codex + Grok code audits (both
"accept"/"do not block"). Substantive audit catches, all adopted:

1. **check-fast wiring claim false** (Codex + Grok): the check-fast foreach omits `physics`;
   fixed with the explicit top-level hook.
2. **Include-path provenance false** (Grok): `Add_Test.cmake`'s `${SSF_SRC_DIR}` paths are junk
   (variable set nowhere); the working `-I`s leak from the banner executable's top-level
   `include_directories`. The gate sets its own includes and was proven self-contained by
   compiling both TUs with only the three project `-I`s.
3. **`remove_definitions` over directory-property surgery** (Grok): adopted, verified in a
   scratch mock.
4. **`ARMA_DONT_USE_SUPERLU`** (Codex wanted it stripped, Grok refused as scope creep): Grok's
   position adopted — it is not a backend selector and is inert here.
5. **Register wording** (Grok): "implemented, compile gate owed" — not "closed".

Post-audit `git status` checked after each Grok round: no writes outside the audited scope.

## Verification status

**Reviewed, not verified.** Executable probes ran (real g++ with the `libbelfem_materials`
`flags.make` flags): gate TU clean with the backend define stripped; `#error` fires with it
present; both TUs compile with only the three explicit `-I`s; `remove_definitions` semantics
confirmed in a scratch CMake project. What has NOT run is the real-tree configure + `make` —
that is Christian's gate and is recorded as owed in the DR-37 row.

## Files Updated

- `tests/physics/backendfree/test_MaterialBackendFree.cpp` (new)
- `tests/physics/backendfree/CMakeLists.txt` (new)
- `tests/physics/CMakeLists.txt`
- `CMakeLists.txt`
- `todo/debt_register.md` (DR-37 status cell)
