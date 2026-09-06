# User Material Backend Boundary Refactor

**Status:** CLOSED 2026-07-03 — Goals 1+2 implemented as `SplineLookupTable` and verified in the tree: `cl_Material.hpp` is fully backend-decoupled (no `cl_Spline.hpp`/`cl_Vector.hpp`/`cl_SpMatrix.hpp`), and Metal/Alloy/Magnesia are reparented to `SplineLookupTable`. **Deferred, non-blocking follow-ups (separate pass):** R7 (collapse harmless duplicate inline `spline_property`/`l`), R11/R12 (full build + material-value verification — only `-fsyntax-only` done so far), R13 (user-material CMake/`MatData` must import host BELFEM compile defs + TPL include paths), and O3 (linalg-import policy for the full plugin template). See `devlog/dl20260702_spline_lookup_table_refactor.md`.
**Scope:** `src/physics/materials`, focused tests. Build-template / `MatData` CMake work (R8-R9) is deferred to a separate pass.
**Non-goal:** Do not remove linalg access from advanced user materials; make the simple scalar/polynomial path backend-light.

## Confirmed Dependency Edge (2026-07-02 re-trace)

The hard linalg coupling in `cl_Material.hpp` is narrow:

- The **only** linalg-dragging include is `#include "cl_Spline.hpp"` (line 59). `cl_JcFunction.hpp` (line 58) is clean (`typedefs` + `Bitset`); `powerlaws.hpp` (line 1712) only re-includes `cl_Material.hpp`.
- `Vector<real>` is **already forward-declared** (line 63); every base signature naming it is a declaration and compiles backend-free. The backend only bites when a user *constructs* a `Vector` — which is why Goal 2 (polynomial API) is separable.
- `Cell<Spline*> mSplines` (line 283) needs only a **forward declaration** of `Spline`; pointer storage does not need the full type.
- What actually forces `cl_Spline.hpp` into the header is exactly **two inline bodies that dereference a spline**: `Material::spline_property()` (line 1689, `->eval`) and `Material::l()` (line 1698, `->integrate`). `Material::spline()` (line 1681) only returns the pointer, so it is fine.

**Spline API users** (`set_spline`/`create_spline`/`mSplines`/`spline_property`): `Metal` (+ its 7 children), `Alloy`, and **`Magnesia`** (a direct `Material` subclass — easy to miss). Only `UserDefinedMaterial` is spline-free. So `Material_LookupTable` must be the parent of **Metal, Alloy, and Magnesia**.

**Diligence hazards:**
1. Member-function-pointer contravariance: `mFunctionCp`/`mFunctionRho`/… are `real (Material::*)(real) const` and **private**. `&Material_LookupTable::cp_spline` is not assignable to them. Keep the `*_spline` dispatch stubs in `Material` (they only call virtual `spline_property`, so they stay backend-neutral) and add a protected `Material::set_spline_dispatch(property)` to perform the private pointer assignment (R3).
2. `Material::density_custom` (inline) calls `this->l(T)`, so `l` must stay declared+virtual in `Material` (body moves down; base version errors).
3. Spline ownership/destruction moves from `~Material` to `~Material_LookupTable`; the spline-touching bodies in `cl_Material.cpp` move to `cl_Material_LookupTable.cpp`.
4. `Alloy` and `Metal` re-declare `spline_property`/`l`; collapse the duplicates into `Material_LookupTable` where identical (Metal keeps any genuinely different override).

## Problem

Standalone user-defined material libraries currently include `cl_Material.hpp`. That header pulls in spline/linalg details, so the material compilation must know whether the host BELFEM build used `BELFEM_ARMADILLO` or `BELFEM_BLAZE`.

Observed dependency chain:

- `cl_Material.hpp` includes `cl_Spline.hpp`, stores `Cell<Spline*>`, and has inline methods that dereference `Spline`.
- `cl_Spline.hpp` includes `cl_Vector.hpp` and `cl_SpMatrix.hpp`, and exposes `Vector<real>`, `Matrix<real>`, and `SpMatrix`.
- `cl_SpMatrix.hpp` also exposes `Vector<real>` in multiply APIs/operators.
- `cl_Vector.hpp` defines `Vector` only under `BELFEM_ARMADILLO` or `BELFEM_BLAZE`.
- `set_user_defined_polynomial()` currently takes `const Vector<real>&`, so even a simple polynomial user material imports linalg.

The build-template side is also incomplete: the MatData material target imports include paths but not the host BELFEM compile definitions or third-party include paths.

## Working Decision

There are two valid user-material tiers:

- **Simple scalar material:** constants, scalar property functions, and polynomial coefficients should not require linalg.
- **Full BELFEM C++ material:** may use `Vector`, `Matrix`, splines, and other BELFEM internals, but must compile against the exact host BELFEM build configuration.

The immediate refactor should make the simple path backend-light without pretending the full C++ plugin ABI can be backend-free.

## Proposed Shape

Introduce an intermediate lookup-table class:

```cpp
class Material;
class UserDefinedMaterial : public Material;
class Material_LookupTable : public Material;
class Metal : public Material_LookupTable;
class Alloy : public Material_LookupTable;
```

Responsibilities:

- `Material`: scalar property dispatch, constants, dependencies, user-function registration, no spline storage, no `Vector` in the primary user-facing polynomial API.
- `UserDefinedMaterial`: stores polynomial coefficients in `Cell<Cell<real>>`, evaluates them directly with Horner's method.
- `Material_LookupTable`: owns `Cell<Spline*>`, implements `set_spline()`, `spline()`, `spline_property()`, and spline-backed `l(T)` behavior.
- `Metal` and `Alloy`: inherit lookup-table support and keep their spline/table-heavy implementation there.

## `set_spline()` Migration

Do not expose `Material`'s private dispatch fields to subclasses. Split today's `Material::set_spline()` into two parts:

- `Material::set_spline_dispatch(MaterialProperty)` - protected, backend-light; clears constants/dependencies, marks property available, and sets `mFunctionCp = &Material::cp_spline`, etc.
- `Material_LookupTable::set_spline(MaterialProperty, Spline*)` - protected; owns/deletes/replaces the spline pointer, then calls `set_spline_dispatch()`.

Keep `Material::cp_spline()`, `rho_spline()`, etc. as thin dispatch stubs that call a virtual hook:

```cpp
virtual real spline_property(MaterialProperty aProperty, real aX) const;
```

Base `Material::spline_property()` should error. `Material_LookupTable` overrides it and evaluates the owned spline.

## Polynomial API

Prefer a backend-light primary API:

```cpp
virtual void
set_user_defined_polynomial(
    MaterialProperty aProperty,
    const Cell<real> & aCoefficients );
```

`UserDefinedMaterial::evaluate_polynomial()` should use Horner's method directly on `Cell<real>`, avoiding `fn_polyval.hpp` and `Vector`.

Compatibility options:

- Add a `Vector<real>` overload only in linalg-aware code, implemented by copying into `Cell<real>`.
- Or leave `Vector<real>` as an advanced API and document that it requires the host BELFEM build configuration.

Implementation note: `Cell<T>(n)` reserves but does not resize; use `set_size()` when creating indexed `Cell<Cell<real>>` storage.

## Work Plan (ordered to keep the build green at each phase)

**Implemented 2026-07-02** as `material::SplineLookupTable` (Christian's initial cut + Claude edge-polish). All material sources and downstream users (mesh, fem/kernel, fem/maxwell) pass `g++ -fsyntax-only` with real build flags; a user-material TU compiles with **no** backend macro.

### Phase A — introduce the class as a pass-through (no behavior change)

- [x] R1. Create `cl_Material_SplineLookupTable.{hpp,cpp}` with `class SplineLookupTable : public Material`, forwarding the constructor to `Material(type, isotropic)`. `#include "cl_Spline.hpp"` lives here, not in `cl_Material.hpp`.
- [x] R2. Reparent `Metal`, `Alloy`, and `Magnesia` from `: public Material` to `: public SplineLookupTable` (headers **and** ctor init lists: `cl_Material_Metal.cpp:31`, `cl_Material_Alloy.cpp:38,53`, `cl_Material_Magnesia.cpp:24`).

### Phase B — move spline storage + logic down (the real surgery)

- [x] ~~R3. Add protected `Material::set_spline_dispatch(MaterialProperty)`~~ — superseded: `Material` declares `friend class material::SplineLookupTable`, so `SplineLookupTable::set_spline` assigns the private dispatch pointers directly. Note: friendship is not inherited, so the `mFunctionRhoKohler…mNFunction` group stays `protected:` for Metal/UserDefined.
- [x] R4. Move `mSplines`, `set_spline`, `spline`, spline destruction into `SplineLookupTable`; `create_spline(property, dYdX0, dYdX1)` stays in `Material` (protected) and dispatches to the virtual `create_spline(fptr, property, BCs, dYdX0, dYdX1)` implemented in `SplineLookupTable`. A `using Material::create_spline ;` in `SplineLookupTable` un-hides the convenience overload. `extend_alpha_to_zero` stays in `Material` (declaration-only `Vector<real>` use, backend-free).
- [x] R5. `Material::spline_property` / `dspline_property` / `reset_spline` / `l` are virtual-and-erroring in `cl_Material.cpp`; real impls in `SplineLookupTable` (inline in its header). `Material::l` stays declared so inline `density_custom` compiles.
- [x] R6. `cl_Material.hpp` no longer includes `cl_Spline.hpp` (only `Spline_Enums.hpp`, which is backend-clean). **Checkpoint passed:** parses with no backend macro.
- [ ] R7. Collapse the duplicated inline `spline_property`/`l` in `Alloy.hpp` (and `Metal`'s `spline_property`) into `SplineLookupTable`; keep only genuinely different overrides. *(Deferred cleanup — duplicates are harmless, they compile against inherited `mSplines`.)*

### Phase C — backend-light polynomial API (Goal 2)

- [x] R8. `set_user_defined_polynomial` now has `Cell<real>` and `std::vector<real>` overloads (the `Vector<real>` overload was removed; no external callers existed). `std::vector` is the documented user-facing form in `example_user_material.cpp`.
- [x] ~~R9. Convert `mPolynomials` to `Cell<Cell<real>>` + direct Horner~~ — superseded: storage stays `Cell<Vector<real>>` + `polyval` *inside* `UserDefinedMaterial` (BELFEM-internal, always compiled with a backend); the backend-neutral conversion happens at the API boundary via `Vector<real>( Cell/std::vector )` ctors (exist in both backends: `cl_AR_Vector.hpp:94,102`, `cl_BZ_Vector.hpp:108,116`).
- [x] ~~R10. Remove `cl_Vector.hpp` from `cl_Material_UserDefined.hpp`~~ — not needed: users include `cl_Material.hpp` only; `UserDefinedMaterial` is internal to BELFEM.

### Phase D — verify

- [ ] R11. Build serial and run existing material tests; spot-check Metal/Alloy/Magnesia spline values are unchanged. *(Syntax-only verified so far — full build/tests pending, Christian runs builds.)*
- [ ] R12. Add a focused test to the repo: a simple polynomial user material that compiles against `cl_Material.hpp` with **no** `BELFEM_ARMADILLO`/`BELFEM_BLAZE` defined. *(Verified ad-hoc 2026-07-02 with a scratch TU + real include flags minus the backend define; not yet a checked-in test.)*

### Deferred (separate pass — was R8-R9 of the original plan)

- [ ] R13. Fix user-material CMake/`MatData` to import the host BELFEM compile definitions and TPL include paths (header decoupling alone does not give a working standalone build — see `dl20260701_user_material_backend_boundary.md`).

## Open Questions

- [x] O1. Final class name: **`SplineLookupTable`** (in `namespace material`).
- [x] O2. `Material::l(T)` stays in base, virtual and erroring; `SplineLookupTable` provides the spline-integrated implementation.
- [ ] O3. Should the full user material template always import host BELFEM linalg settings, even after the simple scalar path is backend-light?

## Definition of Done

- `cl_Material.hpp` can be included by a simple user-defined material without including `cl_Spline.hpp`, `cl_Vector.hpp`, or `cl_SpMatrix.hpp`.
- Existing built-in metal/alloy spline behavior is unchanged.
- Simple polynomial user materials can use `Cell<real>` coefficients without selecting Armadillo or Blaze in their own source.
- Advanced user materials still have a documented path to use BELFEM linalg, compiled against the host BELFEM configuration.
