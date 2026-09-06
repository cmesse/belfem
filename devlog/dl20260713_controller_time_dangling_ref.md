# Controller::time() dangling-reference fix

**Date:** 2026-07-13
**Purpose:** Fix a dangling reference behind a clang warning at `cl_FEM_Calculator.cpp:66`
**Module:** fem/kernel
**AIs involved:** Claude (lead), Codex (audit, high), Grok (third voice, high ~95%)

## Symptom

clang warned on the `MaxwellData` constructor initializer
(`src/fem/kernel/cl_FEM_Calculator.cpp:66`):

```cpp
mTime( aCalculator->group()->parent()->parent()->controller()->time() ),
```

## Root cause

`MaxwellData::mTime` is a **reference member** (`const real & mTime;`,
`cl_FEM_Calculator.hpp:133`), but `Controller::time()` returned `real`
**by value** (`cl_FEM_Controller.hpp:200-201`, `cl_FEM_Controller.cpp:69-73`,
body `return mTime;` — a copy of the live member `real mTime`,
`cl_FEM_Controller.hpp:42`).

A reference member bound to a by-value return aliases a **temporary** that dies
at the end of the constructor's full-expression → dangling reference; every
later read of `mTime` is UB.

The reference is *intended*: `mTime` must alias the controller's advancing time
so the calculator sees it change across timesteps. The sibling reference member
`mTimestep` (`cl_FEM_Calculator.hpp:341`) is bound correctly because its sources
return lvalue references to live members — `IWG::delta_time()` → `real &`
(`cl_IWG.hpp:736`) and `Mesh::time_stamp()` → `real &` (`cl_Mesh.hpp:1535`).

The user's initial hypothesis — that a `shared_ptr<Controller>` was being
converted to raw via `.get()` — was **refuted**. Kernel's raw back-pointer
`mController` (set via `set_controller(this)`) is a legitimate non-owning
back-reference; the `.get()` at `cl_MaxwellFactory.cpp:723` is on the *Kernel*
smart pointer, not the Controller. The warning is purely the dangling temporary.

## Fix

Make the accessors return a reference to their live member, matching the
`delta_time()`/`time_stamp()` pattern:

- `cl_FEM_Controller.hpp:200-204` — `time()` and `time_thermal()` now return
  `const real &`.
- `cl_FEM_Controller.cpp:69-79` — bodies unchanged (`return mTime;` /
  `return mTime2;`), now yielding an lvalue reference to the member.

`simulation_time()` was intentionally left by-value (a loop-bound getter, no
aliasing needed).

## Caller impact

All 68 `time()` call sites (`mt_maxwell_h.cpp` 48, `mt_thermal_h.cpp` 10,
`hphiTrun.cpp` 5, `hphirun.cpp` 2, `cl_MaxwellPostprocessor.cpp` 2, plus the
one intentional bind) and the single `time_thermal()` use
(`hphiTrun.cpp:163`) consume the result as a value; `const real &` decays
cleanly, so no call site breaks. Verified independently by Codex and Grok.

## Status

- [x] Root cause confirmed (Claude + Codex + Grok, all high)
- [x] Fix applied to `cl_FEM_Controller.{hpp,cpp}`
- [ ] Build/verify — handed off to user
