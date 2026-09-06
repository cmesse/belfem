# Devlog 2026-04-08 — Powerlaws Review

**Date:** 2026-04-08
**Topic:** Read-only review of `src/physics/materials/powerlaws.hpp`
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Reviewed `src/physics/materials/powerlaws.hpp` under the repository's read-only review rule. The main correctness problem is a contract mismatch between the documented "constant `jc`/`n`" setup and the inline implementation: the constant-value overloads still dereference `mJcFunction` and `mNFunction`, even though `set_constant()` leaves those pointers unset. I also found two follow-on API robustness issues: the `jc()` / `n()` helper accessors assume function-backed properties, and the piecewise model has no internal guard for mathematically singular parameter sets such as `n <= 1`.

## Key Findings

- `rho_powerlaw(const real normJ)` and `rho_piecewise(const real normJ)` dereference `mJcFunction` / `mNFunction` in the constant-property path at [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L26) and [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L195). That contradicts the documented and implemented storage model for constant HTS properties: `set_constant()` returns early for `jc` / `n` without creating helper functions in [cl_Material.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.cpp#L279), `JcFunction` explicitly documents that constants are stored directly in `Material` in [cl_JcFunction.hpp](/home/christian/codes/belfem/src/physics/materials/cl_JcFunction.hpp#L55), and the usage guide recommends `set_constant(MaterialProperty::jc, ...)` / `set_constant(MaterialProperty::n, ...)` in [materials_usage_guide.md](/home/christian/codes/belfem/src/physics/materials/doc/materials_usage_guide.md#L130) and [materials_usage_guide.md](/home/christian/codes/belfem/src/physics/materials/doc/materials_usage_guide.md#L893). As written, the documented configuration can null-dereference before any assertion fires.
- The public helper accessors `Material::n(...)` and `Material::jc(...)` make the same hidden assumption that the property is function-backed. `n()` unconditionally calls `mNFunction->eval(...)` at [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L687), and `jc()` only checks `have(MaterialProperty::jc)` before dereferencing `mJcFunction` at [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L694). Because `set_constant(MaterialProperty::jc, ...)` sets the `have(jc)` bit in [cl_Material.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.cpp#L235), `jc(B, angle, T)` can also null-dereference on a documented constant-`jc` setup.
- The piecewise branches rely on invariants that are never enforced locally. Every overload computes `1.0 / (n - 1.0)` and later divides by the Bezier coefficient `a`; representative sites are [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L231), [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L242), [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L915), and [powerlaws.hpp](/home/christian/codes/belfem/src/physics/materials/powerlaws.hpp#L927). There is no `BELFEM_ASSERT` in this file or elsewhere that rejects `n <= 1`, zero / negative defect-scaled `jc`, or degenerate transition geometry. In those cases the implementation will generate `inf` / `nan` values instead of a controlled runtime error.

## Changes Made / Proposed

- No source edits.
- Wrote this devlog to preserve the audit trail.

## Open Questions

- Should the constant-`jc` / constant-`n` support remain a public contract, or should the documentation be tightened so HTS power-law evaluation always requires function-backed `JcFunction` objects?
- For defect-aware overloads, should the API enforce `have_defect()` internally with `BELFEM_ASSERT`, or is the intended contract that callers must branch on `have_defect()` first?

## Files Updated

- /home/christian/codes/belfem/devlog/dl20260408_powerlaws_review.md
