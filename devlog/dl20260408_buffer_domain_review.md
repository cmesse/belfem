# Devlog 2026-04-08 — Buffer Domain Review

**Date:** 2026-04-08
**Topic:** Read-only review of the new `DomainType::Buffer` plumbing in `en_DomainType`, `ThinShellFactory`, and Maxwell dispatch
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Reviewed the new `Buffer`-domain draft as a code review. The overall direction is reasonable: tag true buffer layers separately, avoid ghost interfaces against them, and route them through scalar-`phi` assembly rather than conductor thin-shell kernels. The current implementation is not yet runnable, however. There are compile-breaking syntax errors in touched Maxwell files, and the core Maxwell block plumbing still omits `DomainType::Buffer` in the two places that determine whether a block enters the magnetic equation and receives air-like DOFs.

## Key Findings

- Compile blocker: several touched Maxwell files currently contain invalid syntax, e.g. [cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L35), [cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L248), [cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L578), [cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L678), [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L51), [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L529), and [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L1924).
- Functional blocker: `Buffer` blocks are not selected into the magnetic kernel in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L443), and the block-update switch also omits `DomainType::Buffer` in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L499). That means the new `IWG_Maxwell` `case DomainType::Buffer` in [cl_IWG_Maxwell.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_IWG_Maxwell.cpp#L252) will never be reached for buffer blocks yet.
- Functional blocker: `FieldList::collect_block_dofs()` still has no `DomainType::Buffer` case. The block DOF switch in [cl_Maxwell_FieldList.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_Maxwell_FieldList.cpp#L253) maps `Conductor`, `ThinShell`, `Coil`, `Ferro`, and `Air`, but not `Buffer`. Even if buffer blocks were selected, they would not yet receive the air-like `phi` DOF table.
- Runtime blocker: buffer detection in [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L1689) is based on the literal layer label being exactly `"buffer"`. In this factory, `aMaterials` is just the first token from each layer line in the input file, read in [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L96) and [cl_ThinShellFactory.cpp](/home/christian/codes/belfem/src/mesh/cl_ThinShellFactory.cpp#L105). So labels like `buffer1` or any custom material name will not be tagged as `Buffer`.
- Runtime blocker: `mMaterialBlockAssignment` still assigns the original layer label to every generated thin-shell block in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L896). But `create_and_assign_materials()` only loads kernel materials for `Conductor`, `ThinShell`, and `Ferro` blocks in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L2368). Later it still calls `tBlock->set_material(label)` for any non-`air` label in [cl_MaxwellFactory.cpp](/home/christian/codes/belfem/src/fem/maxwell/cl_MaxwellFactory.cpp#L2406), and `Kernel::material()` hard-errors if the label was never loaded in [cl_FEM_Kernel.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Kernel.cpp#L713). So buffer blocks can currently trip a material-lookup failure unless that path is special-cased.
- Positive finding: the scalar `phi` operator on `PENTA6TS` is viable in principle. `PENTA6TS` already maps to the regular 6-node prism Lagrange interpolation in [cl_IF_InterpolationFunctionFactory.cpp](/home/christian/codes/belfem/src/fem/interpolation/cl_IF_InterpolationFunctionFactory.cpp#L170), and the calculator uses `Bscalar` plus `dV_ts` on the scalar path in [cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L533), [cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L540), and [cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L661). So the draft is blocked by plumbing, not obviously by a missing prism-scalar operator.

## Changes Made / Proposed

- No source edits.
- Wrote this devlog to record the current blockers before any implementation proceeds.

## Open Questions

- Should buffer identification be driven by explicit layer metadata / domain assignment rather than by the literal material label?
- Should `Buffer` blocks be treated exactly like `Air` for material assignment too, or do you still want a lightweight material object attached for postprocessing / labeling?

## Files Updated

- /home/christian/codes/belfem/devlog/dl20260408_buffer_domain_review.md
