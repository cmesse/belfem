# Devlog 2026-07-01 - Circuit Dependency Inversion Audit

**Date:** 2026-07-01
**Topic:** Read-only audit of whether `src/circuit/notes.md` dependency-inversion refactor note is still current
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Checked the circuit refactoring notes against the current FEM, Maxwell, executable, and CMake wiring. No source code was modified.

## Key Findings

- The core code-level dependency inversion has been implemented. `src/numerics/sources/cl_Circuit.hpp` defines the abstract `belfem::Circuit` interface, and `src/circuit/cl_ElectricalCircuit.hpp` implements `belfem::electronics::ElectricalCircuit : public Circuit`.
- FEM and Maxwell headers use the abstract `Circuit` interface, not the concrete `ElectricalCircuit` class. The main use sites are `src/fem/kernel/cl_FEM_Controller.hpp`, `src/fem/kernel/cl_FEM_Controller.cpp`, and `src/fem/maxwell/cl_MaxwellBoundaryConditionFactory.hpp`.
- Concrete circuit construction is isolated in the executable layer: `src/executables/hphirun.cpp` and `src/executables/hphiTrun.cpp` instantiate `electronics::ElectricalCircuitFactory` and pass its `Circuit *` to `Controller::set_circuit()`.
- The CMake optionality described in `src/circuit/notes.md` is still incomplete. `src/CMakeLists.txt` still adds `src/circuit` unconditionally, and `src/executables/CMakeLists.txt` still includes and links the `circuit` library unconditionally for the solver executables.

## Changes Made / Proposed

- Added this devlog entry only.
- Proposed interpretation: mark the source-level inversion note as mostly resolved, but keep a build-system task if the circuit module should become actually optional.

## Open Questions

- Should `USE_CIRCUIT` or a similar CMake option be added, with `hphirun`/`hphiTrun` compiled in a no-circuit mode when disabled?

## Files Updated

- `devlog/dl20260701_circuit_dependency_inversion_audit.md`
- `devlog/README.md`
- `tmp/ai_exchange/circuit_dependency_inversion.md`
