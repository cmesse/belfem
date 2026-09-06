# FEM Module — Cross-Cutting Documentation {#fem_index}

**Date:** 2026-08-09
**Purpose:** Documentation that spans the fem submodules (`kernel`, `iwg`, `maxwell`, ...).
Submodule-specific documentation lives in `src/fem/<submodule>/doc/`.

| Document | Content |
|---|---|
| [timestepping_strategy.md](timestepping_strategy.md) | Time integration schemes (BDF1–5, ratio clamps, startup ramp), the PID step-size controller and its legacy fallback, failure/floor policy, warm-start persistence of Δt AND the BDF order/step-size history (2026-08-15), literature (PI/PID step control) |

Related:
- `src/fem/kernel/doc/nonlinear_controller_theory.md` — the inner nonlinear loop
  (Picard/Newton hybrid, relaxation adaptation, guards)
- `doc/input_file_reference.md` — deck keys
