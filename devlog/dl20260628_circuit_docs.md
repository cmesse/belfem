# Devlog 2026-06-28 — Circuit module usage guide + ngspice parser plan

**Date:** 2026-06-28
**Topic:** Document the `src/circuit` module (driven from `hphirun` via
`ElectricalCircuitFactory`); propose an ngspice/SPICE netlist importer and a file-driven
standalone circuit runner.
**AIs involved:** Claude, Codex, Grok
**Claude Confidence:** high (post-audit)
**Codex Audit Confidence:** high (A–G), medium-high (SPICE completeness)
**Grok Audit Confidence:** high (A–E, G, H core), medium (F, parser completeness)
**Literature References:** none (infrastructure/documentation task); `tmp/ngspice-manual.pdf`
for SPICE netlist syntax.

## Summary

Traced the circuit module from `ElectricalCircuitFactory`
(`src/executables/hphirun.cpp:60,78`). Produced two AI+human documents and had both Codex
and Grok independently audit them; applied every confirmed correction.

- **Usage guide:** `src/circuit/doc/circuit_usage_guide.md` (+ `src/circuit/doc/README.md`).
- **Parser plan:** `todo/ngspice_parser_plan.md` (+ entry in `todo/README.md`).

No source code was modified (read-only investigation + docs only).

## Key Findings (circuit module)

- **MNA solver.** Unknowns are node voltages (N−1) plus branch currents of
  "unknown-current" elements. The **last** node index is ground (excluded from the system;
  voltage implicit 0). `cl_ElectricalCircuit.cpp:31-46,536-537,856`.
- **Unknown-current elements are only voltage sources and switches**
  (`cl_ElectricalCircuit.cpp:290-334`). Voltage sources stamp the constant MNA; switches
  stamp in `compute_jacobian_and_rhs()` (`:743`). Current sources and terminal pairs are
  RHS-only (`:728`).
- **Nine components** (`Component_Enums.hpp:21`): resistor, capacitor, inductor, voltage
  source, current source, switch, diode, superconductor, terminal pair.
- **Terminal pair = FEM coupling.** `FEMTwoTerminals` relays the FEM terminal current; it
  keeps I/V/h `ShiftRegister`s but runs **no** BDF companion model and its `order` key is
  **never read** (factory parses `order` only for L/C). Factory always tags these BCs
  `CircuitVoltage` (`cl_ElectricalCircuitFactory.cpp:481`); the `Voltage`/`length` 2-D
  branch (`:523`) is dead code. `mCircuitCurrentBCs` is therefore never populated, so the
  controller's `current()`-fix loop is inert for the standard path.
- **Curve linkage.** Terminal-pair `input curves : [..]` reference the logical curve ids in
  the FEM `topology{curves{}}` block (`tmp/examples/Circuit_Coupling/input.conf`); 2-D uses
  `input terminals` against sideset ids. When output groups are absent the factory fills
  them *from the inputs* and concatenates (inputs duplicated; `set_domains` does not dedup).
- **Rank-0 only.** The circuit solves serially on rank 0; only the success flag is
  broadcast (`cl_FEM_Controller.cpp:131-161`). `set_circuit()` is at `:1749`.
- **A FEM-free runner already exists:** `src/executables/electricalCircuit.cpp` (hardcoded
  circuit, own Newton + adaptive-Δt loop). It recomputes the MNA matrix only when Δt changes
  and `shift()`s at the *end* of the step — opposite phasing to the controller, which
  `shift()`s before solving.

## Audit corrections applied

- **Claim F overstated (both auditors):** `solve_circuit()` is only the Newton loop, not the
  per-step driver. Rewrote the runner section around the real `electricalCircuit.cpp` loop;
  fixed the pseudocode (MNA recompute on Δt change, `shift_back()` on failure, output-header
  before `save_timestep`, deliberate phasing choice).
- **Citation fix:** `set_circuit()` `:332` → `:1749` (verified by grep).
- **2-D terminal-pair path:** corrected "only inputs used" → inputs duplicated.
- **"results broadcast"** → only the OK flag is broadcast.
- **Resistor floor / terminal-pair BDF / `order`-ignored / "time-invariant MNA"** wordings
  corrected to match the code.
- **Parser risk list extended:** node-count (not max-index), implicit/aliased ground,
  `MIL`/`A` suffixes, degrees→radians phase, SIN `Td`≠`time_offset`, `//` comments + `.end`,
  hard-error on `E/F/G/H/B/K/T` and other unsupported elements.

## Parser plan conclusions (the three questions)

1. **Pull circuit from ngspice?** Yes for the lumped subset (R/L/C/V/I/D → existing
   `create_*`, which take raw SI reals). `terminal pair` + `superconductor` need a BELFEM
   extension (recommended: `* belfem:` comment directives so stock ngspice still loads the
   file for cross-checking).
2. **Minimise `input.conf`?** Yes — hybrid `circuit { file : magnet.cir ; }`; keep FEM
   linkage in `input.conf` (recommended) or in `.cir` directives.
3. **Run circuit without FEM?** Already possible via `electricalCircuit.cpp`; the work is to
   make it file-driven (a `circuitrun` reading a netlist/trimmed conf).

## Changes Made / Proposed

- Created `src/circuit/doc/circuit_usage_guide.md`, `src/circuit/doc/README.md`.
- Created `todo/ngspice_parser_plan.md`; indexed in `todo/README.md`.
- No source modified.

## Open Questions

- Refactor `ElectricalCircuit` so ground is index 0 (SPICE-native) vs keep the remap layer.
- `.subckt` flattening needed for real magnet decks?
- Extend `SourceFunction` with `Pulse`/`PWL` types to close the SPICE source gap.
- Factor the duplicated circuit Newton loop (controller / `electricalCircuit.cpp` / new
  runner) into one shared `CircuitSolver`.

## Files Updated

- src/circuit/doc/circuit_usage_guide.md (new)
- src/circuit/doc/README.md (new)
- todo/ngspice_parser_plan.md (new)
- todo/README.md (index entry)
