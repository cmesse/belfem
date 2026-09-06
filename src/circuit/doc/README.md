# Circuit Module Documentation {#circuit_index}

**Module:** `src/circuit` (`belfem::electronics`)
**Date:** 2026-06-28

A lumped-element circuit simulator embedded in BELFEM. It supplies/receives terminal
current and voltage to/from the Maxwell FEM problem (h-φ), and models the external power
supply, dump resistors, charge/dump inductors, filter capacitors, switches, diodes and
lumped HTS branches around a superconducting magnet.

## Documents

| Document | Contents |
|----------|----------|
| [circuit_usage_guide.md](circuit_usage_guide.md) | Full usage guide: driving the circuit from `belfem`, the `circuit{}` input format, every available component, FEM coupling via terminal pairs/curves, source functions, the MNA + Newton + BDF numerics, a worked example, and pitfalls |

## Quick reference

| Topic | Where |
|-------|-------|
| Entry point from the executable | `src/executables/belfem.cpp:164` (`ElectricalCircuitFactory`) |
| Reads the `circuit{}` section | `cl_ElectricalCircuitFactory.cpp` |
| MNA / Jacobian / solve / timestep | `cl_ElectricalCircuit.cpp` |
| Abstract interface used by the FEM controller | `src/numerics/sources/cl_Circuit.hpp` |
| Time-dependent excitation | `src/numerics/sources/cl_SourceFunction.hpp` |
| Component list | `Component_Enums.hpp` |

## Available components

`resistor`, `capacitor`, `inductor`, `voltage source`, `current source`, `switch`,
`diode`, `superconductor`, `terminal pair` (FEM-coupling). See the usage guide §4 for the
per-component key tables.

## Related

- `src/circuit/doc/notes.md` — refactoring notes (namespace, dependency inversion).
- `todo/ngspice_parser_plan.md` — the ngspice/SPICE netlist importer and a proposed
  standalone circuit-only runner. **The parser shipped**: `cl_NetlistParser` and
  `cl_NgspiceCircuitFactory` exist, and a deck reads a netlist through `circuit { file : … }`.
  No `circuitrun` executable exists, and that part remains deferred.
