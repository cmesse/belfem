# Circuit Cut A: Dependency Hygiene and the Terminal-Pair Freeze

**Date:** 2026-08-15
**Purpose:** Record the execution of the audited Cut A plan — circuit/kernel
dependency cleanup plus the DR-75 grammar freeze
**Module:** circuit, fem/kernel, fem/maxwell, executables
**Round:** `tmp/ai_exchange/circuit_cut_a_plan.md` + `circuit_cut_a_impl.patch`

## Origin

Grok's circuit-dependency sweep (relayed by Christian) claimed a CMake link
cycle and a heavyweight include graph. The plan-phase audits (Codex + Grok,
both in the thread) refuted the cycle — `LIBLIST` never feeds
`target_link_libraries` for libraries, only example executables — but
confirmed the hygiene targets and one real defect: the ordinal
circuit↔FEM pairing silently overruns when a terminal pair carries more
than one bracket group (DR-75).

## What landed (A1–A6, A3 before A4 — the audited order)

- **A1** — `circuit` out of the kernel/maxwell `LIBLIST` rows. Zero effect
  on the default tree; the real executables link circuit in
  `src/executables/CMakeLists.txt`.
- **A2** — the dead `cl_Circuit.hpp` includes out of the two Maxwell
  factory headers, and the never-defined
  `MaxwellBoundaryConditionFactory::read_circuit` declaration with them.
  The Controller keeps its include: `Circuit * mCircuit` is the sanctioned
  exchange point.
- **A3** — `cl_FEM_PhysicalBoundaryCondition.hpp` thinned: its includes
  had sat ABOVE the include guard since creation; now guarded, with
  `DofManager` and `Dof` forward-declared (pointer-only members) and the
  fat `cl_FEM_DofManager.hpp` moved into the `.cpp`. This is what takes
  the circuit library off the FEM include graph.
- **A4** — `src/circuit/CMakeLists.txt` include rows cut from 13 to 3
  (`numerics/ode`, `numerics/sources`, and `fem/kernel` until Cut B); the
  rest were either never used or auto-injected by `Add_Library.cmake`.
- **A5** — `electricalCircuit.cpp` builds its circuit by hand and used
  none of its four Maxwell-side includes. Dropping them surfaced two
  hidden transitive dependencies — `using namespace fem;` (nothing fem
  used; deleted) and `constant::pi` (now includes `constants.hpp`
  directly). The gate caught both, which is the A5 story in miniature:
  the includes were dead, but not inert.
- **A6** — the DR-75 freeze: a terminal pair with more than one bracket
  group is refused at setup with a named error (component label, group
  count, and the bracket hint — an unbracketed `1,2` parses as one group
  per id and trips the same error). The now-unreachable group suffix on
  the circuit-side BC label was removed; the Maxwell-side suffix stays,
  it serves the live multi-curve `current` path. All three contract
  documents rewritten in the same session: reference §11 row, schema
  `group_count: 1` constraint, usage-guide grammar paragraph.

## Gates

All static, all green: the complete 15-TU circuit SOURCES list compiled
against the reduced include set (the A4 proof); the PBC fan-out
(kernel/maxwell/thermal factories, Kernel, DofManager); the three
executables; the two named test TUs. Owed at the next rebuild: real
`make` of the four libraries, `make check-fast`, a synthetic two-group
deck aborting with the named message, and the two circuit example decks
still parsing.

Status: **reviewed pending phase-3 audits** — implementation audits
(Codex + Grok) dispatched on the diff; findings land in the exchange
thread. Cut B (relocating the PBC construction out of the circuit
factory) remains future work; `fem/kernel` stays on the circuit include
list until then.
