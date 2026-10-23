# Circuit Module Usage Guide {#circuit_circuit_usage_guide}

**Date:** 2026-06-28
**Purpose:** How the `belfem::electronics` circuit module works, how it is driven from
`belfem` via the `ElectricalCircuitFactory`, what every component is, and how the
circuit couples to the Maxwell FEM problem.
**Module:** `src/circuit`

---

## 1. Overview

The circuit module is a small lumped-element circuit simulator embedded in BELFEM. Its
purpose is to drive (and be driven by) the electromagnetic FEM problem: a coil or tape in
the FEM mesh is exposed to the circuit as a two-terminal "black box", and the circuit
supplies the terminal current while the FEM returns the terminal voltage (and vice
versa). It can model the external power supply, protection resistors, dump/charge
inductors, filter capacitors, switches and diodes that surround a real superconducting
magnet.

Numerically the module solves the circuit with **Modified Nodal Analysis (MNA)** and a
**Newton–Raphson** loop on a sparse Jacobian (SuperLU by default — `SolverType::SUPERLU` in `cl_ElectricalCircuit.hpp`). Energy-storage
elements (L, C) are integrated in time with the BDF scheme from `numerics/ode`. The FEM
coupling element has no companion model: it passes the last FEM current into the circuit
(§4.9).

The module lives entirely in the `belfem::electronics` namespace. The abstract interface
that the FEM controller talks to is `belfem::Circuit`
(`src/numerics/sources/cl_Circuit.hpp`); `electronics::ElectricalCircuit` is the only
concrete implementation.

### Key files

| File | Role |
|------|------|
| `cl_ElectricalCircuitFactory.{hpp,cpp}` | Reads the `circuit{}` section of `input.conf`, builds the circuit, and **emits the FEM boundary conditions** for terminal pairs |
| `cl_ElectricalCircuit.{hpp,cpp}` | The circuit itself: node/component storage, MNA + Jacobian assembly, solve, time stepping, output |
| `cl_Circuit.hpp` (in `numerics/sources`) | Abstract base class — the contract the FEM `Controller` uses |
| `cl_Component.{hpp,cpp}` | Abstract component base (a `graph::Vertex`); holds terminals + label |
| `cl_TwoTerminals.{hpp,cpp}` | Common base for all 2-terminal components (adds `mCurrent`) |
| `cl_ElectricNode.{hpp,cpp}` | A circuit node (a `graph::Vertex`) holding a voltage |
| `cl_Resistor`, `cl_Capacitor`, `cl_Inductor`, `cl_VoltageSource`, `cl_CurrentSource`, `cl_Switch`, `cl_Diode`, `cl_Superconductor` | Concrete components |
| `cl_FEMTwoTerminals.{hpp,cpp}` | The FEM-coupling element (`terminal pair`) |
| `Component_Enums.{hpp,cpp}` | `ComponentType` enum + `to_string` / `component_type` string conversion |
| `cl_SourceFunction.{hpp,cpp}` (in `numerics/sources`) | Time-dependent excitation used by voltage/current sources |

---

## 2. Driving the circuit from `belfem`

The wiring in `src/executables/belfem.cpp` is the canonical example of how to use the
module. (The retired `hphirun.cpp` shows the same sequence and is still readable as a
reference driver, but it is no longer built.)

```cpp
// 1. Build the FEM factory (reads mesh + topology + boundary conditions)
MaxwellFactory tMFactory( "input.conf" );

// 2. Build the circuit. The factory is handed the FEM boundary-condition list by
//    reference so that each 'terminal pair' can PUSH a new PhysicalBoundaryCondition
//    onto it.                                              belfem.cpp:164
electronics::ElectricalCircuitFactory tElFactory(
        "input.conf", tMFactory.boundary_conditions() );

// 3. Create the FEM kernel (copies the boundary-condition list) and controller
auto tKernel  = tMFactory.create_magnetic_kernel();
auto tControl = tMFactory.create_controller();

// 4. Hand the circuit to the controller (non-owning)       belfem.cpp:204
tControl->set_circuit( tElFactory.circuit() );

// 5. Time loop — the controller drives both the FEM and the circuit
while ( tControl->time() < tControl->simulation_time() )
{
    tControl->initialize_timestep();
    ...
    tControl->solve_coupled();   // certified-exit nonlinear loop
    ...
}
```

Important ordering detail: the `ElectricalCircuitFactory` must be constructed **before**
`create_magnetic_kernel()`. Constructing it appends the circuit-derived boundary
conditions to `tMFactory.boundary_conditions()`. Then `create_magnetic_kernel()` copies
that list into the kernel. Later, `set_circuit()` looks for the `CircuitVoltage` /
`CircuitCurrent` BCs in the kernel. If the factory is built after the kernel, those BCs
are missing.

Ownership: `circuit()` returns a non-owning pointer. The factory owns the
`ElectricalCircuit` through a `std::unique_ptr` member, and the controller does not delete
it. Therefore, the factory must outlive every controller use of the circuit. In
`belfem.cpp`, both objects live until the end of `main`.

### Constructor flow (`ElectricalCircuitFactory`)

`cl_ElectricalCircuitFactory.cpp`, constructor

1. Opens the input file (its own `InputFile`).
2. If there is **no** `circuit{}` section, `mCircuit` stays `nullptr` and the whole
   module is inert (`tElFactory.circuit()` returns `nullptr`, the controller runs a pure
   FEM problem). The circuit is therefore fully optional.
3. If the section carries a `file` key, the lumped topology is read from that ngspice
   netlist (`read_circuit_from_netlist()`; `number of nodes` must then be absent and the
   deck's `topology{}` may only add `terminal pair`s). Otherwise it reads
   `number of nodes`, allocates `ElectricalCircuit`, and calls `read_circuit()` which
   parses `topology{}` then `output{}`.

---

## 3. The input file format (`circuit{}` section)

The circuit is described inside `input.conf` next to the FEM problem. There are two
modes. In **classic mode**, described in this section, the whole circuit is written in the
deck and nodes are integer indices. In **hybrid mode** (`circuit { file : … }`, see below),
the lumped circuit comes from an ngspice netlist and every node is referenced by its
netlist **name**. Classic-mode skeleton:

```
circuit
{
    number of nodes : 2 ;        // includes the ground node (last index)

    topology
    {
        <component> { ... }
        <component> { ... }
        ...
    }

    output                       // optional
    {
        file     : CircuitResults.txt ;
        currents : Z1, Z2, L, C, Is ;   // by component label
        voltages : 0, 1 ;               // by node index
    }
}
```

### Nodes and ground

- `number of nodes` counts **all** nodes including ground.
- The **last node index** (`number of nodes - 1`) **is the ground / reference node**
  (`ElectricalCircuit` constructor); its voltage is pinned to 0 and it is excluded from
  the unknowns. This differs from SPICE, where node `0` is ground; the netlist importer performs that remap for you.
- Every component requires `node +` and `node -` (`read_circuit()` in
  `cl_ElectricalCircuitFactory.cpp`). In classic mode both are integers `< number of nodes`.

### Hybrid mode (`circuit { file : … }`)

With a `file` key, `read_circuit_from_netlist()` builds the lumped circuit from the
netlist. `number of nodes` must then be absent, and `topology{}` may contain only
`terminal pair` sections. Every node reference in the deck — terminal-pair `node +` /
`node -` and `output → voltages` — is a netlist node **name** (`0` and `gnd` are ground),
resolved through the netlist's node map; raw indices are not accepted. Labels and the
`currents` list are folded to lower case. Shortened from `examples/tapestack_circuit/input.conf`:

```
circuit
{
    file : tapestack.cir ;
    topology
    {
        terminal pair { label : Z1 ; node + : live ; node - : 0 ;
                        input curves : [1,2] ; output curves : [3,4] ; }
    }
    output { file : CircuitResults.txt ; currents : Z1, L1, Is ; voltages : live, 0 ; }
}
```

The full hybrid grammar is in §11.1 of `doc/input_file_reference.md`.

### Per-component keys

Each subsection of `topology{}` is one component. The subsection **key** selects the type
(matched case-insensitively by `component_type()` in `Component_Enums.cpp`). All
components accept an optional `label` (used for output selection). Values are parsed with
BELFEM's unit machinery (`get_value`/`check_unit`/`unit_to_si`), so units are required and
checked.

See Section 4 for the per-component key tables.

### Output section

`read_output()` in `cl_ElectricalCircuitFactory.cpp`

- `file` (required if `output{}` is present) — output text file name.
- `currents` — comma-separated list of component **labels**; each matching component's
  current is written.
- `voltages` — comma-separated list of nodes; each node voltage is written. In classic
  mode these are **node indices**; in hybrid mode they are netlist node **names**. The
  column header is always `V<internal index>`; in hybrid mode rank 0 prints the
  name-to-column mapping at startup.
- A header row (`t I<label> … V<index> …`) is written by `init_output_file()`. One row per
  accepted timestep is appended by `save_timestep()`, called from the controller's
  `finalize()`; the leading column is time.
- Terminal-pair column timing: `finalize()` calls `compute_circuit_current()`, which
  passes the FEM terminal current for the accepted step to the terminal pair before
  `save_timestep()` writes the row. Every column, including a terminal pair's `I<label>`
  column, therefore belongs to the same accepted step. The circuit solve at the start of
  that step used the terminal current from the preceding accepted step (the initial state
  on a cold start), so the row's `I<label>` value is the current used to start the
  following step's circuit solve.

---

## 4. Available components

The complete component set is the `ComponentType` enum (`Component_Enums.hpp`). The
table below lists the input-file key, the required/optional keys, and the underlying
class.

| Input key | `ComponentType` | Class | Required keys | Optional keys |
|-----------|-----------------|-------|---------------|---------------|
| `resistor` | `RESISTOR` | `Resistor` | `value` [Ohm] | — |
| `capacitor` | `CAPACITOR` | `Capacitor` | `value` [F] | `order` (BDF order 1–6, default 1) |
| `inductor` | `INDUCTOR` | `Inductor` | `value` [H] | `order` (BDF order 1–6, default 1) |
| `voltage source` | `VOLTAGESOURCE` | `VoltageSource` | `type` + source params | (see source functions) |
| `current source` | `CURRENTSOURCE` | `CurrentSource` | `type` + source params | (see source functions) |
| `switch` | `SWITCH` | `Switch` | `initial state` (`open`/`closed`), `switch time` [s] | — |
| `diode` | `DIODE` | `Diode` | — | `Is` [A] (default 1e-4), `Vt` [V] (default 0.026) |
| `superconductor` | `SUPERCONDUCTOR` | `Superconductor` | `Ic` [A], `n` [-], `Ec` [V/m], `length` [m] | — |
| `terminal pair` | `TERMINALPAIR` | `FEMTwoTerminals` | `input terminals`/`input curves` (see §5) | `output terminals`/`output curves`; `length` required on the 2-D path (no output key; see the BC-type note in §5). **`order` is accepted but ignored** — see §4.9 |

All components additionally take `label`, `node +`, `node -`.

`order` must be written with a colon (`order : 2 ;`). The classic path passes it to the
`ode::BDF` constructor unchecked, which stops with an error for an order of 7 or more. A
line without a colon, such as `order 1 ;`, defines a different flag key named `order 1`;
it is never read, so the element silently keeps order 1.

### 4.1 Resistor

`cl_Resistor.cpp` — Ohm's law `I = (V+ − V−)/R`. Enters the MNA matrix as the conductance
stamp `1/R`. The constructor stores the **raw** value; only `set_value()` clamps to a
`1e-10` Ω floor, and the normal `create_resistor` path does **not** call it — so a
zero/near-zero resistor is *not* auto-protected on the input path (give a sane value).

```
resistor { label : R1 ; node + : 1 ; node - : 2 ; value : 0.0001 Ohm ; }
```

### 4.2 Capacitor

`cl_Capacitor.cpp` — companion-model (BDF) discretization. Per timestep it computes a
discretized resistance `mRC = h/(C·β0)` and a history current source `miCh`; both feed the
MNA/RHS each step. `order` selects the BDF order via a `ShiftRegister`. `I = C dV/dt`.

```
capacitor { label : C1 ; node + : 2 ; node - : 3 ; value : 100 F ; order : 1 ; }
```

### 4.3 Inductor

`cl_Inductor.cpp` — companion-model (BDF) discretization. Discretized resistance
`mRL = L·β0/h` plus history source `miLh`. `V = L dI/dt`.

```
inductor { label : L ; node + : 0 ; node - : 1 ; value : 2.5e-6 H ; order : 2 ; }
```

### 4.4 Voltage source

`cl_VoltageSource.cpp` — an "unknown current" element: it adds a row/column to the MNA
system (its branch current is a solved unknown) and enforces `V+ − V− = f(t)`. The
time-dependent value comes from a `SourceFunction` (Section 6).

```
voltage source
{
    label : Vs ; node + : 0 ; node - : 3 ;
    type : sine ; amplitude : 5 kV ; frequency : 50 Hz ; phase : 0 deg ;
}
```

### 4.5 Current source

`cl_CurrentSource.cpp` — injects `I = f(t)` from a `SourceFunction`. It is **not** an
unknown-current element; it stamps directly onto the RHS. Note: current sources are
excluded from the nodal adjacency graph (`compute_adjacency()` in `cl_ElectricalCircuit.cpp`) because they do
not couple node voltages.

```
current source
{
    label : Is ; node + : 0 ; node - : 1 ;
    type : sine ; amplitude : 500 A ; frequency : 10 Hz ; phase : 0 deg ;
}
```

### 4.6 Switch

`cl_Switch.cpp` — an unknown-current element that toggles state at `switch time`.
`closed` ⇒ enforce `V+ − V− = 0` (short); `open` ⇒ enforce branch current `= 0`. The
`shift()` call flips the state once when simulated time passes `switch time`.

```
switch { label : S1 ; node + : 1 ; node - : 2 ; initial state : open ; switch time : 20 ms ; }
```

### 4.7 Diode

`cl_Diode.cpp` — Shockley diode `I = Is·(exp((V+−V−)/Vt) − 1)`, contributing a nonlinear
`dI/dV` to the Jacobian. In classic mode the defaults are `Is = 1e-4 A` and
`Vt = 0.026 V` (`read_circuit()`, `cl_ElectricalCircuitFactory.cpp`). A diode imported
from a netlist follows SPICE instead: `Is` defaults to `1e-14 A` (`is=` on the `.model`
card) and `Vt = n · 0.026 V` (`cl_NgspiceCircuitFactory.cpp`).

```
diode { label : D1 ; node + : 2 ; node - : 0 ; Is : 1e-9 A ; Vt : 0.026 V ; }
```

### 4.8 Superconductor

`cl_Superconductor.cpp` — a **lumped** power-law (E–J) resistor, distinct from the FEM
HTS material. `I = sign(ΔV)·Ic·(|ΔV|/(Ec·ℓ))^(1/n)` with a clamped `dI/dV`. Use this when
you want an HTS branch in the circuit without meshing it.

```
superconductor
{
    label : SC1 ; node + : 0 ; node - : 1 ;
    Ic : 200 A ; n : 20 ; Ec : 1e-4 V/m ; length : 1 m ;
}
```

### 4.9 Terminal pair (FEM coupling) — the important one

`cl_FEMTwoTerminals.cpp` + `read_terminal_pair()` in `cl_ElectricalCircuitFactory.cpp`. This is the bridge
between the lumped circuit and the meshed Maxwell problem. In the circuit it behaves like
a current source whose value is the FEM terminal current (in `compute_jacobian_and_rhs`
it shares the RHS-only case with the current source); in the FEM it appears as a boundary
condition whose value is supplied by the circuit. See Section 5 for the full coupling
story and the curve/terminal linkage.

Unlike the inductor/capacitor, `FEMTwoTerminals` does **not** run a BDF companion model:
it hardcodes `ShiftRegister<real>(1)` for its I/V/h history, never instantiates its
`mBDF`, and `compute_current()` just relays the last FEM current `mIn`. Consequently the
`order` key shown in some examples is **never read for a terminal pair** (the factory
parses `order` only in the capacitor and inductor cases of `read_circuit()`)
and is silently ignored — keep it if you like, but it has no effect.

```
terminal pair
{
    label : Z1 ; node + : 0 ; node - : 1 ;
    input curves  : [1,2] ;   // references FEM topology{ curves{ } } ids
    output curves : [3,4] ;
    order : 1 ;               // NOTE: ignored for terminal pairs
}
```

---

## 5. FEM coupling: terminal pairs, terminals, and curves

### What the factory does for a `terminal pair`

`read_terminal_pair()` in `cl_ElectricalCircuitFactory.cpp`

1. Creates an `FEMTwoTerminals` component in the circuit (`create_terminal_pair`).
2. Reads the input geometry key. Four spellings are accepted; the *curve* forms set the
   "thin shell" flag:
   - `input terminal` / `input terminals` → volume/sideset terminals (thin shell = false)
   - `input curve` / `input curves` → thin-shell curves (thin shell = true)
3. Optionally reads matching `output terminal(s)` / `output curve(s)`.
   - If output is present (3-D problem): each `[...]` group of inputs is concatenated with
     the corresponding output group; the number of input and output groups must match
     (hard error otherwise).
   - If output is absent (2-D problem): the output groups are filled
     **from the input groups** and then concatenated with the inputs, so each BC's
     domain list contains the input ids **duplicated** (`PhysicalBoundaryCondition::set_domains`
     in `cl_FEM_PhysicalBoundaryCondition.cpp` does not deduplicate). This path requires
     `length` and scales the BC by `1/length` (see the BC-type note below).
4. It pushes **one** `fem::PhysicalBoundaryCondition` of type
   **`CircuitVoltage`** onto the shared boundary-condition list, with the
   concatenated domain ids and the thin-shell flag.

Within one key, each bracketed list is a *group*: `input curves : [1,2]` is one group
containing curves 1 and 2. **Exactly one group per terminal pair**: a terminal pair is one lumped component, and the controller pairs
boundary conditions with components by position, so several groups on one pair
would overrun that pairing — the factory refuses it with a hard error. Note
that an unbracketed list (`input curves : 1,2`) parses as one group *per id*
and trips the same error; write the bracketed `[1,2]`. For several pairs,
define several `terminal pair` components — in the example, `Z1` and `Z2` are
two separate components, each with one group.

### How `input curves` reference the FEM mesh

The curve ids are **not** mesh entity ids directly — they are the logical curve ids
defined in the FEM side of `input.conf` under `topology{ curves{ } }`. The `topology`
block below is from `examples/tapestack_circuit/input.conf`; the terminal pair is written
in classic mode (the shipped deck uses hybrid mode, with node names, see §3):

```
topology
{
    thinshell : tape { sideset : 5,6,7,8 ; }
    air       { blocks : 1:2 ; }
    curves
    {
        1 : 10 @ 5 ;     // logical curve 1 = mesh curve 10 on sideset 5
        2 : 10 @ 6 ;
        3 : 12 @ 5 ;
        ...
    }
}

circuit
{
    ...
    terminal pair { label : Z1 ; node + : 0 ; node - : 1 ;
                    input curves : [1,2] ; output curves : [3,4] ; }
}
```

Without thin shells you reference sideset/domain ids directly via `input terminals` /
`output terminals`. In `examples/rlc_circuit/input.conf`, the terminal pair `L1` uses
`input terminals : [4]` and `output terminals : [4]`; because it gives an output key, it
takes the 3-D path. A 2-D terminal pair omits the output key and must give `length`.

### The coupled solve (controller side)

`src/fem/kernel/cl_FEM_Controller.cpp`

- `set_circuit()` splits the kernel's boundary conditions into `mCircuitCurrentBCs`
  (`CircuitCurrent`) and `mCircuitVoltageBCs` (`CircuitVoltage`). **Note:** the electrical
  factory only ever emits `CircuitVoltage` BCs, so `mCircuitCurrentBCs` stays empty and
  the controller's `current()`-fix loop is inert for the standard `terminal pair` path.
- The circuit is solved **once per timestep attempt, on rank 0 only** (it is not
  distributed), at the start of the step in `initialize_timestep()` (and in
  `initialize_magnetic()` on the segregated path), before the FEM nonlinear loop:
  `set_timestep → shift → compute_MNA_matrix → solve_circuit()`. `solve_circuit()` is a
  Newton loop calling `compute_jacobian_and_rhs / solve / residual` with adaptive
  relaxation `omega`. The circuit is not re-solved inside the FEM nonlinear iterations.
- The resulting circuit voltages are `fix()`-ed into the FEM BCs. The circuit success flag
  is broadcast so that all ranks take the same reset branch if the rank-0 solve fails.
- After the FEM solve, in `finalize()`, `compute_circuit_current()` reads the terminal
  current from the FEM abstract DOF and writes both current and voltage back into the
  corresponding `FEMTwoTerminals` via `Circuit::set_current_and_voltage()`. The next
  step's circuit solve uses that current.
- `save_timestep()` (circuit output, called in `finalize()` right after
  `compute_circuit_current()`, see §3) and `shift_back()` (on reset) are also called from
  the controller.

> **Note on the BC type.** The factory always tags terminal-pair BCs as `CircuitVoltage`
> (`read_terminal_pair()`). A 2-D terminal pair — one that omits the output
> terminal/curve list — scales its BC by `1/length`, so it **requires** `length` and
> stops with an error without it. `length` must also be positive (zero would give an
> infinite scale, a negative value would flip the sign of the imposed voltage).

---

## 6. Source functions (voltage/current excitation)

Voltage and current sources use a `type` key with the parameters shown below. These
are parsed in `read_circuit()` (`cl_ElectricalCircuitFactory.cpp`) into a
`belfem::SourceFunction` (`src/numerics/sources/cl_SourceFunction.hpp`).

| `type` | Parameters | Formula |
|--------|------------|---------|
| `constant` | `amplitude` | `A` |
| `ramp` | `amplitude`, `period`, `offset` | rises linearly from 0 to A over `period`, starting after `offset`, then stays at A |
| `sigmoid` | `amplitude`, `period`, `offset`, opt. `fuzzyness` (0.01) | logistic ramp; the curve is set by `fuzzyness` |
| `sine` | `amplitude`, `frequency` or `period`, opt. `phase` | `A·sin(ωt + φ)` |
| `square` | `amplitude`, `frequency` or `period`, opt. `phase` | ±A square wave, evaluated at `t + φ/ω` |
| `triangle` | `amplitude`, `frequency` or `period`, opt. `phase` | triangle wave, evaluated at `t + φ/ω` |
| `sawtooth` | `amplitude`, `frequency` or `period`, opt. `phase` | sawtooth wave, evaluated at `t + φ/ω` |
| `userdefined` | `file`, `label`, `units` | function loaded from a shared library (see below) |

A positive `phase` advances all four periodic waveforms. This follows the SPICE
convention and uses the sign already shown in the `sine` formula above. The phase
produces a circular shift, not a delayed start. To delay the start, use `ramp` or
`sigmoid`, which instead take an `offset` in seconds.

A source that carries an `expk` key is refused. **Delete the key**: it never reached the
waveform, so deleting it leaves the results unchanged. Do not rewrite it as
`fuzzyness = 1/(1+expk)` (`expk = 99` would be `fuzzyness = 0.01`); that conversion
would change the waveform.

For `userdefined`, `SourceFunction::read_user_defined()` opens the library named by `file`
and calls `<label>_init`. That function receives the `SourceFunction *` and configures it.
Because the function is looked up by its plain name, declare it `extern "C"`. The `units`
key specifies the unit of the value returned by the plugin. BELFEM checks it against the
source kind and passes its SI factor to `SourceFunction::set_user_scale()`, which
multiplies every value returned by the plugin. A plugin that returns kA with
`units : kA` therefore drives the circuit in A.

Units are checked against the source kind: `V` for a voltage source and `A` for a
current source. Use `s` for `period`, `Hz` for `frequency`, and `rad` for `phase`
(`deg` is converted).

---

## 7. Numerical method (how the solve works)

### Modified Nodal Analysis

`compute_MNA_matrix()` (`cl_ElectricalCircuit.cpp`) builds the **linear** part of the
system (everything except the Newton nonlinearities) once the sparsity is known. It is
not time-*invariant*: it re-stamps the L/C companion conductances, which change whenever
the timestep changes, so it must be recomputed on every Δt change (both the FEM controller and the
test helper in `tests/circuit/test_ElectricalCircuit.cpp` recompute it at the top of every attempt):

- The unknown vector is `[ node voltages (N−1) ; unknown branch currents ]`. Unknown-current
  elements are voltage sources and switches (`mComponentsUnknown`, filled by
  `create_voltage_source()` and `create_switch()`).
- R, L, C stamp conductances (L/C use their *discretized* resistance) into the nodal
  block. Voltage sources stamp the ±1 incidence entries linking their branch-current
  unknown to the two nodes.
- The ground node (last index) is never assembled.

`compute_jacobian_and_rhs()` builds the Newton residual/Jacobian each iteration:
nonlinear stamps (diode, superconductor `dI/dV`), source RHS contributions, switch
logic, and the L/C history sources; then it adds the constant MNA matrix into the Jacobian.

`solve()` does one Newton update `x −= ω·J⁻¹·r` (SuperLU), pushes voltages to
nodes and currents to components, then recomputes a true RHS norm for the residual.
`residual()` returns `‖r‖ / ‖RHS‖`.

### Time integration

- `set_timestep()` updates L and C companion resistances.
- `shift()` advances time, pushes the BDF history registers (L current, C voltage, FEM
  terminal I/V), recomputes companion coefficients, and resets currents.
- `shift_back()` reverts the history (used when the controller rejects a timestep).
- L and C use `ode::BDF` + `ShiftRegister<real>`; the `order` key sets the BDF order per
  element. The FEM coupling element (`FEMTwoTerminals`) keeps I/V/h `ShiftRegister`s for
  history but does **not** run a BDF companion model (it relays the FEM current directly),
  so its `order` is ignored (see §4.9).

### Abstract `Circuit` interface

The FEM controller only ever sees the abstract base (`numerics/sources/cl_Circuit.hpp`):
`set_timestep, set_omega, compute_MNA_matrix, compute_jacobian_and_rhs, solve, residual,
shift, shift_back, current, voltage, set_current_and_voltage, save_timestep, save_state,
load_state`. Anything
that wants to replace the circuit engine (e.g. a SPICE-backed one) implements this
interface.

---

## 8. Worked example: a tape stack driven by a current source

A current source `Is` (sine, 500 A, 10 Hz) is in parallel with two FEM tapes (`Z1`, `Z2`),
a dump inductor `L1` and a filter capacitor `C1`. All are connected between circuit nodes
0 and 1; node 1 is ground. This is the classic-mode form of the circuit topology in
`examples/tapestack_circuit`. The shipped deck writes the lumped part as the netlist
`tapestack.cir` and uses hybrid mode, with the node named `live` (§3).

```
circuit
{
    number of nodes : 2 ;
    topology
    {
        current source { label : Is ; node + : 0 ; node - : 1 ;
                         type : sine ; amplitude : 500 A ; frequency : 10 Hz ; phase : 0 deg ; }
        terminal pair  { label : Z1 ; node + : 0 ; node - : 1 ;
                         input curves : [1,2] ; output curves : [3,4] ; }
        terminal pair  { label : Z2 ; node + : 0 ; node - : 1 ;
                         input curves : [5,6] ; output curves : [7,8] ; }
        inductor       { label : L1 ; node + : 0 ; node - : 1 ; value : 2.5e-6 H ; order : 2 ; }
        capacitor      { label : C1 ; node + : 0 ; node - : 1 ; value : 100 F ; order : 2 ; }
    }
    output { file : CircuitResults.txt ; currents : Z1, Z2, L1, C1, Is ; voltages : 0,1 ; }
}
```

`CircuitResults.txt` then has the header `t IZ1 IZ2 IL1 IC1 IIs V0 V1` and one row for
each accepted timestep. The `IZ1` and `IZ2` columns contain the FEM terminal current from
that same timestep (§3, output section).

---

## 9. Restart: the `/circuit` group in `memdump.hdf5`

When a circuit is attached, rank 0 writes its state to the `/circuit` group in
the memdump, alongside the field data. `load_memdump` restores this state before
the first post-restart timestep.

The dump layout (format version 2) is:

```
/circuit
    time            scalar  real     circuit clock at the dump
    delta_time      scalar  real     Δt of the last completed step
    x               [n]     real     solution vector (node voltages + branch currents)
    prev_x          [n]     real     previous solution vector
    n_components    scalar  uint     component count, guards the walk below
    c<NNN>_type     scalar  uint     component type at creation index NNN
    c<NNN>_...               ...     per-component state, prefix = creation index
```

Here, `n = number_of_nodes - 1 + number_of_unknown_currents`.

The following per-component datasets are present only for stateful components:

| Component | Datasets under its `c<NNN>_` prefix |
|---|---|
| inductor | `h`, `h_cap`, `i`, `i_cap` (time/current histories, newest first, with capacities), `current` |
| capacitor | `h`, `h_cap`, `v`, `v_cap`, `current` |
| terminal pair | `h`, `i`, `v` (+ `_cap` each), `in`, `vn` (last accepted FEM current/voltage) |
| switch | `is_closed`, `is_switched` |
| voltage source | `value_time` |
| current source | `current` |

Resistors, diodes, and superconductors do not store restart state, and `load_state()` does
not recompute their currents. It restores only the node voltages and the unknown branch
currents. Until the first circuit `solve()` after the restart calls `compute_current()` on
them, these components keep their constructed current of 0.

The loader enforces the following rules. If any rule is violated, loading stops
with a hard error. The remedy is to delete the memdump and restart with a cold
circuit.

- The dof count, component count, per-component types, and register capacities
  must match the constructed circuit. A restart continues the same deck; it is
  not a migration path for an edited deck. These checks are best-effort. The
  loader cannot detect a reordering of two components of the same type, so such
  a change silently swaps their histories. For this reason, using the same deck
  is a contract, not merely a recommendation.
- A dump without `n_components` (an older format) is refused. This prevents the
  loader from silently cold-starting the reactive histories.

After `load_state()`, callers may use either the controller sequence
`set_timestep → shift → compute_MNA_matrix` or stamp without first calling
`shift`. Do not call `shift_back()` before the first post-restore `shift()`,
because there is no attempt to revert. A restored switch retains its fired
latch and does not re-fire.

## 10. Common pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| Forgetting that the **last** node is ground | Wrong currents / singular matrix; node you think is ground is floating | Put ground last; size `number of nodes` to include it |
| Using SPICE's node `0` as ground in a hand-written `topology{}` | Node 0 is a normal node here | Put ground last — or import the netlist through `circuit { file : … }`, where `NgspiceCircuitFactory` maps `0`/`gnd` to the last index automatically |
| Missing units on a value | `check_unit` BELFEM_ERROR | Always give units (`Ohm`, `H`, `F`, `A`, `V`, `s`, `Hz`, `deg`) |
| `terminal pair` curve ids not defined in FEM `topology{curves{}}` | Curve/terminal lookup fails downstream | Define logical curves on the FEM side first |
| Mismatched input/output group counts on a 3-D terminal pair | BELFEM_ERROR in `read_terminal_pair()` | Equal number of `[...]` groups in `input` and `output` |
| Expecting the circuit to run in parallel | Circuit solves on rank 0 only | By design. Only the success flag is broadcast (`initialize_timestep()` in `cl_FEM_Controller.cpp`); the circuit-derived BC values are `fix()`-ed on rank 0 — the *FEM* solve is what is distributed |
| Constructing the factory after the kernel | Circuit BCs missing from the kernel | Build `ElectricalCircuitFactory` before `create_magnetic_kernel()` |
| Destroying the factory while the controller still runs | Dangling circuit pointer | Keep `ElectricalCircuitFactory` alive as long as the controller uses `circuit()` |
| Writing `order 2 ;` without a colon | Element silently stays at BDF order 1 | Write `order : 2 ;` |

---

## 11. Related work

- **FEM-free circuit runs** — there is no standalone circuit executable. The template for
  driving an `ElectricalCircuit` without a mesh is the test helper in
  `tests/circuit/test_ElectricalCircuit.cpp`: `set_timestep()` → `shift()` →
  `compute_MNA_matrix()` → Newton iterations until the residual tolerance is met. This is the
  same sequence the FEM controller runs. The transient cases in that file (resistive, switched
  RLC, diode bridge, and a netlist twin) are the circuit-only examples.
- **SPICE/ngspice import** — `NetlistParser` and
  `NgspiceCircuitFactory` build an `ElectricalCircuit` from a `.cir` netlist through the same
  `create_*` calls, and the hybrid `circuit → file` key pulls the lumped topology out of the deck
  while terminal pairs stay in `input.conf` (§11.1 of `doc/input_file_reference.md`;
  `examples/tapestack_circuit` runs from `tapestack.cir`). Netlist sources may be a bare
  value, `DC <value>`, or `SIN(...)` with zero offset, delay and damping; `PULSE`, `PWL`,
  `EXP`, `SFFM` and `AM` are refused, as are nonzero `SIN` offset, delay or damping. The
  full card list and refusals are in `cl_NgspiceCircuitFactory.hpp`. There is no standalone circuit
  executable; a circuit always runs coupled inside `belfem`.
- **Refactoring notes** — `src/circuit/doc/notes.md` (namespace cleanup, dependency
  inversion, factory split).
