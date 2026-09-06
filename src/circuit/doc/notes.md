# Circuit Model Refactoring Notes {#circuit_notes}

## Testing Strategy

Before diving into major refactoring, it's crucial to validate that our current circuit model approach is fundamentally sound. The circuit integration with the electromagnetic solver is complex, and we want to catch any conceptual issues early before we invest significant time in restructuring the code. A simple test case will help us understand the current coupling mechanism and identify potential problems.

* Before continuing with refactoring, implement and test a simple circuit model example to validate the current approach
* This will help identify any fundamental issues before extensive restructuring

## Naming Convention Improvements

The current naming scheme lacks consistency and clarity. We need to establish a coherent naming convention that makes the code more readable and maintainable. The introduction of a proper namespace will help organize the circuit-related functionality and prevent naming conflicts with other parts of BELFEM.

### Namespace Organization

Creating a dedicated namespace for circuit functionality will cleanly separate this domain-specific code from the rest of the FEM framework. This follows good software engineering practices and makes the API more intuitive for users.

* Introduce `belfem::electronics` namespace to clearly separate circuit functionality
* Rename classes for consistency:
  - `ElectricalCircuit` → `belfem::electronics::Circuit`
  - `ElectricNode` → `belfem::electronics::Node`
* This creates a cleaner, more intuitive API structure

### Enum and Class Naming

Consistent naming conventions improve code readability significantly. The current enum naming and class names could be more descriptive and follow established C++ conventions.

* Use PascalCase for enum values: `ComponentType::Resistor` (more readable than ALL_CAPS)
* Rename `TwoTerminals` → `TwoTerminalComponent` for explicit clarity about the class purpose

## Architecture Refactoring

The current architecture has some fundamental design issues that need to be addressed. The main problems are tight coupling between modules and incorrect dependency directions that make the code hard to maintain and test.

### Factory Pattern Implementation

Right now, component creation is scattered throughout the code, making it difficult to manage and test. A factory pattern will centralize this logic and make it easier to add new component types or modify existing ones.

* Move component creation logic to dedicated factories:
  ```cpp
  Component* tComponent = tFactory.create_resistor(...);
  ```
* Create separate factories:
  - `electronics::ComponentFactory` for circuit components
  - `electronics::BoundaryConditionFactory` for electrical boundary conditions

### Dependency Inversion

This is the most critical architectural issue we need to fix. Currently, core FEM modules depend on the circuit code, which is backwards from a software architecture perspective. This creates unnecessary coupling and prevents us from making the circuit module optional.

* **Former problem (fixed)**: FEM modules depended on the Circuit directory; the controller now sees only `numerics/sources/cl_Circuit.hpp`
* **Solution**: Reverse dependency - Electronics depends on FEM, not vice versa
* **Benefits**: 
  - Circuit module can be optionally disabled in CMake
  - Better separation of concerns
  - Cleaner build dependencies

### Integration Strategy

The integration between the circuit model and the electromagnetic solver needs to be carefully managed. We need a clean interface that allows the Maxwell solver to interact with circuits when present, but doesn't break when no circuit is involved.

1. Virtual base class `cl_Circuit.hpp` in numerics provides interface abstraction
2. In the driver (`belfem.cpp:146,164,204`; the retired `hphirun.cpp` shows the same shape):
   - `MaxwellFactory` is built first, `ElectricalCircuitFactory` second (it appends the
     terminal-pair BCs to the factory's list)
   - The circuit is handed to the Controller via `set_circuit()`
   - Handle nullptr case when no circuit is present
3. Controller setter methods (done): `set_circuit()` splits the BCs into
   - `mCircuitCurrentBCs` (current boundary conditions)
   - `mCircuitVoltageBCs` (voltage boundary conditions)

## Time Integration (BDF)

For time-dependent circuit problems, we need robust time integration. The Backward Differentiation Formula (BDF) methods are particularly well-suited for stiff differential equations that often arise in circuit simulations, especially when dealing with different time scales in electromagnetic and circuit dynamics.

### BDF Implementation

The BDF implementation provides the numerical method for solving the time-dependent differential equations. It's specifically designed to handle the stiff nature of coupled electromagnetic-circuit problems where rapid changes in circuit variables must be resolved accurately.

* Implementation available in `numerics/ode/bdf.cpp`
* **Key Features**:
  - Automatic timestep adaptation (configure once, use everywhere)
  - Handles stiff differential equations efficiently
  - Designed specifically for electromagnetic-circuit coupling

### ShiftRegister Data Container

The ShiftRegister is a specialized data structure that supports the BDF implementation by managing the historical data needed for the multi-step BDF schemes. It handles the storage and retrieval of previous timestep values efficiently.

* New data container that works with the BDF implementation
* **Key Features**:
  - Manages historical timestep data for BDF methods
  - `revert()` function allows undoing timestep operations
  - Optimized for the memory access patterns of BDF schemes

### Usage Guidelines

The BDF and ShiftRegister work together to provide consistent time integration across all time-dependent solvers in BELFEM. The BDF handles the numerical method while ShiftRegister manages the data storage efficiently.

* Check implementation example in `numerics->ode/bdf.cpp`
* Timestep management is centralized through the BDF implementation
* ShiftRegister provides the underlying data management