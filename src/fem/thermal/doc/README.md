# Thermal Module Documentation {#fem_thermal_index}

**Module:** `src/fem/thermal`
**Purpose:** Index of documentation for BELFEM's thermal and magneto-thermal coupling module

---

## Overview

The `thermal` module builds the heat-conduction side of a coupled magneto-thermal
problem. It is the counterpart of `src/fem/maxwell`. That module assembles the electromagnetic
h-φ formulation; `thermal` assembles transient heat conduction and connects it to the
electromagnetic solution. The field problem computes ohmic and hysteretic losses, which
become the thermal heat source.

This module backs the coupled h-ɸ/T problem that `belfem` selects when a deck carries an
unlabeled `linear thermal` or `nonlinear thermal` solver section; without one, the run is
magnetic-only and this module is not built into the solve.

**Key capabilities:**
- Builds a thermal `Kernel` from an `input.conf` deck, either standalone on a mesh or
  coupled to an existing magnetic `Kernel`
- Assembles a transient heat conduction weak form specialized for the Maxwell coupling
- Builds the thermal boundary conditions from the input deck
- Provides element matrices for the h-formulation and the φ-formulation regions

---

## Key Classes

| Class | Files | Role |
|-------|-------|------|
| **`ThermalFactory`** | cl_ThermalFactory.{hpp,cpp} | High-level orchestrator; builds the thermal kernel from an input deck |
| **`IWG_MaxwellThermal`** | cl_IWG_MaxwellThermal.{hpp,cpp} | Coupled weak form; derives from `IWG_TransientHeatConduction` |
| **`ThermalBoundaryConditionFactory`** | cl_ThermalBoundaryConditionFactory.{hpp,cpp} | Builds the thermal boundary conditions |

Element matrices for the two formulation regions live in `matrices/`
(`mt_thermal_h.{hpp,cpp}`, `mt_thermal_phi.{hpp,cpp}`).

## Entry Points

`ThermalFactory` has separate constructors for standalone and coupled setup:

```cpp
ThermalFactory( const string & aInputFile, Mesh   * aMesh );            // standalone
ThermalFactory( const string & aInputFile, Kernel * aMagneticKernel );  // coupled
```

`create_thermal_kernel()` then returns the assembled kernel. Its consumer in the
open-source tree is `src/executables/belfem.cpp`, which constructs the `ThermalFactory` and
hands the thermal kernel to the controller only when the deck asks for the coupled problem.
(`src/executables/hphiTrun.cpp` is the retired dedicated driver — the source is still present
but is no longer built.)

---

## See Also

- **Source code** - header files in `src/fem/thermal/`
- [Maxwell module](../../maxwell/doc/README.md) - the electromagnetic half of the coupling
- [IWG module](../../iwg/doc/README.md) - `IWG_TransientHeatConduction`, the base weak form
- [Executables](../../../executables/doc/README.md) - `belfem`, which solves the coupled problem when the deck carries a thermal solver section

---

**Status:** this index was written from the module's headers and its single consumer. The
coupling itself — how the loss term is transferred and how the two kernels are stepped
relative to one another — is not yet documented and deserves a usage guide of its own.
