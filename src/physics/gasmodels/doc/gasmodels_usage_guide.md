# Gas Models Usage Guide {#physics_gasmodels_gasmodels_usage_guide}

**Date:** 2026-01-30
**Purpose:** Comprehensive guide for using the gasmodels module
**Module:** `src/physics/gasmodels`

---

## Table of Contents

1. [Introduction](#introduction)
2. [Basic Usage](#basic-usage)
3. [Equation of State Selection](#equation-of-state-selection)
4. [Gas Mixtures](#gas-mixtures)
5. [Property Evaluation](#property-evaluation)
6. [Chemical Equilibrium](#chemical-equilibrium)
7. [Compressible Flow Applications](#compressible-flow-applications)
8. [Advanced Topics](#advanced-topics)
9. [Performance Optimization](#performance-optimization)

---

## Introduction

The **gasmodels** module provides a unified interface to various **equations of state** (EoS) for gas property calculations. It builds on the `gastables` module to create complete thermodynamic models with:

- **Ideal gas behavior** for low-pressure applications
- **Cubic equations of state** (SRK, Peng-Robinson) for moderate-to-high pressure
- **Helmholtz free energy models** for cryogenic fluids (H2 — para/normal/ortho — plus O2, CH4 and N2)
- **Gas mixtures** with composition-dependent properties
- **Chemical equilibrium** via Gibbs minimization
- **Compressible flow utilities** for aerospace/propulsion applications

**Design philosophy:**

Following BELFEM's HPC-first approach (see `doc/coding_philosophy.md`):
- **Manual memory management** for `Gas` objects (stack or explicit `delete`)
- **Function pointer dispatch** for zero-overhead EoS polymorphism
- **Preallocated work arrays** for mixture property calculations
- **Spline-based evaluation** for repeated property lookups

**Building the module:** `gasmodels` and the `gastables` module it uses are compiled only when
the build is configured with `-DUSE_GASMODELS=ON`. The option defaults to OFF
(`CMakeLists.txt`, `option( USE_GASMODELS ... )`). When it is ON, the configure also adds the
preprocessor define `BELFEM_GASMODELS`. A default build contains no `Gas` class, so the examples
below cannot link against it.

**Ownership:** a `Gas` owns its component `RefGas` objects and its EoS object, and deletes them
in its destructor. Copy construction and copy assignment are deleted. The raw pointers returned
by `eos()`, `components()` and `component( i )` are non-owning borrows. Never `delete` them, and
do not use them after the `Gas` that produced them is destroyed.

---

## Basic Usage

### Creating an Ideal Gas

The simplest use case - ideal gas model for a pure component:

```cpp
#include "cl_Gas.hpp"

using namespace belfem;

int main()
{
    // Create ideal gas (air)
    Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

    // Define state
    real T = 300.0;      // Temperature (K)
    real p = 101325.0;   // Pressure (Pa)

    // Evaluate properties
    real rho = air.rho(T, p);        // Density (kg/m³)
    real cp = air.cp(T, p);          // Specific heat (J/(kg·K))
    real h = air.h(T, p);            // Enthalpy (J/kg)
    real s = air.s(T, p);            // Entropy (J/(kg·K))
    real mu = air.mu(T, p);          // Viscosity (Pa·s)
    real lambda = air.lambda(T, p);  // Thermal conductivity (W/(m·K))

    return 0;
}
```

**Key points:**
- `Gas` constructor creates internal `RefGas` objects from gastables
- `Gas` destructor cleans up all internal resources
- All property methods are `const` — but **`const` here is logical constness, not thread safety.**
  The const evaluators write a shared cache, so two threads must not call them on the same `Gas`
  object even through a const reference (`cl_Gas.hpp:74-78`). That matches the framework-wide
  posture: BELFEM parallelizes with MPI and is not internally thread safe
- Temperature in Kelvin, pressure in Pascal (SI units throughout)

### Creating a Custom Mixture

```cpp
// Define species and composition
Cell<string> species = {"N2", "O2", "Ar", "CO2"};
Vector<real> molar_fractions = {0.78084, 0.20946, 0.00934, 0.00040};

// Create custom air mixture
Gas custom_air(species, molar_fractions, GasModel::IDGAS);

// Evaluate mixture properties
real T = 300.0, p = 101325.0;
real cp_mix = custom_air.cp(T, p);   // Mixture-averaged cp
real mu_mix = custom_air.mu(T, p);   // RP-1311 mixing rule
real lambda_mix = custom_air.lambda(T, p);  // RP-1311 mixing rule
```

### Using a Real Gas (Cubic EoS)

```cpp
// Peng-Robinson model for nitrogen
Gas N2("N2", GasModel::PR);

// Real gas behavior at high pressure
real T = 120.0;   // K (cryogenic)
real p = 10e6;    // Pa (100 bar)

real v = N2.v(T, p);           // Specific volume (m³/kg)
real Z = p * v / (N2.R(T,p) * T);  // Compressibility factor (Z < 1)

real h = N2.h(T, p);           // Enthalpy with departure function
real s = N2.s(T, p);           // Entropy with departure function
real cp = N2.cp(T, p);         // Heat capacity with departure

// Access departure functions directly (borrowed pointer, owned by N2)
gasmodels::EoS * eos = N2.eos();
real h_dep = eos->hdep(T, p);  // h(T, p) = h_ideal(T) + h_dep(T, p) - h_dep(T, 1 bar); see README, Reference Pressure Convention
```

### Using a Helmholtz EoS (Cryogenic)

```cpp
// High-accuracy hydrogen model
Gas hydrogen(HelmholtzModel::NormalHydrogen);

// Valid from triple point to supercritical
real T = 20.0;    // K (liquid H2)
real p = 101325.0;

real v = hydrogen.v(T, p);     // Accurate for liquid/vapor/supercritical
real h = hydrogen.h(T, p);
real cp = hydrogen.cp(T, p);
real w = hydrogen.c(T, p);     // Speed of sound

// Check phase of the state evaluated last
if (hydrogen.is_liquid()) {
    // Liquid phase at these conditions
}
```

`is_liquid()` takes no state. It reports a flag that the Helmholtz EoS sets during each volume
solve, and the calls to `v`, `h`, `cp` and `c` above all trigger such a solve. Before any
evaluation, the flag is false. For ideal-gas and cubic (SRK/PR) gases, no code in the module sets
it, so it stays false.

---

## Equation of State Selection

### Decision Tree

```
Need gas properties?
│
├─ Cryogenic fluid (H2, O2, CH4, N2)?
│  └─ Use HELMHOLTZ
│     - Accurate across phase boundaries
│     - Wide pressure/temperature range
│     - Handles liquid, vapor, supercritical
│
├─ Pressure > 10 bar?
│  ├─ Non-polar or weakly polar?
│  │  └─ Use SRK
│  │     - Good for hydrocarbons
│  │     - Moderate accuracy, fast
│  │
│  └─ Polar or near critical point?
│     └─ Use PR (Peng-Robinson)
│        - Better for polar molecules
│        - Vapour root only (see note below)
│
└─ Low pressure (< 10 bar)?
   └─ Use IDGAS
      - Simple, fast
      - Accurate for most gases at ambient conditions
```

**Root selection in the cubic models:** when Cardano's formula returns three real roots,
`EoS_Cubic::v` takes the smallest root above the covolume if `Gas::is_liquid()` is true and the
largest (vapour) root otherwise. Only the Helmholtz EoS writes that flag. `set_liquid_flag` is
public, but nothing in the module calls it for SRK or PR, so a cubic-EoS `Gas` always returns the
vapour root.

### Comparison Table

| EoS | Accuracy | Speed | Pressure Range | Typical Applications |
|-----|----------|-------|----------------|---------------------|
| **IDGAS** | Good (low p) | Fast | < 10 bar | Combustion, HVAC, ambient conditions |
| **SRK** | Good (moderate p) | Fast | 1-100 bar | Natural gas, petrochemical processing |
| **PR** | Excellent (high p) | Moderate | 1-1000 bar | Supercritical extraction, refrigeration |
| **HELMHOLTZ** | Excellent (within its validity range) | Slow | Up to the model's `mPmax` (see below) | Cryogenic storage, liquefaction, reference data |

Each Helmholtz model has an upper pressure limit, and `Helmholtz::v` aborts with an always-active
`BELFEM_ERROR` ("Pressure ... exceeds the validity limit ...") above it: 81.8 MPa for oxygen,
1000 MPa for methane, 2000 MPa for the hydrogen models, and 2200 MPa for nitrogen (`mPmax` in
`cl_GM_EoS_Oxygen.cpp`, `cl_GM_EoS_Methane.cpp`, `cl_GM_EoS_Hydrogen.cpp`,
`cl_GM_EoS_Nitrogen.cpp`).

### When to Use Each Model

**IDGAS:**
```cpp
// Good for:
Gas air_combustion;   // default ctor: the 14-species air mixture         // Combustion at ambient p
Gas exhaust_gas("CO2", GasModel::IDGAS);            // Exhaust at low p
Gas ammonia("NH3", GasModel::IDGAS);                // If p < 5 bar

// Not recommended for:
Gas liquid_nitrogen("N2", GasModel::IDGAS);  // Cryogenic liquid - use HELMHOLTZ
Gas supercritical_CO2("CO2", GasModel::IDGAS);  // High p - use PR
```

**SRK (Soave-Redlich-Kwong):**
```cpp
// Good for:
Gas natural_gas("CH4", GasModel::SRK);       // Natural gas pipelines (< 100 bar)
Gas ethane("C2H6", GasModel::SRK);           // LPG-range light hydrocarbon
Gas ethylene("C2H4", GasModel::SRK);         // petrochemical processing

// Less accurate for:
Gas CO2_near_crit("CO2", GasModel::SRK);     // Near critical - use PR instead
Gas polar_molecule("NH3", GasModel::SRK);    // Polar - use PR instead
```

**PR (Peng-Robinson):**
```cpp
// Good for:
Gas CO2_supercritical("CO2", GasModel::PR);  // Supercritical CO2 extraction
Gas ammonia("NH3", GasModel::PR);            // Refrigeration (polar)
Gas water("H2O", GasModel::PR);              // Steam (polar)

// Also good for:
Gas nitrogen_high_p("N2", GasModel::PR);     // High pressure (> 100 bar)
Gas oxygen_storage("O2", GasModel::PR);      // Compressed gas storage
```

**HELMHOLTZ:**
```cpp
// Best for:
Gas liquid_H2(HelmholtzModel::NormalHydrogen);     // Cryogenic H2 storage
Gas LNG(HelmholtzModel::Methane);            // Liquefied natural gas
Gas liquid_O2(HelmholtzModel::Oxygen);       // Rocket propellant (LOX)

// Available fluids: ParaHydrogen, NormalHydrogen, OrthoHydrogen, Oxygen, Methane, Nitrogen
// For other cryogenic fluids, use PR
```

---

## Gas Mixtures

### Creating Mixtures

**Molar fraction basis (most common):**
```cpp
Cell<string> species = {"H2", "N2", "O2"};
Vector<real> molar_fractions = {0.5, 0.3, 0.2};  // Sum = 1.0

Gas mixture(species, molar_fractions, GasModel::IDGAS);

// Properties use mixture rules
real T = 300.0, p = 101325.0;
real M_mix = mixture.M(T, p);        // Mixture molar mass in kg/mol (not g/mol)
real R_mix = mixture.R(T, p);        // Mixture gas constant
real cp_mix = mixture.cp(T, p);      // Mixture-averaged cp
```

**Mass fraction basis:**
```cpp
Vector<real> mass_fractions = {0.1, 0.7, 0.2};  // Sum = 1.0

// Create dummy, then remix by mass
Gas mixture(species, {1.0, 0.0, 0.0}, GasModel::IDGAS);
mixture.remix_mass(mass_fractions);

// Or use molar fractions derived from mass fractions manually
```

### Mixture Rules

**Thermodynamic properties** (additive by molar fraction):

```cpp
// Molar heat capacity (mixture rule)
// Cp_mix = Σ χᵢ·Cp_i

real T = 500.0, p = 101325.0;
real cp_mix = mixture.cp(T, p);

// Ideal-gas part, stored in the mixture heat spline by remix():
// cp_mix = Σ molar_fraction[i] * component[i]->Cp(T) / M_mix
```

That formula is the whole of `cp` only for IDGAS. SRK and PR add the departure
`cpdep(T, p) - cpdep0(T)` to the spline value, and a Helmholtz gas returns the EoS's own `cp`.

**Viscosity** (NASA RP-1311, Gordon & McBride, Eqs. 5.3, 5.5, 5.7):

```cpp
// μ_mix = Σ χᵢ·μᵢ / Σ χⱼ·Φᵢⱼ
// where Φᵢⱼ = interaction parameter

real mu_mix = mixture.mu(T, p);

// Φᵢⱼ is Wilke's form from molecular weight and viscosity ratios ( 5.5 ),
// replaced by the tabulated interaction viscosity where one exists ( 5.7 )
// Computationally expensive for large mixtures
```

**Thermal conductivity** (NASA RP-1311, Eqs. 5.4, 5.6):

```cpp
// λ_mix = Σ χᵢ·λᵢ / Σ χⱼ·Aᵢⱼ
// where Aᵢⱼ = Φᵢⱼ times the RP-1311 molar-mass correction ( 5.6 )

real lambda_mix = mixture.lambda(T, p);
```

### Remixing

**Change molar composition:**
```cpp
Gas mixture(species, initial_fractions, GasModel::IDGAS);

// After chemical reaction, update composition
Vector<real> new_fractions = compute_products(/* ... */);

mixture.remix(new_fractions);  // Updates internal state

// Properties automatically reflect new composition
real cp_new = mixture.cp(T, p);
```

**Change mass composition:**
```cpp
Vector<real> mass_fractions = {0.2, 0.5, 0.3};
mixture.remix_mass(mass_fractions);
```

**Reset to initial composition:**
```cpp
mixture.reset_mixture();  // Returns to composition at construction
```

**Performance note:**
```cpp
// Remix updates:
// - Molar mass M
// - Gas constant R
// - Molar/mass fraction arrays
// - Spline tables (if used)

// Expensive operations during remix:
// - Spline reconstruction (if remix_heat=true, remix_transport=true)
// - Mixture entropy update

// Skip the rebuilds:
mixture.remix(new_fractions,
              false,  // aRemixHeat: skip EoS remix, heat spline, critical point
              false); // aRemixTransport: skip transport splines
```

`remix()` always updates the fractions, `M` and `R`. With `aRemixHeat = false` it skips
`mEoS->remix()`, the heat-spline rebuild and the critical-point update. With
`aRemixTransport = false` it skips the transport-spline rebuild. Those rebuild methods are
private, so the only way to bring the gas back in step is another `remix( x, true, true )`.
Until then, `cp`, `h`, `s`, `mu` and `lambda` use the splines and EoS state of the old
composition alongside the new fractions, `M` and `R`. Skip the rebuilds only when no property
evaluation occurs before the next full `remix()`.

### Accessing Mixture Composition

```cpp
// Current molar fractions
const Vector<real> & chi = mixture.molar_fractions();

// Current mass fractions
const Vector<real> & zeta = mixture.mass_fractions();

// Individual fraction
real chi_N2 = mixture.molar_fraction(0);  // Component 0

// Number of components
uint n = mixture.number_of_components();

// Access component RefGas objects (borrowed pointers, owned by mixture)
Cell<gastables::RefGas*> & components = mixture.components();
gastables::RefGas * N2 = components(0);
```

---

## Property Evaluation

### State Variables

**Finding state from different input pairs:**

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Given T and p, find v
real T = 300.0, p = 101325.0;
real v = air.v(T, p);  // Specific volume (m³/kg)

// Given T and v, find p
real p_calc = air.p(T, v);  // Should equal 101325.0

// Given p and v, find T
real T_calc = air.T(p, v);  // Should equal 300.0

// Density (inverse of specific volume)
real rho = air.rho(T, p);  // = 1/v
```

**Note for real gas:**
- `v(T, p)` is closed-form (Cardano) for SRK/PR and a Newton iteration for the Helmholtz EoS
- `p(T, v)` is direct from EoS
- `T(p, v)` iterates for the Helmholtz EoS and is not implemented for SRK/PR (`BELFEM_ERROR`)

### Thermodynamic Properties

**Caloric properties:**

```cpp
real T = 400.0, p = 200000.0;

// Heat capacities
real cp = air.cp(T, p);      // Constant pressure (J/(kg·K))
real cv = air.cv(T, p);      // Constant volume (J/(kg·K))
real gamma = air.gamma(T, p);// cp/cv (dimensionless)

// Enthalpy and internal energy
real h = air.h(T, p);        // Specific enthalpy (J/kg)
real u = air.u(T, p);        // Specific internal energy (J/kg)
                             // u = h - p·v

// Entropy
real s = air.s(T, p);        // Specific entropy (J/(kg·K))

// Speed of sound
real c = air.c(T, p);        // m/s
```

**Ideal gas relations:**
```cpp
// For IDGAS:
// cp - cv = R
// γ = cp/cv
// c = √(γ·R·T)
// h = ∫cp dT + h_ref
// s = ∫(cp/T) dT - R·ln(p/p_ref) + s_ref

Gas idgas("N2", GasModel::IDGAS);
real cp = idgas.cp(300.0, 101325.0);
real cv = idgas.cv(300.0, 101325.0);
real R = idgas.R(300.0, 101325.0);

BELFEM_ASSERT(std::abs((cp - cv) - R) < 1e-6, "cp - cv = R");
```

**Real gas departure functions:**
```cpp
Gas realgas("N2", GasModel::PR);

real T = 120.0, p = 10e6;

// Total property = ideal + departure
real h = realgas.h(T, p);

// Access departure explicitly
gasmodels::EoS * eos = realgas.eos();
real h_dep = eos->hdep(T, p);
real h_ideal = /* compute from RefGas */;
real h_total = h_ideal + h_dep - eos->hdep0(T);  // Same as h: the 1 bar departure is subtracted (see README)

// Similarly for entropy, cp
real s_dep = eos->sdep(T, p);
real cp_dep = eos->cpdep(T, p);
```

### Transport Properties

```cpp
real T = 300.0, p = 101325.0;

// Dynamic viscosity
real mu = air.mu(T, p);  // Pa·s = kg/(m·s)

// Thermal conductivity
real lambda = air.lambda(T, p);  // W/(m·K)

// Prandtl number
real Pr = air.Pr(T, p);  // cp·μ/λ (dimensionless)
```

**Pressure dependence by model:**
- IDGAS: none. `mu` and `lambda` come from the mixture transport splines, which depend on `T` only
- SRK/PR: `mu` adds a departure term from a Lucas-type correlation in reduced temperature
  `T/T_crit` and reduced pressure `p/p_crit`, capped at 10 times the low-pressure value
  (`Gas::mu_dep`). `lambda` adds a departure term that is a polynomial in `v_crit/v`
  (`Gas::lambda_dep`)
- HELMHOLTZ: `mu` and `lambda` come from the fluid's own transport model, not from the
  functions above

### Thermodynamic Coefficients

```cpp
real T = 300.0, p = 101325.0;

// Thermal expansion coefficient: α = (1/v)(∂v/∂T)_p
real alpha = air.alpha(T, p);  // 1/K

// Isochoric stress coefficient: β = (1/p)(∂p/∂T)_v
real beta = air.beta(T, p);    // 1/K

// Isothermal compressibility: κ = -(1/v)(∂v/∂p)_T
real kappa = air.kappa(T, p);  // 1/Pa
```

**Ideal gas values:**
```cpp
// For IDGAS:
// α = 1/T
// β = 1/T
// κ = 1/p

Gas idgas;   // default ctor: the 14-species air mixture
real T = 300.0, p = 101325.0;
real alpha = idgas.alpha(T, p);
BELFEM_ASSERT(std::abs(alpha - 1.0/T) < 1e-9, "α = 1/T for ideal gas");
```

### Derivatives

```cpp
real T = 400.0, p = 101325.0;

// Temperature derivatives
real dcpdT = air.dcpdT(T, p);    // ∂cp/∂T (J/(kg·K²))
real dsdT = air.dsdT(T, p);      // ∂s/∂T (J/(kg·K²))

// Pressure derivatives
real dsdp = air.dsdp(T, p);      // ∂s/∂p (J/(kg·K·Pa))
real dhdp = air.dhdp(T, p);      // ∂h/∂p (J/(kg·Pa))

// Useful for sensitivity analysis and Newton solvers
```

---

## Chemical Equilibrium

### Gibbs Minimization

The `Gas` class can compute equilibrium composition by minimizing Gibbs free energy subject to elemental mass balance.

**Theory:**

At equilibrium, the total Gibbs energy is minimized:
```
G_total = Σ nᵢ·μᵢ → minimum
```
subject to:
```
Σ aᵢⱼ·nᵢ = bⱼ  (elemental mass balance)
```

where:
- `nᵢ` = moles of species i
- `μᵢ` = chemical potential of species i
- `aᵢⱼ` = formation-table entry: atoms of element j in species i, normalized by the element's
  reference species (see [Elemental Balance Check](#elemental-balance-check))
- `bⱼ` = total amount of element j, in the same normalization

### Computing Equilibrium

**Example: Hydrogen combustion**

```cpp
#include <cstdio>

Cell<string> species = {"H2", "O2", "H2O", "H", "O", "OH", "H2O2"};
Vector<real> initial_guess = {0.2, 0.1, 0.5, 0.05, 0.05, 0.05, 0.05};

Gas combustion(species, initial_guess, GasModel::IDGAS);

// Compute equilibrium at high temperature
real T = 2500.0;  // K
real p = 101325.0;

// In/out: the input fixes the element inventory, the output is the equilibrium
Vector<real> equilibrium_fractions = initial_guess;
combustion.compute_equilibrium(T, p, equilibrium_fractions);

// Apply equilibrium composition
combustion.remix(equilibrium_fractions);

// Analyze results
for (index_t i = 0; i < species.size(); ++i) {
    std::printf("%s: chi = %.6f\n",
                species(i).c_str(),
                equilibrium_fractions(i));
}
```

`compute_equilibrium( T, p, aX )` reads and writes `aX`. On entry, `aX` defines the elemental
inventory that the solver conserves (`trans( formation_table() ) * aX`) and provides the starting
point. After the inventory is taken, entries below `1e-9` are raised to that floor and the
iteration begins.

On exit, `aX` contains the equilibrium molar fractions, normalized to sum to 1. The solver does
not read the composition stored in the `Gas`, so a freshly sized vector would equilibrate an
arbitrary inventory. `remix_to_equilibrium( T, p )` passes the gas's current molar fractions to
the solver and then calls `remix()`.

The example prints with `std::printf`. BELFEM's `message( ... )` from `cl_Logger.hpp` also works,
but `cl_Gas.hpp` does not include that header. In addition, `message` writes through `gLog`,
which every executable must define.

### Remix to Equilibrium (In-Place)

```cpp
// Direct remix to equilibrium
combustion.remix_to_equilibrium(T, p);

// With control over spline updates
combustion.remix_to_equilibrium(T, p,
                                 false,  // aRemixHeat
                                 false); // aRemixTransport
```

The two flags are passed to `remix()`, so the same caveat applies. Disabling them leaves the
EoS, heat spline, critical point and transport splines on the old composition until the next
full `remix()`. With a single component, `remix_to_equilibrium` does nothing.

### Gibbs Energy and Formation Enthalpies

```cpp
// Molar Gibbs energy for each species (output vectors must be presized)
uint num_species = combustion.number_of_components();
Vector<real> gibbs(num_species);
combustion.Gibbs(T, gibbs);  // J/mol

// Temperature derivative of Gibbs energy
Vector<real> dgibbs_dT(num_species);
combustion.dGibbsdT(T, dgibbs_dT);  // J/(mol·K)

// Formation enthalpies
Vector<real> Hf(num_species);
combustion.Hf(T, Hf);  // J/mol
```

### Elemental Balance Check

```cpp
// Access formation table (species × elements)
const Matrix<real> & formation = combustion.formation_table();

// formation(i, j) = atoms of element j in species i, divided by the
//                   atoms of element j in that element's reference species

// Check total elemental mass before and after equilibrium
// Σ χᵢ·formation(i, j) should be constant for element j
```

The entries are not raw atom counts. Each column is normalized by the element's reference
species, such as `H2` for hydrogen and `O2` for oxygen (`Gas::create_formation_table`). The H
column is therefore 1 for both H2 and H2O, and the O column is 0.5 for H2O.

`Σ χᵢ·formation(i, j)` counts moles of the reference species (for example, moles of H2), not
moles of atoms. The columns follow the sorted element names.

### Common Equilibrium Scenarios

**1. Combustion products:**
```cpp
Cell<string> products = {"CO2", "H2O", "N2", "O2", "CO", "H2", "NO", "OH"};
Vector<real> guess = {0.1, 0.2, 0.6, 0.05, 0.01, 0.01, 0.01, 0.01};

Gas exhaust(products, guess, GasModel::IDGAS);
exhaust.remix_to_equilibrium(1800.0, 101325.0);
```

**2. Dissociation at high temperature:**
```cpp
Cell<string> dissoc = {"H2", "H", "O2", "O", "N2", "N"};
Vector<real> initial = {0.4, 0.1, 0.2, 0.1, 0.15, 0.05};

Gas hot_gas(dissoc, initial, GasModel::IDGAS);

// Dissociation increases with temperature
for (real T = 1000.0; T <= 5000.0; T += 500.0) {
    hot_gas.remix_to_equilibrium(T, 101325.0);
    real chi_H = hot_gas.molar_fraction(1);  // Atomic hydrogen fraction
    // chi_H increases with T
}
```

**3. Cryogenic equilibrium (ortho/para hydrogen):**
```cpp
// For specialized applications
// (Currently not implemented in BELFEM - would require custom species)
```

---

## Compressible Flow Applications

### Isentropic Relations

**Finding state after isentropic compression/expansion:**

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Initial state
real T1 = 300.0, p1 = 101325.0;

// Find T after isentropic compression to p2
real p2 = 500000.0;
real T2 = air.isen_T(T1, p1, p2);

// Constant-cp estimate: T2/T1 = (p2/p1)^((γ-1)/γ)
real gamma = air.gamma(T1, p1);
real T2_estimate = T1 * std::pow(p2/p1, (gamma-1)/gamma);

// Similarly, find p after isentropic expansion to T2
real T2_target = 250.0;
real p2_calc = air.isen_p(T1, p1, T2_target);
```

The closed form assumes constant `cp` and is only an estimate. `isen_T` uses
`T0 * (p1/p0)^(R/cp)` as its starting value, then iterates until `s(T1, p1) = s(T0, p0)`. Because
`cp` depends on temperature, even for the ideal-gas air mixture, `T2` and `T2_estimate` differ.

**Isentropic efficiency:**
```cpp
// Compressor with non-ideal efficiency
real T1 = 300.0, p1 = 101325.0, p2 = 500000.0;

real T2_isentropic = air.isen_T(T1, p1, p2);  // Ideal

real eta_c = 0.85;  // Compressor efficiency
real T2_actual = T1 + (T2_isentropic - T1) / eta_c;

// Turbine with non-ideal efficiency
real T3 = 1500.0, p3 = 500000.0, p4 = 101325.0;

real T4_isentropic = air.isen_T(T3, p3, p4);  // Ideal

real eta_t = 0.90;  // Turbine efficiency
real T4_actual = T3 - eta_t * (T3 - T4_isentropic);
```

### Total (Stagnation) Conditions

**Computing stagnation temperature and pressure:**

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Static conditions
real T_static = 288.15;  // K
real p_static = 101325.0;  // Pa
real u = 250.0;  // m/s, flow velocity (subsonic)

real T_total, p_total;
air.total(T_static, p_static, u, T_total, p_total);

// Constant-cp estimates (starting values only):
// T_total ≈ T_static + u²/(2·cp)
// p_total/p_static ≈ (T_total/T_static)^(γ/(γ-1))
```

`total()` starts from those perfect-gas estimates and then solves
`h(T_total, p_total) = h(T_static, p_static) + u²/2` together with
`s(T_total, p_total) = s(T_static, p_static)`. With temperature-dependent `cp`, the returned values
differ from the closed forms.

**Mach number from static and total conditions (perfect-gas estimate):**
```cpp
real gamma = air.gamma(T_static, p_static);

// T_total/T_static = 1 + (γ-1)/2 · Ma²
real Ma2 = (T_total/T_static - 1.0) * 2.0 / (gamma - 1.0);
real Ma = std::sqrt(Ma2);

// Exact for the gas model: Ma = u / c
real Ma_exact = u / air.c(T_static, p_static);
```

### Normal Shock Relations

**Analyzing a normal shock:**

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Upstream (before shock) - supersonic
real T1 = 300.0, p1 = 101325.0, u1 = 600.0;  // M ≈ 1.7

// Downstream (after shock)
real T2, p2, u2;
air.shock(T1, p1, u1, T2, p2, u2);

// Verify shock relations:
// - Entropy increase: s2 > s1
// - Temperature jump: T2 > T1
// - Pressure jump: p2 > p1
// - Velocity decrease: u2 < u1 (now subsonic)

real s1 = air.s(T1, p1);
real s2 = air.s(T2, p2);
BELFEM_ASSERT(s2 > s1, "Entropy must increase across shock");

// Mach number calculation
real c1 = air.c(T1, p1);
real Ma1 = u1 / c1;  // Supersonic

real c2 = air.c(T2, p2);
real Ma2 = u2 / c2;  // Subsonic

std::printf("Normal shock: Ma1 = %.3f -> Ma2 = %.3f, p2/p1 = %.3f\n",
            Ma1, Ma2, p2/p1);
```

### Oblique Shock

**Computing oblique shock with flow deflection:**

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Upstream conditions
real T1 = 300.0, p1 = 101325.0, u1 = 700.0;  // M ≈ 2.0

// Flow deflection angle
real alpha = 20.0 * constant::deg;  // 20 degrees

// Solve for shock angle and downstream conditions
real T2, p2, u2, beta;
air.shock(T1, p1, u1, alpha, T2, p2, u2, beta);

// beta = shock angle (relative to upstream flow)
// Flow deflects by alpha
// Tangential velocity component unchanged
// Normal component undergoes normal shock relations

std::printf("Oblique shock: alpha = %.1f deg, beta = %.1f deg\n",
            alpha / constant::deg,
            beta / constant::deg);
```

### Prandtl-Meyer Expansion

**Supersonic expansion around a corner:**

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Upstream conditions (supersonic)
real T1 = 300.0, p1 = 101325.0, u1 = 500.0;

// Turning angle
real deflection = 15.0 * constant::deg;  // expansion ( positive ); a negative angle is a smooth isentropic compression

real T2, p2, u2;
real Ma2 = air.prandtl_meyer(T1, p1, u1, deflection, T2, p2, u2);

// Expansion fan:
// - T2 < T1 (cooling)
// - p2 < p1 (pressure drop)
// - u2 > u1 (acceleration)
// - Process is isentropic

real s1 = air.s(T1, p1);
real s2 = air.s(T2, p2);
BELFEM_ASSERT(std::abs(s2 - s1) < 1e-6, "Isentropic process");
```

The closed-form perfect-gas angle is a private helper that only seeds the
iteration; `prandtl_meyer()` returns the downstream Mach number.

`prandtl_meyer()` aborts with a `BELFEM_ERROR` unless all three preconditions hold: the gas uses
the ideal-gas model (`is_idgas()`), `T1` is above the caloric table limit `gastables::gTmin`, and
the upstream flow is supersonic (`u1 > c(T1, p1)`). The example above, default air at 500 m/s
and 300 K, meets all three.

### Duct Area Changes

`expand()` and `compress()` compute the state after a sudden change of duct cross section from
`A1` to `A2`:

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

real A1 = 1.0, A2 = 2.0;                        // m², A2 >= A1 for expand
real T1 = 300.0, p1 = 101325.0, u1 = 100.0;
real T2, p2, u2;
air.expand(A1, T1, p1, u1, A2, T2, p2, u2);

real A3 = 1.5;                                  // m², A3 <= A2 for compress
real T3, p3, u3;
air.compress(A2, T2, p2, u2, A3, T3, p3, u3);
```

- `expand()`: a subsonic flow separates at the step and loses total pressure through a
  Borda-Carnot shock, so the routine solves for mass, momentum and energy conservation. A
  supersonic flow expands around the corner instead, so the routine solves for mass, entropy
  and energy conservation.
- `compress()`: the flow stays isentropic on both the subsonic and supersonic branches. It aborts
  if the duct is choked, that is, if `A2` is smaller than the sonic area.

`T_from_h( h, p )` inverts the enthalpy and returns the temperature at which `h(T, p)` equals
the given `h`.

### Area-Mach Number Relation (Isentropic Flow)

```cpp
// Compute area ratio for isentropic nozzle
real M1 = 0.5;  // Subsonic
real gamma = air.gamma(300.0, 101325.0);

// A/A* = (1/M) * ((2/(γ+1)) * (1 + (γ-1)/2 * M²))^((γ+1)/(2(γ-1)))
real term = (2.0/(gamma+1.0)) * (1.0 + (gamma-1.0)/2.0 * M1*M1);
real A_over_Astar = (1.0/M1) * std::pow(term, (gamma+1.0)/(2.0*(gamma-1.0)));

// For M = 1 (throat): A/A* = 1.0
// For M < 1 (subsonic): A/A* > 1 (converging)
// For M > 1 (supersonic): A/A* > 1 (diverging)
```

---

## Advanced Topics

### Component-Specific Properties

For mixtures, access properties of individual components:

```cpp
Cell<string> species = {"N2", "O2", "Ar"};
Vector<real> fractions = {0.78, 0.21, 0.01};
Gas air(species, fractions, GasModel::IDGAS);

real T = 400.0, p = 101325.0;

// Component 0 (N2) properties
index_t i = 0;
real v_N2 = air.v(i, T, p);     // Specific volume of N2
real h_N2 = air.h(i, T, p);     // Specific enthalpy of N2
real cp_N2 = air.cp(i, T, p);   // Specific heat of N2

// Useful for:
// - Partial molar properties
// - Fugacity calculations
// - Activity coefficient models
```

### Direct EoS Access

```cpp
Gas nitrogen("N2", GasModel::PR);

// Get underlying EoS object
gasmodels::EoS * eos = nitrogen.eos();

// Use EoS methods directly
real T = 120.0;
real v = 0.001;  // m³/kg

// State equation
real p = eos->p(T, v);  // Pressure from T, v

// Derivatives
real dpdT = eos->dpdT(T, v);    // (∂p/∂T)_v
real dpdv = eos->dpdv(T, v);    // (∂p/∂v)_T
real dvdT = eos->dvdT(T, v);    // (∂v/∂T)_p = -(∂p/∂T)_v / (∂p/∂v)_T

// Departure functions
p = 10e6;
real h_dep = eos->hdep(T, p);   // h_real - h_ideal
real s_dep = eos->sdep(T, p);   // s_real - s_ideal
real cp_dep = eos->cpdep(T, p); // cp_real - cp_ideal

// Critical point
real T_crit, p_crit, v_crit;
eos->eval_critical_point(T_crit, p_crit, v_crit);
```

### Helmholtz Derivatives

For Helmholtz EoS, access fundamental equation derivatives:

```cpp
#include "cl_Gas.hpp"
#include "cl_GM_Helmholtz.hpp"   // cl_Gas.hpp does not declare gasmodels::Helmholtz

Gas methane(HelmholtzModel::Methane);

// The EoS is actually a Helmholtz object (borrowed pointer, owned by methane)
gasmodels::Helmholtz * helm = dynamic_cast<gasmodels::Helmholtz*>(methane.eos());

if (helm != nullptr) {
    real T = 120.0, p = 5e6;

    // High-level properties
    real h = helm->h(T, p);
    real s = helm->s(T, p);
    real cp = helm->cp(T, p);
    real w = helm->w(T, p);  // Speed of sound

    // Reduced Helmholtz energy phi = phi0 + phir and its derivatives
    // (tau = T_crit/T, delta = v_crit/v): read from the state set by the last
    // update_Tp / update_Tv, so set the state first
    helm->update_Tp(T, p);
    real phir_d  = helm->phir_d();   // ∂phir/∂delta
    real phir_t  = helm->phir_t();   // ∂phir/∂tau
    real phir_dt = helm->phir_dt();  // ∂²phir/∂delta∂tau
    real phi0_t  = helm->phi0_t();   // ∂phi0/∂tau
}
```

The public reduced-energy accessors are `phi0()`, `phi0_t()`, `phi0_tt()`, `phir()`, `phir_d()`,
`phir_dd()`, `phir_t()`, `phir_tt()` and `phir_dt()`. They take no arguments. `helm->alpha( T, p )`
is the thermal expansion coefficient, not the reduced Helmholtz energy.

### Vapor Pressure Curves

For Helmholtz models, vapor pressure correlation:

```cpp
Gas methane(HelmholtzModel::Methane);
gasmodels::EoS * eos = methane.eos();

// Vapor pressure at given temperature
real T_sat = 111.0;  // K (< T_crit)
real p_vap = eos->p_vap(T_sat);

// Saturation temperature at given pressure
real p_sat = 101325.0;  // Pa
real T_vap = eos->T_vap(p_sat);

// Check: p_vap(T_vap(p)) ≈ p
```

---

## Performance Optimization

### Minimize Object Creation

Following BELFEM's manual memory philosophy:

```cpp
// BAD - creates Gas repeatedly
for (int i = 0; i < 10000; ++i) {
    Gas air;   // default ctor builds the 14-species air mixture as an ideal gas  // Expensive!
    real cp = air.cp(300.0, 101325.0);
}

// GOOD - reuse Gas object
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas  // Create once
for (int i = 0; i < 10000; ++i) {
    real cp = air.cp(300.0 + i*0.1, 101325.0);  // Just evaluate
}
```

### Spline Evaluation

For IDGAS, SRK and PR, `Gas` evaluates its own mixture splines (`mHeatSpline`, `mViscositySpline`, `mConductivitySpline`), which `remix()` rebuilds. There is no mode switch, and the `RefGas` mode of the components does not affect these splines. SRK and PR add EoS departure terms to the spline values.

A HELMHOLTZ gas does not use these splines for `cp`, `mu` or `lambda`. Its `cp` comes from the EoS, while `mu` and `lambda` come from the fluid's transport model.

### Preallocate Work Arrays

For mixture calculations, Gas preallocates work arrays:

```cpp
// Internally in Gas constructor:
// mWorkMatrix(num_species, num_species);
// mWorkVector(num_species);
// mWorkMu( gNumberOfSplinePoints );
// mWorkLambda( gNumberOfSplinePoints );   ( one entry per spline knot, filled by remix_transport )

// These are reused across property evaluations
// No allocation during property calls
```

### Avoid Unnecessary Remixing

```cpp
// BAD - remixes unnecessarily
for (int i = 0; i < 1000; ++i) {
    mixture.remix(same_fractions);  // Wasteful
    real cp = mixture.cp(T, p);
}

// GOOD - remix only when composition changes
mixture.remix(new_fractions);
for (int i = 0; i < 1000; ++i) {
    real cp = mixture.cp(T + i*0.1, p);  // No remix
}
```

### Use Appropriate EoS for Speed

Performance ranking (fastest to slowest):
1. **IDGAS** - Mixture spline evaluation
2. **SRK** - Cubic solution (closed-form Cardano root)
3. **PR** - Cubic solution (closed-form Cardano root, same cost as SRK)
4. **HELMHOLTZ** - Complex derivatives (10-100× slower)

```cpp
// For performance-critical loop
Gas fast_air;   // default ctor: the 14-species air mixture  // Fastest

// For accuracy-critical calculation
Gas accurate_H2(HelmholtzModel::NormalHydrogen);  // Most accurate
```

### Benchmark Example

```cpp
#include <cstdio>
#include "cl_Timer.hpp"

Timer timer;   // the constructor starts the clock

// Benchmark IDGAS
Gas idgas("N2", GasModel::IDGAS);
timer.reset();
for (int i = 0; i < 100000; ++i) {
    real cp = idgas.cp(300.0 + i*0.001, 101325.0);
}
uint64_t time_idgas = timer.next();   // ms, and restarts the clock

// Benchmark PR
Gas pr("N2", GasModel::PR);
for (int i = 0; i < 100000; ++i) {
    real cp = pr.cp(300.0 + i*0.001, 101325.0);
}
uint64_t time_pr = timer.stop();      // ms

std::printf("IDGAS: %llu ms, PR: %llu ms, Ratio: %.2f\n",
            (unsigned long long)time_idgas, (unsigned long long)time_pr,
            (real)time_pr / time_idgas);

// Typical result: PR is 3-5× slower than IDGAS
```

---

## Troubleshooting

### Issue: Convergence Failure in v(T,p) (Helmholtz EoS only)

```
BELFEM_ERROR: Too many iterations for T=... K, p=... bar, rho=... kg/m^3, relax=...
```

**Cause:** Near critical point or phase boundary, EoS may have multiple solutions or poor conditioning.

**Solution:**
```cpp
// Provide better initial guess
// (the Helmholtz Newton starts from the SRK volume, or a liquid-volume polynomial below the vapor curve)

// Or switch to direct iteration using p(T, v)
real v_guess = /* estimate */;
for (int iter = 0; iter < 100; ++iter) {
    real p_calc = eos->p(T, v_guess);
    real dpdv = eos->dpdv(T, v_guess);
    v_guess -= (p_calc - p) / dpdv;  // Newton step
    if (std::abs(p_calc - p) < 1e-6) break;
}
```

### Issue: Negative Heat Capacity

No check exists; a negative cp from a cubic EoS inside the two-phase dome is returned as is.

**Cause:** Unphysical state (inside two-phase region or extrapolation beyond valid range).

**Solution:**
```cpp
real T = 100.0, p = 1e6;

Gas gas("N2", GasModel::PR);
gasmodels::EoS * eos = gas.eos();   // borrowed pointer, owned by gas

real T_crit, p_crit, v_crit;
eos->eval_critical_point(T_crit, p_crit, v_crit);

if (T < T_crit) {
    // Subcritical: the cubic may have three real roots here, and gas.v(T, p)
    // returns the largest (vapour) one regardless of p
}
```

Below `T_crit`, a cubic-EoS `Gas` always returns the vapour root, because nothing in the module
sets its liquid flag. The module has no saturation test for SRK or PR. `p_vap` and `T_vap` are
implemented only by the Helmholtz models; the `EoS` base class aborts with
`p_vap not implemented for this EoS`. A Helmholtz model classifies a state as liquid when `T < T_triple`, or
when `T < T_crit` and `p > p_vap(T)` (`Helmholtz::is_liquid`). For N2, O2, CH4 and H2, the
Helmholtz models are the ones that distinguish liquid from vapour.

### Note: Mixture Composition Need Not Sum to 1.0

Molar and mass fractions are normalized to unity inside `remix()` /
`remix_mass()`; a non-normalised input is accepted.

### Issue: Equilibrium Calculation Fails

```
BELFEM_ERROR: To many iterations while trying to find chemical equilibrium.
```

**Cause:** Poor initial guess or elemental imbalance.

**Solution:**
```cpp
// Check elemental balance
Cell<string> species = {"H2", "O2", "H2O"};
Vector<real> fractions = {0.5, 0.25, 0.25};

Gas mixture(species, fractions, GasModel::IDGAS);

const Matrix<real> & formation = mixture.formation_table();
// Column 0 is H (reference species H2), column 1 is O (reference species O2):
// formation(0, 0) = 1   (H2:  2 H atoms / 2 in H2)
// formation(1, 0) = 0   (O2:  no H)
// formation(2, 0) = 1   (H2O: 2 H atoms / 2 in H2)
// formation(2, 1) = 0.5 (H2O: 1 O atom  / 2 in O2)

// Hydrogen inventory in moles of H2: Σ χᵢ·formation(i, 0)
// Oxygen inventory in moles of O2:   Σ χᵢ·formation(i, 1)

// Ensure these remain constant during equilibrium
```

---

## Further Reading

### Module Documentation
- **`README.md`** - Quick reference for gasmodels module

### Related Modules
- **`src/physics/gastables/doc/gastables_usage_guide.md`** - Reference gas data (foundation)
- **`src/numerics/spline/doc/spline_usage_guide.md`** - Spline interpolation

### Literature References

**Equations of State:**
- Soave (1972), "Equilibrium Constants from a Modified Redlich-Kwong Equation of State", *Chemical Engineering Science*, 27(6):1197-1203
- Peng & Robinson (1976), "A New Two-Constant Equation of State", *Industrial & Engineering Chemistry Fundamentals*, 15(1):59-64
- Reid, Prausnitz & Poling (1987), "The Properties of Gases and Liquids" (4th ed.), McGraw-Hill

**Helmholtz Models:**
- Leachman et al. (2009), "Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen", *J. Phys. Chem. Ref. Data*, 38(3):721-748
- Setzmann & Wagner (1991), "A New Equation of State and Tables of Thermodynamic Properties for Methane", *J. Phys. Chem. Ref. Data*, 20(6):1061-1155
- Schmidt & Wagner (1985), "A New Form of the Equation of State for Pure Substances", *Fluid Phase Equilibria*, 19(3):175-200
- Span et al. (2000), "A Reference Equation of State for the Thermodynamic Properties of Nitrogen", *J. Phys. Chem. Ref. Data*, 29(6):1361-1433

**Mixture Rules:**
- Gordon & McBride (1994), NASA RP-1311, Eqs. (5.3)-(5.7)
- Wilke (1950), "A Viscosity Equation for Gas Mixtures", *J. Chem. Phys.*, 18:517-519
- Bird, Stewart & Lightfoot (2007), "Transport Phenomena" (2nd ed.), Wiley

**Chemical Equilibrium:**
- Smith & Missen (1982), "Chemical Reaction Equilibrium Analysis", Wiley-Interscience
- Gordon & McBride (1994), "Computer Program for Calculation of Complex Chemical Equilibrium Compositions", NASA RP-1311

**Compressible Flow:**
- Anderson (2003), "Modern Compressible Flow" (3rd ed.), McGraw-Hill
- Shapiro (1953), "The Dynamics and Thermodynamics of Compressible Fluid Flow", Ronald Press

---

**Last Updated:** 2026-10-22
**Maintainer:** BELFEM development team
