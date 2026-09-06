# Gas Models Module Documentation {#physics_gasmodels_index}

**Date:** 2026-01-30
**Purpose:** Quick reference and navigation guide for the gasmodels module
**Module:** `src/physics/gasmodels`

---

## Overview

The **gasmodels** module provides comprehensive **fluid property models** for gases, combining thermodynamic and transport properties with various **equations of state** (EoS). It builds on the `gastables` module to create complete gas models suitable for:

- Compressible flow analysis
- Chemical equilibrium calculations
- Cryogenic fluid simulations
- Real gas thermodynamics
- Gas mixture modeling

**Key features:**
- **Multiple EoS options:** Ideal gas, SRK, Peng-Robinson, Helmholtz free energy
- **Gas mixtures** with arbitrary composition and remixing capability
- **Chemical equilibrium** via Gibbs minimization
- **Compressible flow utilities:** Shock analysis, Prandtl-Meyer expansion, isentropic relations
- **Transport properties** with the NASA RP-1311 (Gordon & McBride) mixing rules, Eqs. (5.3)-(5.7): Wilke-type Φᵢⱼ with tabulated interaction viscosities where available, and the RP-1311 conductivity correction
- **Cryogenic models:** specialized Helmholtz EoS for H2 (para/normal/ortho), O2, CH4 and N2

---

## Quick Reference

### Key Classes

| Class | Purpose | Location |
|-------|---------|----------|
| `Gas` | Main interface for fluid properties and EoS | `cl_Gas.hpp:80` |
| `EoS` | Abstract base class for equations of state | `cl_GM_EoS.hpp:33` |
| `EoS_Idgas` | Ideal gas equation of state | `cl_GM_EoS_Idgas.hpp:32` |
| `EoS_Cubic` | Cubic equations of state (SRK, PR) | `cl_GM_EoS_Cubic.hpp:37` |
| `Helmholtz` | Helmholtz free energy EoS (cryogenic) | `cl_GM_Helmholtz.hpp:59` |
| `EoS_Hydrogen` | Specialized Helmholtz EoS for H2 | `cl_GM_EoS_Hydrogen.hpp:43` |
| `EoS_Methane` | Specialized Helmholtz EoS for CH4 | `cl_GM_EoS_Methane.hpp:29` |
| `EoS_Oxygen` | Specialized Helmholtz EoS for O2 | `cl_GM_EoS_Oxygen.hpp:27` |
| `Statevals` | Container for state variables | `cl_GM_Statevals.hpp:86` |

### Key Enums

| Enum | Values | Purpose | Location |
|------|--------|---------|----------|
| `GasModel` | `IDGAS`, `SRK`, `PR`, `HELMHOLTZ` | Select EoS type | `en_GM_GasModel.hpp:17` |
| `HelmholtzModel` | `ParaHydrogen`, `NormalHydrogen`, `OrthoHydrogen`, `Oxygen`, `Methane`, `Nitrogen` | Select Helmholtz fluid | `en_Helmholtz.hpp:17-26` |

### Key Files

| File | Purpose |
|------|---------|
| `cl_Gas.{hpp,cpp}` | Main Gas class with all property methods |
| `cl_GM_EoS.{hpp,cpp}` | Abstract EoS base class |
| `cl_GM_EoS_Idgas.{hpp,cpp}` | Ideal gas implementation |
| `cl_GM_EoS_Cubic.{hpp,cpp}` | SRK and PR cubic EoS |
| `cl_GM_Helmholtz.{hpp,cpp}` | Helmholtz free energy EoS |
| `cl_GM_EoS_Hydrogen.{hpp,cpp}` | H2 reference EoS (cryogenic) |
| `cl_GM_EoS_Methane.{hpp,cpp}` | CH4 reference EoS (cryogenic) |
| `cl_GM_EoS_Oxygen.{hpp,cpp}` | O2 reference EoS (cryogenic) |
| `cl_GM_EoS_AlphaFunction.{hpp,cpp}` | Alpha functions for cubic EoS temperature correction |
| `cl_GM_Statevals.hpp` | State variable container ( header only ) |

---

## Common Operations

### Create an Ideal Gas

```cpp
#include "cl_Gas.hpp"

using namespace belfem;

// Pure gas, ideal gas model
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Evaluate properties at T=300 K, p=101325 Pa
real T = 300.0;
real p = 101325.0;

real cp = air.cp(T, p);        // Specific heat (J/(kg·K))
real cv = air.cv(T, p);        // Constant volume heat capacity
real gamma = air.gamma(T, p);  // Heat capacity ratio
real h = air.h(T, p);          // Specific enthalpy (J/kg)
real s = air.s(T, p);          // Specific entropy (J/(kg·K))
real rho = air.rho(T, p);      // Density (kg/m³)
real mu = air.mu(T, p);        // Dynamic viscosity (Pa·s)
real lambda = air.lambda(T, p);// Thermal conductivity (W/(m·K))
```

### Create a Gas Mixture

```cpp
// Define composition
Cell<string> species = {"N2", "O2", "Ar", "CO2"};
Vector<real> molar_fractions = {0.78, 0.21, 0.01, 0.0003};

// Create mixture
Gas mixture(species, molar_fractions, GasModel::IDGAS);

// Properties use mixture rules
real T = 400.0, p = 200000.0;
real cp_mix = mixture.cp(T, p);  // Weighted average of component cp values
real mu_mix = mixture.mu(T, p);  // RP-1311 mixing rule
```

### Use Real Gas (Cubic EoS)

```cpp
// Peng-Robinson equation of state for nitrogen
Gas nitrogen("N2", GasModel::PR);

real T = 100.0;   // Low temperature (cryogenic)
real p = 5e6;     // High pressure (50 bar)

// Real gas behavior
real v = nitrogen.v(T, p);       // Specific volume (m³/kg)
real Z = p * v / (nitrogen.R(T,p) * T);  // Compressibility factor (Z ≠ 1)

real cp = nitrogen.cp(T, p);     // Includes departure function
real h = nitrogen.h(T, p);       // Enthalpy with real gas correction
real s = nitrogen.s(T, p);       // Entropy with real gas correction
```

### Use Helmholtz EoS (Cryogenic)

```cpp
// High-accuracy methane model for cryogenic applications
Gas methane(HelmholtzModel::Methane);

// Valid over wide range including liquid phase
real T = 120.0;  // K (below critical point)
real p = 5e6;    // Pa

real v = methane.v(T, p);        // May be liquid or vapor
real h = methane.h(T, p);        // Accurate across phase boundary
real cp = methane.cp(T, p);      // From Helmholtz derivatives

// Check phase
if (methane.is_liquid()) {
    // In liquid phase
}
```

---

## Gas Models (Equations of State)

### IDGAS - Ideal Gas

**Equation of state:**
```
p·v = R·T
```

**Use when:**
- Pressure < 10 bar
- Temperature well above critical point
- Fast approximate calculations

**Properties:**
- Thermodynamic: From `RefGas` polynomials (NASA CEA)
- Transport: From `RefGas` polynomials with mixture rules
- No departure functions (ideal behavior)

**Example:**
```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas
real v = air.v(300.0, 101325.0);  // v = R·T/p
```

### SRK - Soave-Redlich-Kwong

**Equation of state:**
```
p = R·T/(v-b) - a(T)/(v(v+b))
```

**Use when:**
- Moderate pressure (< 100 bar)
- Non-polar or weakly polar molecules
- Good accuracy/speed trade-off

**Properties:**
- Departure functions for h, s, cp from EoS
- Requires critical point data and acentric factor

**Example:**
```cpp
Gas ethane("C2H6", GasModel::SRK);
real v = ethane.v(300.0, 2e6);   // Real gas volume
```

### PR - Peng-Robinson

**Equation of state:**
```
p = R·T/(v-b) - a(T)/(v(v+b) + b(v-b))
```

**Use when:**
- High pressure (up to 1000 bar)
- Polar molecules
- Better accuracy than SRK near critical point

**Properties:**
- More accurate density prediction than SRK
- Especially good for liquids and near-critical fluids

**Example:**
```cpp
Gas CO2("CO2", GasModel::PR);
real v_subcritical = CO2.v(250.0, 5e6);  // Near critical point
```

### HELMHOLTZ - Fundamental Equation

**Equation of state:**
```
a(ρ, T) = a_ideal(ρ, T) + a_residual(ρ, T)
```
where `a` is Helmholtz free energy.

**Use when:**
- Cryogenic applications (H2, CH4, O2, N2)
- Maximum accuracy required
- Two-phase calculations
- Wide pressure/temperature range

**Available fluids:**
- `HelmholtzModel::NormalHydrogen`, `ParaHydrogen`, `OrthoHydrogen` - Leachman et al. (2009)
- `HelmholtzModel::Methane` - Setzmann & Wagner (1991)
- `HelmholtzModel::Oxygen` - Schmidt & Wagner (1985)
- `HelmholtzModel::Nitrogen` - Span et al. (2000)

**Properties:**
- Highly accurate reference equations
- Valid from triple point to high temperatures
- Can handle liquid, vapor, supercritical regions

**Example:**
```cpp
Gas hydrogen(HelmholtzModel::NormalHydrogen);

// Valid over huge range
real cp_cryo = hydrogen.cp(20.0, 101325.0);      // Liquid H2
real cp_high = hydrogen.cp(300.0, 100e6);        // High pressure
```

---

## Reference Pressure Convention

The heat splines hold the CEA standard state: the **ideal gas at 1 bar**. The ideal gas model evaluates them directly and never calls a departure function.

The cubic models add the departure at the requested state and subtract the departure at 1 bar:

```text
cp( T, p ) = cp_spline( T ) + cpdep( T, p ) - cpdep( T, 1 bar )
```

At 1 bar the two departure terms cancel, so **SRK and PR reduce exactly to the ideal gas model at the reference pressure**. The subtraction enforces that continuity, and the departure splines are rebuilt whenever the gas model changes. Helmholtz does not use this path because it takes caloric properties from the equation of state.

Two consequences matter:

- At pressures other than 1 bar, the result differs from the true real gas property by the departure at 1 bar. For nitrogen that is about 6 % in cp near 80 K, falling below 0.2 % above 250 K, so it matters for cryogenic work and is negligible at ambient conditions.
- The per-component functions `Gas::cp( aIndex, ... )` and `Gas::h( aIndex, ... )` apply the same convention, carried inside the indexed departure functions. The finite rate combustion solver depends on this because it pairs those enthalpies with the ideal gas Gibbs energy.

BELFEM has not settled whether to keep this continuity convention or instead store real gas values at 1 bar and drop the subtraction. Both choices are self consistent. The current code combines the ideal gas spline with the reference pressure subtraction. See the note on `Gas::realgas_cp` in `cl_Gas.cpp`.

## Property Evaluation API

### Thermodynamic State

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `p(T, v)` | Pressure from T and specific volume | Pa | From EoS |
| `v(T, p)` | Specific volume from T and pressure | m³/kg | Iterative for real gas |
| `T(p, v)` | Temperature from p and specific volume | K | Iterative |
| `rho(T, p)` | Density | kg/m³ | `1/v(T,p)` |

### Caloric Properties

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `cp(T, p)` | Specific heat at constant pressure | J/(kg·K) | From RefGas + departure |
| `cv(T, p)` | Specific heat at constant volume | J/(kg·K) | From relations or EoS |
| `gamma(T, p)` | Heat capacity ratio (cp/cv) | - | Dimensionless |
| `h(T, p)` | Specific enthalpy | J/kg | From RefGas + departure |
| `u(T, p)` | Specific internal energy | J/kg | `h - p·v` |
| `s(T, p)` | Specific entropy | J/(kg·K) | From RefGas + departure |
| `c(T, p)` | Speed of sound | m/s | From thermodynamic relations |

### Transport Properties

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `mu(T, p)` | Dynamic viscosity | Pa·s | RP-1311 Eqs. (5.3), (5.5), (5.7) for mixtures |
| `lambda(T, p)` | Thermal conductivity | W/(m·K) | RP-1311 Eqs. (5.4), (5.6) for mixtures |
| `Pr(T, p)` | Prandtl number | - | `cp·μ/λ` |

### Thermodynamic Coefficients

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `alpha(T, p)` | Thermal expansion coefficient | 1/K | `(1/v)(∂v/∂T)_p` |
| `beta(T, p)` | Isochoric stress coefficient | 1/K | `(1/p)(∂p/∂T)_v` |
| `kappa(T, p)` | Isothermal compressibility | 1/Pa | `-(1/v)(∂v/∂p)_T` |

### Derivatives

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `dcpdT(T, p)` | Temperature derivative of cp | J/(kg·K²) | For sensitivity |
| `dsdT(T, p)` | Temperature derivative of s | J/(kg·K²) | |
| `dsdp(T, p)` | Pressure derivative of s | J/(kg·K·Pa) | |
| `dhdp(T, p)` | Pressure derivative of h | J/(kg·Pa) | For total temperature |

---

## Gas Mixture Operations

### Remixing

```cpp
// Create initial mixture
Cell<string> species = {"H2", "O2", "N2"};
Vector<real> initial_fractions = {0.5, 0.25, 0.25};

Gas mixture(species, initial_fractions, GasModel::IDGAS);

// Change composition (e.g., after chemical reaction)
Vector<real> new_fractions = {0.3, 0.1, 0.6};
mixture.remix(new_fractions);  // Updates molar fractions

// Properties automatically recalculated
real cp_new = mixture.cp(300.0, 101325.0);
```

### Remix by Mass Fractions

```cpp
// Specify composition by mass instead of moles
Vector<real> mass_fractions = {0.1, 0.8, 0.1};  // Mass basis
mixture.remix_mass(mass_fractions);

// Access current composition
const Vector<real> & molar = mixture.molar_fractions();
const Vector<real> & mass = mixture.mass_fractions();
```

### Reset to Initial Composition

```cpp
// Return to composition at construction
mixture.reset_mixture();
```

---

## Chemical Equilibrium

### Compute Equilibrium Composition

```cpp
Cell<string> species = {"H2", "O2", "H2O", "H", "O", "OH"};
Vector<real> initial_guess = {0.3, 0.15, 0.5, 0.02, 0.01, 0.02};

Gas combustion(species, initial_guess, GasModel::IDGAS);

real T = 2500.0;  // High temperature
real p = 101325.0;

// Compute equilibrium composition at T, p
Vector<real> equilibrium_fractions(species.size());
combustion.compute_equilibrium(T, p, equilibrium_fractions);

// Apply equilibrium composition
combustion.remix(equilibrium_fractions);
```

### Remix to Equilibrium (In-Place)

```cpp
// Direct remix to equilibrium at given T, p
combustion.remix_to_equilibrium(T, p,
                                 true,   // Remix heat properties
                                 true);  // Remix transport properties
```

### Gibbs Energy

```cpp
Vector<real> gibbs(num_species);
combustion.Gibbs(T, gibbs);  // Molar Gibbs energy for each species (J/mol)

Vector<real> dgibbs_dT(num_species);
combustion.dGibbsdT(T, dgibbs_dT);  // Temperature derivative

Vector<real> formation_enthalpies(num_species);
combustion.Hf(T, formation_enthalpies);  // Formation enthalpies (J/mol)
```

---

## Compressible Flow Utilities

### Isentropic Relations

```cpp
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// State 1
real T1 = 300.0, p1 = 101325.0;

// Find T2 for isentropic compression to p2
real p2 = 500000.0;
real T2 = air.isen_T(T1, p1, p2);

// Or find p2 for isentropic expansion to T2
real T2_target = 250.0;
real p2_calculated = air.isen_p(T1, p1, T2_target);
```

### Total (Stagnation) Conditions

```cpp
real T_static = 300.0;
real p_static = 101325.0;
real velocity = 200.0;  // m/s

real T_total, p_total;
air.total(T_static, p_static, velocity, T_total, p_total);

// T_total = T + U²/(2·cp)
// p_total from isentropic relation
```

### Normal Shock

```cpp
// Upstream conditions
real T1 = 300.0, p1 = 101325.0, U1 = 500.0;  // Supersonic

// Downstream conditions (after shock)
real T2, p2, U2;
air.shock(T1, p1, U1, T2, p2, U2);

// T2 > T1 (entropy increase)
// p2 > p1 (pressure jump)
// U2 < U1 (velocity decrease)
```

### Oblique Shock

```cpp
real alpha = 30.0 * constant::deg;  // Flow deflection angle

real T1 = 300.0, p1 = 101325.0, U1 = 600.0;
real T2, p2, U2, beta;

air.shock(T1, p1, U1, alpha, T2, p2, U2, beta);

// beta = shock angle
// Flow deflects by alpha
```

### Prandtl-Meyer Expansion

```cpp
// Supersonic expansion fan
real T1 = 300.0, p1 = 101325.0, U1 = 500.0;
real deflection = 15.0 * constant::deg;  // expansion ( positive ); a negative angle is a smooth isentropic compression

real T2, p2, U2;
air.prandtl_meyer(T1, p1, U1, deflection, T2, p2, U2);

// T2 < T1 (cooling)
// p2 < p1 (pressure drop)
// U2 > U1 (acceleration)
```

The closed-form perfect-gas angle is a private helper that only seeds the
iteration; `prandtl_meyer()` returns the downstream Mach number.

---

## Advanced Features

### Access Component Properties

For gas mixtures, access individual component properties:

```cpp
Cell<string> species = {"N2", "O2"};
Vector<real> fractions = {0.79, 0.21};
Gas air(species, fractions, GasModel::IDGAS);

// Access component 0 (N2)
index_t i = 0;
real cp_N2 = air.cp(i, 300.0, 101325.0);  // Component-specific cp
real h_N2 = air.h(i, 300.0, 101325.0);    // Component-specific enthalpy
real v_N2 = air.v(i, 300.0, 101325.0);    // Component-specific volume
```

### Access Underlying RefGas Objects

```cpp
// Get RefGas for component i
gastables::RefGas * N2_refgas = air.component(0);

// Use RefGas methods directly
real M = N2_refgas->M();  // Molar mass
real cp_poly = N2_refgas->cp(300.0);  // Specific heat (no pressure dependence)
```

### Access GasData

```cpp
gastables::GasData * data = air.data(0);  // GasData for component 0

real T_crit = data->T_crit();
real p_crit = data->p_crit();
real omega = data->acentric();
```

### Access Equation of State Directly

```cpp
Gas nitrogen("N2", GasModel::PR);

// Get EoS object
gasmodels::EoS * eos = nitrogen.eos();

// Use EoS methods
real p = eos->p(300.0, 0.001);  // Pressure from T, v
real v = eos->v(300.0, 101325.0);  // Volume from T, p

// EoS-specific departure functions
real h_dep = eos->hdep(300.0, 101325.0);  // Enthalpy departure
real s_dep = eos->sdep(300.0, 101325.0);  // Entropy departure
```

### State Variables Container

```cpp
// Access state variables object
gasmodels::Statevals & state = air.statevals();

// State variables are automatically updated during property calls
real M_current = state.get(BELFEM_STATEVAL_M);  // Current mixture molar mass
real R_current = state.get(BELFEM_STATEVAL_R);  // Current mixture gas constant
```

---

## Memory Management

**Ownership model:**
- `Gas` object owns its internal `RefGas*` components
- `Gas` destructor deletes all `RefGas*` objects
- `Gas` owns its `EoS*` object

**Best practice:**
```cpp
// Stack allocation (automatic cleanup)
{
    Gas air;   // default ctor builds the 14-species air mixture as an ideal gas
    real cp = air.cp(300.0, 101325.0);
}  // air destructor called automatically

// Heap allocation (manual cleanup)
Gas * air = new Gas();   // default ctor: the 14-species air mixture
real cp = air->cp(300.0, 101325.0);
delete air;  // Must delete manually
```

**For long-lived objects:**
```cpp
// Smart pointer (automatic cleanup)
std::unique_ptr<Gas> air( new Gas() );
real cp = air->cp(300.0, 101325.0);
// Automatic cleanup when unique_ptr goes out of scope
```

---

## Transport Properties of the Helmholtz Fluids

Viscosity and thermal conductivity for the cryogenic fluids come from published
correlations. Each fluid family has its own source. These correlations are
separate from the equation of state and were fitted independently.

| fluid | viscosity | thermal conductivity |
|---|---|---|
| nitrogen | Lemmon & Jacobsen 2004 | Lemmon & Jacobsen 2004 |
| oxygen | Lemmon & Jacobsen 2004 | Lemmon & Jacobsen 2004 |
| methane | Friend et al. 1989 | Friend et al. 1989 |
| parahydrogen | Muzny et al. 2013 + 2022 erratum | Assael et al. 2011, para tables |
| normal hydrogen | Muzny et al. 2013 + 2022 erratum | Assael et al. 2011, normal tables |
| orthohydrogen | Muzny et al. 2013 + 2022 erratum | Assael et al. 2011, normal tables |

Every fluid now answers both `mu()` and `lambda()`. A fluid with no correlation
would raise `BELFEM_ERROR` from the `HelmholtzTransport` base. It would not
return a wrong number. No current fluid takes that path.

### Where a correlation does not match the isomer

Hydrogen needs extra explanation because the code cannot express this caveat by
itself.

The **viscosity** correlation of Muzny et al. is for *normal* hydrogen only. No
published correlation exists for the spin isomers, so BELFEM installs it for all
three. REFPROP does the same.

The isomers can still return different values. The correlation uses the density
from the gas's equation of state, and those equations of state are isomer
specific. The correlation itself does not distinguish the isomers.

The **thermal conductivity** correlation of Assael et al. is different. It has
separate coefficient tables for normal and parahydrogen, and BELFEM uses both.
Orthohydrogen has no table of its own, so BELFEM uses the normal table for it.

So `lambda()` is isomer aware where `mu()` is not.

### The near-critical caveat for hydrogen

Hydrogen thermal conductivity has a critical enhancement term. Its magnitude is
inversely proportional to the background viscosity.

Assael et al. fitted that term in 2011 against the viscosity correlation that
REFPROP carried at the time. That correlation predates Muzny et al. 2013.
BELFEM uses the newer viscosity, so the enhancement is evaluated with a
background it was not fitted against.

Away from the critical point, this makes no difference. The paper's own
verification points reproduce to better than 0.01 %. Close to the critical point,
the enhancement becomes a real fraction of the total, and the mismatch shows:

| state | contribution of the enhancement | resulting error in `lambda()` |
|---|---|---|
| 35 K, 30 kg/m^3, normal | 20 % of the total | +1.9 % |
| 35 K, 30 kg/m^3, para | 16 % of the total | +5.9 % |

Parahydrogen is the worse of the two for the same reason as above. Its background
viscosity is the normal-hydrogen correlation, so the enhancement inherits that
error on top.

Within a few kelvin of the critical point, treat hydrogen thermal conductivity as
good to a few percent, not to the correlation's nominal uncertainty.

---

## Const Correctness and Thread Safety

Property accessors are `const`. Calling one evaluates a property. It does not
change what the gas is:

```cpp
void report( const Gas & aGas )        // read-only handle is enough
{
    std::cout << aGas.cp( 300.0, 101325.0 ) << std::endl ;
    std::cout << aGas.mu( 300.0, 101325.0 ) << std::endl ;
}
```

This rule applies across the module: `Gas`, `EoS`, `EoS_Idgas`, `EoS_Cubic`,
`AlphaFunction`, `Helmholtz`, the species equations of state, and the transport
models:

| Group | `const` | Examples |
|---|---|---|
| Evaluates a property of a fixed mixture | yes | `cp`, `h`, `s`, `mu`, `lambda`, `alpha`, `shock`, `prandtl_meyer`, `Gibbs` |
| Changes what the gas is | no | `remix`, `remix_mass`, `reset_mixture`, `remix_to_equilibrium` |
| Builds or rebuilds internal tables | no | the spline builders, `init_departure_splines`, `set_reference_point` |

The signature tells you which kind of call it is. A routine that only reads
properties can take `const Gas &`. That keeps accidental remixing out of that
routine.

### What `mutable` covers

The evaluators are *logically* const, not bitwise const. They may update internal
storage, but they do not change the fluid. That storage has two uses:

- **memoization caches**: `Statevals` in `Gas`, `mHelmholtzVals`/`mHelmholtzBits`
  in `Helmholtz`, `mCubicStatevals` in `EoS_Cubic`. A repeated call at the same
  state is served from the cache instead of recomputed.
- **preallocated scratch**: the mixture work matrices, the Cardano work vectors,
  and the `delta^d` and `tau^t` tables of the Helmholtz residual terms. These
  exist so the hot paths allocate nothing (see `CLAUDE.md`, "No temporary
  `Vector`/`Matrix` in frequently-called member functions").

Those members are not part of the identity of the fluid. They are `mutable` for
that reason. A second call with the same arguments still returns the same number.

Scratch for a **non**-const path is deliberately *not* `mutable`. This includes
the RAND equilibrium block, `mFormationTable`, and the viscosity/conductivity
spline build vectors. Leaving them plain keeps the compiler proving that the
composition-changing code stays out of the const paths. Prefer adding a
`mutable` only when the compiler demands it.

### Paired accessors and returning by value

Accessors that expose an internal handle are paired. `Gas::data( index )`
returns `GasData *` for a writable gas. For a const gas, it returns
`const GasData *`.

`RefGas` follows the same pattern for `data()` and for its three `*_spline()`
getters. The caller does not choose the overload directly. C++ picks it from the
constness of the handle.

Pointers are stored writable. `mComponents` is a `Cell< RefGas * >`, not a
`Cell< const RefGas * >`. This keeps the build paths simple.

When a const method only reads a component, it binds a read-only handle at the
point of use:

```cpp
real
Gas::cp( const uint aIndex, const real T, const real p ) const
{
    // read-only handle: the component is only evaluated here
    const gastables::RefGas * tComponent = mComponents( aIndex );

    return tComponent->heat_spline()->deval( T )
        / tComponent->data()->M()
        + mEoS->cpdep( aIndex, T, p );
}
```

Inside this function, `tComponent` is read-only. The compiler enforces that local
rule, so `heat_spline()` resolves to the const overload. The guarantee is narrow:
it holds inside this function, not across the class.

Scalar accessors return by value. `M( T, p )` and `R( T, p )` used to return
references into the state cache. That made each result a live view. A later remix
could change what the reference observed.

That behavior was useful if you expected it. It was a trap if you did not. A
`real` is a register return, so the copy is free. The accessor now returns a
value, not a live view.

### Not a thread-safety guarantee

> **Warning:** `const` here means logically const. It does **not** mean
> reentrant. Because the const evaluators write the shared cache, two threads
> must not call them on the same `Gas` object, even through a `const Gas &`.

This follows BELFEM's framework-wide policy. It is not a gap in this module.
BELFEM parallelizes with MPI and is deliberately not internally thread safe. See
`doc/coding_philosophy.md`.

Readers arriving from modern C++, where `const` usually implies a safe concurrent
read, should note the difference. Give each thread its own `Gas` object, or
serialize the calls.

---

## Common Pitfalls

### 1. Confusing Molar and Mass Fractions

```cpp
// BAD - using mass fractions as molar fractions
Vector<real> mass_fractions = {0.233, 0.767};  // Mass basis
Gas bad(species, mass_fractions, GasModel::IDGAS);  // Wrong!

// GOOD - convert first or use remix_mass
Gas good(species, {1.0, 0.0}, GasModel::IDGAS);  // Dummy initialization
good.remix_mass(mass_fractions);  // Correct mass-based remix
```

### 2. Real Gas Without Critical Data

```cpp
// Some gases lack critical point data in database
Gas rare("RareGas", GasModel::PR);  // May fail!

// Solution: Check if critical data exists
gastables::GasData * data = rare.data(0);
if (!data->has_crit()) {
    // Use IDGAS instead
    Gas rare_idgas("RareGas", GasModel::IDGAS);
}
```

### 3. Forgetting Pressure Dependence

```cpp
// BAD - comparing RefGas with Gas properties
gastables::RefGas * ref_N2 = /* from factory */;
Gas gas_N2("N2", GasModel::IDGAS);

real cp_ref = ref_N2->cp(300.0);           // No pressure argument
real cp_gas = gas_N2.cp(300.0, 101325.0);  // Requires pressure

// For IDGAS, results are identical (no pressure effect)
// For real gas, cp_gas includes departure function
```

### 4. Equilibrium on Non-Reactive Mixture

```cpp
Cell<string> inert = {"N2", "Ar"};  // No reactions possible
Gas mix(inert, {0.5, 0.5}, GasModel::IDGAS);

// Equilibrium calculation will not change composition
mix.remix_to_equilibrium(2000.0, 101325.0);
// Composition unchanged (no reaction pathways)
```

---

## Further Reading

### Module Documentation
- **`gasmodels_usage_guide.md`** - Comprehensive usage guide with examples

### Related Modules
- **`src/physics/gastables/doc/README.md`** - Reference gas data (foundation for gasmodels)
- **`src/numerics/spline/doc/spline_usage_guide.md`** - Spline interpolation used in property evaluation

### External References
**Equations of State:**
- Soave (1972), "Equilibrium Constants from a Modified Redlich-Kwong Equation of State", *Chemical Engineering Science*, 27(6):1197-1203
- Peng & Robinson (1976), "A New Two-Constant Equation of State", *Industrial & Engineering Chemistry Fundamentals*, 15(1):59-64

**Helmholtz Models:**
- Leachman et al. (2009), "Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen", *J. Phys. Chem. Ref. Data*, 38(3):721-748
- Setzmann & Wagner (1991), "A New Equation of State and Tables of Thermodynamic Properties for Methane", *J. Phys. Chem. Ref. Data*, 20(6):1061-1155
- Schmidt & Wagner (1985), "A New Form of the Equation of State for Pure Substances", *Fluid Phase Equilibria*, 19(3):175-200

**Mixture Rules:**
- Gordon & McBride (1994), "Computer Program for Calculation of Complex Chemical Equilibrium Compositions", NASA RP-1311, Eqs. (5.3)-(5.7)
- Wilke (1950), "A Viscosity Equation for Gas Mixtures", *J. Chem. Phys.*, 18:517-519

---

**Last Updated:** 2026-08-25
**Maintainer:** BELFEM development team
