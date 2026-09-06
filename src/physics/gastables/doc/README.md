# Gas Tables Module Documentation {#physics_gastables_index}

**Date:** 2026-01-30
**Purpose:** Quick reference and navigation guide for the gastables module
**Module:** `src/physics/gastables`

---

## Overview

The **gastables** module provides temperature-dependent thermodynamic and transport property data for pure gases using **NASA CEA polynomial format**. It serves as the foundation for the `gasmodels` module, which builds complete fluid models for gas mixtures and real gas equations of state.

**Key features:**
- NASA CEA 9-coefficient polynomial format for thermodynamic properties
- 4-coefficient polynomial format for transport properties (viscosity, thermal conductivity)
- Factory-synthesized extrapolation below and above the tabulated intervals (no supplemental files)
- Spline-based evaluation (faster alternative to polynomial evaluation)
- Gas composition tracking (elemental composition)
- Critical point data for real gas models
- Factory pattern for creating reference gases from data files

**Data sources:**
- NASA CEA (Chemical Equilibrium with Applications) thermodynamic database
- Supplemental critical point and acentric factor data for cubic equations of state

---

## Quick Reference

### Key Classes

| Class | Purpose | Location |
|-------|---------|----------|
| `RefGas` | Reference gas with temperature-dependent properties | `cl_GT_RefGas.hpp:74` |
| `RefGasFactory` | Factory for creating reference gases from data files | `cl_GT_RefGasFactory.hpp:58` |
| `GasData` | Container for gas metadata (M, R, critical properties, composition) | `cl_GT_GasData.hpp:57` |
| `HeatPoly` | Thermodynamic property polynomial (Cp, H, S) | `cl_GT_HeatPoly.hpp:46` |
| `TransportPoly` | Transport property polynomial (μ, λ) | `cl_GT_TransportPoly.hpp:61` |

### Key Files

| File | Purpose |
|------|---------|
| `cl_GT_RefGas.{hpp,cpp}` | Main reference gas class |
| `cl_GT_RefGasFactory.{hpp,cpp}` | Factory for creating reference gases |
| `cl_GT_GasData.{hpp,cpp}` | Gas metadata container |
| `cl_GT_HeatPoly.{hpp,cpp}` | Thermodynamic polynomial base class |
| `cl_GT_TransportPoly.{hpp,cpp}` | Transport polynomial base class |
| `cl_GT_InputThermo.{hpp,cpp}` | Parser for NASA CEA thermo.inp |
| `cl_GT_InputTransport.{hpp,cpp}` | Parser for NASA CEA trans.inp |
| `cl_GT_InputData.{hpp,cpp}` | Parser for critical point and composition data |
| `cl_GT_InputAlpha.{hpp,cpp}` | Parser for alpha function coefficients (cubic EoS) |

---

## Common Operations

### Create a Reference Gas

```cpp
#include "cl_GT_RefGasFactory.hpp"

using namespace belfem;
using namespace gastables;

// Create factory (loads data files from default path)
RefGasFactory factory;

// Create reference gas
RefGas * nitrogen = factory.create_refgas("N2");

// Query properties
real M = nitrogen->M();                      // Molar mass in kg/mol
real cp = nitrogen->cp(300.0);               // Specific heat at 300 K (J/(kg·K))
real h = nitrogen->h(300.0);                 // Specific enthalpy (J/kg)
real s = nitrogen->s(300.0);                 // Specific entropy (J/(kg·K))
real mu = nitrogen->mu(300.0);               // Dynamic viscosity (Pa·s)
real lambda = nitrogen->lambda(300.0);       // Thermal conductivity (W/(m·K))

// Clean up
delete nitrogen;
```

### Switch Between Polynomial and Spline Mode

```cpp
// Default: the factory leaves every gas in SPLINE mode
RefGas * gas = factory.create_refgas("O2");
real cp_spline = gas->cp(500.0);

// set_mode() only rebinds function pointers - the splines were built once
// by the factory; switch to POLY for the exact piecewise polynomials
gas->set_mode(RefGasMode::POLY);
real cp_poly = gas->cp(500.0);  // Same result at interpolation accuracy, slower (interval search)
```

### Access Gas Data

```cpp
RefGas * air_component = factory.create_refgas("N2");
const GasData * data = air_component->data();

// Query metadata
real M = data->M();                   // Molar mass (kg/mol)
real R = data->R();                   // Specific gas constant (J/(kg·K))
real T_crit = data->T_crit();         // Critical temperature (K)
real p_crit = data->p_crit();         // Critical pressure (Pa)
real omega = data->acentric();        // Acentric factor (dimensionless)

// Query elemental composition
const Map<string, real> & composition = data->composition();
// For N2: {"N": 2.0}

// Check data availability
bool has_crit = data->has_crit();           // Critical point data available?
bool has_cubic = data->has_cubic();         // Cubic EoS coefficients available?
```

---

## Property Evaluation API

### Thermodynamic Properties (Molar Basis)

**Capital letter functions return molar properties (J/mol):**

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `Cp(T)` | Molar heat capacity at constant pressure | J/(mol·K) | From polynomial |
| `H(T)` | Molar enthalpy | J/mol | Includes the formation enthalpy: H(0 K) = ΔHf(298.15 K), see the class header |
| `S(T)` | Molar entropy | J/(mol·K) | Absolute, referenced to 0 K |
| `dCpdT(T)` | Temperature derivative of Cp | J/(mol·K²) | For sensitivity analysis |
| `d2CpdT2(T)` | Second temperature derivative of Cp | J/(mol·K³) | For Newton solvers |
| `dSdT(T)` | Temperature derivative of S | J/(mol·K²) | For equilibrium |

### Thermodynamic Properties (Specific Basis)

**Lowercase letter functions return specific properties (J/kg):**

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `cp(T)` | Specific heat capacity at constant pressure | J/(kg·K) | `Cp(T) / M` |
| `h(T)` | Specific enthalpy | J/kg | `H(T) / M` |
| `s(T)` | Specific entropy | J/(kg·K) | `S(T) / M` |
| `dcpdT(T)` | Temperature derivative of cp | J/(kg·K²) | `dCpdT(T) / M` |
| `d2cpdT2(T)` | Second temperature derivative of cp | J/(kg·K³) | `d2CpdT2(T) / M` |

### Transport Properties

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `mu(T)` | Dynamic viscosity | Pa·s | From polynomial or spline |
| `lambda(T)` | Thermal conductivity | W/(m·K) | From polynomial or spline |
| `dmudT(T)` | Temperature derivative of μ | Pa·s/K | For sensitivity analysis |
| `dlambdadT(T)` | Temperature derivative of λ | W/(m·K²) | For sensitivity analysis |
| `d2mudT2(T)` | Second temperature derivative of μ | Pa·s/K² | For Newton solvers |
| `d2lambdadT2(T)` | Second temperature derivative of λ | W/(m·K³) | For Newton solvers |

### Reference State Properties

| Function | Returns | Units | Notes |
|----------|---------|-------|-------|
| `M()` | Molar mass | kg/mol | From GasData |
| `reference_formation_enthalpy()` | Formation enthalpy at 298.15 K | J/mol | For chemical equilibrium |
| `h_ref()` | Specific enthalpy at 298.15 K | J/kg | Reference state |
| `H_ref()` | Molar enthalpy at 298.15 K | J/mol | Reference state |

---

## Data Availability Checks

Before using properties, check if data is available:

```cpp
RefGas * gas = factory.create_refgas("He");

// Check what data is available
if (gas->has_thermo()) {
    real cp = gas->cp(300.0);  // Safe to use
}

if (gas->has_viscosity()) {
    real mu = gas->mu(300.0);  // Safe to use
}

if (gas->has_conductivity()) {
    real lambda = gas->lambda(300.0);  // Safe to use
}

// Below the lowest tabulated interval the factory has synthesized an
// extrapolation, so there is no separate cryo flag -- just check the record.
if (gas->has_thermo()) {
    real cp_cryo = gas->cp(10.0);
}

if (gas->has_viscosity()) {
    real mu_cryo = gas->mu(10.0);
}

// Check if elemental composition is available
if (gas->has_components()) {
    real multiplicity = gas->component_multiplicity("He");  // Returns 1.0
}

// is_noble() only reports the species class ( He, Ne, Ar, Kr, Xe, Rn );
// no transport path in the tree branches on it
if (gas->is_noble()) {
    // monatomic, no rotational/vibrational modes
}
```

---

## Temperature Ranges

### Standard NASA CEA Data

**Typical ranges:**
- Low temperature: 200-300 K
- High temperature: 6000 K

**Most gases have multiple polynomials:**
- Low-T polynomial: 200-1000 K
- Mid-T polynomial: 1000-6000 K

### Outside the tabulated intervals

There is **no per-species cryogenic dataset**. Whatever a record's lowest interval is,
`RefGasFactory` synthesizes an extrapolation below it at construction — together with glue
polynomials across interval junctions, a hot extrapolation above the top interval, and, for
species that have critical data but no transport record, viscosity and conductivity from the
Lucas and Chung correlations
(`cl_GT_RefGas.hpp:55-61`; `create_cryo_poly_heat()` at `cl_GT_RefGas.cpp:243`).

So there is nothing to check for availability beyond whether the species has a record at all:

```cpp
if (gas->has_thermo()) {
    // evaluation below the table works, but the value is extrapolated
}
```

Nothing signals that you have left the tabulated range. Validate extrapolated values against
measurement before relying on them.

---

## Polynomial Formats

### NASA CEA 9-Coefficient Format (Thermodynamic)

For a temperature range [T_min, T_max], the molar properties are:

These are the **nine-coefficient** forms of NASA RP-1311, Eqs. (4.9)–(4.11) — seven polynomial
coefficients `a1…a7` plus the two integration constants `b1`, `b2`.

**Heat capacity:**
```
Cp/R = a1·T⁻² + a2·T⁻¹ + a3 + a4·T + a5·T² + a6·T³ + a7·T⁴
```

**Enthalpy:**
```
H/(R·T) = −a1·T⁻² + a2·T⁻¹·ln(T) + a3 + a4·T/2 + a5·T²/3 + a6·T³/4 + a7·T⁴/5 + b1/T
```

**Entropy:**
```
S/R = −a1·T⁻²/2 − a2·T⁻¹ + a3·ln(T) + a4·T + a5·T²/2 + a6·T³/3 + a7·T⁴/4 + b2
```

where:
- `a1-a7`: temperature-dependent coefficients
- `b1`: integration constant for enthalpy (`mEnthalpyConstant`)
- `b2`: integration constant for entropy (`mEntropyConstant`)
- `R`: universal gas constant (8.314462618 J/(mol·K))

Implemented in `HeatPoly::Cp` / `H` / `S` (`cl_GT_HeatPoly.cpp:36-78`); the shipped `thermo.inp`
records carry nine coefficients per interval.

### 4-Coefficient Format (Transport)

**Dynamic viscosity and thermal conductivity** — NASA RP-1311 Eq. (5.1), in the **natural**
logarithm:
```
ln(X) = A·ln(T) + B/T + C/T² + D
```

The coefficients in `trans.inp` are for this form (`TransportPoly::rawpoly`,
`cl_GT_TransportPoly.cpp:38-44`). Writing the same correlation in base 10 would require every
coefficient divided by `ln 10`; feeding the shipped coefficients into a `10^f` form gives an
error that varies with temperature, not a constant factor.

---

## Memory Management

**Ownership model:**
- `RefGasFactory` creates `RefGas*` with `new`
- **Caller owns** the returned pointer and must `delete` it
- `RefGas` internally manages `HeatPoly*` and `TransportPoly*` containers
- `RefGas` destructor cleans up internal polynomials

**Best practice:**
```cpp
RefGasFactory factory;

// Create reference gas
RefGas * gas = factory.create_refgas("Ar");

// Use gas...
real cp = gas->cp(300.0);

// Clean up when done
delete gas;
```

**For long-lived reference gases:** Consider `std::unique_ptr`:
```cpp
std::unique_ptr<RefGas> gas(factory.create_refgas("Kr"));
real cp = gas->cp(300.0);
// Automatic cleanup
```

---

## Common Patterns

### Creating Multiple Reference Gases

```cpp
RefGasFactory factory;

Cell<string> species = {"N2", "O2", "Ar", "CO2"};
Cell<RefGas*> gases(species.size());

for (index_t i = 0; i < species.size(); ++i) {
    gases(i) = factory.create_refgas(species(i));
}

// Use gases...

// Clean up
for (RefGas * gas : gases) {
    delete gas;
}
```

### Evaluating Properties at Multiple Temperatures

```cpp
RefGas * gas = factory.create_refgas("H2");

// For many temperature evaluations, splines are faster
gas->set_mode(RefGasMode::SPLINE);

Vector<real> temperatures(100);
Vector<real> cp_values(100);

for (index_t i = 0; i < temperatures.length(); ++i) {
    real T = 300.0 + i * 10.0;  // 300 K to 1290 K
    temperatures(i) = T;
    cp_values(i) = gas->cp(T);
}

delete gas;
```

### Accessing Spline Objects Directly

```cpp
RefGas * gas = factory.create_refgas("O2");
gas->set_mode(RefGasMode::SPLINE);

// Get spline objects for advanced use
Spline * heat_spline = gas->heat_spline();
Spline * viscosity_spline = gas->viscosity_spline();
Spline * conductivity_spline = gas->conductivity_spline();

// Use spline methods directly
real h      = heat_spline->eval(300.0);    // enthalpy, not cp
real cp     = heat_spline->deval(300.0);   // the heat spline samples H
real dcp_dT = heat_spline->ddeval(300.0);
```

---

## Integration with Gas Models

The gastables module is used by `gasmodels::Gas` for ideal gas mixtures:

```cpp
#include "cl_Gas.hpp"

// Gas class internally uses RefGas objects
Gas air;   // default ctor builds the 14-species air mixture as an ideal gas

// Properties are computed from RefGas polynomials
real cp = air.cp(300.0, 101325.0);  // Uses N2, O2, Ar RefGas objects
```

See `src/physics/gasmodels/doc/gasmodels_usage_guide.md` for details.

---

## Const Correctness

Property evaluation is `const`. Calls such as `RefGas::Cp( T )`, `h( T )`,
`mu( T )`, `lambda( T )`, and their derivatives do not modify the `RefGas`.
A read-only handle is enough:

```cpp
void report( const gastables::RefGas & aGas )
{
    std::cout << aGas.cp( 300.0 ) << std::endl ;
    std::cout << aGas.mu( 300.0 ) << std::endl ;
}
```

`HeatPoly`, `TransportPoly`, and `GasData` were already const-correct. `RefGas`
now follows the same rule.

This module also needs **no** `mutable`. Unlike `gasmodels`, `RefGas` does not
memoize property values. For each call, it finds the polynomial or spline
interval for the requested temperature and evaluates it immediately. A property
call has no class state to update.

The methods that stay non-const are the build-time methods:
`add_heat_poly`, `add_transport_poly`, `create_splines`, `create_*_poly_*`,
`finalize`, `finalize_thermo`, `finalize_transport`, `fix_switches`,
`fix_reference_points`, `set_mode`, the `set_*`/`unset_*` flags, and the
`delete_*` cleanups. The factory calls these methods. Normal consumers do not.

`data()` and the three `*_spline()` getters come in pairs. The writable form
hands out `GasData *` / `Spline *` and is what the factory and the mixture
builder use. The `const` form returns `const GasData *` / `const Spline *` and
is what a const `RefGas` resolves to. Overload resolution picks the right one
from the constness of the handle, so neither caller has to think about it.

---

## Common Pitfalls

### 1. Temperature Out of Range

```cpp
RefGas * gas = factory.create_refgas("N2");

// Standard range: ~200-6000 K
real cp_ok = gas->cp(300.0);     // OK

real cp_low = gas->cp(50.0);     // synthesized extrapolation below the table
real cp_high = gas->cp(10000.0); // synthesized extrapolation above it

// There is nothing to "enable" -- both extrapolations are built at
// construction. Check that the species has a record, and treat values
// outside the tabulated range as extrapolated.
if (gas->has_thermo()) {
    real cp_cryo = gas->cp(50.0);
}
```

### 2. Forgetting to Delete

```cpp
// BAD - memory leak
void compute_property() {
    RefGasFactory factory;
    RefGas * gas = factory.create_refgas("He");
    real cp = gas->cp(300.0);
    // Forgot to delete gas!
}

// GOOD - clean up
void compute_property() {
    RefGasFactory factory;
    RefGas * gas = factory.create_refgas("He");
    real cp = gas->cp(300.0);
    delete gas;  // Clean up
}
```

### 3. Mode Switching

```cpp
RefGas * gas = factory.create_refgas("Ar");

// The factory has already built the splines and left the gas in SPLINE mode;
// set_mode() only rebinds function pointers, so switching is cheap
gas->set_mode(RefGasMode::POLY);    // exact piecewise polynomials, interval search per call
gas->set_mode(RefGasMode::SPLINE);  // back to the splines, no rebuild
```

### 4. Missing Data

```cpp
RefGas * gas = factory.create_refgas("SomeRareGas");

// Check before using
if (!gas->has_viscosity()) {
    // Viscosity data not available - handle gracefully
    BELFEM_ERROR(false, "Gas %s has no viscosity data", gas->label().c_str());
}
```

---

## Development Notes

### Adding New Gas Data

1. Add thermodynamic data to `share/fluid/thermo.inp` (NASA CEA 9-coefficient format)
2. Add transport data to `share/fluid/trans.inp` (NASA CEA, natural-log correlation)
3. Add the critical-point row to `share/fluid/gasdata.inp`
4. If the species needs a cubic EOS, add its alpha coefficients to `share/fluid/cubicalpha.inp`
5. Test with the `gas` tool (`src/physics/gasmodels/main.cpp`), or with `gastable` if the tree was configured with `-DUSE_EXAMPLES=ON`

### Data File Locations

The directory is resolved at run time, not fixed at build time
(`fn_GT_data_path.cpp`). In order:

1. `$BELFEM_DATA` + `/fluid`, if that global is set (`gBelfemDataPath`, populated from the
   environment by `Communicator::set_globals()`, and settable by a caller);
2. otherwise a walk up from the working directory — `share/fluid`, `../share/fluid`,
   `../../share/fluid`, …

The **fallback** candidates only count if they contain **`gasdata.inp`**, the marker file. An
explicitly set `$BELFEM_DATA` is returned **unchecked** (`fn_GT_data_path.cpp`) and a wrong value
surfaces later, as a failure to open a file.

Files — exactly these four:
- `thermo.inp`: NASA CEA nine-coefficient caloric polynomials
- `trans.inp`: NASA CEA transport coefficients, `ln(X) = A ln T + B/T + C/T² + D`
- `gasdata.inp`: per-species constants — molar mass, critical point, `Z_crit`, acentric factor,
  dipole moment. Also the marker the path resolver looks for
- `cubicalpha.inp`: alpha-function coefficients for the cubic equations of state

---

## Further Reading

### Module-Specific Documentation

- **`gastables_usage_guide.md`** - Comprehensive usage guide with examples

### Related Modules

- **`src/physics/gasmodels/doc/README.md`** - Gas models using gastables
- **`src/numerics/spline/doc/spline_usage_guide.md`** - Spline interpolation

### External References

**NASA CEA Documentation:**
- Gordon & McBride (1994), "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications"
- NASA Reference Publication 1311

**Thermodynamic Property Methods:**
- NIST Chemistry WebBook: https://webbook.nist.gov/chemistry/
- Poling, Prausnitz & O'Connell (2001), "The Properties of Gases and Liquids" (5th ed.), Appendix A

---

**Last Updated:** 2026-08-25
**Maintainer:** BELFEM development team
