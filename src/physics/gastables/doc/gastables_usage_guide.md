# Gas Tables Usage Guide {#physics_gastables_gastables_usage_guide}

**Date:** 2026-01-30
**Purpose:** Comprehensive guide for using the gastables module
**Module:** `src/physics/gastables`

---

## Table of Contents

1. [Introduction](#introduction)
2. [Basic Usage](#basic-usage)
3. [Property Evaluation](#property-evaluation)
4. [Performance Considerations](#performance-considerations)
5. [Advanced Topics](#advanced-topics)
6. [Integration with Gas Models](#integration-with-gas-models)
7. [Data File Format](#data-file-format)
8. [The `gastable` Command-Line Tool](#the-gastable-command-line-tool)

---

## Introduction

The **gastables** module provides a C++ interface to **NASA CEA thermodynamic and transport property data** for pure gases. It is the foundation for the `gasmodels` module, which uses these reference gases to build complete fluid models with:

- Gas mixtures with arbitrary composition
- Real gas equations of state (cubic EoS, Helmholtz EoS)
- Chemical equilibrium calculations
- Compressible flow analysis

**Design philosophy:**

Following BELFEM's HPC-first approach (see `doc/coding_philosophy.md`):
- **Manual memory management**: Factory returns raw pointers, caller deletes
- **Performance-critical paths**: Spline mode for repeated evaluations
- **Zero overhead**: Polynomial evaluation compiles to tight loops
- **Cache-friendly**: Polynomial coefficients stored contiguously

**Build switch:** the module is compiled only when the tree is configured with
`-DUSE_GASMODELS=ON` (`option( USE_GASMODELS … OFF )` in the top-level `CMakeLists.txt`;
`src/physics/CMakeLists.txt` adds `gastables` and `gasmodels` only under this option). In a
default tree, none of the classes below exist.

---

## Basic Usage

### Creating a Reference Gas

The `RefGasFactory` loads data files and creates `RefGas` objects:

```cpp
#include "cl_GT_RefGasFactory.hpp"

using namespace belfem;
using namespace gastables;

int main()
{
    // Create factory (loads data from default path)
    RefGasFactory factory;

    // Create a reference gas
    RefGas * nitrogen = factory.create_refgas("N2");

    // Evaluate properties
    real T = 300.0;  // Temperature in Kelvin
    real cp = nitrogen->cp(T);           // Specific heat (J/(kg·K))
    real h = nitrogen->h(T);             // Specific enthalpy (J/kg)
    real s = nitrogen->s(T);             // Specific entropy (J/(kg·K))
    real mu = nitrogen->mu(T);           // Dynamic viscosity (Pa·s)
    real lambda = nitrogen->lambda(T);   // Thermal conductivity (W/(m·K))

    // Access metadata
    real M = nitrogen->M();  // Molar mass (kg/mol)

    // Clean up
    delete nitrogen;

    return 0;
}
```

**Key points:**
- Factory loads data files once (amortized cost)
- `create_refgas()` returns raw pointer - **caller must delete**
- Property evaluation is thread-safe (const methods on read-only data)
- Factory can be reused to create multiple gases

### Custom Data Path

```cpp
// Use custom data directory
string custom_path = "/path/to/gastables/data/";
RefGasFactory factory(custom_path);

RefGas * gas = factory.create_refgas("O2");
```

The path is resolved at run time: `$BELFEM_DATA/fluid` if that global is set, otherwise a walk up from the working directory (`share/fluid`, `../share/fluid`, …), accepting the first candidate that contains the marker file `gasdata.inp` (`fn_GT_data_path.cpp`).

---

## Property Evaluation

### Thermodynamic Properties

**Molar properties** (capital letters) - units per mole:

```cpp
RefGas * gas = factory.create_refgas("Ar");

real T = 500.0;  // Kelvin

// Molar properties
real Cp_molar = gas->Cp(T);         // J/(mol·K)
real H_molar = gas->H(T);           // J/mol, on the factory's enthalpy scale (below)
real S_molar = gas->S(T);           // J/(mol·K), standard state p = 1 bar

// Derivatives (for sensitivity analysis, Newton solvers)
real dCp_dT = gas->dCpdT(T);        // J/(mol·K²)
real d2Cp_dT2 = gas->d2CpdT2(T);    // J/(mol·K³)
real dS_dT = gas->dSdT(T);          // J/(mol·K²)
```

**Specific properties** (lowercase letters) - units per kg:

```cpp
// Specific properties (automatically divided by molar mass)
real cp_specific = gas->cp(T);       // J/(kg·K)
real h_specific = gas->h(T);         // J/kg
real s_specific = gas->s(T);         // J/(kg·K)

// Derivatives
real dcp_dT = gas->dcpdT(T);         // J/(kg·K²)
real d2cp_dT2 = gas->d2cpdT2(T);     // J/(kg·K³)

// Relationship: cp = Cp / M
BELFEM_ASSERT(std::abs(cp_specific - Cp_molar/gas->M()) < 1e-10, "Inconsistent!");
```

`S( T )` and `s( T )` are the ideal gas's absolute entropy at the standard pressure
`gastables::gPref` = 1 bar (`GT_globals.hpp`, `BELFEM_PREF`). `RefGas` has no pressure argument;
the `gasmodels` layer (`Gas`) adds the pressure term `ln( gPref / p )`.

**Enthalpy scale.** `RefGasFactory` discards the record's integration constant and shifts each
species so that `H( 298.15 K ) = ΔHf( 298.15 K ) + [ H( 298.15 K ) − H( 0 K ) ]`. Both terms
come from `thermo.inp` (the `RefGas` class documentation in `cl_GT_RefGas.hpp` states this
contract). Each species therefore has its own shifted scale. Differences within one species,
such as `H( T2 ) − H( T1 )`, are unaffected.

**Reference state properties:**

```cpp
// Properties at 298.15 K
real H_ref = gas->H_ref();                         // Molar enthalpy at 298.15 K (J/mol), shifted scale
real h_ref = gas->h_ref();                         // Specific enthalpy at 298.15 K (J/kg), shifted scale
real H_formation = gas->reference_formation_enthalpy();  // Formation enthalpy at 298.15 K (J/mol)

delete gas;
```

### Transport Properties

```cpp
RefGas * gas = factory.create_refgas("He");

real T = 300.0;

// Viscosity
real mu = gas->mu(T);                // Dynamic viscosity (Pa·s)
real dmu_dT = gas->dmudT(T);         // ∂μ/∂T (Pa·s/K)
real d2mu_dT2 = gas->d2mudT2(T);     // ∂²μ/∂T² (Pa·s/K²)

// Thermal conductivity
real lambda = gas->lambda(T);        // W/(m·K)
real dlambda_dT = gas->dlambdadT(T); // ∂λ/∂T (W/(m·K²))
real d2lambda_dT2 = gas->d2lambdadT2(T);  // ∂²λ/∂T² (W/(m·K³))

delete gas;
```

A species with no `trans.inp` record gets synthesized transport if its critical data are complete.
Otherwise these accessors return zero rather than raising (see Handling Missing Data). Check
`has_viscosity()` and `has_conductivity()` separately; the two flags are independent.

**Common error:**
```cpp
// Don't confuse with kinematic viscosity
real mu = gas->mu(T);       // Dynamic viscosity (Pa·s)
real rho = /* compute density */;
real nu = mu / rho;         // Kinematic viscosity (m²/s) - you compute this
```

### Gas Metadata

```cpp
RefGas * gas = factory.create_refgas("CO2");
const GasData * data = gas->data();

// Basic properties
string label = data->label();        // Chemical symbol ("CO2")
string name = data->name();          // Full name ("carbon dioxide")
string cas = data->cas();            // CAS number ("124-38-9")
real M = data->M();                  // Molar mass (kg/mol)
real R = data->R();                  // Specific gas constant (J/(kg·K)) = Ru/M

// Critical point (if available)
if (data->has_crit()) {
    real T_crit = data->T_crit();    // Critical temperature (K)
    real p_crit = data->p_crit();    // Critical pressure (Pa)
    real rho_crit = data->rho_crit();// Critical density (kg/m³)
    real Z_crit = data->Z_crit();    // Critical compressibility factor
    real omega = data->acentric();   // Acentric factor (for cubic EoS)
    real dipole = data->dipole();    // Dipole moment (Debye)
}

// Cubic EoS coefficients (if available)
if (data->has_cubic()) {
    const Vector<real> & srk_coeffs = data->srk();  // SRK alpha function coeffs
    const Vector<real> & pr_coeffs = data->pr();    // PR alpha function coeffs
}

// Elemental composition
const Map<string, real> & composition = data->composition();
// For CO2: {"C": 1.0, "O": 2.0}

const Cell<string> & elements = data->elements();
// For CO2: {"C", "O"}

// Query specific element
real oxygen_count = data->component_multiplicity("O");  // Returns 2.0 for CO2
real nitrogen_count = data->component_multiplicity("N"); // Returns 0.0 (not present)

delete gas;
```

`RefGas` also reports `has_components()` (whether an elemental composition was read) and
`is_liquid()` (whether the CEA record describes a condensed phase), in addition to
`has_thermo()`, `has_viscosity()`, `has_conductivity()` and `is_noble()` (`cl_GT_RefGas.hpp`).

---

## Performance Considerations

### Polynomial vs Spline Mode

**Spline mode (default: the factory leaves every gas in this mode)**
- Cubic splines are sampled from the NASA CEA polynomials. `RefGasFactory::create_refgas()`
  builds them once on a fixed grid from 0 K to `gTmax` = 13000 K in 5 K steps (`GT_globals.hpp`)
- There is no interval search, so this mode is expected to be faster than polynomial mode. No
  measured speedup is documented
- Interpolation accuracy
- `d2CpdT2()` / `d2cpdT2()` still evaluate the polynomials in this mode (`RefGas::set_mode`)

**Polynomial mode**
- Evaluates the polynomials directly, with a linear interval search on every call. Inside the
  tabulated intervals these are the NASA CEA 9-coefficient polynomials; outside them they are the
  glue and extrapolation polynomials the factory synthesized
- `set_mode()` only rebinds pointers in either direction; nothing is rebuilt
- **Synthesized transport always reads the spline.** If the factory synthesized a species' viscosity
  or conductivity from critical-point data because `trans.inp` has no record for it, those values
  exist only in the splines; the polynomial list contains an empty placeholder. `set_mode( POLY )`
  therefore binds the matching accessor (`mu()` or `lambda()`) and its derivatives to the spline, so
  both modes return the same synthesized values (`RefGas::set_mode`, `mViscositySynthesized`,
  `mConductivitySynthesized`).

To compare the two modes on your machine (`Timer::stop()` returns milliseconds):

```cpp
RefGas * gas = factory.create_refgas("N2");   // arrives in SPLINE mode
real tSum = 0.0;                              // keeps the loops from being optimized away

// Time polynomial mode
gas->set_mode(RefGasMode::POLY);
Timer tTimerPoly;
for (int i = 0; i < 100000; ++i) {
    tSum += gas->cp(300.0 + i * 0.01);
}
uint64_t tMillisecondsPoly = tTimerPoly.stop();

// Switch back to spline mode (pointer rebind, nothing is rebuilt)
gas->set_mode(RefGasMode::SPLINE);

// Time spline mode
Timer tTimerSpline;
for (int i = 0; i < 100000; ++i) {
    tSum += gas->cp(300.0 + i * 0.01);
}
uint64_t tMillisecondsSpline = tTimerSpline.stop();

delete gas;
```

**When to use splines:**
- Evaluating properties at >100 different temperatures
- Inside time-stepping loops
- During iterative solvers

**When to use polynomials:**
- When the polynomial values themselves are wanted instead of their spline interpolation. The
  splines exist either way, so a small number of evaluations is no reason to switch

### Memory Management Best Practices

Following BELFEM's manual memory philosophy:

```cpp
void bad_example()
{
    RefGasFactory factory;

    // BAD - creates factory repeatedly
    for (int i = 0; i < 1000; ++i) {
        RefGasFactory factory_inner;  // Reloads data files every iteration!
        RefGas * gas = factory_inner.create_refgas("Ar");
        real cp = gas->cp(300.0);
        delete gas;
    }
}

void good_example()
{
    // GOOD - reuse factory
    RefGasFactory factory;

    for (int i = 0; i < 1000; ++i) {
        RefGas * gas = factory.create_refgas("Ar");  // Fast - data already loaded
        real cp = gas->cp(300.0);
        delete gas;
    }
}

void better_example()
{
    // BETTER - reuse RefGas object
    RefGasFactory factory;
    RefGas * gas = factory.create_refgas("Ar");  // Create once

    for (int i = 0; i < 1000; ++i) {
        real cp = gas->cp(300.0 + i * 0.1);  // Just evaluate
    }

    delete gas;  // Clean up once
}
```

### Cache-Friendly Access Patterns

For evaluating multiple gases at same temperature:

```cpp
// tGases holds RefGas pointers created by the factory (see Pattern 3)
index_t tNumSpecies = tGases.size();
Cell<real> tCp(tNumSpecies, 0.0);    // two-argument constructor creates the elements

// GOOD - same temperature for every species
real T = 500.0;
for (index_t i = 0; i < tNumSpecies; ++i) {
    tCp(i) = tGases(i)->cp(T);
}

// LESS GOOD - random access if T varies unpredictably
for (index_t i = 0; i < tNumSpecies; ++i) {
    real tTrandom = get_random_temperature();
    tCp(i) = tGases(i)->cp(tTrandom);
}
```

---

## Advanced Topics

### Accessing Internal Polynomials

For debugging or custom property evaluations:

The polynomial lookup `find_heat_poly( T )` is `protected` (it returns `const HeatPoly *`) and
cannot be called from user code. The splines are public:

```cpp
RefGas * gas = factory.create_refgas("H2");
real T = 1500.0;

const Spline * heat_spline = gas->heat_spline();
const Spline * visc_spline = gas->viscosity_spline();
const Spline * cond_spline = gas->conductivity_spline();

// The heat spline samples MOLAR ENTHALPY, so the derivative ladder is shifted
// by one and every value is per mole: eval is H, deval is Cp, ddeval is dCp/dT.
real H      = heat_spline->eval(T);    // J/mol
real Cp     = heat_spline->deval(T);   // J/(mol·K)
real dCp_dT = heat_spline->ddeval(T);  // J/(mol·K²)

delete gas;
```

Prefer the public accessors `H( T )`, `Cp( T )` and `dCpdT( T )`, or their per-kg forms `h`, `cp`
and `dcpdT`. In SPLINE mode, they dispatch to the spline and already account for the derivative
ladder.

### Handling Missing Data

A label that is not in `thermo.inp` never produces a `RefGas`: `create_refgas()` stops the run
with an always-active `BELFEM_ERROR` ("Species … not found in thermo.inp ( labels are case
sensitive )", `RefGasFactory::create_refgas`). There is no thermo fallback to write.

Transport data may be missing. A species without a `trans.inp` record gets viscosity and
conductivity synthesized from the Lucas and Chung correlations if its critical data are complete
(`data()->has_crit()`); otherwise `mu()` and `lambda()` return 0. Check each property separately:

```cpp
RefGas * gas = factory.create_refgas("O2");

if (gas->has_viscosity()) {
    real mu = gas->mu(300.0);
}

if (gas->has_conductivity()) {
    real lambda = gas->lambda(300.0);
}

delete gas;
```

Synthesized transport is evaluated from its spline in both modes (see Polynomial vs Spline Mode).

### Cryogenic Temperature Handling

**There is no separate cryogenic dataset and no `has_cryo_*` predicate.** The shipped tables stop
at their lowest interval, and `RefGasFactory` *synthesizes* the range below it as an extrapolation
at construction — together with glue polynomials across interval junctions and, for species with
critical data but no transport record, viscosity and conductivity from the Lucas and Chung
correlations
(`cl_GT_RefGas.hpp:55-61`; `create_cryo_poly_heat()` at `cl_GT_RefGas.cpp:243`,
`create_cryo_poly_transport()` at `:290,306`).

So the only question to ask is whether the species has a record at all:

```cpp
RefGas * helium = factory.create_refgas("He");

if (helium->has_thermo()) {
    // the cryogenic extrapolation is already in place
    real cp_cryo = helium->cp(10.0);
    real h_cryo  = helium->h(4.2);
}

if (helium->has_viscosity()) {
    real mu_cryo = helium->mu(10.0);
}

delete helium;
```

Treat values far below the lowest tabulated interval as an extrapolation and sanity-check them
against measurement; the accessor will not warn you.

### Interaction Parameters for Mixtures

For gas mixtures, viscosity interaction parameters may be available:

```cpp
RefGasFactory factory;

// Check if interaction parameter exists
bool has_interaction = factory.interaction_viscosity_exists("N2", "O2");

if (has_interaction) {
    // Create interaction RefGas object
    RefGas * interaction = factory.create_interaction_viscosity("N2", "O2");

    // Use in mixture viscosity calculation (see gasmodels documentation)
    real mu_interaction = interaction->mu(300.0);

    delete interaction;
}

// Note: This is typically handled internally by Gas class
```

### Spline Discretization

The spline grid is fixed and cannot be changed through the public API. `RefGasFactory` builds
the grid once in its constructor: `create_temperature_steps()` fills
`linspace( 0.0, gTmax, gNumberOfSplinePoints )`, and `create_helpmatrix()` builds the matching
help matrix for a 5 K step (`gDeltaT`, `GT_globals.hpp`). Both functions are public, but neither
accepts a custom grid. `RefGas::create_splines()`, which consumes them, is `protected` and called
only by `create_refgas()`. The resolution is set at compile time by the constants `gTmax` and `gDeltaT`
in `GT_globals.hpp`.

---

## Integration with Gas Models

The gastables module is the foundation for the `gasmodels` module:

### Example: Gas Mixture

```cpp
#include "cl_Gas.hpp"

// Create a gas mixture (internally uses RefGas objects)
Cell<string> species = {"N2", "O2", "Ar"};
Vector<real> molar_fractions = {0.78, 0.21, 0.01};

Gas air(species, molar_fractions, GasModel::IDGAS);

// Properties are computed from weighted RefGas contributions
real T = 300.0;
real p = 101325.0;
real cp_mix = air.cp(T, p);  // Mixture rule applied to RefGas cp values
```

### Example: Ideal Gas vs Real Gas

```cpp
// Ideal gas - uses RefGas for thermodynamic properties only
Gas idgas_N2("N2", GasModel::IDGAS);
real cp_idgas = idgas_N2.cp(300.0, 101325.0);  // From RefGas polynomial
real v_idgas = idgas_N2.v(300.0, 101325.0);    // From ideal gas law: v = RT/p

// Real gas (cubic EoS) - uses RefGas data + cubic EoS corrections
Gas realgas_N2("N2", GasModel::PR);  // Peng-Robinson EoS
real cp_real = realgas_N2.cp(300.0, 101325.0); // RefGas cp + departure function
real v_real = realgas_N2.v(300.0, 101325.0);   // From cubic EoS (v ≠ RT/p)

// RefGas data used:
// - Molar mass M
// - Critical point (T_crit, p_crit)
// - Acentric factor ω
// - Thermodynamic polynomials for ideal gas reference
```

See `src/physics/gasmodels/doc/gasmodels_usage_guide.md` for complete details.

---

## Data File Format

### Thermodynamic Data (`thermo.inp`)

NASA CEA nine-coefficient format. Each species has a header, then one block per temperature
interval. The **exponent row** is part of the record and states the powers explicitly, which is
the quickest confirmation that this is the 9-coefficient form and not the 7-coefficient CHEMKIN
one:

```
<name>            <reference text>                                     @<n>
 <n_int> <date> <composition>  <phase>  <molecular weight>  <heat of formation>
    <T_lo>    <T_hi><n_coeff> -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0     <H(298)-H(0)>
 a1              a2              a3              a4              a5
 a6              a7                              b1              b2
    ... one T-range line + two coefficient lines per further interval ...
```

**Example — the first interval of N2, as shipped** (`share/fluid/thermo.inp:44-48`):
```
N2                Ref-Elm. Gurvich,1978 pt1 p280 pt2 p207.  @1
 4 tpis78 N   2.00    0.00    0.00    0.00    0.00 0   28.0134000          0.000
     63.651    250.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8670.104
 3.317622110D+02-1.698335042D+01 3.846858399D+00-3.609755503D-03 2.025052710D-05
-5.794621947D-08 6.632584479D-11                -9.842625425D+02 1.629503524D+00
```

Note the `D` exponent marker — Fortran double-precision notation, not `E`.

### Transport Data (`trans.inp`)

Four coefficients `A B C D` per interval, for the **natural**-log correlation
`ln(X) = A·ln(T) + B/T + C/T² + D`. Each line begins with the property selector — `V` for
viscosity, `C` for thermal conductivity — followed by the interval bounds:

```
<name>            <reference text>                                     @<n>
 <V|C>   <T_lo>   <T_hi>   A               B               C               D
```

**Example — N2 viscosity, as shipped** (`share/fluid/trans.inp:51-55`):
```
N2                                V4C5  BOUSHEHRI ET AL (1987)  SVEHLA (1994) @1
 V     63.7    250.0 6.18933491E-01-5.43168794E 01 9.78471406E 02 1.82870322E 00
 V    250.0 1000.0   0.62526577E 00-0.31779652E 02-0.16407983E 04 0.17454992E 01
 V 1000.0   5000.0   0.87395209E 00 0.56152222E 03-0.17394809E 06-0.39335958E 00
 V 5000.0  15000.0   0.88503551E+00 0.90902171E+03-0.73129061E+06-0.53503838E+00
```

Note `E 01` with a space where a sign would normally sit — the readers parse fixed columns
(`cl_GT_InputTransport.cpp`), so the spacing is significant.

### Critical Point Data (`gasdata.inp`)

Fixed-column table of per-species constants. **The file's units are not the stored units:** the
reader converts the molar-mass column from g/mol to kg/mol and the critical-pressure column from
**bar** to Pa (`InputData::read_data` in `cl_GT_InputData.cpp`).

```
<symbol>  <name>  <CAS>  <M, g/mol>  <T_crit, K>  <p_crit, bar>  <Z_crit>  <omega>  <dipole>  <sources...>  <symmetry flag>
```

**Example — helium, as shipped** (`share/fluid/gasdata.inp`):
```
He          helium            7440-59-7   4.003   5.195   2.2832  0.3040 -0.3835  0.0   CoolProp-8.0.0 ... symmetry
```

(Helium's negative acentric factor is correct, not a sign error.)

### Cubic-EOS Alpha Coefficients (`cubicalpha.inp`)

The fourth shipped file contains the per-species alpha-function coefficients used by the cubic
equations of state. Entries are keyed by symbol and CAS and include a critical temperature and
pressure:

```
<symbol>  <CAS>  <T_crit, K>  <p_crit, bar>  <alpha coefficients ...>
```

For a species listed here, the reader **overwrites** the `T_crit` and `p_crit` read from
`gasdata.inp` with these values (`InputAlpha::read_data` in `cl_GT_InputAlpha.cpp`), so the
critical point matches the alpha coefficients. `cubicalpha.inp` is read after `gasdata.inp`
(`RefGasFactory::create_refgas`).

There is **no** `crthermo.inp` or `crtrans.inp`. The factory opens exactly four files —
`thermo.inp`, `trans.inp`, `gasdata.inp` and `cubicalpha.inp`
(`cl_GT_RefGasFactory.cpp:39-49`) — and `share/fluid/` contains exactly those four.

---

## Common Usage Patterns

### Pattern 1: Evaluate Multiple Properties

```cpp
RefGas * gas = factory.create_refgas("O2");
real T = 400.0;

// Evaluate all thermodynamic properties
real cp = gas->cp(T);
real h = gas->h(T);
real s = gas->s(T);
real mu = gas->mu(T);
real lambda = gas->lambda(T);

// For efficiency, properties are computed independently
// (no caching between calls - each call recomputes)

delete gas;
```

### Pattern 2: Temperature Sweep

```cpp
RefGas * gas = factory.create_refgas("Ar");   // already in SPLINE mode

Cell<real> tTemperatures(100, 0.0);
Cell<real> tCp(100, 0.0);

for (index_t i = 0; i < tTemperatures.size(); ++i) {
    tTemperatures(i) = 300.0 + i * 10.0;  // 300 K to 1290 K
    tCp(i) = gas->cp(tTemperatures(i));
}

delete gas;
```

### Pattern 3: Multi-Species Evaluation

```cpp
RefGasFactory factory;

Cell<string> tSpecies = {"H2", "O2", "N2", "H2O", "CO2"};

// Cell(n) only reserves capacity; push() appends
Cell<RefGas*> tGases;
tGases.reserve(tSpecies.size());

// Create all reference gases
for (const string & tLabel : tSpecies) {
    tGases.push(factory.create_refgas(tLabel));
}

// Evaluate properties
real T = 500.0;
for (index_t i = 0; i < tGases.size(); ++i) {
    message(InfoLevel::Default,
            "%s: cp = %.3f J/(kg·K), mu = %.3e Pa·s",
            tSpecies(i).c_str(),
            tGases(i)->cp(T),
            tGases(i)->mu(T));
}

// Clean up: Cell<RefGas*> does not own its pointees
for (RefGas * tGas : tGases) {
    delete tGas;
}
```

### Pattern 4: Noble Gas Detection

```cpp
RefGas * gas = factory.create_refgas("He");

if (gas->is_noble()) {
    // Noble gases: He, Ne, Ar, Kr, Xe, Rn
    // (Monatomic, no rotational/vibrational modes)
    // is_noble() only reports the class; no transport path in the tree branches on it
}
```

### Pattern 5: Enthalpy Scale and Reaction Enthalpies

`H( T )` does **not** return a CEA absolute enthalpy. Each species has its own shifted scale:
`H( 298.15 K ) = ΔHf( 298.15 K ) + [ H( 298.15 K ) − H( 0 K ) ]` (see Thermodynamic Properties).
The sensible term `H( 298.15 K ) − H( 0 K )` differs by species (8670.104 J/mol for N2,
`share/fluid/thermo.inp:46`). A stoichiometric sum of `H( T )` over products minus reactants
therefore includes the stoichiometric sum of those terms. They do not cancel, so the result is
not a reaction enthalpy.

The `RefGas` class documentation (`cl_GT_RefGas.hpp`, Caloric properties) explains how the
in-tree combustion code returns to the CEA convention: it subtracts each species' enthalpy at
298.15 K and adds the formation enthalpy back. `H_ref()` and `reference_formation_enthalpy()`
provide these two values.

---

## Troubleshooting

### Issue: Gas Not Found

```
Species XYZ not found in thermo.inp ( labels are case sensitive )
```

Raised by `create_refgas()` before any data is read; the run stops.

**Solution:**
1. Check spelling (case-sensitive: "N2", not "n2")
2. Check if species exists in `thermo.inp`
3. Add custom data if needed

### Issue: Property Evaluation Fails

```
Temperature T=<T> out of bounds for gas <label>.
```

Raised when a property is evaluated from the polynomials at a `T` outside the species' polynomial
coverage. After the factory's extrapolations, that coverage runs from 0 K to the larger of
`gTmax` = 13000 K and the upper bound of the last tabulated interval. This applies separately to
the heat, viscosity and conductivity polynomials (`find_heat_poly`, `find_viscosity_poly`,
`find_conductivity_poly` in `cl_GT_RefGas.cpp`). For N2 the heat polynomials reach 20000 K
(`share/fluid/thermo.inp`) and the viscosity polynomials reach 15000 K (`share/fluid/trans.inp`).
Negative temperatures always raise.

This applies to every polynomial-backed accessor in POLY mode and also to `d2CpdT2()` /
`d2cpdT2()` in SPLINE mode, because those methods always evaluate the polynomials. The other
SPLINE-mode accessors do not raise: beyond the spline grid (0 K to 13000 K), they extrapolate
the edge cubic.

**Solution:**
1. Keep `T` within the coverage above, or use SPLINE mode for the first-order accessors
2. Note where the tabulated range ends (~200–6000 K for most gases); outside it the value comes
   from the extrapolations the factory synthesized, not from the table
3. Sanity-check extrapolated values against measurement — nothing warns you that you have left
   the tabulated range

### Issue: Missing Transport Data

```
mu(T) or lambda(T) returned 0
```

The species has no record for that property in `trans.inp` and no complete critical point to
synthesize one from. The mode does not matter: synthesized transport is evaluated from its spline
in POLY mode as well.

**Solution:**
1. Check with `gas->has_viscosity()` before calling `mu(T)` and with `gas->has_conductivity()`
   before calling `lambda(T)`
2. Use estimation methods (Chapman-Enskog, etc.) for missing data
3. Add transport data to `trans.inp` if available from literature

### Issue: Slow Performance

**Symptom:** Property evaluation is slower than expected

**Solution:**
1. Switch to spline mode: `gas->set_mode(RefGasMode::SPLINE)`
2. Reuse RefGas objects instead of recreating
3. Avoid creating RefGasFactory repeatedly
4. Profile with BELFEM Profiler to identify bottlenecks

---

## Performance Notes

No measured timings are documented for this module; time your own case with the snippet in
Polynomial vs Spline Mode. The structural costs are:

| Operation | Cost |
|-----------|------|
| `RefGasFactory` construction | parses the four data files once |
| `create_refgas()` | reads one record, synthesizes the extra intervals, and samples the heat, viscosity and conductivity splines on the 2601-point grid |
| `set_mode()` | pointer rebind; nothing is rebuilt |
| POLY-mode evaluation | linear search over the species' intervals, then the polynomial |
| SPLINE-mode evaluation | direct index into the uniform grid, then one cubic |

Choosing a mode has no build cost to amortize: `create_refgas()` builds the splines in either
case, and `set_mode()` only rebinds pointers.

---

## The `gastable` Command-Line Tool

`main.cpp` in this module builds an executable `gastable` that prints a property table for one
species. It is built only when the tree is configured with both `-DUSE_GASMODELS=ON` and
`-DUSE_EXAMPLES=ON` (`USE_EXAMPLES` defaults to OFF; `src/physics/gastables/CMakeLists.txt`).

```
gastable --gas N2 --Tmin 200 --Tmax 2000 --deltaT 50
```

| Option | Short | Meaning | Default |
|--------|-------|---------|---------|
| `--gas` | `-g` | species label (case-sensitive) | — |
| `--Tmin` | `-a` | lowest temperature in K | 200 |
| `--Tmax` | `-b` | highest temperature in K | 2000 |
| `--deltaT` | `-d` | temperature step in K | 50 |
| `--molar` | `-m` | print molar heats (J/mol) instead of per-kg | off |
| `--version` | `-V` | print the version banner | — |
| `--help` | `-h` | list the options | — |

The table columns are `T`, `cp`, `h`, `s`, `mu`, `lambda` (or `Cp`, `H`, `S` with `--molar`).
`h` is in kJ/kg (kJ/mol), `mu` in µPa·s, and `lambda` in mW/(m·K). `--Tmin` is raised to at
least `gDeltaT` = 5 K, and `--Tmax` is capped at `gTmax` = 13000 K.

---

## Further Reading

### Module Documentation
- `README.md` - Quick reference and API overview
- `src/physics/gasmodels/doc/gasmodels_usage_guide.md` - Using gastables in Gas class

### External References
- **NASA CEA:** Gordon & McBride (1994), "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications", NASA RP-1311
- **NIST Chemistry WebBook:** https://webbook.nist.gov/chemistry/
- **Poling et al.:** "The Properties of Gases and Liquids" (5th ed., 2001), Appendix A

---

**Last Updated:** 2026-10-22
**Maintainer:** BELFEM development team
