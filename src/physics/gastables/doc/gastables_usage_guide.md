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
real H_molar = gas->H(T);           // J/mol
real S_molar = gas->S(T);           // J/(mol·K)

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

**Reference state properties:**

```cpp
// Properties at 298.15 K
real H_ref = gas->H_ref();                         // Molar enthalpy at 298.15 K (J/mol)
real h_ref = gas->h_ref();                         // Specific enthalpy at 298.15 K (J/kg)
real H_formation = gas->reference_formation_enthalpy();  // Formation enthalpy (J/mol)
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
```

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
```

---

## Performance Considerations

### Polynomial vs Spline Mode

**Spline mode (default: the factory leaves every gas in it)**
- Cubic splines sampled off the NASA CEA polynomials, built once by `RefGasFactory::create_refgas()`
- ~2-3× faster evaluation than polynomial (no interval search)
- Interpolation accuracy

**Polynomial mode**
- Evaluates NASA CEA 9-coefficient polynomial directly; exact, with a linear interval search per call
- `set_mode()` is a pointer rebind in either direction; nothing is rebuilt

```cpp
RefGas * gas = factory.create_refgas("N2");   // arrives in SPLINE mode

// Benchmark polynomial mode
gas->set_mode(RefGasMode::POLY);
Timer t_poly;
for (int i = 0; i < 100000; ++i) {
    real T = 300.0 + i * 0.01;
    real cp = gas->cp(T);
}
real time_poly = t_poly.stop();

// Switch back to spline mode (pointer rebind, nothing is rebuilt)
gas->set_mode(RefGasMode::SPLINE);

// Benchmark spline mode
Timer t_spline;
for (int i = 0; i < 100000; ++i) {
    real T = 300.0 + i * 0.01;
    real cp = gas->cp(T);
}
real time_spline = t_spline.stop();

// Typical result: time_spline ≈ 0.4 × time_poly (2.5× faster)
```

**When to use splines:**
- Evaluating properties at >100 different temperatures
- Inside time-stepping loops
- During iterative solvers

**When to use polynomials:**
- One-off property evaluations
- Small number of temperature points (<10)

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
// GOOD - sequential access, cache-friendly
Cell<RefGas*> gases(n_species);
Vector<real> cp_values(n_species);

real T = 500.0;
for (index_t i = 0; i < n_species; ++i) {
    cp_values(i) = gases(i)->cp(T);
}

// LESS GOOD - random access if T varies unpredictably
for (index_t i = 0; i < n_species; ++i) {
    real T_random = get_random_temperature();
    cp_values(i) = gases(i)->cp(T_random);
}
```

---

## Advanced Topics

### Accessing Internal Polynomials

For debugging or custom property evaluations:

```cpp
RefGas * gas = factory.create_refgas("H2");

// Get thermodynamic polynomial for a specific temperature
real T = 1500.0;
HeatPoly * poly = gas->find_heat_poly(T);  // Protected - friend access only

// Alternative: work with exposed splines
gas->set_mode(RefGasMode::SPLINE);
Spline * heat_spline = gas->heat_spline();
Spline * visc_spline = gas->viscosity_spline();
Spline * cond_spline = gas->conductivity_spline();

// Evaluate the spline directly. The heat spline samples ENTHALPY, so the
// derivative ladder is shifted by one: eval is H, deval is cp, ddeval is dcp/dT.
real h       = heat_spline->eval(T);    // enthalpy
real cp      = heat_spline->deval(T);   // specific heat
real dcp_dT  = heat_spline->ddeval(T);  // its temperature derivative
```

Prefer the public accessors `H( T )`, `Cp( T )` and `dCpdT( T )`: in SPLINE mode they
dispatch to the spline with the derivative ladder already accounted for.

### Handling Missing Data

Some gases lack complete data:

```cpp
RefGas * rare_gas = factory.create_refgas("SomeRareGas");

// Check before using
if (rare_gas->has_thermo()) {
    real cp = rare_gas->cp(300.0);
} else {
    // Fallback: estimate or use default
    real cp_estimated = estimate_cp_from_structure(rare_gas);
}

if (rare_gas->has_viscosity()) {
    real mu = rare_gas->mu(300.0);
} else {
    // Use correlation (e.g., Chapman-Enskog)
    real mu_estimated = chapman_enskog_viscosity(rare_gas, 300.0);
}
```

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

### Creating Splines with Custom Temperature Steps

For advanced users who want to control spline discretization:

```cpp
RefGasFactory factory;
RefGas * gas = factory.create_refgas("Ar");

// Create custom temperature steps
Vector<real> T_steps(100);
for (index_t i = 0; i < T_steps.length(); ++i) {
    T_steps(i) = 200.0 + i * 50.0;  // 200 K to 5150 K in 50 K steps
}

// Create help matrix for spline construction
SpMatrix help_matrix;
factory.create_helpmatrix(help_matrix);

// Build splines (protected method - requires friend access or subclassing)
// Normally done once by RefGasFactory::create_refgas(); set_mode() is a pointer rebind
```

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
```

Note `E 01` with a space where a sign would normally sit — the readers parse fixed columns
(`cl_GT_InputTransport.cpp`), so the spacing is significant.

### Critical Point Data (`gasdata.inp`)

Fixed-column table of per-species constants. **The file's units are not the stored units:** the
reader converts the molar-mass column from g/mol to kg/mol and the critical-pressure column from
**bar** to Pa (`cl_GT_InputData.cpp:110,115`).

```
<symbol>  <name>  <CAS>  <M, g/mol>  <T_crit, K>  <p_crit, bar>  <Z_crit>  <omega>  <dipole>  <sources...>  <symmetry flag>
```

**Example — helium, as shipped** (`share/fluid/gasdata.inp`):
```
He          helium            7440-59-7   4.003   5.195   2.2832  0.3040 -0.3835  0.0   CoolProp-8.0.0 ... symmetry
```

(Helium's negative acentric factor is correct, not a sign error.)

### Cubic-EOS Alpha Coefficients (`cubicalpha.inp`)

The fourth shipped file. It carries the per-species alpha-function coefficients used by the cubic
equations of state, keyed by symbol and CAS with the critical temperature and pressure repeated:

```
<symbol>  <CAS>  <T_crit, K>  <p_crit, bar>  <alpha coefficients ...>
```

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
```

### Pattern 2: Temperature Sweep

```cpp
RefGas * gas = factory.create_refgas("Ar");
gas->set_mode(RefGasMode::SPLINE);  // Faster for sweep

Vector<real> T_range(100);
Vector<real> cp_range(100);

for (index_t i = 0; i < T_range.length(); ++i) {
    T_range(i) = 300.0 + i * 10.0;  // 300 K to 1290 K
    cp_range(i) = gas->cp(T_range(i));
}

delete gas;
```

### Pattern 3: Multi-Species Evaluation

```cpp
RefGasFactory factory;

Cell<string> species = {"H2", "O2", "N2", "H2O", "CO2"};
Cell<RefGas*> gases(species.length());

// Create all reference gases
for (index_t i = 0; i < species.length(); ++i) {
    gases(i) = factory.create_refgas(species(i));
}

// Evaluate properties
real T = 500.0;
for (index_t i = 0; i < gases.length(); ++i) {
    message(InfoLevel::Default,
            "%s: cp = %.3f J/(kg·K), mu = %.3e Pa·s",
            species(i).c_str(),
            gases(i)->cp(T),
            gases(i)->mu(T));
}

// Clean up
for (RefGas * gas : gases) {
    delete gas;
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

### Pattern 5: Formation Enthalpy for Equilibrium

```cpp
Cell<RefGas*> reactants, products;
// ... populate with species ...

real T = 1500.0;  // Combustion temperature

// Compute reaction enthalpy
real delta_H = 0.0;
for (RefGas * product : products) {
    delta_H += product->H(T);  // Molar enthalpy
}
for (RefGas * reactant : reactants) {
    delta_H -= reactant->H(T);
}

// delta_H includes formation enthalpies automatically
// (NASA CEA polynomials are absolute, not relative to elements)
```

---

## Troubleshooting

### Issue: Gas Not Found

```
BELFEM_ERROR: Gas species "XYZ" not found in database
```

**Solution:**
1. Check spelling (case-sensitive: "N2", not "n2")
2. Check if species exists in `thermo.inp`
3. Add custom data if needed

### Issue: Property Evaluation Fails

```
BELFEM_ERROR: Temperature T=20000.0 out of bounds for gas N2.
```

Raised in POLY mode only when `T` leaves `[0, gTmax]`; SPLINE mode clamps to the edge interval.

**Solution:**
1. Check that the species has a record at all: `gas->has_thermo()`
2. Note where the tabulated range ends (~200–6000 K for most gases); outside it the value comes
   from the extrapolations the factory synthesized, not from the table
3. Sanity-check extrapolated values against measurement — nothing warns you that you have left
   the tabulated range

### Issue: Missing Transport Data

```
mu(T) returned 0 - the species has no transport record and no critical point to synthesize one from
```

**Solution:**
1. Check with `gas->has_viscosity()` before calling `mu(T)`
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

## Performance Benchmarks

Typical timings on Intel Xeon (single core, 2.5 GHz):

| Operation | Polynomial Mode | Spline Mode | Notes |
|-----------|----------------|-------------|-------|
| Create RefGas | ~100 μs | ~100 μs | Data already loaded |
| `set_mode()` | ~0 | ~0 | Pointer rebind; splines are built once by the factory |
| Evaluate cp(T) | ~50 ns | ~20 ns | 2.5× faster |
| Evaluate h(T) | ~60 ns | ~20 ns | 3× faster |
| Evaluate mu(T) | ~40 ns | ~20 ns | 2× faster |

There is no break-even point to consider: the splines are built once by the factory, and `set_mode()` costs nothing.

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

**Last Updated:** 2026-01-30
**Maintainer:** BELFEM development team
