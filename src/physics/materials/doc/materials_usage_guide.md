# BELFEM Materials Module - Usage Guide {#physics_materials_materials_usage_guide}

**Module:** `src/physics/materials`
**Purpose:** How to create materials, query their properties, extend the roster, and understand the physical models behind every curve
**Date:** 2026-08-25
**Revision:** 2.0

---

## Revision History

| Date | Version | Changes |
|------|---------|---------|
| 2026-01-16 | 1.0 | Initial documentation (reorganized from README.md) |
| 2026-07-14 | 1.1 | Assembly contract: full-signature policy for jc/n dependency routing; undeformed-mesh density rule |
| 2026-08-25 | 2.0 | Rewritten against the current source: real roster and API, RRR-at-construction and the cached lookup tables, cryogenic α, Wachtman moduli, per-metal Kohler curves, formula alloys; the invented `Material_Metal` sketch, `populate_rho_lambda_databases`, `debye_cp` and the textbook-only model formulas removed |

---

## 1. Common Pitfalls

Read this section first.

### 1.1 Ownership moves on assignment

```cpp
material::JcFunction * tJc = tFactory.create_jc_function( 1e9, 5.0, 0.5, 2.0 );
tYbco->set_jc_function( tJc );   // the material now owns tJc and deletes it in its destructor
delete tJc ;                     // WRONG: double free
```

`set_jc_function()`, `set_n_function()` and `load_bh_curve()` / `set_bh_curve()` all transfer ownership. The base
`Material` destructor deletes the J_c and n functions; `Metal`'s destructor deletes the B-H
curve. The factory keeps no ownership of anything it returns: **the caller deletes the
material**.

### 1.2 A superconductor's resistivity is `rho_powerlaw()`, and its arguments are `(J, T, B, β)`

```cpp
real tRho = tHts->rho( T );                                          // normal-state resistivity — wrong below T_c
real tRho = tHts->rho_powerlaw( tNormJ, tT, tNormB, tAngleNxB );     // E–J power law
```

Because every argument is a `real`, swapped arguments still compile and run. The full overload is
`rho_powerlaw( normJ, T, normB, angleNxB )`; see §7.3 for the assembly contract.

### 1.3 The thermal mass matrix takes the density of the undeformed mesh

```cpp
real tDensity = tMaterial->density( gTroom );   // or ref_density()
real tDensity = tMaterial->density( T );        // WRONG for the mass matrix: the mesh does not expand
```

§7.4 has the reasoning.

### 1.4 Give a metal its RRR at construction

```cpp
Material * tCopper = tFactory.create_material( "copper", 100.0 );
```

`set_RRR()` fixes the residual resistivity ρ₀. Every material constant except `mu` (μ₀), `q` (1)
and `density_correction` (1) starts as NaN, so a metal without an RRR has no ρ(T) and no λ(T). `Metal::set_RRR` also builds the field-dependent
lookup table when tables are enabled (the factory's third argument, default `true`). A formula
alloy takes its RRR through the same factory argument; it does not override `set_RRR`, so the
inherited `Material::set_RRR` aborts with `BELFEM_ERROR`.

### 1.5 The field-dependent accessors are `( T, B, beta )`

```cpp
real tRho    = tMetal->rho( T, B, beta );      // |B| in tesla, beta = angle between B and J in radians
real tLambda = tMetal->lambda( T, B, beta );
```

Temperature comes first. For a `Metal`, these evaluate Kohler's rule directly when no table
exists; for an `Alloy`, they read the table and assert its presence.

### 1.6 Two different angles

`MaterialDependency::angleBxJ` is the field-current angle used by the metals' Kohler channel.
`MaterialDependency::angleNxB` is the tape-normal-field angle used by the HTS J_c channel,
unfolded to [0, π] since 2026-08-16. The FEM calculator routes on `depends()`, so declare the
one you actually consume.

### 1.7 Polynomial coefficients are in descending order

```cpp
// lambda( T ) = 0.01 T + 50
mat->set_user_defined_polynomial( MaterialProperty::lambda, std::vector< real >{ 0.01, 50.0 } );
```

Highest power first (`fn_polyval.hpp`); the overloads take `Cell< real >` or
`std::vector< real >`.

### 1.8 Check `have()` before querying a user-defined material

```cpp
if ( tCustom->have( MaterialProperty::lambda ) ) { real tK = tCustom->lambda( 300.0 ); }
```

A user-defined material has only the properties its init function set. Querying an unset
property dispatches through a null routing pointer.

### 1.9 The label `buffer` is reserved

`buffer` names the thin-shell φ-formulation buffer layer; the factory maps it to `Magnesia` and
rejects a user-defined (`usermat`) material registered under that name; a built-in section with
that label resolves to `Magnesia` and is not checked.

---

## 2. Architecture

### 2.1 Class hierarchy

```
Material                                   property routing, constants, dependencies, jc/n, B-H curve
└── SplineLookupTable                      per-property cubic splines on a uniform T grid
    ├── Metal                              Bloch–Grüneisen ρ, Hust λ, Kohler, c_p, cryo α, Wachtman E, lookup tables
    │   ├── Copper, Silver, Aluminum, Chromium, Indium, WhiteTin, Lead
    │   ├── Ferromagnetic                  magnon resistivity, reduced magnetization, Debye from ρ
    │   │   ├── Iron
    │   │   └── Nickel
    │   ├── HastelloyC276                  lookup alloy: its own tabulated fits
    │   └── YBCO                           HTS: normal-state Metal + power law + Callaway λ
    ├── Alloy                              homogenized from Metal components ( formula alloys )
    └── Magnesia                           ceramic buffer
└── UserDefinedMaterial                    plugin from a shared library
```

Supporting classes: `MaterialFactory`; `BhCurve` (interface) and `BhSplineCurve` (HDF5-backed);
`JcFunction` with `JcFunctionModifiedKim`, `JcFunctionDatabase`, `JcFunctionUserDefined`;
`Abundance` (molar masses, isotope mass variance).

### 2.2 Property routing

Every material property — `E`, `nu`, `cp`, `lambda`, `rho`, `alpha`, `density`, `debye`,
`mu`, … — is served through a function pointer. During construction, the material points it at
one of three sources: a **constant** (`set_constant`), a **spline** sampled on a uniform grid
(`create_spline` / `set_spline`), or a **custom** analytic routine (`set_custom`, dispatching to
the virtual `*_custom` methods). Derivatives (`dcpdT`, `drhodT`, `dlambdadT`, …) follow the same
routing. `have( property )` reports whether any source was set; `depends( property, dependency )`
reports what the source varies with.

### 2.3 `MaterialType`

`UserDefined`, `Ferro`, `HTS`, `PureMetal`, `LookupAlloy`, `CompositeAlloy`, `NonMetal`. In
practice: the nine metals *and* the formula `Alloy` are `PureMetal`; `HastelloyC276` is
`LookupAlloy`; `YBCO` is `HTS`; `Magnesia` is `NonMetal`. `Ferro` is not used by any built-in
class — Iron and Nickel are `PureMetal` through `Ferromagnetic`. The FEM side interrogates
`depends()` rather than the type wherever it can.

---

## 3. Core Classes

### 3.1 `Material`

```cpp
// availability and dependencies
bool  have( MaterialProperty ) const ;
bool  depends( MaterialProperty, MaterialDependency ) const ;
bool  is_constant( MaterialProperty ) const ;
real  constant_property( MaterialProperty ) const ;
void  set_constant( MaterialProperty, real ) ;

// temperature-dependent accessors, all [SI]
real E( T ), nu( T ), K( T ), G( T ), cp( T ), lambda( T ), rho( T ), rho_i( T ),
     alpha( T ), density( T ), debye( T ), Rp02( T ) ;
real dcpdT( T ), d2cpdT2( T ), drhodT( T ), dlambdadT( T ), dEdT( T ) ;
real ref_density() const ;

// field-dependent ( metals and alloys )
real rho( T, B, beta ), lambda( T, B, beta ) ;
real drhodT( T, B, beta ), drhodB( T, B, beta ), drhodbeta( T, B, beta ) ;

// magnetics
real mu( H [, T] ), H( B ) ;  void dmudH( H, real & mu, real & dmudH ) ;   // returns both μ(H) and dμ/dH
void load_bh_curve( const BhCurve * ) ;           // takes ownership AND routes H, mu, dmudH to the curve
void set_bh_curve( const BhCurve * ) ;            // takes ownership only — use load_bh_curve()

// superconductors — see §7.3 for the full overload set
real rho_powerlaw( normJ, T, normB, angleNxB ) ;
real jc( normB, angleNxB, T ), n( normB, angleNxB, T ) ;
void set_jc_function( JcFunction * ), set_n_function( JcFunction * ) ;   // take ownership

// user-defined materials
void set_user_defined_function( MaterialProperty, MaterialDependency, MatFunc1 * ) ;
void set_user_defined_function( MaterialProperty, MaterialDependency, MaterialDependency, MatFunc2 * ) ;
void set_user_defined_function( MaterialProperty, MaterialDependency, MaterialDependency, MaterialDependency, MatFunc3 * ) ;
void set_user_defined_polynomial( MaterialProperty, const Cell< real > & ) ;
void set_user_defined_polynomial( MaterialProperty, const std::vector< real > & ) ;
```

with `MatFunc1 = real ( const Material *, real )` and the two- and three-argument analogs.

### 3.2 `MaterialFactory`

```cpp
Material * create_material( const string & aLabel,
                            const real     aRRR = BELFEM_QUIET_NAN,
                            const bool     aBuildTables = true );      // built-in or formula
Material * create_material( const string & aLibraryPath, const string & aLabel );   // plugin
void       print_material_list( std::ostream & );
BhCurve  * create_bh_curve( const string & aPath, const string & aLabel );          // HDF5 → BhSplineCurve
JcFunction * create_jc_function( real aJc0, real aB0, real aK2, real aAlpha );       // modified Kim
JcFunction * create_jc_function( const string & aPath, const string & aLabel );     // database
```

Labels are case-insensitive. Anything that is not a built-in label but parses as
`<Element><integer>…` becomes an `Alloy` (§4.6). File arguments go through
`material::data_file()`, which searches the run directory first and `$BELFEM_DATA/material`
second.

### 3.3 `BhCurve`

An interface — `B( H )`, `H( B )`, `mu( H )`, `nu( B )` — implemented by `BhSplineCurve`, which
the factory builds from an HDF5 B-H database.

### 3.4 `JcFunction`

`eval( B, angle )` and `eval( B, angle, T )`, each function overriding the arity that matches
the dependencies it declares through `JcParameter`; both base versions raise an error. Three implementations:
`JcFunctionModifiedKim` (analytic, §9.6), `JcFunctionDatabase` (tensor-mesh table over
(T, log₁₀B, angle), log-scaled values; T and B are clamped to the table range, the angle is
wrapped with π periodicity), `JcFunctionUserDefined`
(plugin).

### 3.5 `Abundance`

Molar masses and the isotope/alloy mass-variance parameter Γ from a stoichiometry, used by the
metals (`M`) and by YBCO's Callaway model.

---

## 4. Material Types

### 4.1 Pure metals — `Metal`

Every metal constructor follows a fixed sequence of `create_*()` calls; Copper is the reference
implementation. The order matters because later steps read earlier curves:

1. `set_constants()` — molar mass, melting point, reference density, Sommerfeld γ and Debye β
   (writing β after M derives `debye0K`).
2. `create_cp()` — the Sommerfeld–Debye cubic plus log-log Bézier segments (§9.3).
3. `create_alpha()` — the ΔL/L Bézier, then `create_cryo_expansion()`, which fits the
   cryogenic branch against c_p and builds the α spline and its integral (§9.4).
4. Debye curve and resistivity anchor — from c_p (`compute_debye_from_cv` at 12–14 K: Copper,
   Silver, Lead), from measured ρ(T) (`compute_debye_from_rho`: Iron, Nickel), or a fitted
   curve (Aluminum, Chromium, Indium, WhiteTin); then `set_rho_i_ref()` for the Bloch–Grüneisen
   amplitude.
5. `set_lambda_coefficients()` — the Hust form of λ(T) (§9.2).
6. `create_kohler()` — the longitudinal and transverse magnetoresistance curves (§9.5).
7. `create_mech()` — Wachtman E(T) and the Grüneisen-derived ν(T) (§9.7).
8. `set_RRR()` if an RRR was given — fixes ρ₀ and, with tables enabled, builds
   `<label>_RRR<n>.hdf5`.

```cpp
Material * tCu = tFactory.create_material( "copper", 200.0 );
real tRho77   = tCu->rho( 77.0 );
real tRhoB    = tCu->rho( 77.0, 10.0, 0.5 * constant::pi );   // 10 T transverse
real tLambdaB = tCu->lambda( 77.0, 10.0, 0.5 * constant::pi );
```

All nine metals provide c_p, α, E, ν, ρ(T), λ(T), θ_D(T), density, and a Kohler curve.

### 4.2 Ferromagnetic metals — `Iron`, `Nickel`

Add to `Metal`: a reduced-magnetization curve (Crangle & Goodman 1971), the magnon term of the
resistivity below the Curie point (`rho_mag`), a Debye curve inverted from measured resistivity
with that term removed, and a Bloch–Grüneisen exponent of 4.5 instead of 5. Nickel's E(T)
carries the ΔE dip across the Curie point as two Wachtman branches joined by a Bézier bridge.
Neither class supplies a B-H curve of its own; assign one with `load_bh_curve()` if the magnetic
solver needs it.

### 4.3 Lookup alloy — `HastelloyC276`

A `Metal` of type `LookupAlloy` with its own fitted c_p, α, E, ν, ρ(T), λ(T) and magnetic
susceptibility χ(T). No field dependence: it declares neither Kohler dependencies nor a table,
and only the `( T )` accessors apply.

### 4.4 High-temperature superconductor — `YBCO`

Normal-state c_p, α, E, ν, ρ(T) and a λ(T) from the Callaway phonon model plus a
Wiedemann–Franz electronic term (see `callaway_thermal_conductivity.md` — and its §8 before
trusting the phonon parameters: the electronic term dominates the fit), with `T_crit` = 92.5 K.
The superconducting resistivity is the E–J power law with J_c and n supplied as constants or
`JcFunction` objects:

```cpp
Material * tYbco = tFactory.create_material( "ybco" );
tYbco->set_jc_function( tFactory.create_jc_function( "jc_data.hdf5", "YBCO" ) );
tYbco->set_constant( MaterialProperty::n, 25.0 );
real tRho = tYbco->rho_powerlaw( tNormJ, tT, tNormB, tAngleNxB );
```

### 4.5 Ceramic — `Magnesia`

c_p, α, E, ν and λ(T); no electrical resistivity. This is the material behind the reserved
label `buffer`. Its α, like Hastelloy's and YBCO's, follows c_p below room temperature through
the anchored Grüneisen branch (`thermal_expansion_from_heat_capacity.md`, §7b).

### 4.6 Formula alloys — `Alloy`

`Sn60Pb40`, `Fe71Cr19Ni10`, … in integer mass percent, homogenized from the metal roster with
a caller-supplied RRR. What is mixed, how, and — above all — for which alloys the result is
trustworthy is the subject of `alloy_homogenization.md`. Read its §3 before using one for
anything cryogenic.

### 4.7 User-defined — `UserDefinedMaterial`

See §6.

---

## 5. The `material` Executable

```
material <label> [options]
  -r, --rrr <value>                  RRR ( default 100 for a metal, 10 for a formula alloy )
  -b, --field <tesla>                |B| for the field-dependent columns
  -a, --angle <degrees>              angle between B and J
  -t, --temperatures <min> <max> <step>
  -c, --create                       build the HDF5 lookup table
  -m, --mesh                         also write the table's tensor mesh as an Exodus file
  -l, --list                         print the built-in roster
```

`material copper -r 100 -b 5 -a 90 -t 4 300 4` prints ρ, λ, c_p, E, ν, α over the range as far
as each property exists. Without `-c` the executable constructs the material with tables
disabled, so a `Metal` still evaluates Kohler's rule directly and an `Alloy` cannot be asked
for field dependence.

---

## 6. User-Defined Materials

A user material is a shared library that exports one C-linkage function named `<Label>_init`.
That function receives the `Material` and sets its properties. `example_user_material.cpp` and
`UserMaterialTemplate.cmake` are the templates.

```cpp
#include "cl_Material.hpp"
using namespace belfem ;

real my_rho( const Material * mat, real T )
{
    return 1.7e-8 * ( 1.0 + 0.004 * ( T - 293.15 ) );
}

real my_lambda( const Material * mat, real T )
{
    return 2.44e-8 * T / mat->rho( T );        // other properties of the same material are reachable
}

extern "C" void MyAlloy_init( Material * mat )
{
    mat->set_constant( MaterialProperty::E,             200e9 );
    mat->set_constant( MaterialProperty::nu,            0.3 );
    mat->set_constant( MaterialProperty::ref_density,   8900.0 );
    mat->set_constant( MaterialProperty::T_ref_density, 293.15 );

    mat->set_user_defined_function( MaterialProperty::rho,    MaterialDependency::T, & my_rho );
    mat->set_user_defined_function( MaterialProperty::lambda, MaterialDependency::T, & my_lambda );

    // descending order: cp( T ) = 0.12 T + 385
    mat->set_user_defined_polynomial( MaterialProperty::cp, std::vector< real >{ 0.12, 385.0 } );
}
```

Build with the CMake template (C++17, the project standard), then

```cpp
Material * tAlloy = tFactory.create_material( "./libmyalloy.so", "MyAlloy" );
```

The label must match the function prefix; the factory resolves `<Label>_init` with `dlsym` and
fails loudly otherwise. Two- and three-argument functions declare their dependencies in order
(`normB`, `angleNxB`, `T` for a J_c, for instance) and receive the arguments in that order.
Superconducting user materials supply J_c and n through `JcFunction` objects, not through
`jc_custom` (§7.3).

Rules: use `extern "C"`; make the first parameter of every function `const Material *`; set only
the properties you have — everything else stays `have() == false`; declare dependencies truthfully
because the FEM dispatch trusts them.

---

## 7. Property Queries and Dependencies

### 7.1 Availability and dependencies

```cpp
bool tHasRho = tMat->have( MaterialProperty::rho );
bool tRhoT   = tMat->depends( MaterialProperty::rho, MaterialDependency::T );
bool tRhoB   = tMat->depends( MaterialProperty::rho, MaterialDependency::normB );
bool tConst  = tMat->is_constant( MaterialProperty::cp );
```

`depends()` is declarative: it tells the caller which overload carries information, but it does not
validate the call. `cl_FEM_Calculator` uses this flag to decide whether a block's resistivity
takes the `( T, B, beta )` path.

### 7.2 `MaterialDependency`

`T`, `normB` (in-plane field magnitude), `angleBxJ` (field–current angle, Kohler), `angleNxB`
(normal–field angle, HTS), `normH` (field strength, μ of ferromagnets), `normJ`, `Jc`, `rho`.

`angleNxB` is unfolded to [0, π] since 2026-08-16: β = acos(n·b/|b|), where +n is the
master-side facet normal (the layer-stack direction). β < π/2 means the field has a component
along +n. Measured jc(β) tables consume the full range; Kim-type analytic laws are even in β and
unaffected. Despite the historical name, this is NOT the angle between current and field — that
is the metals' `angleBxJ` Kohler channel.

### 7.3 Assembly contract: the full-signature policy for J_c / n

The FEM assembly path always calls the **full** T-bearing HTS overloads —

```cpp
rho_powerlaw   ( normJ, T, normB, angleNxB [, x, y, z, t] );
rho_piecewise  ( normJ, T, normB, angleNxB [, x, y, z, t] );
drho_powerlaw_dJ / drho_piecewise_dJ  // same argument lists
```

— regardless of what a given material actually depends on. **The material, not the caller,
decides which arguments matter:** internally, `Material::jc_eval()` and `Material::n_eval()`
route on the declared dependencies of the attached `JcFunction` / n-function
(`depends_on( JcParameter::T )` selects the `eval( B, angle, T )` or `eval( B, angle )`
override) and fall back to the plain constants (`MaterialProperty::jc` / `MaterialProperty::n`)
when no function is attached.

Consequences for material implementers:

- Provide jc/n through a `JcFunction` (or as constants via `set_constant`). Declare the
  dependencies truthfully — the routing, and the FEM dispatch above it, trust the `depends()` /
  `dependencies()` flags completely.
- Override the `eval()` arity that matches your declared T-dependence (both defaults raise
  errors).
- There is no T-only callback for J_c or n: registering either with the one-argument
  `set_user_defined_function( …, MaterialDependency::T, f )` is rejected at load time. User
  superconductors go through `JcFunction`.
- The reduced overloads (`(normJ)`, `(normJ, T)`, `(normJ, normB, angleNxB)`, and their defect
  twins) are not used by the assembly path, and nothing in the tree calls them. The `(normJ, T)`
  forms read `jc_custom(T)` / `n_custom(T)`, which no registration path populates any more, so
  they abort on every material; do not call them.
- ρₙ (the normal-matrix branch of the parallel combination) is always evaluated at the local
  temperature `T`; in magnetics-only problems the caller passes `T = gTbulk`, so nothing changes
  there.

### 7.4 Density and the undeformed mesh

Real density changes with temperature — but **all BELFEM computations run on the undeformed
mesh**, and the transport properties are already corrected for thermal expansion. The thermal
mass matrix must therefore use the density at which the mesh is *not* deformed, which is by
default room temperature:

```cpp
// CORRECT — reference density of the undeformed mesh:
real tDensity = tMaterial->have( MaterialProperty::density )
              ? tMaterial->density( gTroom )
              : tMaterial->ref_density();

// WRONG — double-counts thermal expansion:
real tDensity = tMaterial->density( T );
```

`calculator::MaxwellData` caches this value once at construction (`density()` accessor); do not
re-derive it per integration point. `density( T )` itself is ρ_ref / l(T)³ with l the relative
length from the α spline's integral.

---

## 8. Usage Patterns

**A conductor with magnetoresistance**

```cpp
Material * tCopper = tFactory.create_material( "copper", 2000.0 );   // table built here, cached
real tRho    = tCopper->rho( 77.0, 10.0, 0.5 * constant::pi );
real tLambda = tCopper->lambda( 77.0, 10.0, 0.5 * constant::pi );
```

**HTS with an analytic J_c**

```cpp
Material * tYbco = tFactory.create_material( "ybco" );
tYbco->set_jc_function( tFactory.create_jc_function( 1e9, 5.0, 0.5, 2.0 ) );  // Jc0, B0, k², α
tYbco->set_constant( MaterialProperty::n, 25.0 );
real tRho = tYbco->rho_powerlaw( 1e8, 77.0, 2.0, 0.5 );                        // J, T, B, angleNxB
```

**HTS with tabulated J_c and n**

```cpp
tYbco->set_jc_function( tFactory.create_jc_function( "jc_data.hdf5", "YBCO"   ) );
tYbco->set_n_function(  tFactory.create_jc_function( "n_data.hdf5",  "YBCO_n" ) );
```

**A ferromagnet with a B-H curve**

```cpp
Material * tIron = tFactory.create_material( "iron", 30.0 );
tIron->load_bh_curve( tFactory.create_bh_curve( "MatData/bhdata.hdf5", "Iron" ) );   // not set_bh_curve: that stores without routing
real tH  = tIron->H( 1.5 );      // A/m for B = 1.5 T
real tMu = tIron->mu( 1000.0 );  // H/m at H = 1000 A/m
```

**A solder from its composition**

```cpp
Material * tSolder = tFactory.create_material( "Sn60Pb40", 10.0 );
```

**A plugin material**

```cpp
Material * tCustom = tFactory.create_material( "./libmymat.so", "MyMaterial" );
if ( tCustom->have( MaterialProperty::rho ) ) { real tRho = tCustom->rho( 300.0 ); }
```

---

## 9. The Physical Models

Each model is described as the code implements it; the theory documents contain the derivations.

### 9.1 Electrical resistivity of a metal — Bloch–Grüneisen with an effective Debye curve

ρ(T) = ρ₀ + ρ_i(T),  ρ_i(T) = A · J_n(θ(T)/T) / (θ(T)/T)ⁿ

with J_n the Debye integral (tabulated once per exponent by `debye.f90`), n = 5 for the simple
metals and 4.5 for Iron and Nickel (`n_bloch_gruen`), and θ(T) the *effective* Debye curve of
the class — derived from c_p, from measured ρ(T), or fitted — so that the one-parameter formula
reproduces the measured phonon resistivity. `set_rho_i_ref( T_ref, ρ_i,ref [, θ] )` fixes A at
a literature point. ρ₀ follows from the RRR by the regula falsi in `Metal::set_RRR`:

RRR = ( ρ_i(T_ref) + ρ₀ ) / ρ₀  at `T_ref_rho_i` (273.15 K for Copper).

Ferromagnets add a magnon term `rho_mag(T)` below the Curie point.

### 9.2 Thermal conductivity of a metal — Hust

λ(T) = 1 / ( w₀ + w_i + w_i0 ),  w₀ = ρ₀ / (L₀ T),  w_i = hust(P, T),  w_i0 = C·w_i·w₀ / (w_i + w₀)

with the eight Hust coefficients per metal (`set_lambda_coefficients`, NBS Special Publication
260-90). In field, λ(T, B, β) = λ(T) · ρ(T) / ρ(T, B, β) — a Wiedemann–Franz closure with L = L₀.

### 9.3 Specific heat

Sommerfeld–Debye cubic γT + βT³ below 0.02·θ_D, a fifth-order beam polynomial to the first
control point, two or three cubic Béziers in ln c_p over ln T, and a linear tail in T. The
representation and its constraints are in `thermal_expansion_from_heat_capacity.md`, §3.

### 9.4 Thermal expansion

Above a split temperature T\* = min(0.618·θ_D, 273.15 K) the tangent coefficient
α = (1/L) dL/dT comes from a Bézier fitted to ΔL/L referred to 293.15 K. Below it, α = C(T)·c_p(T)
with ln C a cubic matched at T\* — the Grüneisen equation of state, which gives the correct
T³ approach to zero that a ΔL/L fit cannot. `thermal_expansion_from_heat_capacity.md` is the
full account, including the two guards and the Grüneisen diagnostic that reads O(2) for every
metal.

### 9.5 Magnetoresistance — Kohler curves

Δρ/ρ_0T = kohler( B, S, β ) with the similarity parameter S = ρ_ref / ρ_0T, where ρ_0T is the
zero-field resistivity at the current temperature and ρ_ref the one the curve was fitted at.
Each metal carries its own longitudinal and transverse curves in ln(B·S) — Bézier or polynomial
segments with a low-field Hermite cubic, a mid-field fit and a saturated or linear tail — and
combines them by Pippard's angular interpolation A_∥ cos²β + A_⊥ sin²β. Data sources are cited
per class (Lüthi 1960, Fickett 1972, Klaffky & Coleman 1974, Kozlova & Kondorskii 1963, Arajs &
Dunmyre 1965, …). The table `<label>_RRR<n>.hdf5` stores ln ρ over (T, log₁₀B, β) on a tensor
mesh; `Metal` evaluates the curves directly when no table exists, `Alloy` requires the table.

### 9.6 HTS power law and the modified Kim model

E = E_c (J/J_c)ⁿ  →  ρ_powerlaw = (E_c / J_c) · (J/J_c)^(n−1), with E_c = 10⁻⁴ V/m by default.

`JcFunctionModifiedKim( Jc0, B0, k², α )`:

J_c(B, θ) = Jc0 / [ 1 + √(k² sin²θ + cos²θ) · B/B0 ]^α

— the anisotropic Kim–Anderson form; the factory's four arguments map to (Jc0, B0, k², α) in
that order. Tabulated J_c and n come from `JcFunctionDatabase`.

### 9.7 Elastic moduli — Wachtman E and a Grüneisen-derived ν

E(T) = E₀ − b·T·exp(−T₀/T), fitted per metal (mostly against Blanke, *Thermophysikalische
Stoffgrößen*, 1989; Chromium against Armstrong & Brown 1964; Indium against Kim & Ledbetter
1998). ν(T) is not tabulated: `Metal::create_mech( E₀, b, T₀, T₂, ν₂ )` takes one anchor
(T₂, ν₂), forms K_T = E/(3(1 − 2ν)), the adiabatic K_S = K_T / (1 − T α_V² K_T /(ρ c_p)) and the
Grüneisen parameter γ = α_V K_S / (ρ c_p), holds γ constant, and recovers ν(T) = ½ − E/(6K_T(T))
along the spline grid with K_S(T) = γ ρ c_p / α_V. The 0 K value is extrapolated with zero
slope; the result is checked to stay inside (−1, ½). Nickel supplies its own E(T) (two Wachtman
branches and a Bézier bridge across the Curie point) and only the ν construction.

Isothermal versus adiabatic: the tabulated moduli are dynamic (adiabatic); the difference on E is
0.3–0.5 % at room temperature and vanishes at cryogenic temperatures, and is not corrected.

### 9.8 Callaway phonon conductivity

Used by YBCO only; `callaway_thermal_conductivity.md`.

---

## 10. `MaterialProperty` reference

This table is generated from `cl_Material.hpp`. The index is the storage slot, and `UNDEFINED`
sizes the arrays. Constants and non-constant properties share the enumeration; everything from
`T_crit` upward is a constant.

| Property | Index | Unit | Meaning |
|---|---|---|---|
| `density` | 0 | kg/m³ | density |
| `E` | 1 | Pa | Young's modulus |
| `nu` | 2 | - | Poisson's ratio |
| `cp` | 3 | J/(kg·K) | specific heat capacity |
| `lambda` | 4 | W/(m·K) | thermal conductivity |
| `mu` | 5 | H/m | magnetic permeability, defined as ∂B/∂H |
| `rho` | 6 | Ω·m | electric resistivity |
| `alpha` | 7 | 1/K | thermal expansion coefficient, defined as (1/l)·∂l/∂T, not (1/l)·Δl/ΔT! |
| `Rp02` | 8 | Pa | yield stress |
| `debye` | 9 | K | Debye temperature |
| `rho_i` | 10 | Ω·m | inner resistivity of noble metal |
| `kohler_trans` | 11 | - | Kohler parameter for transverse magnetoresistance |
| `kohler_long` | 12 | - | Kohler parameter for longitudinal magnetoresistance |
| `T_crit` | 13 | K | critical temperature |
| `M` | 14 | kg/mol | molar mass |
| `Gamma` | 15 | - | impurity parameter |
| `R` | 16 | J/(kg·K) | specific gas constant |
| `T_max` | 17 | K | maximum temperature |
| `ref_density` | 18 | kg/m³ | reference density |
| `T_ref_density` | 19 | K | temperature at reference density |
| `gamma` | 20 | J/(kg·K²) | linear Debye parameter, cv = γ·T + β·T³ = ∂cp/∂T at T=0 |
| `beta` | 21 | J/(kg·K⁴) | cubic Debye parameter, cv = γ·T + β·T³ |
| `q` | 22 | - | number of atoms per molecule  (default: 1) |
| `debye0K` | 23 | K | Debye temperature at 0 K |
| `rho_i_ref` | 24 | Ω·m | reference inner resistivity |
| `T_ref_rho_i` | 25 | K | temperature at reference inner resistivity |
| `A_bloch_gruen` | 26 | - | A-parameter for Bloch-Grüneisen law |
| `n_bloch_gruen` | 27 | - | n-parameter for Bloch-Grüneisen law |
| `RRR` | 28 | - | residual resistivity ratio |
| `rho_0` | 29 | Ω·m | residual resistivity |
| `layer_thickness` | 30 | m | characteristic length, e.g., layer thickness |
| `ec` | 31 | V/m | critical electric field |
| `jc` | 32 | A/m² | critical current density |
| `n` | 33 | - | exponent for power law |
| `Tcurie` | 34 | K | Curie temperature |
| `A_electron_magnon` | 35 | Ω · m / K² | electron-magnon scattering parameter |
| `A_spin_disorder` | 36 | Ω · m | spin disorder amplitude |
| `grueneisen` | 37 | — | Grüneisen parameter |
| `density_correction` | 38 | - | density scaling factor, to correct solder thickness |
| `rho0_pure` | 39 | — | resistivity from TPCR for pure metals |

Units and meanings are the header's own comments. Two of them are imprecise as of
2026-08-25: `A_bloch_gruen` is marked dimensionless but carries Ω·m, and `rho0_pure` has no unit
and is a resistivity [Ω·m].

The `mu` row's "defined as ∂B/∂H" is imprecise, and worth unpicking because the Newton path
depends on it. `mu` is index 5, i.e. inside the **non-constant** range — a material may hold a
constant permeability, or a B-H curve, or a custom `mu( H )`, and loading a curve repoints
`mFunctionMu` at `mu_bhcurve` (`cl_Material.cpp:958-960`).

- For a **constant** permeability B = μH is linear, so B/H and ∂B/∂H are the same number and the
  label is harmless. `dmudH_const()` returns a zero derivative alongside it
  (`cl_Material.hpp:2112-2116`).
- As soon as μ varies with H they part company. `mu` is then the **secant** permeability B/H
  (`cl_BhCurve.hpp:31`), and `dmudH` is **dμ/dH** (`:32`, `cl_Material.hpp:1043`) — *not* the
  differential permeability.

If you need ∂B/∂H — the tangent a Newton iteration wants — build it from both. For B = μ(H)·H,

```
dB/dH = mu + H * dmudH
```

---

## 11. Thread Safety

BELFEM parallelizes with MPI and is deliberately not thread-safe internally
(`doc/coding_philosophy.md`). Within that policy the materials module behaves as follows:

- Construct and configure every material on the master thread, before any parallel region;
  `MaterialFactory` and every `set_*()` are unsynchronized.
- After construction the property accessors are `const` and touch no mutable state in
  `Material`, `Metal`, `BhSplineCurve` or `JcFunction`; the lookup-table evaluation is `const`
  as well. Concurrent read-only queries from OpenMP threads are therefore expected to be safe,
  but this is a design intent, not a verified contract — protect the calls with
  `#pragma omp critical` if in doubt.
- Never mutate a material after it has been assigned to mesh blocks (see the contracts
  document).

---

## 12. Performance

- A metal with tables enabled builds `<label>_RRR<n>.hdf5` once per (label, RRR) and loads it
  afterwards; delete the file to force a rebuild, and expect an automatic rebuild when a file
  predates the current format. The build samples every node of a tensor mesh over
  (T, log₁₀B, β) on rank 0 and projects collectively — it is a one-off cost at construction, not
  something to trigger in a loop.
- Table queries are shape-function interpolations; direct Kohler evaluation (a `Metal` without
  a table) costs a few polynomial and exponential evaluations per call.
- Constants are cheapest, splines next, custom routines last; polynomial user properties are
  spline-free and fast.

No wall-clock figures are recorded in the tree; measure before relying on any.

---

## 13. See Also

- `materials_contracts_and_invariants.md` — ownership, lifetime, table-build and immutability contracts
- `alloy_homogenization.md` — formula alloys
- `thermal_expansion_from_heat_capacity.md` — c_p and α
- `callaway_thermal_conductivity.md` — YBCO's λ_ph kernel
- `../../database/doc/database_usage_guide.md` — the lookup tables
- `doc/coding_philosophy.md`, `doc/documentation_guidelines.md`
