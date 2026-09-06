# Materials Module Documentation {#physics_materials_index}

**Module:** `src/physics/materials`
**Purpose:** Index of the materials documentation, the current material roster, and the shortest path to a working call
**Date:** 2026-08-25

---

## What the module provides

- Temperature-dependent properties — c_p, α, E, ν, ρ, λ, density, Debye temperature — for a fixed roster: nine pure metals, one lookup alloy, one high-temperature superconductor and one ceramic. Each curve is fitted to a data source cited in its class.
- Field-dependent resistivity and thermal conductivity for the metals, using per-metal Kohler magnetoresistance curves evaluated through a cached lookup table `<label>_RRR<n>.hdf5`.
- Homogenized alloys from a composition formula such as `Sn60Pb40` or `Fe71Cr19Ni10`, mixed from
  the pure-metal roster with a caller-supplied residual resistivity ratio (RRR).
- HTS power-law resistivity with J_c and n as functions of field, angle and temperature.
- B-H curves for ferromagnetic materials, plus a plugin path for user-defined materials.

## Documents

| Document | Use when |
|---|---|
| [materials_usage_guide.md](materials_usage_guide.md) | You call the API: factory, accessors, dependencies, user-defined materials, the physical models behind each curve, thread safety |
| [materials_contracts_and_invariants.md](materials_contracts_and_invariants.md) | You hand objects around: ownership, lifetime, the lookup-table build, immutability during a simulation |
| [resistivity_laws.md](resistivity_laws.md) | You choose or debug an HTS E-J law: what `powerlaw`, `piecewise` and `riva` compute, their tangents, degenerate windows, and which to use |
| [alloy_homogenization.md](alloy_homogenization.md) | Building a material from a formula: what is mixed, how it is mixed, and which alloy results can be trusted |
| [thermal_expansion_from_heat_capacity.md](thermal_expansion_from_heat_capacity.md) | You fit or debug a cryogenic α(T): the Grüneisen branch, its guards, the c_p representation it rests on |
| [callaway_thermal_conductivity.md](callaway_thermal_conductivity.md) | You wire the Callaway phonon-conductivity kernel (YBCO's λ_ph) or read its 20-parameter interface |
| [material_property_sources.md](material_property_sources.md) | You need the literature behind a fit: per-material source table for α, c_p, ρ, E, ν with the full bibliography and the known gaps |
| [../../database/doc/database_usage_guide.md](../../database/doc/database_usage_guide.md) | You touch the lookup tables behind ρ(T, B, β) and λ(T, B, β): projection, collective build, clamping conventions |

## Material roster

The authoritative list is `material --list` (or `MaterialFactory::print_material_list`). Labels are case-insensitive; element symbols are accepted as aliases.

| Label | Class | Fitted properties |
|---|---|---|
| `aluminum` (`aluminium`, `al`) | `Aluminum` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T) |
| `chromium` (`cr`) | `Chromium` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T); Néel/spin-flip anomalies smoothed |
| `copper` (`cu`) | `Copper` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T) |
| `silver` (`ag`) | `Silver` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T) |
| `indium` (`in`) | `Indium` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T) |
| `tin` (`sn`) | `WhiteTin` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T) |
| `lead` (`pb`) | `Lead` | c_p, α, E/ν, ρ(T), λ(T), Kohler, θ_D(T) |
| `iron` (`ferro`, `fe`) | `Iron` (ferromagnetic) | as above plus magnetization and the magnon resistivity term |
| `nickel` (`ni`) | `Nickel` (ferromagnetic) | as above; E carries the ΔE dip across the Curie point |
| `hastelloyc276` (`hastelloy`) | `HastelloyC276` | c_p, α, E/ν, ρ(T), λ(T), magnetic susceptibility χ(T); no field dependence |
| `ybco` | `YBCO` | normal-state c_p, α, E/ν, ρ, λ (Callaway); power-law ρ with J_c(B, β, T) and n |
| `magnesia` (`mgo`, `buffer`) | `Magnesia` | c_p, α, E/ν, λ(T); electrical insulator |
| `<El><pct><El><pct>…` | `Alloy` | homogenized from the metals above, see [alloy_homogenization.md](alloy_homogenization.md) |

Each class states its data sources in the comments of its `create_*()` routines; the
`material` executable prints the resulting tables.

## The shortest working call

```cpp
#include "cl_MaterialFactory.hpp"

MaterialFactory tFactory ;

// label, RRR ( optional ), build the field-dependent lookup tables ( default true )
Material * tCopper = tFactory.create_material( "copper", 100.0 );

real tRho    = tCopper->rho( 77.0 );                               // Ω·m, zero field
real tRhoB   = tCopper->rho( 77.0, 5.0, 0.5 * constant::pi );      // T [K], |B| [T], angle B∠J [rad]
real tLambda = tCopper->lambda( 77.0, 5.0, 0.5 * constant::pi );   // W/(m·K)
real tCp     = tCopper->cp( 77.0 );                                // J/(kg·K)
real tE      = tCopper->E( 77.0 );                                 // Pa

delete tCopper ;   // the caller owns what the factory returns
```

The argument order of the field-dependent accessors is `( T, B, beta )`. When tables are enabled, the first construction of a metal builds and writes `Copper_RRR100.hdf5` (the class label, not the factory label) into the run directory; later constructions load it.

To print the same tables from the command line:

```
material copper --rrr 100 --field 5 --angle 90 --temperatures 4 300 4
material Fe71Cr19Ni10 --rrr 1.13
material --list
```

## Where data files are looked for

Every path handled by the factory — a B-H database, a J_c/n database, or a material or defect plugin — goes through `material::data_file()` (`fn_material_data_path.hpp`), which tries three locations in order:

1. the path as written, relative to the run directory, or absolute;
2. the same relative path below `$BELFEM_DATA/material`;
3. the file name alone below `$BELFEM_DATA/material`.

The run directory wins, so a local copy always overrides the shared database.
Step 3 is what lets `MatData/bhdata.hdf5` in an input file resolve in a run
directory that has no `MatData` of its own. If none of the three exists the
path is passed on unchanged, so the error names the file as written — and a
plugin name still reaches `dlopen` and its `$LD_LIBRARY_PATH` search.
`$BELFEM_DATA` is read by `Communicator::set_globals()` into a global, and on
an installed tree that global falls back to the compiled-in install data
directory when the variable is unset — so an unset `$BELFEM_DATA` leaves only
step 1 in a build-tree run, but not necessarily in an installed one. Since 2026-08-31 the search order itself lives one layer down, in
`belfem::search_data_file()` (`io/filetools.hpp`); `material::data_file()` is a
forwarder that supplies the `material` subdirectory. The move was made so that
a **user-defined source function's** `.so` — loaded from `numerics/sources`,
which sits below `physics` and cannot include this module — resolves the same
way. `gastables::data_path()` deliberately does *not* share it: it is a
directory resolver with marker-file validation and its own relative-path
ladder.

What the shipped files contain — the six critical-current
tables and the B-H curves — is indexed in
`share/material/README.md` at the repository root.

## Three rules that cost the most when broken

```cpp
// 1. Ownership moves on assignment — never delete what you handed over
material::JcFunction * tJc = tFactory.create_jc_function( 1e9, 5.0, 0.5, 2.0 );
tYbco->set_jc_function( tJc );
delete tJc ;                                   // WRONG: double free

// 2. A superconductor's resistivity is rho_powerlaw(), not rho()
real tRhoHts = tYbco->rho_powerlaw( tNormJ, tT, tNormB, tAngleNxB );

// 3. The thermal mass matrix takes the density of the undeformed mesh
real tDensity = tMaterial->density( gTroom );  // not density( T )
```

The usage guide opens with the full list and the reason for each rule.

## Source map

| Part | Files |
|---|---|
| Factory | `cl_MaterialFactory.{hpp,cpp}` |
| Base classes | `cl_Material.{hpp,cpp}` → `cl_Material_SplineLookupTable.{hpp,cpp}` → `cl_Material_Metal.{hpp,cpp}` → `cl_Material_Ferromagnetic.{hpp,cpp}` |
| Materials | `cl_Material_<Name>.{hpp,cpp}` for every roster entry; `cl_Material_Alloy.{hpp,cpp}`; `cl_Material_UserDefined.{hpp,cpp}` |
| Property functions | `cl_BhCurve.{hpp,cpp}`, `cl_Material_BhSplineCurve.{hpp,cpp}`, `cl_JcFunction*.{hpp,cpp}` |
| Composition | `cl_Material_Abundance.{hpp,cpp}` (molar masses, isotope mass-variance parameter Γ) |
| Fortran kernels | `debye.f90` (Debye integral tables, Callaway conductivity), interface in `debye.hpp` |
| Helpers | `fn_hust.hpp` (Hust thermal-conductivity form), `fn_material_data_path.hpp`; polynomial evaluation via `fn_polyval.hpp` from `src/linalg` |
| Executables | `main.cpp` → `material` (tables, lookup-table generation); `mattest.cpp` is a developer scratch driver |

## Related modules

- **Physics/Database** (`src/physics/database/`): the projected lookup tables behind ρ(T, B, β) and λ(T, B, β)
- **FEM Kernel** (`src/fem/kernel/`): where properties are evaluated during assembly; `cl_FEM_Calculator.cpp` decides from `depends()` whether the field-dependent overloads are used
- **Numerics/Bezier, Numerics/Spline** (`src/numerics/`): the curve types every fit is stored in
