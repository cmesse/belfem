/*
 * Example User-Defined Material for BELFEM
 *
 * This file demonstrates how to create a custom material that can be
 * dynamically loaded by BELFEM. Copy and modify this template for your
 * own materials.
 *
 * Compilation:
 *   1. Copy UserMaterialTemplate.cmake to your directory as CMakeLists.txt
 *   2. Edit CMakeLists.txt to set BELFEM_DIR and material name
 *   3. mkdir build && cd build && cmake .. && make
 *
 * Usage in BELFEM:
 *   MaterialFactory factory;
 *   Material* mat = factory.create_material("./libmyalloy.so", "MyAlloy");
 *
 * A complete, running plugin built on this template ships as
 * examples/tape_quench_usermat/src/matlib.cpp -- five materials, measured
 * tables, and a superconductor. Read it when this file runs out of answers.
 */

#include <belfem_user_api>

using namespace belfem;

// =============================================================================
// MATERIAL PROPERTY FUNCTIONS
// =============================================================================
// Define your material property functions here. The first parameter must
// always be "const Material* mat" to allow access to other properties.

/**
 * @brief Custom electrical resistivity function
 * @param mat Pointer to material (allows access to other properties)
 * @param T Temperature [K]
 * @return Electrical resistivity [Ω·m]
 */
real my_resistivity([[maybe_unused]] const Material* mat, real T)
{
    // Example: Linear temperature dependence
    // rho(T) = rho0 * (1 + alpha * (T - T0))

    real rho0 = 1.7e-8;     // Resistivity at reference temperature [Ω·m]
    real alpha = 0.004;     // Temperature coefficient [1/K]
    real T0 = 293.15;       // Reference temperature [K]

    return rho0 * (1.0 + alpha * (T - T0));
}

/**
 * @brief Custom specific heat function
 * @param mat Pointer to material
 * @param T Temperature [K]
 * @return Specific heat capacity [J/(kg·K)]
 */
real my_specific_heat([[maybe_unused]] const Material* mat, real T)
{
    // Example: Temperature-dependent specific heat
    return 385.0 + 0.12 * T;
}

/**
 * @brief Custom thermal conductivity function
 * @param mat Pointer to material
 * @param T Temperature [K]
 * @return Thermal conductivity [W/(m·K)]
 */
real my_thermal_conductivity(const Material* mat, real T)
{
    // Example: You can access other material properties
    real rho = mat->rho(T);  // Get resistivity at this temperature

    // Wiedemann-Franz law: lambda = L * T / rho
    // L = Lorenz number = 2.44e-8 W·Ω/K²
    real L = 2.44e-8;
    return L * T / rho;
}

// =============================================================================
// MATERIAL INITIALIZATION FUNCTION
// =============================================================================
// This function is called when the material is loaded. The function name
// must be: extern "C" void <MaterialName>_init(Material* mat)
// where <MaterialName> matches the second argument in create_material()

extern "C" void MyAlloy_init(Material* mat)
{
    // -------------------------------------------------------------------------
    // Set constant material properties
    // -------------------------------------------------------------------------

    // Mechanical properties
    mat->set_constant(MaterialProperty::E, 200e9);           // Young's modulus [Pa]
    mat->set_constant(MaterialProperty::nu, 0.3);            // Poisson's ratio [-]

    // Density (required for thermal mass matrix)
    mat->set_constant(MaterialProperty::ref_density, 8900.0);   // Density [kg/m³]
    mat->set_constant(MaterialProperty::T_ref_density, 293.15); // Reference temperature [K]

    // Other constant properties (optional)
    // mat->set_constant(MaterialProperty::alpha, 1.2e-5);   // Thermal expansion [1/K]

    // -------------------------------------------------------------------------
    // Set temperature-dependent properties using custom functions
    // -------------------------------------------------------------------------

    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &my_resistivity);

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &my_specific_heat);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &my_thermal_conductivity);

    // -------------------------------------------------------------------------
    // Alternative: Use polynomials for simple temperature dependencies
    // -------------------------------------------------------------------------
    // Coefficients in DESCENDING order (MATLAB style): [c₀, c₁, c₂, ...]
    // f(T) = c₀*T^n + c₁*T^(n-1) + ... + c_n
    // std::vector keeps this material independent of the linalg backend.

    // Example: alpha(T) = 1.2e-5 + 3e-9*T  (coefficients: [3e-9, 1.2e-5])
    // std::vector<real> alpha_coeffs = {3e-9, 1.2e-5};
    // mat->set_user_defined_polynomial(MaterialProperty::alpha, alpha_coeffs);
}

// =============================================================================
// A SUPERCONDUCTOR: jc AND n TAKE A DIFFERENT ROUTE
// =============================================================================
// jc and n are NOT ordinary properties. They are evaluated through a
// JcFunction object, and only the two- and three-dependency overloads of
// set_user_defined_function build one. See note 6 below for the trap.
//
// The three-dependency signature is fixed: ( normB, angleNxB, T ), in that
// order. angleNxB is the UNFOLDED field-to-normal angle in [0, pi].

real my_jc([[maybe_unused]] const Material* mat,
           [[maybe_unused]] real normB,
           [[maybe_unused]] real angleNxB,
                            real T)
{
    // The fit's ROOT (Tc0) sits ABOVE the T_crit registered below. That gap is
    // the whole point: the solver's normal-state gate is `T > T_crit`, so at
    // exactly T == T_crit it still takes the SUPERCONDUCTING branch and must
    // find a strictly positive jc there. A fit whose root coincides with
    // T_crit returns 0 at that temperature -- an assert in a debug build, and
    // a division by zero in rho_powerlaw(), which has no T_crit gate at all.
    const real Tc0 = 92.5;    // where THIS fit reaches zero
    const real Jc0 = 3.0e10;

    return T <= Tc0 ? Jc0 * (1.0 - T / Tc0) : 1.0;
}

real my_n([[maybe_unused]] const Material* mat,
          [[maybe_unused]] real normB,
          [[maybe_unused]] real angleNxB,
                           real T)
{
    // Same reasoning as my_jc: n must be STRICTLY above 1 everywhere the
    // superconducting branch is taken, T == T_crit included. Falling through 1
    // exactly at the cutoff is the failure this offset avoids.
    const real Tc0 = 92.5;

    return T <= Tc0 ? 1.0 + 24.0 * (1.0 - T / Tc0) : 1.0;
}

extern "C" void MySuperconductor_init(Material* mat)
{
    mat->set_constant(MaterialProperty::ref_density, 6390.0);
    mat->set_constant(MaterialProperty::T_ref_density, 293.15);

    // ec is required by every E-J law reachable from here: the two- and
    // three-dependency jc registration bypasses the default-ec path, so
    // constant_property( ec ) asserts in debug and yields NaN in release when
    // it was never set.
    //
    // T_crit is required by the Piecewise and Riva laws, which gate on it. The
    // plain PowerLaw never reads it. Omitting it is NOT a loud failure: the
    // default is a quiet NaN, and `T > NaN` is FALSE, so a piecewise material
    // silently stays on its superconducting branch at every temperature.
    //
    // T_crit must sit strictly INSIDE the validity window of your own fits --
    // below the root of jc and below where n falls through 1 -- because the
    // cutoff is inclusive. See my_jc above.
    mat->set_constant(MaterialProperty::T_crit, 90.0);   // < Tc0 = 92.5 above
    mat->set_constant(MaterialProperty::ec, 1e-4);       // critical field [V/m]

    mat->set_user_defined_function(MaterialProperty::jc,
                                    MaterialDependency::normB,
                                    MaterialDependency::angleNxB,
                                    MaterialDependency::T,
                                    &my_jc);

    mat->set_user_defined_function(MaterialProperty::n,
                                    MaterialDependency::normB,
                                    MaterialDependency::angleNxB,
                                    MaterialDependency::T,
                                    &my_n);

    // The normal-state resistivity the power law falls back to above T_crit.
    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &my_resistivity);

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &my_specific_heat);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &my_thermal_conductivity);
}

// =============================================================================
// NOTES
// =============================================================================
//
// 1. FUNCTION NAMING:
//    - The init function must be extern "C" to prevent name mangling
//    - Name pattern: <MaterialName>_init where MaterialName matches usage
//    - Example: MyAlloy_init for create_material("libmyalloy.so", "MyAlloy")
//    - One library may hold many init functions; the label in the deck picks
//      which one runs
//
// 2. WHERE THE LIBRARY IS LOOKED FOR:
//    The path from the deck is resolved in this order --
//      1. the path as written, relative to the run directory, or absolute
//      2. the same relative path below <data root>/material
//      3. for a path CONTAINING a slash, the file name alone below
//         <data root>/material
//      4. otherwise the name is handed to dlopen unchanged
//    The root is gBelfemDataPath, not $BELFEM_DATA directly. Note that dlopen
//    searches $LD_LIBRARY_PATH only for a name with NO slash, so adding a
//    directory component DISABLES the loader search rather than enabling it.
//    The copy in your build directory works as-is.
//
// 3. THE ONE INCLUDE:
//    <belfem_user_api> is the umbrella for the whole plugin API. Do not
//    include the individual BELFEM headers -- the umbrella is the supported
//    surface, and it is deliberately free of the linear-algebra backend, so
//    your plugin does not have to match the host's Armadillo/Blaze choice.
//
// 4. AVAILABLE PROPERTIES (for set_constant or set_user_defined_function):
//    - Mechanical:    E, nu, alpha, Rp02
//    - Thermal:       cp, lambda, debye, gamma, beta, T_max
//    - Electrical:    rho, rho_0, rho_i, RRR, kohler_trans, kohler_long
//    - Superconductor: jc, n, ec, T_crit, layer_thickness
//    - Magnetic:      mu, Tcurie
//    - Physical:      ref_density, T_ref_density, M, R
//    The complete enumeration is MaterialProperty in cl_Material.hpp.
//
// 5. DEPENDENCIES:
//    - MaterialDependency::T         - Temperature [K]
//    - MaterialDependency::normB     - Magnetic flux density magnitude [T]
//    - MaterialDependency::angleBxJ  - Angle between B and J [rad]
//    - MaterialDependency::angleNxB  - Angle between normal and B [rad]
//    - MaterialDependency::normH     - Magnetic field strength [A/m]
//
// 6. THE ARGUMENT ORDER IS FIXED, AND SO IS THE OVERLOAD:
//    Each multi-parameter property accepts exactly one dependency order, and
//    BELFEM_ERRORs on any other:
//
//      mu     ( normH, T )                    two dependencies
//      jc, n  ( normB, angleNxB, T )          three dependencies
//      rho    ( T, normB, angleBxJ )          three dependencies  [see below]
//      lambda ( T, normB, angleBxJ )          three dependencies  [see below]
//
//    CAUTION on the three-dependency rho and lambda: the registration is
//    accepted and the order guards are real, but the path does not currently
//    work. The registration installs T, normB and angleBxJ and then calls
//    set_custom(), which RESETS the dependency set and restores T alone. Treat
//    those two as temperature-only until that is fixed; jc, n and mu are
//    unaffected.
//
//    A two-dependency jc or n is also accepted and gives a field- and
//    angle-dependent fit with no temperature argument. Be aware that this
//    overload does NOT validate the dependencies you declare -- it forces
//    normB and angleNxB regardless of what you pass.
//
//    Note that rho and lambda are T-FIRST while jc and n are T-LAST. All the
//    arguments are `real`, so getting this wrong compiles cleanly and
//    evaluates the property at the wrong point -- this has bitten BELFEM
//    itself more than once.
//
//    THE TRAP: the ONE-dependency overload accepts jc and n (with
//    MaterialDependency::T -- it rejects any other dependency) and does the
//    wrong thing. It stores the pointer and calls set_custom(), but never
//    builds the JcFunction, so jc_eval() falls back to constant_property( jc )
//    -- which set_custom() has just set to NaN.
//
//    Precisely: the fit is unreachable from the FEM ASSEMBLY route, which
//    always evaluates through the full-signature overloads. It IS still
//    reached by the temperature-only power-law overloads and by
//    postprocessing, so the material is inconsistent rather than inert.
//    In a debug build constant_property() asserts ("Property is not
//    constant"); in a release build it returns NaN and the jc > 0 test then
//    fails too. Always register jc and n through the two- or three-dependency
//    form, as MySuperconductor_init does.
//
// 7. ACCESSING OTHER PROPERTIES:
//    Within your functions, you can call mat->property(args) to access
//    other material properties, enabling coupled property definitions.
//
// 8. REBUILD AFTER A HOST REBUILD:
//    There is no version handshake at the dlopen boundary. A plugin built
//    against older headers loads without complaint and calls into whatever the
//    host's layout is now. Rebuild the plugin whenever BELFEM's headers change,
//    and keep NDEBUG/DEBUG in agreement with the host.
//
// =============================================================================
