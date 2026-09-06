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
 */

#include <belfem_user_api>

#include "user_table.hpp"
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

// -----------------------------------------------------------------------------
// jc and n take the THREE-dependency signature ( normB, angleNxB, T ), which is
// the only one that builds a JcFunctionUserDefined and populates mJcFunction /
// mNFunction. The one-dependency ( T ) overload stores the pointer but leaves
// mJcFunction null, so jc_eval() falls back to constant_property( jc ) -- which
// set_custom() has just set to NaN, i.e. the fit is loaded and unreachable.
//
// These fits are field- and angle-independent, so normB and angleNxB are unused.
// angleNxB is the UNFOLDED field-to-normal angle in [0, pi].
// -----------------------------------------------------------------------------

// Critical temperature of THIS fit, in K.
//
// It is not a free parameter: the jc quartic has a root at T = 90.0371 K and is
// NEGATIVE above it, and the n quintic falls through 1 at T = 91.8784 K and
// through 0 at T = 92.0791 K. Both are hard failures for the piecewise power
// law, which asserts jc > 0 and n > 1. The YBCO builtin's 92.5 K sits past both
// roots -- using it here leaves a live superconducting window in which this
// plugin's own fits are already invalid.
//
// The same number MUST be the cutoff in hts_jc / hts_n below and the value
// handed to MaterialProperty::T_crit in hts_init(), which is why it is written
// once. rho_piecewise() gates on T > T_crit, so the cutoff is inclusive: at
// exactly T = T_crit the solver still takes the superconducting branch, and it
// must find the polynomials there ( jc = 5.86e7 A/m^2, n = 9.48 ) rather than
// the normal-state fallbacks ( jc = 1 A/m^2, n = 1, which asserts ).
const real hts_Tcrit = 90.0 ;

real hts_jc(const Material* mat, real normB, real angleNxB, real T)
{
    
    real Jc0 = 1.6171*1*4.66e10 ;
    real jc = 1.0 ;
    
    if ( T <= hts_Tcrit )
    {
        jc = Jc0*(0.00004016*std::pow(T,4.0)-0.0089*std::pow(T,3.0)+0.7256*std::pow(T,2.0)-35.7*T+1189)/202.7 ;
    }
    
    return jc ;
    
}

real hts_n(const Material* mat, real normB, real angleNxB, real T)
{

    real n = 1.0 ;
    
    if ( T <= hts_Tcrit )
    {
        n = (-2.99e-7*std::pow(T,5.0)+7.546e-5*std::pow(T,4.0)-7.844e-3*std::pow(T,3.0)+0.4259*std::pow(T,2.0) -12.12*T+183.4) ;
    }
    
    return n ;
    
}

real hts_rhon(const Material* mat, real T)
{
    
    return 1e-6 + 4.7e-9*(T-90.0) ;
    
}


real cu_resistivity(const Material* mat, real T)
{
    
    if ( T >= 1.0 && T < 7.0 )
    {
        return 2.0e-11 ;
    }
    if ( T < 40.0 )
    {
        return 1.002557e-11 + T*2.945503e-12 - T*T*2.767806e-13 + T*T*T*8.665115e-15 ;
    }
    if ( T < 100.0 )
    {
        return 1.370786e-9 - T*8.741734e-11 + T*T*1.738251e-12 - T*T*T*6.532611e-15 ;
    }
    if ( T < 1358.0 )
    {
        return -3.514582e-9 + T*7.064722e-11 - T*T*8.917638e-15 - T*T*T*1.026538e-17 ;
    }
    
    return -3.514582e-9 + 1358.0*7.064722e-11 - 1358.0*1358.0*8.917638e-15 - 1358.0*1358.0*1358.0*1.026538e-17 ;

}

real ag_resistivity(const Material* mat, real T)
{
    static usermat::Table tTable( "rho_Ag.txt" ) ;
    return tTable( T ) ;
}

real hast_resistivity(const Material* mat, real T)
{
    static usermat::Table tTable( "rho_Hast.txt" ) ;
    return tTable( T ) ;
}

real rint_resistivity(const Material* mat, real T)
{

    return 1e-3 ;
}

/**
 * @brief Custom specific heat function
 * @param mat Pointer to material
 * @param T Temperature [K]
 * @return Specific heat capacity [J/(kg·K)]
 */
 
real cu_cp(const Material* mat, real T)
{
    static usermat::Table tTable( "cp_Cu.txt" ) ;
    return tTable( T ) ;
}

real ag_cp(const Material* mat, real T)
{
    static usermat::Table tTable( "cp_Ag.txt" ) ;
    return tTable( T ) ;
}

real hast_cp(const Material* mat, real T)
{
    static usermat::Table tTable( "cp_Hast.txt" ) ;
    return tTable( T ) ;
}

real rint_cp(const Material* mat, real T)
{
    static usermat::Table tTable( "cp_HTS.txt" ) ;
    return tTable( T ) ;
}

/**
 * @brief Custom thermal conductivity function
 * @param mat Pointer to material
 * @param T Temperature [K]
 * @return Thermal conductivity [W/(m·K)]
 */
real cu_k(const Material* mat, real T)
{
    static usermat::Table tTable( "k_Cu.txt" ) ;
    return tTable( T ) ;
}

real ag_k(const Material* mat, real T)
{
    static usermat::Table tTable( "k_Ag.txt" ) ;
    return tTable( T ) ;
}

real hast_k(const Material* mat, real T)
{
    static usermat::Table tTable( "k_Hast.txt" ) ;
    return tTable( T ) ;
}

real rint_k(const Material* mat, real T)
{
    static usermat::Table tTable( "k_HTS.txt" ) ;
    return tTable( T ) ;
}

// =============================================================================
// MATERIAL INITIALIZATION FUNCTION
// =============================================================================
// This function is called when the material is loaded. The function name
// must be: extern "C" void <MaterialName>_init(Material* mat)
// where <MaterialName> matches the second argument in create_material()

extern "C" void cu_init(Material* mat)
{
    // -------------------------------------------------------------------------
    // Set constant material properties
    // -------------------------------------------------------------------------

    // Density (required for thermal mass matrix)
    mat->set_constant(MaterialProperty::ref_density, 8700.0);   // Density [kg/m³]
    mat->set_constant(MaterialProperty::T_ref_density, 293.15); // Reference temperature [K]


    // -------------------------------------------------------------------------
    // Set temperature-dependent properties using custom functions
    // -------------------------------------------------------------------------

    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &cu_resistivity);

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &cu_cp);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &cu_k);

}

extern "C" void ag_init(Material* mat)
{
    // -------------------------------------------------------------------------
    // Set constant material properties
    // -------------------------------------------------------------------------

    // Density (required for thermal mass matrix)
    mat->set_constant(MaterialProperty::ref_density, 10490.0);   // Density [kg/m³]
    mat->set_constant(MaterialProperty::T_ref_density, 293.15); // Reference temperature [K]


    // -------------------------------------------------------------------------
    // Set temperature-dependent properties using custom functions
    // -------------------------------------------------------------------------

    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &ag_resistivity);

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &ag_cp);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &ag_k);

}

extern "C" void hast_init(Material* mat)
{
    // -------------------------------------------------------------------------
    // Set constant material properties
    // -------------------------------------------------------------------------

    // Density (required for thermal mass matrix)
    mat->set_constant(MaterialProperty::ref_density, 8890.0);   // Density [kg/m³]
    mat->set_constant(MaterialProperty::T_ref_density, 293.15); // Reference temperature [K]


    // -------------------------------------------------------------------------
    // Set temperature-dependent properties using custom functions
    // -------------------------------------------------------------------------

    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &hast_resistivity);

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &hast_cp);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &hast_k);

}

extern "C" void rint_init(Material* mat)
{
    // -------------------------------------------------------------------------
    // Set constant material properties
    // -------------------------------------------------------------------------

    // Density (required for thermal mass matrix)
    mat->set_constant(MaterialProperty::ref_density, 6390.0);   // Density [kg/m³]
    mat->set_constant(MaterialProperty::T_ref_density, 293.15); // Reference temperature [K]

    // -------------------------------------------------------------------------
    // Set temperature-dependent properties using custom functions
    // -------------------------------------------------------------------------

    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &rint_resistivity);

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &rint_cp);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &rint_k);

}

extern "C" void hts_init(Material* mat)
{
    // -------------------------------------------------------------------------
    // Set constant material properties
    // -------------------------------------------------------------------------

    // Density (required for thermal mass matrix)
    mat->set_constant(MaterialProperty::ref_density, 6390.0);   // Density [kg/m³]
    mat->set_constant(MaterialProperty::T_ref_density, 293.15); // Reference temperature [K]

    // -------------------------------------------------------------------------
    // Set temperature-dependent properties using custom functions
    // -------------------------------------------------------------------------

    // T_crit is REQUIRED by the piecewise power law: rho_piecewise() reads it
    // through constant_property(), which asserts in a debug build and returns
    // NaN in a release build when it was never set. The value is the cutoff of
    // this plugin's own fits -- see hts_Tcrit above for why it is not 92.5 K.
    mat->set_constant(MaterialProperty::T_crit, hts_Tcrit);

    mat->set_user_defined_function(MaterialProperty::jc,
                                    MaterialDependency::normB,
                                    MaterialDependency::angleNxB,
                                    MaterialDependency::T,
                                    &hts_jc);

    mat->set_user_defined_function(MaterialProperty::n,
                                    MaterialDependency::normB,
                                    MaterialDependency::angleNxB,
                                    MaterialDependency::T,
                                    &hts_n);
                                    
    mat->set_user_defined_function(MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &hts_rhon);
                                    
    mat->set_constant(MaterialProperty::ec, 1e-4); 

    mat->set_user_defined_function(MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &rint_cp);

    mat->set_user_defined_function(MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &rint_k);

}

// =============================================================================
// NOTES
// =============================================================================
//
// 1. FUNCTION NAMING:
//    - The init function must be extern "C" to prevent name mangling
//    - Name pattern: <MaterialName>_init where MaterialName matches usage
//    - Example: MyAlloy_init for create_material("libmyalloy.so", "MyAlloy")
//
// 2. AVAILABLE PROPERTIES (for set_constant or set_user_defined_function):
//    - Mechanical: E, nu, alpha, Rp02
//    - Thermal: cp, lambda
//    - Electrical: rho
//    - Physical: ref_density, T_ref_density
//
// 3. DEPENDENCIES:
//    - MaterialDependency::T         - Temperature [K]
//    - MaterialDependency::normB     - Magnetic field magnitude [T]
//    - MaterialDependency::angleBxJ  - Angle between B and J [rad]
//    - MaterialDependency::angleNxB  - Angle between normal and B [rad]
//    - MaterialDependency::normH     - Magnetic field strength [A/m]
//
// 4. MULTI-PARAMETER FUNCTIONS:
//    See README.md for examples of mu(H,T), rho(T,B,angle), jc(B,angle,T)
//
// 5. ACCESSING OTHER PROPERTIES:
//    Within your functions, you can call mat->property(args) to access
//    other material properties, enabling coupled property definitions.
//
// =============================================================================
