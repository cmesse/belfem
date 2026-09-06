/*
 * Example User-Defined Material for BELFEM
 *
 * This file demonstrates how to create a custom material that can be
 * dynamically loaded by BELFEM. Copy and modify this template for your
 * own materials.
 *
 * Compilation:
 *   1. Copy UserMaterialTemplate.cmake to your directory as CMakeLists.txt
 *   2. Edit CMakeLists.txt to set BELFEM_DIR and the library name
 *   3. mkdir build && cd build && cmake .. && make
 *
 * Usage in BELFEM:
 *   MaterialFactory factory;
 *   Material* mat = factory.create_material("./libmyalloy.so", "MyAlloy");
 */

#include <belfem_user_api>

using namespace belfem;

// =============================================================================
// MATERIAL PROPERTY FUNCTIONS
// =============================================================================
// Define your material property functions here. The first parameter must
// always be "const Material* mat" to allow access to other properties.

/**
 * @brief Plugin entry point: sets the constant properties of the bulk material
 * @param mat Pointer to the material being initialised
 */

extern "C" void bulk_init(Material* mat)
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

    mat->set_constant(MaterialProperty::rho, 1e-8);

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
