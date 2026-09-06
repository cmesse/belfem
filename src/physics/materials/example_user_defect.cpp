/*
 * Example User-Defined Defect for BELFEM
 *
 * This file demonstrates how to create a custom defect that can be
 * dynamically loaded by BELFEM. Copy and modify this template for your
 * own defects.
 *
 * A defect is a multiplier on the local critical current density: it returns
 * 1 where the conductor is pristine and something below 1 where it is not.
 * The signature is fixed at ( x, y, z, t ), so a defect may be switched on in
 * time as well as placed in space.
 *
 * Compilation:
 *   1. Copy UserMaterialTemplate.cmake to your directory as CMakeLists.txt
 *      (a defect alone needs only the material API. If you ship a defect and
 *      a source function in ONE library -- as examples/tape_quench_usermat
 *      does -- use UserLibraryTemplate.cmake instead)
 *   2. Edit CMakeLists.txt to set BELFEM_DIR and the library name
 *   3. mkdir build && cd build && cmake .. && make
 *
 * A complete, running defect built on this template ships as
 * examples/tape_quench_usermat/src/defect.cpp, which also shows how to size
 * one so that it actually propagates a quench.
 */

#include <belfem_user_api>

using namespace belfem;

// =============================================================================
// User defined defect functions
// =============================================================================
// Define your defect function here. The signature is fixed by BELFEM:
// real ( const real x, const real y, const real z, const real t ).
// Ignoring an argument you do not need is normal and not a defect in itself.

/**
 * @brief Custom defect
 * @param x,y,z spatial coordinates (m)
 * @param t time (s)
 * @return multiplier on jc: 1 = pristine, 0 = fully blocked
 */
real my_defect(const real x,
               [[maybe_unused]] const real y,
               const real z,
               [[maybe_unused]] const real t)
{
    real factor = 0.99 ;
    real dx = 0.5e-3 ;
    real dz = 0.5e-3 ;
    real x0 = 0e-3 ;
    real z0 = 0e-3 ;

    return 1-factor*exp(-((x-x0)*(x-x0)/(2*dx*dx)) - ((z-z0)*(z-z0)/(2*dz*dz))) ;
}

// =============================================================================
// DEFECT INITIALIZATION FUNCTION
// =============================================================================
// This function is called when the defect is loaded. The function name
// must be: extern "C" void <DefectName>_init(Material* mat)
// where <DefectName> matches the label passed to Material::read_defect().
//
// The library path given to read_defect() is resolved in this order:
//   1. the path as written, relative to the run directory, or absolute
//   2. the same relative path below <data root>/material
//   3. for a path CONTAINING a slash, the file name alone below that directory
//   4. otherwise the name is handed to dlopen unchanged
// The root is gBelfemDataPath, not $BELFEM_DATA directly. dlopen searches
// $LD_LIBRARY_PATH only for a slashless name, so adding a directory component
// DISABLES the loader search rather than enabling it.

extern "C" void MyDefect_init(Material* mat)
{

    mat->set_user_defined_defect(&my_defect);

}
