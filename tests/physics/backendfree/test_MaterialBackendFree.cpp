/*
 * Backend-free compile gate for the user-material API.
 *
 * This translation unit is compiled — never linked, never run — as part of the
 * test_material_backendfree OBJECT library, with the linear-algebra backend
 * macro deliberately stripped (see CMakeLists.txt in this directory). It pins
 * the contract established by the 2026-07-02 SplineLookupTable refactor:
 * cl_Material.hpp and the user-material surface must compile with no
 * BELFEM_ARMADILLO / BELFEM_BLAZE macro defined, so that a simple user
 * material never needs the backend.
 *
 * If this file fails to build, someone re-coupled cl_Material.hpp (or its
 * include closure) to the backend, or the build system started leaking the
 * backend define into this directory. Both are contract breaks, not problems
 * with this gate.
 */

// Guard 1: the build system must not inject the backend macro here. If this
// fires, the define stripping in CMakeLists.txt has been bypassed and the
// gate would otherwise pass vacuously.
#if defined( BELFEM_ARMADILLO ) || defined( BELFEM_BLAZE )
#error "backend define leaked into the backend-free material gate"
#endif

#include <vector>

#include "cl_Material.hpp"

// Guard 2: a backend-coupled header included without the backend macro
// degrades to forward declarations and would NOT fail the compile by itself.
// Its include guard still gets defined, and that is what we test here. The
// native wrapper headers (cl_BZ_*, cl_AR_*) include Blaze/Armadillo without
// any macro gate, so their guards are checked too.
// The umbrella is included HERE, before the guard block, so the guards below
// cover its whole closure ( cl_SourceFunction.hpp, cl_JcFunction*.hpp,
// cl_Material_UserDefined.hpp, powerlaws.hpp ) and not cl_Material.hpp alone.
#include <belfem_user_api>

#if defined( BELFEM_CL_VECTOR_HPP ) || defined( BELFEM_CL_MATRIX_HPP ) \
 || defined( BELFEM_SPLINE_HPP ) || defined( BELFEM_CL_SPMATRIX_HPP ) \
 || defined( BELFEM_CL_BZ_VECTOR_HPP ) || defined( BELFEM_CL_BZ_MATRIX_HPP ) \
 || defined( BELFEM_CL_AR_VECTOR_HPP ) || defined( BELFEM_CL_AR_MATRIX_HPP )
#error "cl_Material.hpp pulled a linalg/spline header into the backend-free material gate"
#endif

using namespace belfem;

// -----------------------------------------------------------------------------
// The user-material surface, as a plugin author writes it
// (mirrors example_user_material.cpp; external linkage on purpose — a static
// unreferenced function would trip -Wunused-function under -Wall -Werror)
// -----------------------------------------------------------------------------

real
gate_resistivity( const Material * aMaterial, const real aT )
{
    // linear temperature dependence, coefficients arbitrary
    return 1.7e-8 * ( 1.0 + 0.004 * ( aT - 293.15 ) );
}

real
gate_thermal_conductivity( const Material * aMaterial, const real aT )
{
    // Wiedemann-Franz form, exercises a property accessor on the base class
    return 2.44e-8 * aT / aMaterial->rho( aT );
}

extern "C" void
GateMaterial_init( Material * aMaterial )
{
    aMaterial->set_constant( MaterialProperty::E, 200e9 );
    aMaterial->set_constant( MaterialProperty::ref_density, 8900.0 );

    aMaterial->set_user_defined_function( MaterialProperty::rho,
                                          MaterialDependency::T,
                                          &gate_resistivity );

    aMaterial->set_user_defined_function( MaterialProperty::lambda,
                                          MaterialDependency::T,
                                          &gate_thermal_conductivity );

    // the backend-neutral polynomial overload (descending coefficients)
    std::vector< real > tCoefficients = { 3e-9, 1.2e-5 };
    aMaterial->set_user_defined_polynomial( MaterialProperty::alpha,
                                            tCoefficients );
}
