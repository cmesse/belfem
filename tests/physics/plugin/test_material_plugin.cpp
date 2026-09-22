/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

//------------------------------------------------------------------------------
// Fixture for the two-dependency jc and n registration guards.
//
// Each label below is one registration a plugin author could write. The
// accepted one declares ( normB, angleNxB ); the others declare a pair the
// host must refuse, so their _init raises out of the UserDefinedMaterial
// constructor. The build places this at
//
//     <build>/tests/physics/plugindata/material/libtest_material_plugin.so
//
// and the test passes that absolute path to the constructor, so no data-path
// search is involved.
//
// It links nothing: every call below is the virtual set_user_defined_function
// on the material the executable hands in, so the plugin needs no BELFEM
// symbol of its own and the executable's -rdynamic link is not what this
// fixture exercises.
//------------------------------------------------------------------------------

#include "cl_Material.hpp"

using belfem::Material ;
using belfem::MaterialDependency ;
using belfem::MaterialProperty ;
using belfem::real ;

namespace
{
    //! The two sentinel scales differ by three decades, so a transposed
    //! delivery is a different number rather than a sign:
    //! jc( 2.5, 0.75 ) = 2500750000 against jc( 0.75, 2.5 ) = 752500000
    real
    two_arg_jc( [[maybe_unused]] const Material * aMaterial,
                const real aNormB,
                const real aAngleNxB )
    {
        return 1.0e9 * aNormB + 1.0e6 * aAngleNxB ;
    }

    //! n( 2.5, 0.75 ) = 22 against n( 0.75, 2.5 ) = 11.5
    real
    two_arg_n( [[maybe_unused]] const Material * aMaterial,
               const real aNormB,
               const real aAngleNxB )
    {
        return 20.0 + 2.0 * aNormB - 4.0 * aAngleNxB ;
    }
}

//------------------------------------------------------------------------------
// the one accepted declaration
//------------------------------------------------------------------------------

extern "C" void
TwoArgHts_init( Material * aMaterial )
{
    aMaterial->set_user_defined_function( MaterialProperty::jc,
                                          MaterialDependency::normB,
                                          MaterialDependency::angleNxB,
                                          & two_arg_jc );

    aMaterial->set_user_defined_function( MaterialProperty::n,
                                          MaterialDependency::normB,
                                          MaterialDependency::angleNxB,
                                          & two_arg_n );
}

//------------------------------------------------------------------------------
// declarations the host must refuse
//------------------------------------------------------------------------------

//! jc with the two arguments the other way round
extern "C" void
ReversedJc_init( Material * aMaterial )
{
    aMaterial->set_user_defined_function( MaterialProperty::jc,
                                          MaterialDependency::angleNxB,
                                          MaterialDependency::normB,
                                          & two_arg_jc );
}

//! jc declared against a dependency the two-argument path never delivers
extern "C" void
UnrelatedJc_init( Material * aMaterial )
{
    aMaterial->set_user_defined_function( MaterialProperty::jc,
                                          MaterialDependency::T,
                                          MaterialDependency::normB,
                                          & two_arg_jc );
}

//! jc with the field first and the second argument wrong
extern "C" void
SecondArgJc_init( Material * aMaterial )
{
    aMaterial->set_user_defined_function( MaterialProperty::jc,
                                          MaterialDependency::normB,
                                          MaterialDependency::T,
                                          & two_arg_jc );
}

//! n with the two arguments the other way round
extern "C" void
ReversedN_init( Material * aMaterial )
{
    aMaterial->set_user_defined_function( MaterialProperty::n,
                                          MaterialDependency::angleNxB,
                                          MaterialDependency::normB,
                                          & two_arg_n );
}

//! n with the field first and the second argument wrong
extern "C" void
WrongSecondN_init( Material * aMaterial )
{
    aMaterial->set_user_defined_function( MaterialProperty::n,
                                          MaterialDependency::normB,
                                          MaterialDependency::T,
                                          & two_arg_n );
}

//------------------------------------------------------------------------------
