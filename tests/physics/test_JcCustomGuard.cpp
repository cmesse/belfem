/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * jc and n have no T-only custom callback.
 *
 * Material::set_custom( jc ) and set_custom( n ) once wrote their callbacks
 * into the rho_i and debye dispatch slots, while the power law evaluates jc
 * and n through a JcFunction and never read them. Both cases now abort, so
 * the resistivity and Debye routing of a material cannot be hijacked.
 */

#include <string>

#include <gtest/gtest.h>

#include "cl_Material.hpp"
#include "cl_Material_Copper.hpp"
#include "cl_Material_UserDefined.hpp"

TEST( JcCustomGuard, SetCustomRejectsJcAndN )
{
    // a built-in metal: tables off, analytic Kohler path
    belfem::material::Copper tCopper( 50.0, false );

    EXPECT_THROW( tCopper.set_custom( belfem::MaterialProperty::jc ), std::exception );
    EXPECT_THROW( tCopper.set_custom( belfem::MaterialProperty::n ),  std::exception );

    // the guard fires before the have-flag or any dispatch slot is touched
    EXPECT_FALSE( tCopper.have( belfem::MaterialProperty::jc ) );
    EXPECT_FALSE( tCopper.have( belfem::MaterialProperty::n ) );
    EXPECT_GT( tCopper.rho_i( 77.0 ), 0.0 );
    EXPECT_GT( tCopper.debye( 77.0 ), 0.0 );
}

//------------------------------------------------------------------------------
// The two-dependency jc and n registration honors the pair it is given.
//
// set_user_defined_function( Property, Dependency1, Dependency2, MatFunc2* )
// always calls the installed function as f( material, normB, angleNxB ), so
// ( normB, angleNxB ) is the only declaration it can accept; any other pair
// would transpose the author's arguments silently. The cases below run
// through a dlopen-loaded fixture because that is the only surface a plugin
// author has. The refusal cases were accepted before the guard landed.
//------------------------------------------------------------------------------

#ifdef BELFEM_TEST_PLUGIN_ROOT

namespace
{
    //! the fixture built by tests/physics/CMakeLists.txt; the constructor
    //! takes this absolute path straight to dlopen, so no data-path search
    //! is involved
    const std::string tPluginPath =
            std::string( BELFEM_TEST_PLUGIN_ROOT ) + "/material/libtest_material_plugin.so" ;

    // runs aCall, requires a std::runtime_error whose what() carries the
    // literal part of the BELFEM_ERROR text
    template< typename F >
    void
    expect_refusal( F && aCall, const std::string & aSubstring )
    {
        bool tThrown = false ;
        try
        {
            aCall();
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;
            EXPECT_NE( std::string( aError.what() ).find( aSubstring ), std::string::npos )
                << "message was: " << aError.what();
        }
        EXPECT_TRUE( tThrown ) << "expected a refusal carrying \"" << aSubstring << "\"";
    }

    // builds the fixture material of a label; the refusal cases raise out of
    // the constructor
    void
    load_plugin_material( const char * aLabel )
    {
        const belfem::material::UserDefinedMaterial tMaterial( tPluginPath, aLabel );
    }
}

TEST( JcCustomGuard, TwoArgumentJcAndNHonorTheDeclaredFieldAndAngle )
{
    const belfem::material::UserDefinedMaterial tMaterial( tPluginPath, "TwoArgHts" );

    EXPECT_TRUE( tMaterial.have( belfem::MaterialProperty::jc ) );
    EXPECT_TRUE( tMaterial.have( belfem::MaterialProperty::n ) );

    // the flags the assembly and post-processing routes read
    EXPECT_TRUE(  tMaterial.depends( belfem::MaterialProperty::jc, belfem::MaterialDependency::normB ) );
    EXPECT_TRUE(  tMaterial.depends( belfem::MaterialProperty::jc, belfem::MaterialDependency::angleNxB ) );
    EXPECT_FALSE( tMaterial.depends( belfem::MaterialProperty::jc, belfem::MaterialDependency::T ) );

    EXPECT_TRUE(  tMaterial.depends( belfem::MaterialProperty::n, belfem::MaterialDependency::normB ) );
    EXPECT_TRUE(  tMaterial.depends( belfem::MaterialProperty::n, belfem::MaterialDependency::angleNxB ) );
    EXPECT_FALSE( tMaterial.depends( belfem::MaterialProperty::n, belfem::MaterialDependency::T ) );

    // the sentinel scales of test_material_plugin.cpp: both values are exactly
    // representable, and a transposed delivery would give 752500000 and 11.5
    EXPECT_DOUBLE_EQ( tMaterial.jc( 2.5, 0.75 ), 2500750000.0 );
    EXPECT_DOUBLE_EQ( tMaterial.n( 2.5, 0.75 ), 22.0 );
}

TEST( JcCustomGuard, TwoArgumentJcAndNRejectReversedDependencies )
{
    expect_refusal( [](){ load_plugin_material( "ReversedJc" ); },
                    "first parameter of function jc" );

    expect_refusal( [](){ load_plugin_material( "UnrelatedJc" ); },
                    "must be magnetic field in T" );

    expect_refusal( [](){ load_plugin_material( "SecondArgJc" ); },
                    "second parameter of function jc" );

    expect_refusal( [](){ load_plugin_material( "ReversedN" ); },
                    "first parameter of function n" );

    expect_refusal( [](){ load_plugin_material( "WrongSecondN" ); },
                    "must be an angle in rad" );
}

#endif
