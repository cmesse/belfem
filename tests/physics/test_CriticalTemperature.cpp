/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * The 'critical temperature' material key. It overrides T_crit on a
 * user-defined superconductor and is refused everywhere else:
 *
 *  - on a builtin material, because a builtin's T_crit is part of a
 *    calibration already consumed in its constructor ( YBCO samples its
 *    lambda spline there ), so moving it afterwards would desynchronize the
 *    resistivity gate from the thermal conductivity
 *  - on a user-defined material that registers no jc ( untested, below )
 *  - with a value that is not finite and positive
 *
 * The superconductor cases load the plugin fixture built by
 * tests/physics/CMakeLists.txt, whose TwoArgHts label registers jc and n.
 * The no-jc refusal is not reachable here: every label of that fixture
 * registers jc or n, and the n-only ones are refused by their own
 * registration guard before the factory reads the key.
 *
 * Every factory here is built with tables off. On a refused deck whatever
 * was built before the throw leaks ( the factory does not own its
 * materials, and a throwing constructor runs no destructor ); one material
 * per refused deck bounds that to a single case.
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_InputFile.hpp"
#include "cl_Material.hpp"
#include "cl_MaterialFactory.hpp"

using belfem::real ;
using belfem::Material ;
using belfem::MaterialFactory ;
using belfem::MaterialProperty ;

namespace
{
    // a 'materials' deck in the test temp directory under an absolute path,
    // so the case does not depend on the working directory; removed on
    // scope exit
    class Deck
    {
        const std::string mPath ;

    public:

        Deck( const std::string & aName, const std::string & aBody ) :
            mPath( ::testing::TempDir() + "tcrit_" + aName + ".conf" )
        {
            std::ofstream tFile( mPath );
            tFile << "materials\n{\n" << aBody << "\n}\n";
        }

        Deck( const Deck & ) = delete ;
        Deck( Deck && ) = delete ;
        Deck & operator=( const Deck & ) = delete ;
        Deck & operator=( Deck && ) = delete ;

        ~Deck()
        {
            std::remove( mPath.c_str() );
        }

        const std::string &
        path() const
        {
            return mPath ;
        }
    };

    void
    delete_materials( MaterialFactory & aFactory )
    {
        for ( auto & tEntry : aFactory.materials() )
        {
            delete tEntry.second ;
        }
    }

    // build the factory on a deck and hand back what its constructor threw;
    // an empty string means the deck was accepted. Only std::runtime_error
    // is caught, so any other exception fails the test
    std::string
    refusal( const std::string & aName, const std::string & aBody )
    {
        Deck tDeck( aName, aBody );
        belfem::InputFile tFile( tDeck.path() );

        try
        {
            MaterialFactory tFactory( tFile.section( "materials" ), false );
            delete_materials( tFactory );
            return "" ;
        }
        catch ( const std::runtime_error & aError )
        {
            return aError.what();
        }
    }

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    const char * gFiniteRefusal = "must be finite and positive" ;
}

//------------------------------------------------------------------------------

TEST( CriticalTemperature, BuiltinOverrideIsRefused )
{
    // a builtin superconductor, where the key would otherwise make sense
    std::string tMessage = refusal( "builtin_ybco",
        "    ybco\n    {\n        builtin : ybco ;\n        jc : 3e10 A/m^2 ;\n"
        "        n : 25 ;\n        critical temperature : 90 K ;\n    }" );
    EXPECT_TRUE( contains( tMessage, "cannot be used on the builtin material" ) ) << tMessage ;

    // a builtin metal: the builtin refusal is the one that names the cause
    tMessage = refusal( "builtin_copper",
        "    cu\n    {\n        builtin : copper ;\n        critical temperature : 90 K ;\n    }" );
    EXPECT_TRUE( contains( tMessage, "cannot be used on the builtin material" ) ) << tMessage ;
}

#ifdef BELFEM_TEST_PLUGIN_ROOT

namespace
{
    // a usermat superconductor section from the plugin fixture, with the
    // critical-temperature line aTcrit ( empty for none )
    std::string
    usermat_hts( const std::string & aTcrit )
    {
        return "    rebco\n    {\n" + aTcrit
             + "        usermat\n        {\n            file  : "
             + std::string( BELFEM_TEST_PLUGIN_ROOT )
             + "/material/libtest_material_plugin.so ;\n"
               "            label : TwoArgHts ;\n        }\n    }" ;
    }
}

//------------------------------------------------------------------------------

TEST( CriticalTemperature, NonFiniteOrNonPositiveIsRefused )
{
    std::string tMessage = refusal( "usermat_zero",
        usermat_hts( "        critical temperature : 0 K ;\n" ) );
    EXPECT_TRUE( contains( tMessage, gFiniteRefusal ) ) << tMessage ;

    tMessage = refusal( "usermat_negative",
        usermat_hts( "        critical temperature : -5 K ;\n" ) );
    EXPECT_TRUE( contains( tMessage, gFiniteRefusal ) ) << tMessage ;

    tMessage = refusal( "usermat_infinite",
        usermat_hts( "        critical temperature : inf K ;\n" ) );
    EXPECT_TRUE( contains( tMessage, gFiniteRefusal ) ) << tMessage ;
}

//------------------------------------------------------------------------------

// the shipped-deck shape ( examples/tape_quench_usermat ): keeps the cases
// above from degenerating into "the key is always refused"
TEST( CriticalTemperature, UsermatOverrideIsAccepted )
{
    Deck tDeck( "usermat_accepted", usermat_hts( "        critical temperature : 90.0 K ;\n" ) );
    belfem::InputFile tFile( tDeck.path() );

    MaterialFactory * tFactory = nullptr ;
    ASSERT_NO_THROW( tFactory = new MaterialFactory( tFile.section( "materials" ), false ) );

    Material * tMaterial = tFactory->materials()( "rebco" );
    ASSERT_NE( tMaterial, nullptr );
    EXPECT_EQ( tMaterial->constant_property( MaterialProperty::T_crit ), 90.0 );

    delete_materials( *tFactory );
    delete tFactory ;
}

#endif
