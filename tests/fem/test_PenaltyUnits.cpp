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
 * The unit rules of the three penalty slots, as far as they are reachable
 * without a Controller. The slots follow different rules on purpose:
 *
 *   eta    dimensionless, zero is the off state ( test_GhostSwitch.cpp )
 *   chi    dimensionless, zero is the off state
 *   k_reg  must carry a resistance unit, and zero is refused
 *
 * chi and k_reg are read only inside Controller::set_params, through
 * Section::get_value( key, "-" ) and get_value( key, "Ohm" ). What is pinned
 * here is the io-layer half of those reads: the unit check in both
 * directions, and that a zero passes it. The k_reg > 0 refusal and the
 * chi >= 0 refusal sit in the controller and are not reachable from here.
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <string>

#include "typedefs.hpp"
#include "cl_InputFile.hpp"

namespace
{
    // the solver block of a deck with one penalty section, in the test temp
    // directory under an absolute path, so the case does not depend on the
    // working directory; removed on scope exit
    class Deck
    {
        const std::string mPath ;

    public:

        Deck( const std::string & aName,
              const std::string & aBlock,
              const std::string & aBody ) :
            mPath( ::testing::TempDir() + aName )
        {
            std::ofstream tFile( mPath );
            tFile << "solver\n{\n    nonlinear magnetic\n    {\n        "
                  << aBlock << "\n        {\n" << aBody
                  << "\n        }\n    }\n}\n";
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

    const belfem::input::Section *
    penalty_section( const belfem::InputFile & aDeck, const std::string & aBlock )
    {
        return aDeck.section( "solver" )
                    ->section( "nonlinear magnetic" )
                    ->section( aBlock );
    }
}

//------------------------------------------------------------------------------

// a bare number is not a resistance: the opposite direction of the eta and
// chi rule, and the one a reader of "dimensioned values are refused" gets
// backwards
TEST( PenaltyUnits, BareKRegIsRefusedByTheUnitCheck )
{
    Deck tGuard( "penalty_kreg_bare.conf", "nitsche ghost penalty",
                 "            eta : 4 ;\n            k_reg : 1e-3 ;" );
    belfem::InputFile tDeck( tGuard.path() );

    EXPECT_ANY_THROW( penalty_section( tDeck, "nitsche ghost penalty" )
                          ->get_value( "k_reg", "Ohm" ) );
}

//------------------------------------------------------------------------------

TEST( PenaltyUnits, DimensionedKRegIsReadInSi )
{
    Deck tGuard( "penalty_kreg_milli.conf", "nitsche ghost penalty",
                 "            eta : 4 ;\n            k_reg : 2 mOhm ;" );
    belfem::InputFile tDeck( tGuard.path() );

    belfem::value tKreg ;
    ASSERT_NO_THROW( tKreg = penalty_section( tDeck, "nitsche ghost penalty" )
                                 ->get_value( "k_reg", "Ohm" ) );
    EXPECT_DOUBLE_EQ( tKreg.first, 2.0e-3 );
}

//------------------------------------------------------------------------------

// zero passes the unit layer; refusing it is the controller's job
// ( k_reg > 0, because the regularized harmonic mean of two zero layer
// stiffnesses is 0/0 )
TEST( PenaltyUnits, ZeroKRegPassesTheUnitLayer )
{
    Deck tGuard( "penalty_kreg_zero.conf", "nitsche ghost penalty",
                 "            eta : 4 ;\n            k_reg : 0 Ohm ;" );
    belfem::InputFile tDeck( tGuard.path() );

    belfem::value tKreg ;
    ASSERT_NO_THROW( tKreg = penalty_section( tDeck, "nitsche ghost penalty" )
                                 ->get_value( "k_reg", "Ohm" ) );
    EXPECT_EQ( tKreg.first, 0.0 );
}

//------------------------------------------------------------------------------

TEST( PenaltyUnits, DimensionedChiIsRefused )
{
    Deck tGuard( "penalty_chi_dim.conf", "coulomb gauge penalty",
                 "            chi : 1e-4 Ohm ;" );
    belfem::InputFile tDeck( tGuard.path() );

    EXPECT_ANY_THROW( penalty_section( tDeck, "coulomb gauge penalty" )
                          ->get_value( "chi", "-" ) );
}

//------------------------------------------------------------------------------

// zero is the off state of chi, as it is of eta
TEST( PenaltyUnits, ZeroChiIsAccepted )
{
    Deck tGuard( "penalty_chi_zero.conf", "coulomb gauge penalty",
                 "            chi : 0 ;" );
    belfem::InputFile tDeck( tGuard.path() );

    belfem::value tChi ;
    ASSERT_NO_THROW( tChi = penalty_section( tDeck, "coulomb gauge penalty" )
                                ->get_value( "chi", "-" ) );
    EXPECT_EQ( tChi.first, 0.0 );
}
