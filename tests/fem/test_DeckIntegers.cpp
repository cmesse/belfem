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
 * Section::get_int, the io-layer half of every integer deck key: it rounds
 * half away from zero and refuses a value an int cannot hold, and it is
 * deliberately permissive about the sign. The sign and budget rules of the
 * iteration keys ( "max iterations", "stall window", ... ) belong to the
 * controller, which validates each value before its uint cast.
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_InputFile.hpp"

namespace
{
    // a one-section deck in the test temp directory, under an absolute path
    // so the case does not depend on the working directory. The file is
    // removed on scope exit
    class Deck
    {
        const std::string mPath ;

    public:

        Deck( const std::string & aName, const std::string & aBody ) :
            mPath( ::testing::TempDir() + aName )
        {
            std::ofstream tFile( mPath );
            tFile << "solver\n{\n" << aBody << "\n}\n";
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
}

//------------------------------------------------------------------------------

TEST( DeckIntegers, GetIntRoundsHalfAwayFromZero )
{
    Deck tGuard( "deckint_round.conf", "    n : 2.6 ;\n    m : -2.6 ;\n    h : 2.5 ;\n    g : -2.5 ;" );
    belfem::InputFile tDeck( tGuard.path() );
    const belfem::input::Section * tSection = tDeck.section( "solver" );

    // std::round, not a truncating cast ( 2, -2, 2, -2 ) and not
    // floor( x + 0.5 ) ( which gives -2 for -2.5 )
    EXPECT_EQ( tSection->get_int( "n" ), 3 );
    EXPECT_EQ( tSection->get_int( "m" ), -3 );
    EXPECT_EQ( tSection->get_int( "h" ), 3 );
    EXPECT_EQ( tSection->get_int( "g" ), -3 );
}

//------------------------------------------------------------------------------

TEST( DeckIntegers, GetIntRefusesOutOfIntRange )
{
    // the bounds of this platform's int, and ten times beyond them
    const long long tMax = std::numeric_limits< int >::max() ;
    const long long tMin = std::numeric_limits< int >::min() ;
    Deck tGuard( "deckint_range.conf",
        "    big : " + std::to_string( 10 * tMax ) + " ;\n"
        "    small : " + std::to_string( 10 * tMin ) + " ;\n"
        "    edge : " + std::to_string( tMax ) + " ;\n"
        "    low : " + std::to_string( tMin ) + " ;" );
    belfem::InputFile tDeck( tGuard.path() );
    const belfem::input::Section * tSection = tDeck.section( "solver" );

    // the refusal names the key
    try
    {
        tSection->get_int( "big" );
        ADD_FAILURE() << "an out-of-range value was cast to int" ;
    }
    catch ( const std::runtime_error & aError )
    {
        EXPECT_NE( std::string( aError.what() ).find( "big" ), std::string::npos )
            << aError.what() ;
    }

    EXPECT_ANY_THROW( tSection->get_int( "small" ) );

    // the bounds themselves are representable and must pass
    EXPECT_EQ( tSection->get_int( "edge" ), std::numeric_limits< int >::max() );
    EXPECT_EQ( tSection->get_int( "low" ), std::numeric_limits< int >::min() );
}

//------------------------------------------------------------------------------

TEST( DeckIntegers, GetIntItselfAcceptsNegatives )
{
    Deck tGuard( "deckint_negative.conf", "    n : -1 ;" );
    belfem::InputFile tDeck( tGuard.path() );

    int tValue = 0 ;
    EXPECT_NO_THROW( tValue = tDeck.section( "solver" )->get_int( "n" ) );
    EXPECT_EQ( tValue, -1 );
}
