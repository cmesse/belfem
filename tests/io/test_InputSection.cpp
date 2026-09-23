/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * The two id-list readers of an input section, and where they disagree.
 *
 * get_id_groups keeps brackets: ids inside one bracket form one group, and
 * outside a bracket EVERY id is its own group, including each member of an
 * a:b range. So "42:48" is seven groups and "[42:48]" is one. The circuit
 * factory and the cohomology turn each group into one condition, so the
 * difference decides how many conditions a deck imposes.
 *
 * get_ids drops the brackets and returns one flat list, so both spellings
 * give the same seven ids.
 *
 * Both keep the word order of the deck and expand a range in ascending order
 * whichever way it is written.
 */

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_InputFile.hpp"

using namespace belfem ;

namespace
{
    // a deck in the test temp directory, removed again on scope exit; the
    // path is absolute, so the case does not depend on the working directory
    struct Deck
    {
        std::string mPath ;

        explicit Deck( const std::string & aBody ) :
            mPath( ::testing::TempDir() + "test_input_section_"
                + ::testing::UnitTest::GetInstance()->current_test_info()->name()
                + ".conf" )
        {
            std::ofstream tFile( mPath );
            tFile << "test\n{\n" << aBody << "}\n" ;
        }

        ~Deck()
        {
            std::remove( mPath.c_str() );
        }
    };

    const std::string gDeckBody =
        "    a : 42:48 ;\n"
        "    b : [42:48] ;\n"
        "    c : 3, 1, [7, 5], 48:42 ;\n"
        "    d : x, 1 ;\n" ;

    void
    expect_groups( const Cell< Cell< id_t > > & aGroups,
                   const Cell< Cell< id_t > > & aExpected )
    {
        ASSERT_EQ( aGroups.size(), aExpected.size() );
        for ( index_t g = 0; g < aGroups.size(); ++g )
        {
            ASSERT_EQ( aGroups( g ).size(), aExpected( g ).size() ) << "group " << g ;
            for ( index_t k = 0; k < aGroups( g ).size(); ++k )
            {
                EXPECT_EQ( aGroups( g )( k ), aExpected( g )( k ) )
                    << "group " << g << ", entry " << k ;
            }
        }
    }

    void
    expect_ids( const Vector< id_t > & aIDs, const Cell< id_t > & aExpected )
    {
        ASSERT_EQ( aIDs.length(), aExpected.size() );
        for ( index_t k = 0; k < aIDs.length(); ++k )
        {
            EXPECT_EQ( aIDs( k ), aExpected( k ) ) << "entry " << k ;
        }
    }

    // runs aCall and requires a std::runtime_error whose what() carries both
    // literal fragments of the BELFEM_ERROR text
    template< typename F >
    void
    expect_refusal( F && aCall,
        const std::string & aFirst, const std::string & aSecond )
    {
        bool tThrown = false ;
        try
        {
            aCall();
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;

            const std::string tMessage( aError.what() );
            EXPECT_NE( tMessage.find( aFirst ), std::string::npos )
                << "message was: " << tMessage ;
            EXPECT_NE( tMessage.find( aSecond ), std::string::npos )
                << "message was: " << tMessage ;
        }
        EXPECT_TRUE( tThrown )
            << "expected a refusal carrying \"" << aFirst
            << "\" and \"" << aSecond << "\"" ;
    }
}

TEST( InputSection, UnbracketedRangeIsOneGroupPerId )
{
    Deck tDeck( gDeckBody );
    InputFile tInput( tDeck.mPath );

    Cell< Cell< id_t > > tGroups ;
    tInput.section( "test" )->get_id_groups( "a", tGroups );

    expect_groups( tGroups, { { 42 }, { 43 }, { 44 }, { 45 }, { 46 }, { 47 }, { 48 } } );
}

TEST( InputSection, BracketedRangeIsOneGroup )
{
    Deck tDeck( gDeckBody );
    InputFile tInput( tDeck.mPath );

    Cell< Cell< id_t > > tGroups ;
    tInput.section( "test" )->get_id_groups( "b", tGroups );

    expect_groups( tGroups, { { 42, 43, 44, 45, 46, 47, 48 } } );
}

TEST( InputSection, GetIdsIgnoresBrackets )
{
    Deck tDeck( gDeckBody );
    InputFile tInput( tDeck.mPath );

    Vector< id_t > tUnbracketed ;
    Vector< id_t > tBracketed ;
    tInput.section( "test" )->get_ids( "a", tUnbracketed );
    tInput.section( "test" )->get_ids( "b", tBracketed );

    expect_ids( tUnbracketed, { 42, 43, 44, 45, 46, 47, 48 } );
    expect_ids( tBracketed,   { 42, 43, 44, 45, 46, 47, 48 } );
}

TEST( InputSection, WordOrderIsKeptAndRangesAreAscending )
{
    Deck tDeck( gDeckBody );
    InputFile tInput( tDeck.mPath );

    Cell< Cell< id_t > > tGroups ;
    tInput.section( "test" )->get_id_groups( "c", tGroups );

    // the bracket keeps its own word order; only a range is sorted
    expect_groups( tGroups,
        { { 3 }, { 1 }, { 7, 5 }, { 42 }, { 43 }, { 44 }, { 45 }, { 46 }, { 47 }, { 48 } } );

    Vector< id_t > tIDs ;
    tInput.section( "test" )->get_ids( "c", tIDs );
    expect_ids( tIDs, { 3, 1, 7, 5, 42, 43, 44, 45, 46, 47, 48 } );
}

TEST( InputSection, NonIntegerIdIsRefused )
{
    Deck tDeck( gDeckBody );
    InputFile tInput( tDeck.mPath );

    const input::Section * tSection = tInput.section( "test" );

    // the bad word comes first: "1, x" would reach the unit parser when the
    // deck is read, as a number followed by the unit x, and refuse there.
    // get_ids names the section type
    expect_refusal( [ & ]()
        {
            Vector< id_t > tIDs ;
            tSection->get_ids( "d", tIDs );
        },
        "Error reading id list in section test", ": x" );

    expect_refusal( [ & ]()
        {
            Cell< Cell< id_t > > tGroups ;
            tSection->get_id_groups( "d", tGroups );
        },
        "Error reading id list", ": x" );
}
