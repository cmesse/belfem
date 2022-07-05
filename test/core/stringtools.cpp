//
// Created by Christian Messe on 2018-12-26.
//

//
#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "stringtools.hpp"

using namespace belfem;

TEST( Stringtools, clean_string )
{
    string tString = " Alpha  Bravo     Charlie Delta ";
    string tClean = clean_string( tString );
    string tExpect = "Alpha Bravo Charlie Delta";
    EXPECT_EQ( tExpect, tClean );

    Cell< string > tWords = string_to_words( tString );

    EXPECT_EQ( tWords.size(), (uint) 4 );
    EXPECT_EQ( tWords( 0 ), "Alpha" );
    EXPECT_EQ( tWords( 1 ), "Bravo" );
    EXPECT_EQ( tWords( 2 ), "Charlie" );
    EXPECT_EQ( tWords( 3 ), "Delta");

    // also test if function detects only one word
    Cell< string > tWord = string_to_words( "Echo" );

    EXPECT_EQ( tWord.size(), (uint) 1 );
    EXPECT_EQ( tWord( 0 ), "Echo" );

    // test empty string
    Cell< string > tEmpty = string_to_words( " " );
    EXPECT_EQ( tEmpty.size(),(uint) 0 );

}
