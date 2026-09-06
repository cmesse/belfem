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
 * Unit tests for src/core/stringtools.{hpp,cpp}. These functions carry the
 * input-deck parser, so several of the expectations below pin down quirks of
 * the current behavior rather than an idealized contract; those are flagged
 * QUIRK and derived from the implementation, not from a specification.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "stringtools.hpp"
#include "fn_check_unit.hpp"

using namespace belfem;

//------------------------------------------------------------------------------
// clean_string
//------------------------------------------------------------------------------

TEST( StringTools, CleanStringCollapsesAndTrimsSpaces )
{
    EXPECT_EQ( clean_string( " Alpha  Bravo     Charlie Delta " ),
               "Alpha Bravo Charlie Delta" );
}

TEST( StringTools, CleanStringCollapsesTabsAndNewlines )
{
    // tab, carriage return and line feed all count as whitespace
    EXPECT_EQ( clean_string( "Alpha\tBravo\r\nCharlie" ), "Alpha Bravo Charlie" );
}

TEST( StringTools, CleanStringLeavesTidyStringAlone )
{
    // fast path: nothing to clean
    EXPECT_EQ( clean_string( "Alpha" ), "Alpha" );
    EXPECT_EQ( clean_string( "Alpha Bravo" ), "Alpha Bravo" );
}

TEST( StringTools, CleanStringTruncatesAtHash )
{
    // everything from the first '#' onwards is a comment and is dropped.
    // QUIRK: the space in front of the '#' survives, because the trailing-space
    // erase only fires when the loop runs off the end of the string.
    EXPECT_EQ( clean_string( "value 42 # this is a comment" ), "value 42 " );

    // without the leading space there is nothing left to keep
    EXPECT_EQ( clean_string( "#comment only" ), "" );
}

TEST( StringTools, CleanStringPreservesSpacesInsideQuotes )
{
    // quoted runs are copied verbatim, including the quotes themselves
    EXPECT_EQ( clean_string( "key \"two  spaces\" tail" ),
               "key \"two  spaces\" tail" );
}

TEST( StringTools, CleanStringSingleCharacterPassesThrough )
{
    // strings of length <= 1 are returned unchanged
    EXPECT_EQ( clean_string( " " ), " " );
    EXPECT_EQ( clean_string( "x" ), "x" );
}

//------------------------------------------------------------------------------
// first_word
//------------------------------------------------------------------------------

TEST( StringTools, FirstWordCleansBeforeSplitting )
{
    EXPECT_EQ( first_word( "  Alpha   Bravo Charlie " ), "Alpha" );
}

TEST( StringTools, FirstWordOfSingleWordIsWholeString )
{
    EXPECT_EQ( first_word( "Alpha" ), "Alpha" );
}

TEST( StringTools, FirstWordHonorsDelimiter )
{
    EXPECT_EQ( first_word( "Alpha,Bravo,Charlie", ',' ), "Alpha" );
    EXPECT_EQ( first_word( "Alpha:Bravo", ':' ), "Alpha" );

    // the default delimiter does not split on a comma
    EXPECT_EQ( first_word( "Alpha,Bravo Charlie" ), "Alpha,Bravo" );
}

//------------------------------------------------------------------------------
// string_to_words
//------------------------------------------------------------------------------

TEST( StringTools, StringToWordsSplitsOnSpace )
{
    Cell< string > tWords = string_to_words( " Alpha  Bravo     Charlie Delta " );

    ASSERT_EQ( tWords.size(), (uint) 4 );
    EXPECT_EQ( tWords( 0 ), "Alpha" );
    EXPECT_EQ( tWords( 1 ), "Bravo" );
    EXPECT_EQ( tWords( 2 ), "Charlie" );
    EXPECT_EQ( tWords( 3 ), "Delta" );
}

TEST( StringTools, StringToWordsSingleWord )
{
    Cell< string > tWords = string_to_words( "Echo" );

    ASSERT_EQ( tWords.size(), (uint) 1 );
    EXPECT_EQ( tWords( 0 ), "Echo" );
}

TEST( StringTools, StringToWordsWhitespaceOnlyYieldsNothing )
{
    EXPECT_EQ( string_to_words( " " ).size(), (uint) 0 );
}

TEST( StringTools, StringToWordsCustomDelimiterSkipsEmptyFields )
{
    // with a non-space delimiter the string is NOT cleaned first,
    // and empty fields between two delimiters are dropped
    Cell< string > tWords = string_to_words( "Alpha,Bravo,,Charlie", ',' );

    ASSERT_EQ( tWords.size(), (uint) 3 );
    EXPECT_EQ( tWords( 0 ), "Alpha" );
    EXPECT_EQ( tWords( 1 ), "Bravo" );
    EXPECT_EQ( tWords( 2 ), "Charlie" );
}

TEST( StringTools, StringToWordsCustomDelimiterKeepsSpaces )
{
    Cell< string > tWords = string_to_words( "Alpha Bravo;Charlie", ';' );

    ASSERT_EQ( tWords.size(), (uint) 2 );
    EXPECT_EQ( tWords( 0 ), "Alpha Bravo" );
    EXPECT_EQ( tWords( 1 ), "Charlie" );
}

//------------------------------------------------------------------------------
// search_and_replace
//------------------------------------------------------------------------------

TEST( StringTools, SearchAndReplaceReplacesAllOccurrences )
{
    EXPECT_EQ( search_and_replace( "a-b-c", "-", "+" ), "a+b+c" );
    EXPECT_EQ( search_and_replace( "aaa", "a", "bb" ), "bbbbbb" );
}

TEST( StringTools, SearchAndReplaceHandlesMultiCharacterNeedle )
{
    EXPECT_EQ( search_and_replace( "hello world", "world", "there" ),
               "hello there" );
}

TEST( StringTools, SearchAndReplaceCanDelete )
{
    EXPECT_EQ( search_and_replace( "{1,2,3}", "{", "" ), "1,2,3}" );
}

TEST( StringTools, SearchAndReplaceNoMatchReturnsInput )
{
    EXPECT_EQ( search_and_replace( "Alpha", "z", "y" ), "Alpha" );
}

TEST( StringTools, SearchAndReplaceEmptyNeedleReturnsInput )
{
    // guards against an infinite loop
    EXPECT_EQ( search_and_replace( "Alpha", "", "y" ), "Alpha" );
}

//------------------------------------------------------------------------------
// string_to_lower / string_to_upper
//------------------------------------------------------------------------------

TEST( StringTools, StringToLower )
{
    EXPECT_EQ( string_to_lower( "AlPhA 42" ), "alpha 42" );
}

TEST( StringTools, StringToUpper )
{
    EXPECT_EQ( string_to_upper( "AlPhA 42" ), "ALPHA 42" );
}

//------------------------------------------------------------------------------
// string_to_bool
//------------------------------------------------------------------------------

TEST( StringTools, StringToBoolAcceptsTrueSpellings )
{
    EXPECT_TRUE( string_to_bool( "true" ) );
    EXPECT_TRUE( string_to_bool( "TRUE" ) );
    EXPECT_TRUE( string_to_bool( "On" ) );
    EXPECT_TRUE( string_to_bool( "yes" ) );
    EXPECT_TRUE( string_to_bool( "1" ) );
}

TEST( StringTools, StringToBoolEverythingElseIsFalse )
{
    EXPECT_FALSE( string_to_bool( "false" ) );
    EXPECT_FALSE( string_to_bool( "off" ) );
    EXPECT_FALSE( string_to_bool( "no" ) );
    EXPECT_FALSE( string_to_bool( "0" ) );
    EXPECT_FALSE( string_to_bool( "maybe" ) );

    // QUIRK: not trimmed -- a stray space makes it false
    EXPECT_FALSE( string_to_bool( " true" ) );
}

//------------------------------------------------------------------------------
// to_real
//------------------------------------------------------------------------------

TEST( StringTools, ToRealParsesDecimalAndExponent )
{
    EXPECT_NEAR( to_real( "3.5" ), 3.5, 1e-12 );
    EXPECT_NEAR( to_real( "-1.25" ), -1.25, 1e-12 );
    EXPECT_NEAR( to_real( "1e-3" ), 1.0e-3, 1e-15 );
}

TEST( StringTools, ToRealRejectsNonNumeric )
{
    EXPECT_TRUE( std::isnan( to_real( "Alpha" ) ) );
}

TEST( StringTools, ToRealStopsAtTrailingGarbage )
{
    // QUIRK: strtod-style parsing -- a numeric prefix is accepted and the
    // remainder is silently ignored
    EXPECT_NEAR( to_real( "3.5Alpha" ), 3.5, 1e-12 );
}

TEST( StringTools, ToRealOfEmptyStringIsNan )
{
    // nothing to convert is not a real, as the header promises
    EXPECT_TRUE( std::isnan( to_real( "" ) ) );
}

//------------------------------------------------------------------------------
// is_integer
//------------------------------------------------------------------------------

TEST( StringTools, IsIntegerAcceptsSignedDigits )
{
    EXPECT_TRUE( is_integer( "42" ) );
    EXPECT_TRUE( is_integer( "+42" ) );
    EXPECT_TRUE( is_integer( "-42" ) );
    EXPECT_TRUE( is_integer( "0" ) );
}

TEST( StringTools, IsIntegerRejectsEverythingElse )
{
    EXPECT_FALSE( is_integer( "" ) );
    EXPECT_FALSE( is_integer( "-" ) );
    EXPECT_FALSE( is_integer( "+" ) );
    EXPECT_FALSE( is_integer( "4.2" ) );
    EXPECT_FALSE( is_integer( "42a" ) );
    EXPECT_FALSE( is_integer( "Alpha" ) );
    EXPECT_FALSE( is_integer( " 42" ) );
}

//------------------------------------------------------------------------------
// path helpers
//------------------------------------------------------------------------------

TEST( StringTools, BasenameKeepsExtension )
{
    // belfem:: is required here: glibc also declares a char* basename( const
    // char* ) in <string.h>, and for a string literal that C overload is the
    // exact match, so an unqualified call compares pointers instead of text.
    // Callers passing a std::string get belfem::basename either way.
    EXPECT_EQ( belfem::basename( "/home/user/mesh.msh" ), "mesh.msh" );
    EXPECT_EQ( belfem::basename( "mesh.msh" ), "mesh.msh" );

    // note: unlike POSIX basename(1), no extension is stripped
    EXPECT_EQ( belfem::basename( string( "/home/user/mesh.msh" ) ), "mesh.msh" );
}

TEST( StringTools, Dirname )
{
    EXPECT_EQ( belfem::dirname( "/home/user/mesh.msh" ), "/home/user" );

    // no separator at all -> empty directory
    EXPECT_EQ( belfem::dirname( "mesh.msh" ), "" );
}

TEST( StringTools, Filename )
{
    EXPECT_EQ( belfem::filename( "/home/user/mesh.msh" ), "mesh.msh" );
}

TEST( StringTools, Filetype )
{
    EXPECT_EQ( belfem::filetype( "/home/user/mesh.msh" ), "msh" );
    EXPECT_EQ( belfem::filetype( "archive.tar.gz" ), "gz" );
}

TEST( StringTools, FiletypeOfExtensionlessPathReturnsWholeString )
{
    // QUIRK: with no '.' present, find_last_of returns npos and npos + 1
    // wraps to 0, so the entire path is returned
    EXPECT_EQ( belfem::filetype( "mesh" ), "mesh" );
}

//------------------------------------------------------------------------------
// format_with_leading_zeros
//------------------------------------------------------------------------------

TEST( StringTools, FormatWithLeadingZerosReturnsAFormatString )
{
    // this returns a printf format, not a formatted number
    EXPECT_EQ( format_with_leading_zeros( 7 ), "%01u" );
    EXPECT_EQ( format_with_leading_zeros( 42 ), "%02u" );
    EXPECT_EQ( format_with_leading_zeros( 999 ), "%03u" );
    EXPECT_EQ( format_with_leading_zeros( 1234 ), "%04u" );
    EXPECT_EQ( format_with_leading_zeros( 123456789 ), "%09u" );
}

//------------------------------------------------------------------------------
// utf8_character_count
//------------------------------------------------------------------------------

TEST( StringTools, Utf8CharacterCountAscii )
{
    EXPECT_EQ( utf8_character_count( "Alpha" ), (size_t) 5 );
    EXPECT_EQ( utf8_character_count( "" ), (size_t) 0 );
}

TEST( StringTools, Utf8CharacterCountMultiByte )
{
    // "\xC2\xB5" is U+00B5 MICRO SIGN, a two-byte sequence
    EXPECT_EQ( utf8_character_count( "\xC2\xB5" "m" ), (size_t) 2 );

    // "\xC2\xB0" is U+00B0 DEGREE SIGN
    EXPECT_EQ( utf8_character_count( "0 \xC2\xB0" "C" ), (size_t) 4 );
}

//------------------------------------------------------------------------------
// unit_to_si
//------------------------------------------------------------------------------

TEST( StringTools, UnitToSiDimensionless )
{
    // both the placeholder "-" and the empty string mean "no unit"
    for ( const string & tSpec : { string( "-" ), string( "" ) } )
    {
        value tValue = unit_to_si( tSpec );

        EXPECT_NEAR( tValue.first, 1.0, 1e-12 );

        for ( uint k = 0; k < 7; ++k )
        {
            EXPECT_NEAR( tValue.second[ k ], 0.0, 1e-12 );
        }
    }
}

//------------------------------------------------------------------------------
// unit_to_si -- magnetic flux density
//------------------------------------------------------------------------------

// Tesla is V*s/m^2 = kg/( s^2 * A ). Until 2026-08-30 every tesla token carried
// the VOLT exponents, because the block was copied from "V" with only the string
// changed. check_unit compares nothing but the exponent array, so the two units
// were indistinguishable. These cases fail against that table and pass against
// the corrected one.

TEST( StringTools, UnitToSiTeslaHasFluxDensityDimension )
{
    const value tTesla = unit_to_si( "T" );

    EXPECT_NEAR( tTesla.first, 1.0, 1e-12 );

    // L, M, T, I, theta, N, J
    EXPECT_NEAR( tTesla.second[ 0 ],  0.0, 1e-12 );   // NOT 2 -- the volt bug
    EXPECT_NEAR( tTesla.second[ 1 ],  1.0, 1e-12 );
    EXPECT_NEAR( tTesla.second[ 2 ], -2.0, 1e-12 );   // NOT -3 -- the volt bug
    EXPECT_NEAR( tTesla.second[ 3 ], -1.0, 1e-12 );

    for ( uint k = 4; k < 7; ++k )
    {
        EXPECT_NEAR( tTesla.second[ k ], 0.0, 1e-12 );
    }
}

TEST( StringTools, UnitToSiTeslaIsNotVolt )
{
    // the discriminating case: these two compared EQUAL before the fix
    EXPECT_FALSE( check_unit( unit_to_si( "T" ), "V" ) );
    EXPECT_FALSE( check_unit( unit_to_si( "V" ), "T" ) );
}

TEST( StringTools, UnitToSiTeslaEqualsCompoundForm )
{
    // the other half of the defect: the correctly written compound reduces to
    // the true tesla vector, and used to be REFUSED where "T" was accepted
    EXPECT_TRUE( check_unit( unit_to_si( "V*s/m^2" ), "T" ) );
    EXPECT_TRUE( check_unit( unit_to_si( "kg/(s^2*A)" ), "T" ) );
}

TEST( StringTools, UnitToSiTeslaFamilyScalesAndAgrees )
{
    // the scale factors were never wrong; keep them pinned, and keep the whole
    // family on one dimension so conversions within it stay legal
    const std::pair< string, real > tCases[] = {
            { "G",   1e-4 }, { "muT", 1e-6 }, { "mT", 1e-3 },
            { "T",   1.0  }, { "kT",  1e3  }, { "MT", 1e6  } };

    for ( const auto & tCase : tCases )
    {
        const value tValue = unit_to_si( tCase.first );

        EXPECT_NEAR( tValue.first, tCase.second, tCase.second * 1e-12 );
        EXPECT_TRUE( check_unit( tValue, "T" ) );
    }
}

TEST( StringTools, UnitToSiForceFamilyScalesAndAgrees )
{
    // MN carried kN's factor until 2026-09-04
    const std::pair< string, real > tCases[] = {
            { "muN", 1e-6 }, { "mN", 1e-3 }, { "N", 1.0 },
            { "kN",  1e3  }, { "MN", 1e6  } };

    for ( const auto & tCase : tCases )
    {
        const value tValue = unit_to_si( tCase.first );

        EXPECT_NEAR( tValue.first, tCase.second, tCase.second * 1e-12 );
        EXPECT_TRUE( check_unit( tValue, "N" ) );
    }
}

TEST( StringTools, UnitToSiVoltFamilyScalesAndAgrees )
{
    // muV carried mV's factor until 2026-09-04
    const std::pair< string, real > tCases[] = {
            { "muV", 1e-6 }, { "mV", 1e-3 }, { "V", 1.0 }, { "kV", 1e3 } };

    for ( const auto & tCase : tCases )
    {
        const value tValue = unit_to_si( tCase.first );

        EXPECT_NEAR( tValue.first, tCase.second, tCase.second * 1e-12 );
        EXPECT_TRUE( check_unit( tValue, "V" ) );
    }
}

//------------------------------------------------------------------------------
// string_to_cell
//------------------------------------------------------------------------------

TEST( StringTools, StringToCellReadsPlainList )
{
    Cell< int > tValues;
    string_to_cell< int >( "{ 1, 2, 5 }", tValues );

    ASSERT_EQ( tValues.size(), (uint) 3 );
    EXPECT_EQ( tValues( 0 ), 1 );
    EXPECT_EQ( tValues( 1 ), 2 );
    EXPECT_EQ( tValues( 2 ), 5 );
}

TEST( StringTools, StringToCellExpandsAscendingRange )
{
    Cell< int > tValues;
    string_to_cell< int >( "{ 3:6 }", tValues );

    ASSERT_EQ( tValues.size(), (uint) 4 );
    EXPECT_EQ( tValues( 0 ), 3 );
    EXPECT_EQ( tValues( 1 ), 4 );
    EXPECT_EQ( tValues( 2 ), 5 );
    EXPECT_EQ( tValues( 3 ), 6 );
}

TEST( StringTools, StringToCellExpandsDescendingRange )
{
    Cell< int > tValues;
    string_to_cell< int >( "{ 6:3 }", tValues );

    ASSERT_EQ( tValues.size(), (uint) 4 );
    EXPECT_EQ( tValues( 0 ), 6 );
    EXPECT_EQ( tValues( 1 ), 5 );
    EXPECT_EQ( tValues( 2 ), 4 );
    EXPECT_EQ( tValues( 3 ), 3 );
}

TEST( StringTools, StringToCellStopsAtSemicolon )
{
    Cell< int > tValues;
    string_to_cell< int >( "{ 1, 2 ; 3, 4 }", tValues );

    ASSERT_EQ( tValues.size(), (uint) 2 );
    EXPECT_EQ( tValues( 0 ), 1 );
    EXPECT_EQ( tValues( 1 ), 2 );
}

//------------------------------------------------------------------------------
// to_pair
//------------------------------------------------------------------------------

TEST( StringTools, ToPairSplitsSymbolNumberGroups )
{
    Cell< std::pair< string, int > > tResult;
    to_pair< int >( "Pb37Sn63", tResult );

    ASSERT_EQ( tResult.size(), (uint) 2 );
    EXPECT_EQ( tResult( 0 ).first, "Pb" );
    EXPECT_EQ( tResult( 0 ).second, 37 );
    EXPECT_EQ( tResult( 1 ).first, "Sn" );
    EXPECT_EQ( tResult( 1 ).second, 63 );
}

TEST( StringTools, ToPairClearsPreviousContent )
{
    Cell< std::pair< string, int > > tResult;
    to_pair< int >( "Cu100", tResult );
    to_pair< int >( "RRR50", tResult );

    ASSERT_EQ( tResult.size(), (uint) 1 );
    EXPECT_EQ( tResult( 0 ).first, "RRR" );
    EXPECT_EQ( tResult( 0 ).second, 50 );
}

//------------------------------------------------------------------------------
// datatype_string
//------------------------------------------------------------------------------

TEST( StringTools, DatatypeString )
{
    EXPECT_EQ( datatype_string< bool >(), "bool" );
    EXPECT_EQ( datatype_string< int >(), "int" );
    EXPECT_EQ( datatype_string< double >(), "double" );
    EXPECT_EQ( datatype_string< std::string >(), "string" );

    // unspecialized fallback
    EXPECT_EQ( datatype_string< char >(), "unknown" );
}
