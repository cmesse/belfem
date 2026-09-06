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
 * Unit tests for src/circuit/fn_spice_number.{hpp,cpp}, the single point
 * where netlist number tokens become SI values. Expectations follow the
 * ngspice manual ("Ngspice scale factors" table, "Letters following a
 * number" rule); QUIRK-flagged cases pin ngspice-compatible behavior that
 * is easy to get wrong rather than an idealized contract.
 */

#include <gtest/gtest.h>
#include <clocale>
#include <cmath>
#include <stdexcept>

#include "typedefs.hpp"
#include "fn_spice_number.hpp"

using belfem::electronics::spice_number_to_si;

//------------------------------------------------------------------------------
// plain numbers
//------------------------------------------------------------------------------

TEST( SpiceNumber, PlainIntegersAndFloats )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( "12" ),      12.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "-44" ),    -44.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "+5" ),       5.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "3.14159" ),  3.14159 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "100" ),    100.0 );
}

TEST( SpiceNumber, BareAndTrailingDecimalPoint )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( ".5" ),  0.5 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1." ),  1.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "-.25" ), -0.25 );
}

TEST( SpiceNumber, IntegerExponents )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1e-14" ),  1.0e-14 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "2.65e3" ), 2.65e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1E3" ),    1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1e+3" ),   1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "2.5e-6" ), 2.5e-6 );
}

//------------------------------------------------------------------------------
// scale factors
//------------------------------------------------------------------------------

TEST( SpiceNumber, AllScaleFactors )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1T" ),   1.0e12 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1G" ),   1.0e9 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1MEG" ), 1.0e6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1k" ),   1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1mil" ), 25.4e-6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1m" ),   1.0e-3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1u" ),   1.0e-6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1n" ),   1.0e-9 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1p" ),   1.0e-12 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1f" ),   1.0e-15 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1a" ),   1.0e-18 );
}

TEST( SpiceNumber, MIsMilliAndMegIsMega )
{
    // the classic trap: M = milli, MEG = mega, case-insensitively
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1M" ),    1.0e-3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1Meg" ),  1.0e6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1meg" ),  1.0e6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "10MEG" ), 1.0e7 );
}

TEST( SpiceNumber, FIsFemtoNotFarad )
{
    // the showcase footgun: a 100 F capacitor is written "100"; "100F"
    // is 100 femto
    EXPECT_DOUBLE_EQ( spice_number_to_si( "100F" ), 1.0e-13 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "100" ),  100.0 );
}

TEST( SpiceNumber, AIsAttoNotAmpere )
{
    // the twin trap for current sources: a 1 A source is written "1"
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1A" ),  1.0e-18 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "10A" ), 1.0e-17 );
}

TEST( SpiceNumber, CaseInsensitiveFactors )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1K" ),   1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1MIL" ), 25.4e-6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1Mil" ), 25.4e-6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "2.2NF" ), 2.2e-9 );  // n + trailing F
}

TEST( SpiceNumber, SignedValueWithFactorAndExponent )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( "-1.5k" ),    -1.5e3 );
    // exponent and scale factor compose
    EXPECT_DOUBLE_EQ( spice_number_to_si( "-1.5e-3k" ), -1.5 );
}

//------------------------------------------------------------------------------
// trailing letters (units) are ignored
//------------------------------------------------------------------------------

TEST( SpiceNumber, TrailingUnitLettersIgnored )
{
    // manual: 10, 10V, 10Volts and 10Hz all represent the same number
    EXPECT_DOUBLE_EQ( spice_number_to_si( "10V" ),     10.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "10Volts" ), 10.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "10Hz" ),    10.0 );
}

TEST( SpiceNumber, LettersAfterScaleFactorIgnored )
{
    // manual: M, MA, MSec and MMhos all represent the same scale factor
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1kOhm" ),   1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1kHz" ),    1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1MSec" ),   1.0e-3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1MMhos" ),  1.0e-3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1megohm" ), 1.0e6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "2.5kOhm" ), 2.5e3 );
}

TEST( SpiceNumber, LongestMatchWinsOverMilli )
{
    // QUIRK (BELFEM policy, forced by the manual's table): "mil" is a
    // 3-letter scale factor, so longest match is required -- otherwise
    // "1mil" would read as milli + "il" and the table entry would be dead.
    // "1milli" = mil + trailing "li" is the consequence; not traced
    // against an ngspice binary
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1milli" ), 25.4e-6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1mils" ),  25.4e-6 );
}

TEST( SpiceNumber, LoneEIsAUnitLetterNotAnExponent )
{
    // "1e" has no exponent digits, so per the manual's rules the 'e' is
    // not an integer exponent and not a scale factor -- it is an ignored
    // unit letter ( BELFEM reading; not traced against an ngspice binary )
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1e" ),     1.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1exact" ), 1.0 );
}

TEST( SpiceNumber, BoundaryAndCompositionCases )
{
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1.e3" ),  1.0e3 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "+.5" ),   0.5 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( ".5k" ),   500.0 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1E3K" ),  1.0e6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1MEGV" ), 1.0e6 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "1e308" ),  1.0e308 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "-1e308" ), -1.0e308 );
    EXPECT_DOUBLE_EQ( spice_number_to_si( "0" ),  0.0 );
    // signed zero survives the scale multiplication
    EXPECT_TRUE( std::signbit( spice_number_to_si( "-0" ) ) );
}

TEST( SpiceNumber, LocaleIndependentConversion )
{
    // regression for the strtod hazard: conversion must not depend on
    // LC_NUMERIC ( std::from_chars is locale-independent by contract )
    const char * tOld = std::setlocale( LC_NUMERIC, "de_DE.UTF-8" );
    if ( tOld == nullptr )
    {
        GTEST_SKIP() << "de_DE.UTF-8 locale not installed";
    }
    const double tValue = spice_number_to_si( "3.14" );
    std::setlocale( LC_NUMERIC, "C" );
    EXPECT_DOUBLE_EQ( tValue, 3.14 );
}

//------------------------------------------------------------------------------
// malformed tokens hard-error
//------------------------------------------------------------------------------

TEST( SpiceNumber, RejectsNonNumbers )
{
    EXPECT_THROW( spice_number_to_si( "" ),    std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "abc" ), std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "k" ),   std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "-" ),   std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "." ),   std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "+." ),  std::runtime_error );
}

TEST( SpiceNumber, RejectsStrtodExtras )
{
    // strtod would accept these; the SPICE grammar must not
    EXPECT_THROW( spice_number_to_si( "inf" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "nan" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "0x10" ), std::runtime_error );
}

TEST( SpiceNumber, RejectsOutOfRange )
{
    EXPECT_THROW( spice_number_to_si( "1e999" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "-1e999" ), std::runtime_error );
    // pinned policy: accepted magnitudes are zero or the normal double
    // range -- underflow AND subnormals reject loudly, on every toolchain
    EXPECT_THROW( spice_number_to_si( "1e-999" ), std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1e-320" ), std::runtime_error );
    // the same policy applies after scaling
    EXPECT_THROW( spice_number_to_si( "1e-300f" ), std::runtime_error );
}

TEST( SpiceNumber, RejectsScaledOverflow )
{
    // a finite mantissa can still overflow through the scale factor
    EXPECT_THROW( spice_number_to_si( "1e300T" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "-1e300T" ), std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "9e307k" ),  std::runtime_error );
}

TEST( SpiceNumber, RejectsEmbeddedSuffixForm )
{
    // old-style "2k5" ( = 2.5k in some pre-ngspice dialects ) must not
    // silently parse as 2000
    EXPECT_THROW( spice_number_to_si( "2k5" ), std::runtime_error );
}

TEST( SpiceNumber, RejectsWhitespaceAndGarbage )
{
    EXPECT_THROW( spice_number_to_si( " 1" ),   std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1 " ),   std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1 k" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1..2" ), std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1k%" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "--1" ),  std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "e3" ),   std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1.2.3" ), std::runtime_error );
    EXPECT_THROW( spice_number_to_si( "1-2" ),  std::runtime_error );
    // incomplete exponent: the 'e' becomes a unit letter, the '-' is garbage
    EXPECT_THROW( spice_number_to_si( "1e-" ),  std::runtime_error );
}
