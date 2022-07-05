//
// Created by Christian Messe on 19.09.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Communicator.hpp"


#include "fn_r2.hpp"
#include "fn_linspace.hpp"

#include "cl_Gas.hpp"
#include "GT_globals.hpp"
#include "en_GM_GasModel.hpp"



using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

/**
 * This test calculates the entropy of air at 250 bar and compares
 * the values with literature data. Moreover, the entropy departure
 * derivative to p is testet
 */
TEST( GASMODELS, Cubic_Entropy )
{

//------------------------------------------------------------------------------

    Cell<string> tSpecies = { "N2", "O2", "Ar" };
    Vector<real> tMolarFractions = { 0.7812, 0.2096, 0.0092 };

    Gas tAirID( tSpecies, tMolarFractions, GasModel::IDGAS );
    Gas tAirPR( tSpecies, tMolarFractions, GasModel::PR );
    Gas tAirSRK( tSpecies, tMolarFractions, GasModel::SRK );

//------------------------------------------------------------------------------
// REFERENCE DATA, VDI HEAT ATLAS, D2.2 Table 7
//------------------------------------------------------------------------------

    // temperatures in °C
    Vector< real > tT = { -150., -125., -100., -75., -50., -25., 0., 25., 50.,
                          75., 100., 125., 150., 200., 300., 400., 500., 600.,
                          700., 800., 900., 1000 };

    // change unit to K
    tT += 273.15;

    real tP = 250e5;

    // Entropy at 250 bar
    Vector< real > tS = { -3.2325, -2.9055, -2.6306, -2.3979, -2.2026,
                          -2.0388, -1.8999, -1.7801, -1.6749, -1.581, -1.496,
                          -1.4183, -1.3467, -1.2177, -1.0009, -0.8204, -0.6641,
                          -0.5255, -0.4005, -0.2863, -0.1812, -0.0836 };

    // change unit to J/(kg*K)
    tS *= 1000;

    // add offset at reference point
    tS += tAirID.s( 298.15, 1e5 );

    uint tNumSamples = tT.length();

//------------------------------------------------------------------------------
// Data containers
//------------------------------------------------------------------------------

    Vector< real > tValuesPR( tNumSamples );
    Vector< real > tExpectPR( tNumSamples );
    Vector< real > tValuesSRK( tNumSamples );
    Vector< real > tExpectSRK( tNumSamples );

//------------------------------------------------------------------------------
// Entropy Test
//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tAirPR.s( tT( k ), tP );
        tValuesSRK( k ) = tAirSRK.s( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tValuesPR, tS ), 1.0, 0.005 );
    EXPECT_NEAR( r2( tValuesSRK, tS ), 1.0, 0.005 );

//------------------------------------------------------------------------------
// Test for pressure derivative of departure function
//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tAirPR.eos()->dsdepdp( tT( k ), tP );

        tValuesSRK( k ) = tAirSRK.eos()->dsdepdp( tT( k ), tP );

        tExpectPR( k ) = ( tAirPR.eos()->sdep( tT( k ), 1.001 * tP )
                           - tAirPR.eos()->sdep( tT( k ), 0.999 * tP ) ) / ( 0.002 * tP );
        tExpectSRK( k ) = ( tAirSRK.eos()->sdep( tT( k ), 1.001 * tP )
                            - tAirSRK.eos()->sdep( tT( k ), 0.999 * tP ) ) / ( 0.002 * tP );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 1e-6 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 1e-6 );

//------------------------------------------------------------------------------
// Test for pressure derivative of enthalpy departure function
//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tAirPR.eos()->dhdepdp( tT( k ), tP );

        tValuesSRK( k ) = tAirSRK.eos()->dhdepdp( tT( k ), tP );

        tExpectPR( k ) = ( tAirPR.eos()->hdep( tT( k ), 1.001 * tP )
                           - tAirPR.eos()->hdep( tT( k ), 0.999 * tP ) ) / ( 0.002 * tP );
        tExpectSRK( k ) = ( tAirSRK.eos()->hdep( tT( k ), 1.001 * tP )
                            - tAirSRK.eos()->hdep( tT( k ), 0.999 * tP ) ) / ( 0.002 * tP );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 1e-6 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 1e-6 );

//------------------------------------------------------------------------------
}
