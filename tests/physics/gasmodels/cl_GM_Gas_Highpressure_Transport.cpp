//
// Created by Christian Messe on 25.09.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Communicator.hpp"


#include "fn_r2.hpp"
#include "fn_linspace.hpp"

#include "cl_Gas.hpp"
#include "en_GM_GasModel.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

/**
 * This test creates a cubic gas model of air and compares mu and lambda
 * at 250 bar with reference data from VDI heat atlas.
 */
TEST( GASMODELS, HighpressureTransport )
{
//------------------------------------------------------------------------------
// VDI Heat Atlas D2.2 Tables 9 and 10
//------------------------------------------------------------------------------
    Vector< real > tT = { -150., -125., -100., -75., -50., -25., 0., 25., 50.,
                          75., 100., 125., 150., 200., 300., 400., 500., 600.,
                          700., 800., 900., 1000 };

    tT += 273.15;

    Vector< real > tLambda_1 = { 11.679, 13.984, 16.205, 18.347, 20.416, 22.418,
                                 24.360, 26.247, 28.082, 29.872, 31.620, 33.328,
                                 35.000, 38.248, 44.417, 50.240, 55.795, 61.139,
                                 66.312, 71.348, 76.271, 81.099 };

    tLambda_1 *= 0.001;

    Vector< real > tLambda_250 = { 106.410, 81.886, 64.314, 53.295, 46.985, 43.506,
                                   41.545, 41.025, 41.059, 41.445, 42.064, 42.843,
                                   43.735, 45.742, 50.222, 54.946, 59.730, 64.505,
                                   69.243, 73.936, 78.581, 83.183 };
    tLambda_250 *= 0.001;

    Vector< real > tMu_1 = { 8.664, 10.261, 11.780, 13.229, 14.614, 15.942, 17.218,
                             18.448, 19.635, 20.783, 21.896, 23.024, 24.072, 26.087,
                             29.845, 33.314, 36.557, 39.621, 42.538, 45.337, 48.036,
                             50.651};

    tMu_1 *= 1e-6;

    Vector< real > tMu_250 = { 84.346, 57.676, 41.945, 33.032, 28.507, 26.419, 25.590,
                               25.423, 25.624, 26.041, 26.590, 28.248, 28.830, 30.133,
                               32.978, 35.881, 38.736, 41.515, 44.216, 46.842, 49.401,
                               51.900};

    tMu_250 *= 1e-6;

    uint tNumSamples = tT.length();

//------------------------------------------------------------------------------

    Cell<string> tSpecies = { "N2", "O2", "Ar" };
    Vector<real> tMolarFractions = { 0.7812, 0.2096, 0.0092 };

    Gas tAirID( tSpecies, tMolarFractions, GasModel::IDGAS );

    // no need to create SRK since is the same
    Gas tAirPR( tSpecies, tMolarFractions, GasModel::PR );

//------------------------------------------------------------------------------
// Data containers
//------------------------------------------------------------------------------

    Vector< real > tValuesPR( tNumSamples );
    Vector< real > tValuesIDGAS( tNumSamples );


    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesIDGAS( k ) = tAirID.lambda( tT( k ), 1e5 );
        tValuesPR( k ) = tAirPR.lambda( tT( k ), 1e5 );
    }

    // lambda IDGAS p=1bar
    EXPECT_NEAR( r2( tValuesIDGAS, tLambda_1 ),  1.0,  0.01 );

    // lambda PR p=1bar
    EXPECT_NEAR( r2( tValuesPR, tLambda_1 ),  1.0,  0.01 );

    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tAirPR.lambda( tT( k ), 250e5 );
    }

    // lambda PR p=250bar
    EXPECT_NEAR( r2( tValuesPR, tLambda_250 ) ,  1.0,  0.02 );

//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesIDGAS( k ) = tAirID.mu( tT( k ), 1e5 );
        tValuesPR( k ) = tAirPR.mu( tT( k ), 1e5 );
    }

    // mu IDGAS p=1bar
    EXPECT_NEAR( r2( tValuesIDGAS, tMu_1 ),  1.0,  0.01 );

    // mu PR p=1bar
    EXPECT_NEAR( r2( tValuesPR, tMu_1 ),  1.0,  0.01 );

    for( uint k=0; k<5; ++k )
    {
        tValuesPR( k ) = 0;
        tMu_250( k ) = 0;
    }

    for( uint k=5; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tAirPR.mu( tT( k ), 250e5 );
    }

    // mu PR p=250bar
    EXPECT_NEAR( r2( tValuesPR, tMu_250  ),  1.0,  0.01 );
}