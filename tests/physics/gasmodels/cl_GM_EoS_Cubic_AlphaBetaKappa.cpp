//
// Created by Christian Messe on 18.09.19.
//

//
// Created by Christian Messe on 15.09.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Communicator.hpp"
#include "cl_SpMatrix.hpp"

#include "fn_r2.hpp"
#include "fn_linspace.hpp"
#include "fn_cardano.hpp"
#include "fn_max.hpp"
#include "GT_globals.hpp"
#include "cl_Gas.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

TEST( GASMODELS, Cubic_AlphaBetaKappa )
{


/**
 * This test checks if the thermodynamic coefficients are correct
 */
//------------------------------------------------------------------------------

    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    Gas tPR( tSpecies, tMolarFractions, GasModel::PR );
    Gas tSRK( tSpecies, tMolarFractions, GasModel::SRK );

    uint tNumSamples = 501;

    Vector< real > tT;
    tT = linspace( 100.0, 2000.0, tNumSamples );

    real tP = 200e5;

//------------------------------------------------------------------------------

    Vector< real > tValuesPR( tNumSamples );
    Vector< real > tExpectPR( tNumSamples );

    Vector< real > tValuesSRK( tNumSamples );
    Vector< real > tExpectSRK( tNumSamples );

//------------------------------------------------------------------------------
// Test Beta
//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        real tV =   tPR.v( tT(k), tP );

        real tP1 =  tPR.p( tT(k)*0.999, tV );
        real tP2 =  tPR.p( tT(k)*1.001, tV );

        tValuesPR( k ) = tPR.beta( tT(k), tP );
        tExpectPR( k ) = ( tP2 - tP1 ) / ( 0.002 * tT(k) * tP );


        tV = tSRK.v( tT( k ), tP );
        tValuesSRK( k ) = tSRK.beta( tT(k), tP );
        tP1 =  tSRK.p( tT(k)*0.999, tV );
        tP2 =  tSRK.p( tT(k)*1.001, tV );
        tExpectSRK( k ) = ( tP2 - tP1 ) / ( 0.002 * tT(k) * tP );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 1e-6 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 1e-6 );

//------------------------------------------------------------------------------
// Test kappa
//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        real tV =   tPR.v( tT(k), tP );

        real tV1 =  tPR.v( tT(k), tP*0.999 );
        real tV2 =  tPR.v( tT(k), tP*1.001 );

        tValuesPR( k ) = tPR.kappa( tT(k), tP );
        tExpectPR( k ) = -( tV2 - tV1 ) / ( 0.002 *tV * tP );

        tV =   tSRK.v( tT(k), tP );
        tV1 =  tSRK.v( tT(k), tP*0.999 );
        tV2 =  tSRK.v( tT(k), tP*1.001 );

        tValuesSRK( k ) = tSRK.kappa( tT(k), tP );
        tExpectSRK( k ) = -( tV2 - tV1 ) / ( 0.002 *tV * tP );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 1e-6 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 1e-6 );

//------------------------------------------------------------------------------
// Test alpha
//------------------------------------------------------------------------------

    for( uint k=0; k<tNumSamples; ++k )
    {
        real tV =   tPR.v( tT(k), tP );

        real tV1 =  tPR.v( tT(k)*0.999, tP );
        real tV2 =  tPR.v( tT(k)*1.001, tP );

        tValuesPR( k ) = tPR.alpha( tT(k), tP );
        tExpectPR( k ) = ( tV2 - tV1 ) / ( 0.002 *tV * tT( k ) );

        tV =   tSRK.v( tT(k), tP );

        tV1 =  tSRK.v( tT(k)*0.999, tP );
        tV2 =  tSRK.v( tT(k)*1.001, tP );

        tValuesSRK( k ) = tSRK.alpha( tT(k), tP );
        tExpectSRK( k ) = ( tV2 - tV1 ) / ( 0.002 *tV * tT( k ) );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 1e-6 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 1e-6 );

//------------------------------------------------------------------------------
}