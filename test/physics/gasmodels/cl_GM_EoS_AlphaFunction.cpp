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

#include "cl_GT_RefGasFactory.hpp"
#include "cl_GM_EoS_AlphaFunction.hpp"
#include "cl_GM_EoS_AlphaFunctionFactory.hpp"
#include "cl_Spline.hpp"

using namespace belfem;
using namespace gastables;
using namespace gasmodels;

TEST( GASMODELS, AlphaFunction )
{

//------------------------------------------------------------------------------

    // crate a factory
    RefGasFactory tRefgasFactory;

    AlphaFunctionFactory tAlphaFactory;

    // crate two reference gases
    const RefGas * tEthane = tRefgasFactory.create_refgas( "C2H6" );
    const RefGas * tHydrazine = tRefgasFactory.create_refgas( "N2H4" );

    // number of steps
    uint tN = 1001;
    real tTmax = 6000 ;


    Vector< real > tT;
    tT = linspace( 50.0, tTmax, tN );
    real tDeltaT = tT( 1 ) - tT( 0 );

    // help matrix for spline
    SpMatrix tHelpMatrix;
    spline::create_helpmatrix( tN, tDeltaT, tHelpMatrix );

    Vector< real > tValues( tN );
    Vector< real > tExpect( tN );

    // standard soave redlich kwong
    real tTcrit = tHydrazine->data()->T_crit();
    real tOmega = tHydrazine->data()->acentric();

    EXPECT_NEAR( tTcrit,  653,  5.0 );
    EXPECT_NEAR( tOmega,  0.316,  0.1 );

//------------------------------------------------------------------------------
// SRK
//------------------------------------------------------------------------------
    // soave redlich kwong
    real tM = 0.480 + ( 1.574 - 0.176*tOmega ) * tOmega;

    AlphaFunction * tSRK = tAlphaFactory.create_srk( tHydrazine->data() );

    // function test
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tSRK->alpha( tT( k ) );

        // manual calculation
        tExpect( k ) = std::pow( 1.0 + tM * ( 1.0 - std::sqrt( tT(k)/tTcrit ) ) , 2 );
    }

    Spline tSpline( tT, tExpect, tHelpMatrix );
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );

    // check first derivative
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tSRK->dalphadT( tT( k ) );

        // manual calculation
        tExpect( k ) = tSpline.deval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );


    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tSRK->d2alphadT2( tT( k ) );

        // manual calculation
        tExpect( k ) = tSpline.ddeval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.01 );

//------------------------------------------------------------------------------
// MC
//------------------------------------------------------------------------------

    AlphaFunction * tMC = tAlphaFactory.create_ccr_mc_srk( tHydrazine->data() );

    real tC1 = ( 0.1654 - 0.1094 * tOmega ) * tOmega + 0.5178;
    real tC2 = 0.3279 - 0.4291 * tOmega;
    real tC3 = 1.3506 * tOmega + 0.4866;

    // find cutoff temperature
    real tTcut = tTcrit*std::pow(1.0/tC1 + 1.0, 2 );

    // reset matrices
    tValues.fill( 0.0 );
    tExpect.fill( 0.0 );

    // function test
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tMC->alpha( tT( k ) );

        real tX = 1.0 - std::sqrt( tT( k ) / tTcrit );

        if ( tT( k ) < tTcrit )
        {
            tExpect( k ) = std::pow( 1.0 + (( tC3 * tX + tC2 ) * tX + tC1 ) * tX, 2 );
        }
        else if ( tT( k ) < tTcut )
        {
            tExpect( k ) = std::pow( 1.0 + tC1 * tX, 2 );
        }
        else
        {
            break;
        }
    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );

    // update spline data
    tSpline.update_data( tHelpMatrix, tExpect );

    // first derivative
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tMC->dalphadT( tT( k ) );

        tExpect( k ) = tSpline.deval( tT( k ) );

    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );

    // second  derivative
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tMC->d2alphadT2( tT( k ) );

        tExpect( k ) = tSpline.ddeval( tT( k ) );

    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.01 );

//------------------------------------------------------------------------------
// CCR
//------------------------------------------------------------------------------

    AlphaFunction * tCCR = tAlphaFactory.create_ccr_pr( tHydrazine->data() );

    tC1 = ( 1.3569 * tOmega + 0.9957 ) * tOmega + 0.4077;
    tC2 = ( 3.559 - 11.2986 * tOmega ) * tOmega - 0.1146;
    tC3 = ( 11.7802 * tOmega - 3.8901 ) * tOmega + 0.5033;

    // function test
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tCCR->alpha( tT( k ) );

        // manual calculation
        if(  tT( k ) < tTcrit )
        {
            real tX = 1.0 - std::sqrt( tT( k ) / tTcrit );

            tExpect( k ) = std::exp( tC1 * ( 1.0 - tT( k ) / tTcrit ) )
                           * std::pow( ( 1.0  + tX * tX * ( tC2 + tC3 * tX ) ), 2 );
        }
        else
        {
            tExpect( k ) = std::exp( tC1 * ( 1.0 - tT( k ) / tTcrit ) );
        }

    }

    // update spline data
    tSpline.update_data( tHelpMatrix, tExpect );

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );

    // check first derivative
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tCCR->dalphadT( tT( k ) );

        // calculate using spline
        tExpect( k ) = tSpline.deval( tT( k ) );

        // std::cout << tT( k ) << " " << tValues( k ) << " " << tSpline.deval( tT( k ) ) << std::endl;
    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );

    // check second derivative
    for( uint k=0; k<tN; ++k )
    {
        tValues( k ) = tCCR->d2alphadT2( tT( k ) );

        // calculate using spline
        tExpect( k ) = tSpline.ddeval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.01 );

//------------------------------------------------------------------------------
// PM
//------------------------------------------------------------------------------

    AlphaFunction * tPM = tAlphaFactory.create_pm_srk( tEthane->data() );

    tTcrit = tEthane->data()->T_crit();

    tC1 = 0.64281;
    tC2 = 0.75568;
    tC3 = 0.80351;

    // function test
    for( uint k=0; k<tN; ++k )
    {
        real tX = std::sqrt( tT( k ) / tTcrit );

        tValues( k ) = tPM->alpha( tT( k ) );

        // Equation ( 14 )
        tExpect( k ) = std::exp( 2.0 * tC1 * ( 1.0 - tX )
                                 - std::pow( tC2 * ( 1.0 - tX ), 2 )
                                 + 2.0/3.0 * std::pow( tC3 * ( 1.0 - tX ), 3 ) );

    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );

    // update spline data
    tSpline.update_data( tHelpMatrix, tExpect );

    // test first derivative
    for( uint k=0; k<tN; ++k )
    {

        tValues( k ) = tPM->dalphadT( tT( k ) );

        // Equation ( 14 )
        tExpect( k ) = tSpline.deval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.000005 );


    // test second derivative
    for( uint k=0; k<tN; ++k )
    {

        tValues( k ) = tPM->d2alphadT2( tT( k ) );

        // Equation ( 14 )
        tExpect( k ) = tSpline.ddeval( tT( k ) );
    }
    EXPECT_NEAR( r2( tValues, tExpect ), 1.0, 0.01 );

//------------------------------------------------------------------------------
//  TIDY UP
//------------------------------------------------------------------------------

    delete tSRK;
    delete tMC;
    delete tCCR;
    delete tPM;

    delete tEthane;
    delete tHydrazine;
}