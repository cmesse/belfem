//
// Created by Christian Messe on 18.09.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Spline.hpp"
#include "fn_r2.hpp"
#include "fn_linspace.hpp"
#define protected public
#define private   public
#include "cl_Gas.hpp"
#include "GT_globals.hpp"
#include "en_GM_GasModel.hpp"
#include "cl_GM_EoS_Cubic.hpp"
#undef protected
#undef private


using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

/**
 * this test checks the equation of state and tests if the derivatives of p
 * are correct. The derivatives are compared against spline interpolations.
 */
TEST( GASMODELS, Cubic_Departure )
{
//------------------------------------------------------------------------------
// Gasmodels
//------------------------------------------------------------------------------

    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    Gas tPR( tSpecies, tMolarFractions, GasModel::PR );
    Gas tSRK( tSpecies, tMolarFractions, GasModel::SRK );

//------------------------------------------------------------------------------
// Literature data, calculated from  NIST Chemistry WebBook, SRD 69
//------------------------------------------------------------------------------

    Vector< real > tHDEP = { -360.2986558, -351.49876, -342.4490438,
                             -333.1520441, -323.6189744, -313.8719003,
                             -303.9457592, -293.8896137, -283.7661293,
                             -273.6483953, -263.6142392, -253.7395406,
                             -244.0925097, -234.7300875, -225.6964423,
                             -217.022931, -208.7289369, -200.823242,
                             -193.3057303, -186.1692439, -179.4014125,
                             -172.9863108, -166.9058505, -161.1408781,
                             -155.6719916, -150.4801143, -145.5468739,
                             -140.8548341, -136.3876187, -132.1299599,
                             -128.0676971, -124.1877423, -120.4780272,
                             -116.9274383, -113.5257492, -110.2635503,
                             -107.1321824, -104.1236721, -101.2306719,
                             -98.44640523, -95.76461516, -93.17951845,
                             -90.68576317, -88.2783906, -85.95280061,
                             -83.70472045, -81.53017649, -79.42546873,
                             -77.38714773, -75.41199373, -73.49699788,
                             -71.63934511, -69.83639864, -68.08568601,
                             -66.38488626, -64.73181839, -63.1244308,
                             -61.56079166, -60.03908018, -58.55757857,
                             -57.11466473, -55.70880553, -54.33855067,
                             -53.00252707, -51.69943364, -50.42803653,
                             -49.18716475, -47.97570612, -46.79260354,
                             -45.63685152, -44.50749301, -43.40361646,
                             -42.32435303, -41.26887411, -40.23638892,
                             -39.22614237, -38.23741296, -37.26951093,
                             -36.32177645, -35.39357804, -34.48431093,
                             -33.59339575, -32.72027705, -31.86442214,
                             -31.02531986, -30.20247951};
    tHDEP *= 1000;

    // number of steps
    uint tNumSamples = tHDEP.length();

    Vector< real > tT;
    tT = linspace( 200.0, 625.0, tNumSamples );
    real tP = 200e5;

//------------------------------------------------------------------------------
// Spline Help Matrices and data vectors
//------------------------------------------------------------------------------

    // help matrices for spline
    SpMatrix tHelpMatrix;
    spline::create_helpmatrix( tNumSamples, tT(1)-tT(0), tHelpMatrix );

    Vector< real > tValuesPR( tNumSamples );
    Vector< real > tExpectPR( tNumSamples );
    Vector< real > tValuesSRK( tNumSamples );
    Vector< real > tExpectSRK( tNumSamples );

//------------------------------------------------------------------------------
// Important Constants
//------------------------------------------------------------------------------

    // specific gas constant in J/(kg*K)
    const real tR = constant::Rm / tPR.data( 0 )->M();

    // critical temperature in K
    const real tTcrit = tPR.data( 0 )->T_crit();

    // critical pressure in Pa
    const real tPcrit = tPR.data( 0 )->p_crit();

    // SRK Constant
    real tB_srk = 0.086640349964958 * tR * tTcrit / tPcrit;

    // Peng Robinson Constant
    real tB_pr  = 0.077796073903888 * tR * tTcrit / tPcrit;

//------------------------------------------------------------------------------
// Enthalpy Departure test
//------------------------------------------------------------------------------
    /**
     *  this test compares HDEP to an alternative formulas that are found
     *  in standard literature ( it is unlikely that both literature and
     *  implementatuion are wrong )
     */
    EoS_Cubic * tEoS_SRK = reinterpret_cast< EoS_Cubic* >( tSRK.eos() );
    EoS_Cubic * tEoS_PR  = reinterpret_cast< EoS_Cubic* >( tPR.eos() );


    for( uint k=0; k<tNumSamples; ++k )
    {
        real tZ = tP * tPR.v( tT( k ), tP ) / ( tR * tT( k ) );

        real tB = tP * tB_pr / ( tR * tT( k ) );

        tValuesPR( k ) = tEoS_PR->hdep( tT( k ), tP );

        tExpectPR( k ) =      tR * tT( k ) * ( tZ - 1.0 ) +
                              ( tT( k ) * tEoS_PR->dadT( tT( k ) ) - tEoS_PR->a( tT( k ) ) ) *
                              std::log(  ( tZ +  2.414213562373095 * tB )
                              / ( tZ  - 0.414213562373095 * tB ) ) / ( 2.0 * std::sqrt( 2 ) * tB_pr );


        tValuesSRK( k ) = tEoS_SRK->hdep( tT( k ), tP );

        tZ = tP * tSRK.v( tT( k ), tP ) / ( tR * tT( k ) );
        tB = tP * tB_srk / ( tR * tT( k ) );

        tExpectSRK( k ) =  tR * tT( k ) * ( tZ - 1.0 ) +
                           ( tT( k ) * tEoS_SRK->dadT( tT( k ) ) - tEoS_SRK->a( tT( k ) ) ) *
                           std::log(  ( tZ + tB) / ( tZ  ) ) / tB_srk ;

    }


    // check correctness of function
    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );

    // check departure against literature ( this is less perfect )
    EXPECT_NEAR( r2( tValuesPR, tHDEP ), 1.0, 0.01 );
    EXPECT_NEAR( r2( tValuesSRK, tHDEP ), 1.0, 0.01);

//------------------------------------------------------------------------------
//  Departure of specific heat capacity
//------------------------------------------------------------------------------

    /**
     * the implemented CPDEP function is compared with splines derivatives
     *  from the HDEP function
     */
    Spline tSplineSRK( tT, tValuesSRK , tHelpMatrix );
    Spline tSplinePR( tT,  tValuesPR, tHelpMatrix );

    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tEoS_PR->cpdep( tT( k ), tP );
        tExpectPR( k ) = tSplinePR.deval( tT( k ) );

        tValuesSRK( k ) = tEoS_SRK->cpdep( tT( k ), tP );
        tExpectSRK( k ) = tSplineSRK.deval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );

//------------------------------------------------------------------------------

    comm_barrier();

//------------------------------------------------------------------------------
// Entropy departure
//------------------------------------------------------------------------------

    /**
     * like in HDEP, SDEP is compared against alternative formulas from
     * literature. Note that SDEP is the departure from an ideal gas at that
     * pressure.
     */
    for( uint k=0; k<tNumSamples; ++k )
    {
        real tZ = tP * tPR.v( tT( k ), tP ) / ( tR * tT( k ) );

        real tB = tP * tB_pr / ( tR * tT( k ) );

        tValuesPR( k ) = tEoS_PR->sdep( tT( k ), tP );

        //  tR * std::log( tZ - tB ) +
        tExpectPR( k ) = tR * std::log( tZ - tB ) +
                         tEoS_PR->dadT( tT( k ) ) * std::log(  ( tZ +  2.414213562373095 * tB )
                       / ( tZ  - 0.414213562373095 * tB ) ) / ( 2.0 * std::sqrt( 2 ) * tB_pr );


        tZ = tP * tSRK.v( tT( k ), tP ) / ( tR * tT( k ) );
        tB = tP * tB_srk / ( tR * tT( k ) );

        tValuesSRK( k ) = tEoS_SRK->sdep( tT( k ), tP );

        tExpectSRK( k ) = tR * std::log( tZ - tB ) +
                          tEoS_SRK->dadT( tT( k ) ) *
                          std::log(  ( tZ + tB) / ( tZ  ) ) / tB_srk ;
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );

//------------------------------------------------------------------------------
// Departure of Entropy departure of derivative
//------------------------------------------------------------------------------

    /**
     * again, data are compared against spline derivatives.
     */
    // reference data
    tSplineSRK.update_data( tHelpMatrix, tValuesSRK );
    tSplinePR.update_data( tHelpMatrix, tValuesPR );


    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesPR( k ) = tEoS_PR->dsdepdT( tT( k ), tP );
        tExpectPR( k ) = tSplinePR.deval( tT( k ) );
        tValuesSRK( k ) = tEoS_SRK->dsdepdT( tT( k ), tP );
        tExpectSRK( k ) = tSplineSRK.deval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );

//------------------------------------------------------------------------------
}
