//
// Created by Christian Messe on 16.09.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Communicator.hpp"
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
TEST( GASMODELS, Cubic_State )
{
//------------------------------------------------------------------------------
// Gasmodels
//------------------------------------------------------------------------------

    Cell<string> tSpecies = { "CH4" };
    Vector<real> tMolarFractions = { 1.0 };

    // create a Soave-Redlich-Kwong model for Methane
    Gas tSRK( tSpecies, tMolarFractions, GasModel::SRK );

    // create a Peng-Robinson model of Methane
    Gas tPR( tSpecies, tMolarFractions, GasModel::PR );

//------------------------------------------------------------------------------
// Reference Data
//------------------------------------------------------------------------------

/**
 * the reference data are taken from the  NIST Chemistry WebBook, SRD 69
 * for an isochoric state with v= 0.01 kg/m^3 and an isothermal case where T=300K
 */

    // number of samples
    uint tNumSamples = 87;

    // isochoric state
    Vector< real > tT;
    tT = linspace( 195.0, 625.0, tNumSamples );

    // Volume for this state
    real tIsochoricVolume = 0.01;

    // Pressure as function of T in bar
    Vector< real > tIsochoricPressures = { 48.42444675, 52.52370187, 56.55420275,
                                           60.53844565, 64.48902534, 68.41331485,
                                           72.31589534, 76.19982451, 80.06729755,
                                           83.91999543, 87.75927205, 91.58625778,
                                           95.40191929, 99.20709602, 103.0025243,
                                           106.7888541, 110.566662, 114.3364615,
                                           118.0987111, 121.8538217, 125.6021627,
                                           129.344067, 133.0798358, 136.8097424,
                                           140.5340354, 144.2529421, 147.9666707,
                                           151.6754126, 155.3793442, 159.0786288,
                                           162.773418, 166.4638527, 170.1500645,
                                           173.8321766, 177.5103043, 181.1845564,
                                           184.8550352, 188.5218375, 192.1850549,
                                           195.8447742, 199.501078, 203.1540449,
                                           206.8037499, 210.4502645, 214.0936575,
                                           217.7339942, 221.3713378, 225.0057487,
                                           228.637285, 232.2660028, 235.891956,
                                           239.5151967, 243.1357751, 246.7537398,
                                           250.3691379, 253.9820148, 257.5924146,
                                           261.20038, 264.8059526, 268.4091726,
                                           272.010079, 275.60871, 279.2051025,
                                           282.7992924, 286.3913149, 289.9812039,
                                           293.5689929, 297.1547142, 300.7383994,
                                           304.3200793, 307.8997842, 311.4775434,
                                           315.0533857, 318.6273391, 322.1994311,
                                           325.7696885, 329.3381376, 332.904804,
                                           336.4697129, 340.0328889, 343.5943561,
                                           347.1541381, 350.712258, 354.2687384,
                                           357.8236017, 361.3768695, 364.9285632 };

    // change unit to Pa
    tIsochoricPressures *= 1e5;

    // isothermal state
    real tIsothermalTemperature = 300.0;

    // specific volumes in kg/m^3
    Vector< real > tV;
    tV = linspace( 0.004, 0.3, tNumSamples );

    // pressures for this state in bar
    Vector< real > tIsothermalPressures = { 412.09315323, 171.40958114, 119.55970372,
                                            93.65533522, 77.38506956, 66.03790287,
                                            57.62020726, 51.10868905, 45.91423409,
                                            41.67057262, 38.13688575, 35.14788345,
                                            32.58615276, 30.36590303, 28.42294673,
                                            26.7082788, 25.18382548, 23.81955003,
                                            22.59143401, 21.48003924, 20.46946405,
                                            19.54657291, 18.70041938, 17.92180781,
                                            17.20295648, 16.53723568, 15.91896199,
                                            15.34323512, 14.80580747, 14.30297891,
                                            13.83151125, 13.38855832, 12.97160838,
                                            12.57843637, 12.20706411, 11.85572702,
                                            11.522846, 11.20700367, 10.90692422,
                                            10.62145613, 10.3495574, 10.09028281,
                                            9.84277299, 9.60624483, 9.37998321,
                                            9.16333382, 8.95569679, 8.75652116,
                                            8.56529997, 8.38156599, 8.20488789,
                                            8.03486682, 7.87113347, 7.71334534,
                                            7.56118437, 7.41435478, 7.27258117,
                                            7.1356068, 7.00319201, 6.87511285,
                                            6.75115979, 6.63113659, 6.51485927,
                                            6.40215515, 6.292862, 6.18682728,
                                            6.0839074, 5.98396709, 5.88687879,
                                            5.79252214, 5.70078345, 5.61155525,
                                            5.52473586, 5.44022905, 5.35794363,
                                            5.27779313, 5.19969554, 5.12357297,
                                            5.04935144, 4.97696061, 4.90633357,
                                            4.83740664, 4.77011914, 4.70441329,
                                            4.64023395, 4.57752854, 4.51624686, };

    // change unit to Pa
    tIsothermalPressures *= 1e5;

//------------------------------------------------------------------------------
// Spline Help Matrices
//------------------------------------------------------------------------------

    // help matrices for spline
    SpMatrix tHelpMatrixT;
    spline::create_helpmatrix( tNumSamples, tT(1)-tT(0), tHelpMatrixT );

    SpMatrix tHelpMatrixV;
    spline::create_helpmatrix( tNumSamples, tV(1)-tV(0), tHelpMatrixV );

//------------------------------------------------------------------------------
// Initialize solution vectors
//------------------------------------------------------------------------------

    // calculated values for SRK model
    Vector< real > tValuesSRK( tNumSamples );

    // calculated values for PR model
    Vector< real > tValuesPR( tNumSamples );

    // reference values for SRK model
    Vector< real > tExpectSRK( tNumSamples );

    // reference values for PR Model
    Vector< real > tExpectPR( tNumSamples );

//------------------------------------------------------------------------------
// Important data for the gas ( are the same for PR and SRK )
//------------------------------------------------------------------------------

    // specific gas constant in J/(kg*K)
    const real tR = constant::Rm / tSRK.data( 0 )->M();

    // critical temperature in K
    const real tTcrit = tSRK.data( 0 )->T_crit();

    // critical pressure in Pa
    const real tPcrit = tSRK.data( 0 )->p_crit();

//------------------------------------------------------------------------------
// Pressure function
//------------------------------------------------------------------------------

    /**
     * first we test if p( T, v ) has been implemented correctly
     */

    // B-Constant for SRK
    real tB_srk = 0.086640349964958 * tR * tTcrit / tPcrit;

    // B-Constant for PR
    real tB_pr  = 0.077796073903888 * tR * tTcrit / tPcrit;

    // grab the Equations of state from the model
    EoS_Cubic * tEoS_SRK = reinterpret_cast< EoS_Cubic* >( tSRK.eos() );
    EoS_Cubic * tEoS_PR  = reinterpret_cast< EoS_Cubic* >( tPR.eos() );

    for( uint k=0; k<tNumSamples; ++k )
    {
        // see 10.1016/0009-2509(72)80096-4
        tExpectSRK( k ) =
                tR * tT( k ) / ( tV( k ) - tB_srk )
                    - tEoS_SRK->a( tT( k ) ) /
                      ( tV( k ) * ( tV( k ) + tB_srk ) );

        // see 10.1021/i160057a011
        tExpectPR( k ) = tR * tT( k ) / ( tV( k ) - tB_pr )
                    - tEoS_PR->a( tT( k ) ) /
                   ( tV( k ) * ( tV( k ) + tB_pr )
                    + tB_pr * ( tV( k ) - tB_pr ) );

        tValuesSRK( k ) = tSRK.p( tT( k ), tV( k ) );
        tValuesPR( k )  =  tPR.p( tT( k ), tV( k ) );

    }

    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );

//------------------------------------------------------------------------------
// Isochoric Testcase
//------------------------------------------------------------------------------

    // compare function against literature
    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesSRK( k ) = tSRK.p( tT( k ), tIsochoricVolume );
        tValuesPR( k ) = tPR.p( tT( k ), tIsochoricVolume );
    }

    EXPECT_NEAR( r2( tValuesSRK, tIsochoricPressures ), 1.0, 0.015 );
    EXPECT_NEAR( r2( tValuesPR, tIsochoricPressures ), 1.0, 0.005 );

    // create the splines to test derivatives
    Spline tSrkSplineT( tT, tValuesSRK, tHelpMatrixT );
    Spline tPrSplineT( tT, tValuesPR, tHelpMatrixT );

    // test dpdT
    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesSRK( k ) = tSRK.eos()->dpdT( tT( k ), tIsochoricVolume );
        tValuesPR( k )  = tPR.eos()->dpdT( tT( k ), tIsochoricVolume );
        tExpectSRK( k ) = tSrkSplineT.deval( tT( k ) );
        tExpectPR( k )  = tPrSplineT.deval( tT( k ) );
    }

    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );

    // test d2pdT2
    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesSRK( k ) = tSRK.eos()->d2pdT2( tT( k ), tIsochoricVolume );
        tValuesPR( k )  = tPR.eos()->d2pdT2( tT( k ), tIsochoricVolume );
        tExpectSRK( k ) = tSrkSplineT.ddeval( tT( k ) );
        tExpectPR( k )  = tPrSplineT.ddeval( tT( k ) );
    }
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.001 );
    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.001 );

//------------------------------------------------------------------------------
// Isothermal Testcase
//------------------------------------------------------------------------------

    // compare function against literature
    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesSRK( k ) = tSRK.p( tIsothermalTemperature, tV( k ) );
        tValuesPR( k ) = tPR.p( tIsothermalTemperature, tV( k ) );
    }

    EXPECT_NEAR( r2( tValuesSRK, tIsothermalPressures ), 1.0, 0.01 );
    EXPECT_NEAR( r2( tValuesPR, tIsothermalPressures ), 1.0, 0.01 );

    // create the splines to test derivatives
    Spline tSrkSplineV( tV, tValuesSRK, tHelpMatrixV );
    Spline tPrSplineV( tV, tValuesPR, tHelpMatrixV );

    // test dpdv
    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesSRK( k ) = 1.0/tSRK.eos()->dpdv( tIsothermalTemperature, tV( k ) );
        tValuesPR( k )  = 1.0/tPR.eos()->dpdv( tIsothermalTemperature, tV( k ));
        tExpectSRK( k ) = 1.0/tSrkSplineV.deval( tV( k ) );
        tExpectPR( k )  = 1.0/tPrSplineV.deval( tV( k ) );
    }

    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.000001 );
    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.000001 );

    // test d2pdv2
    for( uint k=0; k<tNumSamples; ++k )
    {
        tValuesSRK( k ) = 1.0/tSRK.eos()->d2pdv2( tIsothermalTemperature, tV( k ) );
        tValuesPR( k )  = 1.0/tPR.eos()->d2pdv2( tIsothermalTemperature, tV( k ));
        tExpectSRK( k ) = 1.0/tSrkSplineV.ddeval( tV( k ) );
        tExpectPR( k )  = 1.0/tPrSplineV.ddeval( tV( k ) );
    }
    EXPECT_NEAR( r2( tValuesSRK, tExpectSRK ), 1.0, 0.001 );
    EXPECT_NEAR( r2( tValuesPR, tExpectPR ), 1.0, 0.001 );
}
