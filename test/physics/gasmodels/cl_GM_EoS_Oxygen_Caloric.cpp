//
// Created by Christian Messe on 24.08.20.
//


#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Gas.hpp"

#include "cl_GM_EoS_Oxygen.hpp"
#include "fn_GM_Helmholz_DerivTest.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

TEST( GASMODELS, Oxygen_Caloric )
{
//----------------------------------------------------------------------------

    // crate reference mixture
    Gas tRef( "O2" );

    // create EoS
    EoS_Oxygen tGas( tRef );

//----------------------------------------------------------------------------
// check correct impplementation of derivatives
//----------------------------------------------------------------------------

    Vector< real > tR2( 7 );

    deriv_test( tGas, tR2 );

    // phi0_t
    EXPECT_NEAR( tR2( 0 ), 1.0, 1e-6 );

    // phi0_tt
    EXPECT_NEAR( tR2( 1 ), 1.0, 1e-6 );

    // phir_t
    EXPECT_NEAR( tR2( 2 ), 1.0, 1e-6 );

    // phir_d
    EXPECT_NEAR( tR2( 3 ), 1.0, 1e-6 );

    // phir_tt
    EXPECT_NEAR( tR2( 4 ), 1.0, 1e-6 );

    // phir_dt
    EXPECT_NEAR( tR2( 5 ), 1.0, 1e-6 );

    // phir_dd
    EXPECT_NEAR( tR2( 6 ), 1.0, 1e-6 );

//----------------------------------------------------------------------------
// Physical Properties
//----------------------------------------------------------------------------
    // temperature in K
    Vector< real > tT = { 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.,
                          105., 110., 115., 120., 140., 160., 180., 200., 220.,
                          240., 260., 280., 300., 320., 340., 360., 380., 400.,
                          420., 440., 460., 480., 500., 520., 540., 560., 580.,
                          600., 620., 640., 660., 680., 700., 720., 740., 760.,
                          780., 800., 820., 840., 860., 880., 900., 920., 940.,
                          960., 980., 1000.0 };

//----------------------------------------------------------------------------

    // pressure in Pa
    real tP = 1e6 ;

//----------------------------------------------------------------------------

    // Density in kg/m^3
    Vector< real > tRho = {
            1304.60, 1283.20, 1261.00, 1238.40, 1215.50, 1192.20,
            1168.40, 1144.10, 1119.00, 1093.00, 1065.80, 1037.20,
            1006.70, 38.25, 30.39, 25.65, 22.33, 19.84, 17.88,
            16.30, 14.98, 13.86, 12.91, 12.08, 11.35, 10.71,
            10.14, 9.62, 9.16, 8.74, 8.36, 8.01, 7.68, 7.39,
            7.11, 6.86, 6.62, 6.40, 6.19, 6.00, 5.82, 5.64,
            5.48, 5.33, 5.19, 5.05, 4.92, 4.80, 4.68, 4.57,
            4.46, 4.36, 4.26, 4.17, 4.08, 4.00, 3.92, 3.84};

//----------------------------------------------------------------------------

    Vector< real > tH = { -191.91, -183.57, -175.20, -166.82, -158.44, -150.06,
                          -141.65, -133.21, -124.70, -116.10, -107.37, -98.48,
                          -89.35, 94.18, 116.82, 137.46, 157.22, 176.52, 195.54,
                          214.39, 233.14, 251.84, 270.52, 289.21, 307.93, 326.71,
                          345.55, 364.48, 383.50, 402.63, 421.86, 441.21, 460.68,
                          480.27, 499.97, 519.80, 539.74, 559.80, 579.98, 600.27,
                          620.66, 641.17, 661.78, 682.48, 703.29, 724.19, 745.17,
                          766.25, 787.41, 808.64, 829.96, 851.35, 872.81, 894.35,
                          915.95, 937.61, 959.33, 981.12 };
    tH *= 1000.0 ;

    // offset correction
    tH += tGas.h( tT( 0 ), tP ) -tH( 0 );

//----------------------------------------------------------------------------

    Vector< real > tS = { 2.1092, 2.2544, 2.3884, 2.5126, 2.6282, 2.7365,
                          2.8384, 2.9349, 3.0269, 3.1151, 3.2002, 3.2830,
                          3.3642, 4.8998, 5.0747, 5.2126, 5.3290, 5.4307,
                          5.5213, 5.6033, 5.6784, 5.7477, 5.8121, 5.8724,
                          5.9292, 5.9828, 6.0338, 6.0823, 6.1287, 6.1732,
                          6.2160, 6.2571, 6.2969, 6.3353, 6.3725, 6.4085,
                          6.4435, 6.4775, 6.5106, 6.5428, 6.5742, 6.6048,
                          6.6346, 6.6638, 6.6923, 6.7202, 6.7474, 6.7741,
                          6.8002, 6.8258, 6.8509, 6.8755, 6.8996, 6.9233,
                          6.9465, 6.9693, 6.9917, 7.0137 };

    tS *= 1000.0 ;
    tS += tGas.s( tT( 0 ), tP ) -tS( 0 );

//---------------------------------------------------------------------------

    Vector< real > tCv = { 1.174661, 1.088987, 1.047272, 1.018394, 0.993953,
                           0.971643, 0.950847, 0.931401, 0.913254, 0.896388,
                           0.880815, 0.866600, 0.853896, 0.735520, 0.685068,
                           0.670807, 0.663070, 0.659044, 0.657138, 0.656619,
                           0.657162, 0.658618, 0.660905, 0.663955, 0.667696,
                           0.672051, 0.676934, 0.682257, 0.687932, 0.693877,
                           0.700015, 0.706277, 0.712602, 0.718937, 0.725240,
                           0.731474, 0.737608, 0.743620, 0.749490, 0.755207,
                           0.760759, 0.766141, 0.771348, 0.776379, 0.781234,
                           0.785915, 0.790424, 0.794766, 0.798944, 0.802964,
                           0.806831, 0.810550, 0.814127, 0.817567, 0.820875,
                           0.824059, 0.827122, 0.830071 };
    tCv *= 1000.0 ;

//---------------------------------------------------------------------------

    Vector< real > tCp = { 1.668041, 1.671192, 1.675251, 1.675989, 1.676266,
                           1.678501, 1.684188, 1.694295, 1.709660, 1.731291,
                           1.760658, 1.800091, 1.853486, 1.256768, 1.065498,
                           1.005188, 0.974293, 0.956753, 0.946166, 0.939686,
                           0.935934, 0.934183, 0.934004, 0.935115, 0.937304,
                           0.940396, 0.944239, 0.948694, 0.953638, 0.958963,
                           0.964569, 0.970372, 0.976299, 0.982288, 0.988286,
                           0.994250, 1.000146, 1.005946, 1.011628, 1.017174,
                           1.022574, 1.027818, 1.032900, 1.037818, 1.042570,
                           1.047157, 1.051581, 1.055844, 1.059950, 1.063904,
                           1.067709, 1.071371, 1.074895, 1.078287, 1.081550,
                           1.084692, 1.087716, 1.090628 };

    tCp *= 1000.0 ;

//---------------------------------------------------------------------------

    Vector< real > tW = { 1129.961, 1130.192, 1104.559, 1069.349, 1030.805,
                          991.105, 950.892, 910.234, 868.973, 826.850, 783.528,
                          738.580, 691.433, 189.999, 214.235, 233.742, 250.845,
                          266.332, 280.631, 293.997, 306.593, 318.534, 329.904,
                          340.769, 351.181, 361.187, 370.826, 380.133, 389.141,
                          397.875, 406.363, 414.624, 422.679, 430.545, 438.238,
                          445.771, 453.156, 460.403, 467.522, 474.521, 481.408,
                          488.189, 494.870, 501.456, 507.953, 514.364, 520.694,
                          526.946, 533.123, 539.228, 545.265, 551.235, 557.142,
                          562.987, 568.772, 574.499, 580.171, 585.788};


//---------------------------------------------------------------------------

    // allocate data containers
    uint tNumSamples = tT.length() ;
    Vector< real > tY( tNumSamples ) ;
    Vector< real > tY0( tNumSamples );

//----------------------------------------------------------------------------

    // density
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k ) = 1.0 / tGas.v( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tY, tRho ), 1.0, 1e-3 );

//----------------------------------------------------------------------------


    // enthalpy
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k ) = tGas.h( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tY, tH ), 1.0, 1e-3 );

//----------------------------------------------------------------------------

    // entropy
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k ) = tGas.s( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tY, tS ), 1.0, 1e-3 );

//----------------------------------------------------------------------------

    // specific heat capacity at constant volume
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k ) = tGas.cv( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tY, tCv ), 1.0, 1e-2 );

//----------------------------------------------------------------------------

    // specific heat capacity at constant pressure
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k ) = tGas.cp( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tY, tCp ), 1.0, 1e-2 );

//----------------------------------------------------------------------------

    // speed of sound
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k ) = tGas.w( tT( k ), tP );
    }

    EXPECT_NEAR( r2( tY, tW ), 1.0, 1e-3 );
//----------------------------------------------------------------------------


    // entropy derivative against temperature
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k )   = tGas.dsdT( tT( k ), tP ) ;
        real tS1 = tGas.s( tT( k )*1.001, tP ) ;
        real tS0 = tGas.s( tT( k )*0.999, tP ) ;

        tY0( k  ) = ( tS1 - tS0 ) / ( 0.002 * tT( k ) );

    }
    EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-3 );

//----------------------------------------------------------------------------

    // entropy derivative against pressure
    for( uint k=0; k<tNumSamples; ++k )
    {
        tY( k )   = tGas.dsdp( tT( k ), tP ) ;
        real tS1 = tGas.s( tT( k ), tP * 1.0001) ;
        real tS0 = tGas.s( tT( k ), tP * 0.9999 ) ;
        tY0( k  ) = ( tS1 - tS0 ) / ( 0.0002 * tP );

    }
    EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-3 );

//----------------------------------------------------------------------------
}