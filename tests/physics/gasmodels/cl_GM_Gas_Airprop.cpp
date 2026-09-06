//
// Created by Christian Messe on 12.09.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Gas.hpp"
using namespace belfem;
using namespace belfem::gastables;

TEST( GASMODELS, AirProp )
{
//------------------------------------------------------------------------------
/**
 * This test creates a reference mixture of N2, O2 and Ar, and compares
 * the result with literature data from the VDI Heat Atlas
 */
//------------------------------------------------------------------------------

// VDI HEAT ATLAS D2.2
    Cell<string> tSpecies = { "N2", "O2", "Ar" };
    Vector<real> tMolarFractions = { 0.7812, 0.2096, 0.0092 };

// create the mixture as ideal gas
    Gas tAir( tSpecies, tMolarFractions );

//------------------------------------------------------------------------------
// REFERENCE DATA, VDI HEAT ATLAS, D2.2 Table 1
//------------------------------------------------------------------------------

// Reference temperature and pressure for h and s data
    const real tTref = 298.15;
    const real tPref = 1e5;

// Temperatures in Celsius
    Vector<real> tT = { -70., -60., -50., -40., -30., -20., -10., 0., 10., 20.,
                        30., 40., 50., 60., 70., 80., 90., 100., 120., 140.,
                        160., 180., 200., 250., 300., 350., 400., 450., 500.,
                        550., 600., 650., 700., 750., 800., 850., 900.,
                        950., 1000 };

// shift temperatures to K
    tT += 273.15;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// specific heat capacity in kJ/(kg*K)
    Vector<real> tCp = { 1.0068, 1.0064, 1.0061, 1.0059, 1.0058, 1.0057, 1.0058,
                         1.0059, 1.0061, 1.0064, 1.0067, 1.0071, 1.0077, 1.0082,
                         1.0089, 1.0097, 1.0105, 1.0115, 1.0136, 1.016, 1.0188,
                         1.0218, 1.0252, 1.0347, 1.0454, 1.0568, 1.0688, 1.0808,
                         1.0927, 1.1043, 1.1154, 1.126, 1.1361, 1.1455, 1.1544,
                         1.1628, 1.1706, 1.1778, 1.1846 };

// change unit to J/(kg*K)
    tCp *= 1000.0;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// enthalpy in kJ/kg
    Vector<real> tH = { -95.57, -85.51, -75.44, -65.38, -55.33, -45.27,
                        -35.21, -25.15, -15.09, -5.03, 5.04, 15.11, 25.18,
                        35.26, 45.34, 55.44, 65.54, 75.65, 95.9, 116.19, 136.54,
                        156.95, 177.42, 228.91, 280.9, 333.46, 386.6, 440.33,
                        494.67, 549.6, 605.09, 661.13, 717.68, 774.72, 832.22,
                        890.16, 948.49, 1007.2, 1066.3 };

// change unit to J/kg
    tH *= 1000.0;

// add offset
    tH += tAir.h( tTref, tPref );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// entropy in kJ/(kg*K )
    Vector<real> tS = { -0.3822, -0.3338, -0.2877, -0.2436, -0.2013, -0.1608,
                        -0.1218, -0.0843, -0.0481, -0.0132, 0.0205, 0.0532,
                        0.0849, 0.1156, 0.1454, 0.1744, 0.2026, 0.2301, 0.283,
                        0.3333, 0.3814, 0.4275, 0.4717, 0.5751, 0.67, 0.7579,
                        0.8399, 0.917, 0.9896, 1.0584, 1.1239, 1.1863, 1.2459,
                        1.3031, 1.358, 1.4107, 1.4615, 1.5106, 1.5579 };
    tS *= 1000.0;

// add entropy offset to entropy
    tS += tAir.s( tTref, tPref );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// dynamic viscosity in µPa * s
    Vector<real> tMu = { 13.511, 14.067, 14.614, 15.152, 15.68, 16.201,
                         16.714, 17.218, 17.715, 18.205, 18.689, 19.165,
                         19.635, 20.099, 20.557, 21.009, 21.455, 21.896,
                         22.763, 23.61, 24.439, 25.251, 26.046, 27.97,
                         29.811, 31.579, 33.284, 34.932, 36.53, 38.084,
                         39.597, 41.073, 42.517, 43.931, 45.317, 46.679,
                         48.018, 49.336, 50.635 };

// change unit to Pa*s
    tMu *= 1e-6;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// thermal conductivity in mW/(m*K )
    Vector<real> tLambda = { 18.766, 19.596, 20.416, 21.224, 22.023, 22.811,
                             23.59, 24.36, 25.121, 25.873, 26.618, 27.354,
                             28.082, 28.804, 29.518, 30.225, 30.925, 31.62,
                             32.989, 34.336, 35.66, 36.964, 38.248, 41.382,
                             44.417, 47.367, 50.24, 53.047, 55.795, 58.49,
                             61.139, 63.745, 66.312, 68.846, 71.348, 73.822,
                             76.271, 78.695, 81.099 };

// change unit to W/(m*K)
    tLambda *= 0.001;

//------------------------------------------------------------------------------
// Run the Tests
//------------------------------------------------------------------------------


    EXPECT_NEAR(( tAir.M( gastables::gTref, gastables::gPref ) * 1000 - 28.9583 ) / 28.9583, 0, 1e-4 );

    uint tSteps = tT.length();

    // vector with calculated values
    Vector<real> tValues( tSteps );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// test cp
    for ( uint k = 0; k < tSteps; ++k )
    {
        tValues( k ) = tAir.cp( tT( k ), tPref );
    }

    EXPECT_NEAR( r2( tValues, tCp ), 1.0, 1e-3 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// test h
    for ( uint k = 0; k < tSteps; ++k )
    {
        tValues( k ) = tAir.h( tT( k ), tPref );
    }

    EXPECT_NEAR( r2( tValues, tH ), 1.0, 1e-6 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// test s
    for ( uint k = 0; k < tSteps; ++k )
    {
        tValues( k ) = tAir.s( tT( k ), tPref );

    }

    EXPECT_NEAR( r2( tValues, tS ), 1.0, 1e-4 );


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// test mu
    for ( uint k = 0; k < tSteps; ++k )
    {
        tValues( k ) = tAir.mu( tT( k ), tPref );

    }

    EXPECT_NEAR( r2( tValues, tMu ), 1.0, 1e-3 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// test lambda
    for ( uint k = 0; k < tSteps; ++k )
    {
        tValues( k ) = tAir.lambda( tT( k ), tPref );
    }

    EXPECT_NEAR( r2( tValues, tLambda ), 1.0, 5e-3 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}