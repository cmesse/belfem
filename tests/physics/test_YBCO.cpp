/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * YBCO material regression tests.
 */

#include <algorithm>
#include <array>
#include <cmath>

#include <gtest/gtest.h>

#include "cl_Material_YBCO.hpp"
#include "debye.hpp"

namespace
{
    void
    expect_finite_positive( const belfem::real aTemperature,
                            const belfem::real aLambda )
    {
        EXPECT_TRUE( std::isfinite( aLambda ) ) << "T = " << aTemperature;
        EXPECT_GT( aLambda, 0.0 ) << "T = " << aTemperature;
    }

    double
    callaway_phonon_lambda( const double * aParams )
    {
        double tPhononLambda = 0.0;
        int tStatus = 0;
        callaway_conductivity( aParams, &tPhononLambda, &tStatus );
        EXPECT_EQ( tStatus, 0 );
        return tPhononLambda;
    }
}

TEST( CallawayConductivityParameterLayout, ReadsDFromParams13NotOmegaOptSlot )
{
    double tParams[ 20 ] = {};

    tParams[  0 ] = 100.0;  // T [K]
    tParams[  1 ] = 300.0;  // theta [K]
    tParams[  2 ] = 3000.0; // vg [m/s]
    tParams[  3 ] = 6000.0; // rho [kg/m^3]
    tParams[  4 ] = 50e9;   // G [Pa]
    tParams[  5 ] = 2.0;    // Grueneisen parameter
    tParams[  6 ] = 0.1;    // molar mass [kg/mol]
    tParams[  7 ] = 0.0;    // no impurity scattering
    tParams[  8 ] = 0.0;    // no superconducting reduction
    tParams[  9 ] = 1e9;    // suppress boundary scattering
    tParams[ 10 ] = 4.0;    // integrand exponent
    tParams[ 11 ] = 2.0;    // b, Umklapp parameter
    tParams[ 12 ] = 0.0;    // d, Umklapp temperature correction
    tParams[ 13 ] = 0.0;    // deformation potential [eV]
    tParams[ 14 ] = 0.0;    // effective electron mass factor
    tParams[ 15 ] = 0.0;    // no electron-phonon fallback
    tParams[ 16 ] = 0.0;    // no gap
    tParams[ 17 ] = 0.0;    // acoustic e-ph off
    tParams[ 18 ] = 0.0;    // optical e-ph off
    tParams[ 19 ] = 501.0;  // omega_opt slot, intentionally unused here

    const double tLowD = callaway_phonon_lambda( tParams );

    tParams[ 12 ] = 10.0;
    const double tHighD = callaway_phonon_lambda( tParams );

    EXPECT_TRUE( std::isfinite( tLowD ) );
    EXPECT_TRUE( std::isfinite( tHighD ) );
    EXPECT_GT( tLowD, 0.0 );
    EXPECT_GT( tHighD, 0.0 );
    EXPECT_GT( std::abs( tLowD - tHighD ), 0.25 * tHighD );

    tParams[ 12 ] = 0.64;
    tParams[ 19 ] = 100.0;
    const double tOmegaLow = callaway_phonon_lambda( tParams );

    tParams[ 19 ] = 1000.0;
    const double tOmegaHigh = callaway_phonon_lambda( tParams );

    EXPECT_DOUBLE_EQ( tOmegaLow, tOmegaHigh );
}

TEST( YBCOThermalConductivity, RepresentativeTemperaturesAreFinitePositive )
{
    belfem::material::YBCO tYBCO;

    const std::array< belfem::real, 4 > tTemperatures = { 20.0, 77.0, 100.0, 300.0 };

    // Golden values of the built class. Re-baked 2026-08-25 after the optical
    // channel of debye.f90 was repaired, alpha( T ) below the split was tied to
    // cp ( which fixed the Grueneisen input of the kernel ), and the Callaway
    // parameters were refitted; the previous baseline was
    // { 37.912336735507921, 17.555180059713223, 14.441408654447251, 10.819449811287697 }.
    const std::array< belfem::real, 4 > tExpected =
            { 38.081989818097284, 17.744274544991541, 14.598387629409466, 11.010842982887585 };

    for ( uint k = 0; k < tTemperatures.size(); ++k )
    {
        const belfem::real tTemperature = tTemperatures[ k ];
        const belfem::real tLambda = tYBCO.lambda( tTemperature );
        expect_finite_positive( tTemperature, tLambda );
        EXPECT_NEAR( tExpected[ k ], tLambda, 1e-9 * std::max( 1.0, tExpected[ k ] ) );
    }
}

TEST( YBCOThermalConductivity, CurveIsFinitePositiveAndHasNoLargeSpikes )
{
    belfem::material::YBCO tYBCO;

    belfem::real tPrevious = tYBCO.lambda( 1.0 );
    expect_finite_positive( 1.0, tPrevious );

    for ( uint tK = 2; tK <= 400; ++tK )
    {
        const belfem::real tTemperature = static_cast< belfem::real >( tK );
        const belfem::real tLambda = tYBCO.lambda( tTemperature );
        expect_finite_positive( tTemperature, tLambda );

        if ( tK >= 20 )
        {
            const belfem::real tScale = std::max( std::abs( tLambda ), std::abs( tPrevious ) );
            EXPECT_LT( std::abs( tLambda - tPrevious ) / tScale, 0.08 )
                    << "T = " << tTemperature;
        }

        tPrevious = tLambda;
    }
}

TEST( YBCOThermalConductivity, CriticalRegionIsSmooth )
{
    belfem::material::YBCO tYBCO;

    const std::array< belfem::real, 7 > tTemperatures =
            { 90.0, 92.0, 92.5, 93.0, 95.0, 100.0, 105.0 };

    belfem::real tPrevious = 0.0;
    bool tHavePrevious = false;

    for ( const belfem::real tTemperature : tTemperatures )
    {
        const belfem::real tLambda = tYBCO.lambda( tTemperature );
        expect_finite_positive( tTemperature, tLambda );

        if ( tHavePrevious )
        {
            const belfem::real tScale = std::max( std::abs( tLambda ), std::abs( tPrevious ) );
            EXPECT_LT( std::abs( tLambda - tPrevious ) / tScale, 0.08 )
                    << "T = " << tTemperature;
        }

        tPrevious = tLambda;
        tHavePrevious = true;
    }
}
