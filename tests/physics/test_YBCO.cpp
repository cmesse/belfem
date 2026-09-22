/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * YBCO material regression tests.
 */

#include <algorithm>
#include <array>
#include <cmath>

#include <string>

#include <gtest/gtest.h>

#include "constants.hpp"
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

    // CODATA 2018, written out rather than read from the kernel so that a
    // constant wrong on both sides cannot cancel
    const double tPi     = 3.141592653589793238462643383279502884;
    const double tKB     = 1.380649e-23;
    const double tPlanck = 6.62607015e-34;
    const double tLight  = 299792458.0;
    const double tZeta3  = 1.2020569031595943;

    // The kernel reduced to a single scattering rate. With the Umklapp,
    // impurity and gap channels off, tau = hbar / ( K5 kB T z ) and the
    // Callaway integral collapses to
    //   kappa = kB ( kB T )^2 J3( theta / T ) / ( 2 pi^2 v_g hbar^2 K5 ),
    //   J3( x ) = int_0^x z^3 e^z / ( e^z - 1 )^2 dz,  J3( inf ) = 6 zeta( 3 )
    // ( 10.1103/PhysRev.113.1046 ). The boundary rate v_g / L0 is not zero
    // but 3e-9 1/s against an optical rate of 4e10 1/s, and the tail beyond
    // theta / T = 60 is below 1e-20, so both sit far below the band asserted
    // against this value.
    double
    optical_channel_closed_form( const double aT,
                                 const double aVg,
                                 const double aRho,
                                 const double aNe,
                                 const double aLambdaOpt,
                                 const double aWaveNumber )
    {
        const double tHbar = 0.5 * tPlanck / tPi;

        // the Raman shift is passed in 1/cm, the phonon energy is h c nu
        const double tEnergy = 100.0 * aWaveNumber * tPlanck * tLight;

        const double tU = tEnergy / ( tKB * aT );
        const double tV = std::exp( tU );

        const double tK5 = aLambdaOpt * aNe * tKB * aT / ( aRho * aVg * aVg )
                           * tU * tU * tV / ( tV - 1.0 );

        return tKB * ( tKB * aT ) * ( tKB * aT ) * 6.0 * tZeta3
               / ( 2.0 * tPi * tPi * aVg * tHbar * tHbar * tK5 );
    }

    // the single-rate parameter set of the closed form above: every channel
    // but the optical one is switched off through its own parameter
    void
    fill_optical_only_params( double * aParams,
                              const double aT,
                              const double aLambdaOpt )
    {
        for ( uint k = 0; k < 20; ++k ) aParams[ k ] = 0.0;

        aParams[  0 ] = aT;      // T [K]
        aParams[  1 ] = 3000.0;  // theta [K], so that theta / T truncates nothing
        aParams[  2 ] = 3000.0;  // vg [m/s]
        aParams[  3 ] = 6000.0;  // rho [kg/m^3]
        aParams[  4 ] = 50e9;    // G [Pa], only a divisor of the disabled tau_U
        aParams[  5 ] = 0.0;     // gamma = 0 switches Umklapp off
        aParams[  6 ] = 0.1;     // molar mass [kg/mol], same divisor
        aParams[  7 ] = 0.0;     // Gamma = 0 switches impurity scattering off
        aParams[  8 ] = 0.0;     // T_crit = 0 skips the d-wave reduction
        aParams[  9 ] = 1e12;    // layer thickness [m], boundary rate 3e-9 1/s
        aParams[ 10 ] = 4.0;     // integrand exponent
        aParams[ 15 ] = 2e27;    // conduction electron concentration [1/m^3]
        aParams[ 18 ] = aLambdaOpt;
        aParams[ 19 ] = 501.0;   // Raman shift [1/cm]
    }

    //! YBCO with its raw law made public; the default constructor is the one
    //! the factory calls, so the probe is the object a solver would get
    struct RawYBCO : public belfem::material::YBCO
    {
        using belfem::material::YBCO::lambda_custom;
    };
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

TEST( CallawayConductivityOpticalChannel, SingleRateMatchesClosedForm )
{
    double tParams[ 20 ];

    fill_optical_only_params( tParams, 50.0, 1.06 );
    const double tKernel = callaway_phonon_lambda( tParams );
    const double tClosedForm = optical_channel_closed_form( 50.0, 3000.0, 6000.0, 2e27, 1.06, 501.0 );

    EXPECT_NEAR( tKernel, tClosedForm, 1e-7 * tClosedForm );

    // the same rate at another temperature pins U( T ), which is where the
    // unit of the Raman shift enters
    fill_optical_only_params( tParams, 20.0, 1.06 );
    const double tKernelCold = callaway_phonon_lambda( tParams );
    const double tClosedFormCold = optical_channel_closed_form( 20.0, 3000.0, 6000.0, 2e27, 1.06, 501.0 );

    EXPECT_NEAR( tKernelCold, tClosedFormCold, 1e-7 * tClosedFormCold );

    // K5 is linear in the coupling and kappa is its reciprocal, so doubling
    // the coupling halves the result; this separates the coupling from the
    // phonon energy, which the two cases above share
    fill_optical_only_params( tParams, 50.0, 2.12 );
    const double tKernelStrong = callaway_phonon_lambda( tParams );

    EXPECT_NEAR( tKernelStrong, 0.5 * tKernel, 1e-10 * tKernel );
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

TEST( YBCOThermalConductivity, LambdaSplineOriginSlopeIsLorenzOverRho0 )
{
    belfem::material::YBCO tYBCO;

    const belfem::real tRho0 = tYBCO.constant_property( belfem::MaterialProperty::rho_0 );

    // the lambda spline is clamped at T = 0 to the slope of the electronic
    // term L0 T / ( rho_i + rho_0 ), whose ideal part vanishes as T^5; a
    // Tangent start is an identity row of the slope system, so deval( 0 )
    // returns the imposed value to the roundoff of that solve
    const belfem::real tSlope = belfem::constant::L0 / tRho0;

    EXPECT_NEAR( tYBCO.dlambdadT( 0.0 ), tSlope, 1e-9 * tSlope );

    // the reciprocal is the mutation this case exists for; the two candidates
    // are a factor four apart, so the band above cannot cover both
    EXPECT_GT( std::abs( tSlope - tRho0 / belfem::constant::L0 ), 0.5 * tSlope );
}

TEST( YBCOThermalConductivity, LambdaSplineMatchesRawLawBelowFirstKnots )
{
    RawYBCO tProbe;

    // off-knot ( the knots are 4 K apart from 0 ) and above 1 K, below which
    // the raw law is not finite
    const std::array< belfem::real, 5 > tTemperatures = { 1.5, 2.0, 2.5, 3.0, 3.5 };

    belfem::real tMaxDeviation = 0.0;

    for ( const belfem::real tTemperature : tTemperatures )
    {
        const belfem::real tRaw = tProbe.lambda_custom( tTemperature );
        ASSERT_TRUE( std::isfinite( tRaw ) ) << "T = " << tTemperature;
        ASSERT_GT( tRaw, 0.0 ) << "T = " << tTemperature;

        const belfem::real tDeviation =
                std::abs( tProbe.lambda( tTemperature ) - tRaw ) / tRaw;

        tMaxDeviation = std::max( tMaxDeviation, tDeviation );
    }

    RecordProperty( "max_relative_deviation", std::to_string( tMaxDeviation ) );

    // knot resolution alone: the electronic term is nearly linear here. The
    // measured maximum is 4.8e-4 on both backends, the band is ten times
    // that; a reciprocal start slope gives about -19 % at 2 K
    EXPECT_LT( tMaxDeviation, 0.005 );
}
