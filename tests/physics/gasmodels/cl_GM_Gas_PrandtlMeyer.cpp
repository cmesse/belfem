//
// Created by Christian Messe on 14.10.19.
// exact characteristic method 08.08.26.
//
// validation of the exact characteristic Prandtl-Meyer method against the
// perfect gas closed form, published tables, the second law, path
// independence and the conservation invariants.
//

#ifndef BELFEM_CL_GM_GAS_PRANDTLMEYER_HPP
#define BELFEM_CL_GM_GAS_PRANDTLMEYER_HPP

#include <gtest/gtest.h>

#include <cmath>

#include "typedefs.hpp"
#include "constants.hpp"

#include "cl_Communicator.hpp"
#include "cl_Gas.hpp"

using namespace belfem;
using namespace belfem::gastables;

//------------------------------------------------------------------------------

// closed form nu( Ma ) of the calorically perfect gas
static real
nu_perfect( const real Ma, const real k )
{
    const real G = std::sqrt( ( k + 1. ) / ( k - 1. ) );
    const real L = std::sqrt( Ma * Ma - 1. );

    return G * std::atan( L / G ) - std::atan( L );
}

//------------------------------------------------------------------------------

// cold air case, where the gas is nearly calorically perfect and the
// classic handbook values apply
TEST( GASMODELS, PrandtlMeyerColdCase )
{
    Gas tAir;

    const real tT1    = 216.65;
    const real tP1    = 19395.0;
    const real tU1    = 2.0 * tAir.c( tT1, tP1 );
    const real tAlpha = 20.0 * constant::deg;

    real tT2;
    real tP2;
    real tU2;
    const real tMa2 = tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha,
                                          tT2, tP2, tU2 );

    EXPECT_NEAR( tMa2, 2.83, 0.01 );
    EXPECT_NEAR( tP2 / tP1, 0.271, 0.01 );
    EXPECT_NEAR( tT2 / tT1, 0.689, 0.01 );
}

//------------------------------------------------------------------------------

// the exact method must reproduce the closed form solution when the gas is
// calorically perfect. argon is monatomic, so the closed form turning angle to
// a target Mach number is exact and the solver has to land on the perfect gas
// isentrope state.
//
// the window is narrow, though: only the tabulated interval from 250 K to
// 1000 K carries cp = 2.5 R exactly. below 250 K the NASA table hands over to a
// fitted low temperature interval whose polynomial overshoots 2.5 R by 7e-6 at
// 200 K, and above 1000 K the next interval carries a3 = 2.500069401 with small
// compensating powers of T. the expansion cools the gas, so the upstream
// temperature and the largest target Mach number are chosen such that the whole
// isentrope stays inside that window
TEST( GASMODELS, PrandtlMeyerPerfectGasLimit )
{
    Gas tArgon( "Ar" );

    const real tT1 = 950.0;
    const real tP1 = 1.0e5;

    const real tMaTargets[ 6 ] = { 1.05, 1.2, 1.5, 2.0, 2.5, 3.0 };

    /* the premise: argon must actually be calorically perfect over the whole
     * sweep, coldest state ( Ma = 3 ) to upstream state. the caloric data is a
     * spline through the tabulated enthalpy, so cp is a spline derivative and
     * carries about 1e-9 of interpolation noise even where the underlying
     * polynomial is constant */
    ASSERT_NEAR( tArgon.cp( 318.0, tP1 ) / tArgon.cp( tT1, tP1 ),
                 1.0, 1.0e-8 );

    const real tK   = tArgon.gamma( tT1, tP1 );
    const real tMa1 = 1.01;
    const real tU1  = tMa1 * tArgon.c( tT1, tP1 );

    const real tNu1 = nu_perfect( tMa1, tK );

    for( uint i = 0; i < 6; ++i )
    {
        const real tMaRef = tMaTargets[ i ];

        // exact perfect gas turning angle to the target Mach number
        const real tAlpha = nu_perfect( tMaRef, tK ) - tNu1;

        real tT2;
        real tP2;
        real tU2;
        const real tMa2 = tArgon.prandtl_meyer( tT1, tP1, tU1, tAlpha,
                                                tT2, tP2, tU2 );

        /* the tolerance follows the caloric noise above: the closed form and
         * the integral see cp values that differ in the last few digits, and
         * near the sonic point the inverse map nu -> T amplifies that */
        EXPECT_NEAR( tMa2, tMaRef, 1.0e-6 );

        // perfect gas isentrope temperature
        const real tTRef = tT1
                * ( 1.0 + 0.5 * ( tK - 1. ) * tMa1 * tMa1 )
                / ( 1.0 + 0.5 * ( tK - 1. ) * tMaRef * tMaRef );

        EXPECT_NEAR( tT2 / tTRef, 1.0, 1.0e-6 );
    }
}

//------------------------------------------------------------------------------

// second law test: expanding and compressing back by the same angle must
// recover the upstream state
TEST( GASMODELS, PrandtlMeyerRoundTrip )
{
    Gas tAir;

    const real tT1    = 1400.0;
    const real tP1    = 1.0e5;
    const real tU1    = 2.0 * tAir.c( tT1, tP1 );
    const real tAlpha = 30.0 * constant::deg;

    real tT2;
    real tP2;
    real tU2;
    tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha, tT2, tP2, tU2 );

    real tT3;
    real tP3;
    real tU3;
    tAir.prandtl_meyer( tT2, tP2, tU2, -tAlpha, tT3, tP3, tU3 );

    EXPECT_NEAR( tT3 / tT1, 1.0, 1.0e-8 );
    EXPECT_NEAR( tP3 / tP1, 1.0, 1.0e-8 );
    EXPECT_NEAR( tU3 / tU1, 1.0, 1.0e-8 );
}

//------------------------------------------------------------------------------

// conservation invariants, checked numerically and not by construction
TEST( GASMODELS, PrandtlMeyerInvariants )
{
    Gas tAir;

    const real tT1    = 1550.0;
    const real tP1    = 1.0e5;
    const real tU1    = 2.0 * tAir.c( tT1, tP1 );
    const real tAlpha = 30.0 * constant::deg;

    real tT2;
    real tP2;
    real tU2;
    tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha, tT2, tP2, tU2 );

    const real tHt1 = tAir.h( tT1, tP1 ) + 0.5 * tU1 * tU1;
    const real tHt2 = tAir.h( tT2, tP2 ) + 0.5 * tU2 * tU2;

    EXPECT_NEAR( tHt2 / tHt1, 1.0, 1.0e-12 );

    EXPECT_NEAR( tAir.s( tT2, tP2 ) / tAir.s( tT1, tP1 ), 1.0, 1.0e-10 );
}

//------------------------------------------------------------------------------

// literature cross-check against NACA TN 2125 ( Noyes, R. N.:
// "Prandtl-Meyer Flow for a Diatomic Gas of Variable Specific Heat",
// NACA Technical Note 2125, Lewis Flight Propulsion Laboratory, June 1950,
// Table I ). Noyes integrates the corner flow of a diatomic gas whose
// vibration is a harmonic oscillator, cp / R = 7/2 + tau^2 e^tau /
// ( e^tau - 1 )^2 with tau = theta / t, and recommends the effective
// vibrational temperature theta = 5450 R for air. his tables start at
// Ma = 1 and carry three to four digits, and the harmonic oscillator
// deviates from the NASA fits by a few tenths of a percent, hence the
// percent level tolerances
TEST( GASMODELS, PrandtlMeyerLiterature )
{
    Gas tAir;

    const real tP1    = 1.0e5;
    const real tTheta = 5450.0 / 1.8;   // Rankine to Kelvin

    // Table I (a), T / theta = 0.3: row alpha =  5 deg to row alpha = 30 deg
    // alpha =  5: Ma = 1.250, p/P = 0.394, t/T = 0.783
    // alpha = 30: Ma = 2.092, p/P = 0.110, t/T = 0.553
    {
        const real tT     = 0.3 * tTheta;
        const real tT1    = 0.783 * tT;
        const real tU1    = 1.250 * tAir.c( tT1, tP1 );
        const real tAlpha = 25.0 * constant::deg;

        real tT2;
        real tP2;
        real tU2;
        const real tMa2 = tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha,
                                              tT2, tP2, tU2 );

        EXPECT_NEAR( tMa2, 2.092, 0.021 );
        EXPECT_NEAR( tT2 / tT1, 0.553 / 0.783, 0.007 );
        EXPECT_NEAR( tP2 / tP1, 0.110 / 0.394, 0.006 );
    }

    // Table I (b), T / theta = 0.6: row alpha =  5 deg to row alpha = 30 deg
    // alpha =  5: Ma = 1.248, p/P = 0.402, t/T = 0.806
    // alpha = 30: Ma = 2.047, p/P = 0.119, t/T = 0.601
    {
        const real tT     = 0.6 * tTheta;
        const real tT1    = 0.806 * tT;
        const real tU1    = 1.248 * tAir.c( tT1, tP1 );
        const real tAlpha = 25.0 * constant::deg;

        real tT2;
        real tP2;
        real tU2;
        const real tMa2 = tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha,
                                              tT2, tP2, tU2 );

        EXPECT_NEAR( tMa2, 2.047, 0.021 );
        EXPECT_NEAR( tT2 / tT1, 0.601 / 0.806, 0.008 );
        EXPECT_NEAR( tP2 / tP1, 0.119 / 0.402, 0.006 );
    }
}

//------------------------------------------------------------------------------

// one turn of alpha and N turns of alpha / N along the same isentrope must
// agree. the retired projection method failed this test in the hot regime,
// because the closed form nu at local gamma is not a state function
TEST( GASMODELS, PrandtlMeyerPathIndependence )
{
    Gas tAir;

    // hot case: stagnation temperature about 2500 K at Ma1 = 2
    const real tT1    = 1550.0;
    const real tP1    = 1.0e5;
    const real tU1    = 2.0 * tAir.c( tT1, tP1 );
    const real tAlpha = 30.0 * constant::deg;

    const uint tN = 30;

    // one turn of 30 deg
    real tT2;
    real tP2;
    real tU2;
    const real tMa2 = tAir.prandtl_meyer( tT1, tP1, tU1, tAlpha,
                                          tT2, tP2, tU2 );

    // thirty turns of 1 deg
    real tT = tT1;
    real tP = tP1;
    real tU = tU1;
    real tMa = 0.0;

    for( uint k = 0; k < tN; ++k )
    {
        real tTn;
        real tPn;
        real tUn;
        tMa = tAir.prandtl_meyer( tT, tP, tU, tAlpha / tN, tTn, tPn, tUn );
        tT = tTn;
        tP = tPn;
        tU = tUn;
    }

    EXPECT_NEAR( tT / tT2, 1.0, 1.0e-8 );
    EXPECT_NEAR( tP / tP2, 1.0, 1.0e-8 );
    EXPECT_NEAR( tU / tU2, 1.0, 1.0e-8 );
    EXPECT_NEAR( tMa / tMa2, 1.0, 1.0e-8 );
}

#endif //BELFEM_CL_GM_GAS_PRANDTLMEYER_HPP
