/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the closed-form polynomial solvers fn_cardano.hpp and
 * fn_ferrari.hpp. Promoted from the 2026-07-23 three-AI audit scratch
 * harness (Claude Fable + Codex + Grok; 19k-case random sweep passed).
 *
 * Contract under test (both solvers): distinct real roots, sorted
 * ascending; multiple roots reported once; complex roots never returned.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "fn_cardano.hpp"
#include "fn_ferrari.hpp"

namespace
{
    using belfem::real;
    using belfem::Vector;

    // multiply polynomial ( coefficient list, highest power first ) by ( x - r )
    void
    mul_lin( std::vector< real > & p, const real r )
    {
        std::vector< real > q( p.size() + 1, 0.0 );
        for ( size_t i = 0; i < p.size(); ++i )
        {
            q[ i ]     += p[ i ];
            q[ i + 1 ] -= p[ i ] * r;
        }
        p = q;
    }

    // multiply polynomial by ( x^2 - 2*re*x + re^2 + im^2 ), a complex pair
    void
    mul_quad( std::vector< real > & p, const real re, const real im )
    {
        std::vector< real > q( p.size() + 2, 0.0 );
        for ( size_t i = 0; i < p.size(); ++i )
        {
            q[ i ]     += p[ i ];
            q[ i + 1 ] -= p[ i ] * 2.0 * re;
            q[ i + 2 ] += p[ i ] * ( re * re + im * im );
        }
        p = q;
    }

    // distinct sorted expected roots
    std::vector< real >
    distinct( std::vector< real > roots )
    {
        std::sort( roots.begin(), roots.end() );
        std::vector< real > out;
        for ( real r : roots )
        {
            if ( out.empty() || std::fabs( r - out.back() )
                 > 1e-6 * std::max( 1.0, std::fabs( r ) ) )
            {
                out.push_back( r );
            }
        }
        return out;
    }

    // run ferrari (5 coeffs) or cardano (4 coeffs) and compare root sets
    void
    check( const std::vector< real > & coeffs,
           const std::vector< real > & expected,
           const real tol )
    {
        Vector< real > A( coeffs.size() );
        for ( size_t i = 0; i < coeffs.size(); ++i )
        {
            A( i ) = coeffs[ i ];
        }
        Vector< real > X;

        if ( coeffs.size() == 5 )
        {
            belfem::ferrari( A, X );
        }
        else
        {
            belfem::cardano( A, X );
        }

        const std::vector< real > exp2 = distinct( expected );

        ASSERT_EQ( X.length(), exp2.size() );
        for ( size_t i = 0; i < exp2.size(); ++i )
        {
            EXPECT_NEAR( X( i ), exp2[ i ],
                         tol * std::max( 1.0, std::fabs( exp2[ i ] ) ) );
        }
    }
}

//------------------------------------------------------------------------------
// cardano
//------------------------------------------------------------------------------

TEST( Cardano, ThreeDistinctRoots )
{
    // (x-1)(x-2)(x-3)
    check( { 1.0, -6.0, 11.0, -6.0 }, { 1.0, 2.0, 3.0 }, 1e-8 );
}

TEST( Cardano, DoubleRoot )
{
    // (x-1)^2 (x+2) = x^3 - 3x + 2  ( D == 0 exactly representable )
    check( { 1.0, 0.0, -3.0, 2.0 }, { -2.0, 1.0 }, 1e-6 );
}

TEST( Cardano, OneRealRoot )
{
    // x^3 + x + 1, real root ~ -0.6823278038280193
    check( { 1.0, 0.0, 1.0, 1.0 }, { -0.6823278038280193 }, 1e-9 );
}

TEST( Cardano, DegenerateQuadraticLinear )
{
    check( { 0.0, 1.0, -3.0, 2.0 }, { 1.0, 2.0 }, 1e-9 );  // quadratic
    check( { 0.0, 0.0, 2.0, -4.0 }, { 2.0 }, 1e-12 );      // linear
    check( { 0.0, 1.0, 0.0, 1.0 }, {}, 0.0 );              // no real roots
}

//------------------------------------------------------------------------------
// ferrari
//------------------------------------------------------------------------------

TEST( Ferrari, Biquadratic )
{
    // x^4 - 5x^2 + 4 = (x^2-1)(x^2-4)
    check( { 1.0, 0.0, -5.0, 0.0, 4.0 }, { -2.0, -1.0, 1.0, 2.0 }, 1e-9 );
}

TEST( Ferrari, FourDistinctRoots )
{
    // (x-1)(x-2)(x-3)(x-4)
    check( { 1.0, -10.0, 35.0, -50.0, 24.0 }, { 1.0, 2.0, 3.0, 4.0 }, 1e-8 );
}

TEST( Ferrari, NoRealRoots )
{
    check( { 1.0, 0.0, 5.0, 0.0, 4.0 }, {}, 0.0 );   // (x^2+1)(x^2+4)
    check( { 1.0, 0.0, 0.0, 0.0, 1.0 }, {}, 0.0 );   // x^4 + 1
}

TEST( Ferrari, QuadrupleRoot )
{
    // (x-2)^4 — closed form loses ~half the digits at multiplicity
    check( { 1.0, -8.0, 24.0, -32.0, 16.0 }, { 2.0 }, 1e-3 );
}

TEST( Ferrari, DoubleRootPlusPair )
{
    // (x-1)^2 (x+2)(x+3)
    std::vector< real > p{ 1.0 };
    mul_lin( p, 1.0 );
    mul_lin( p, 1.0 );
    mul_lin( p, -2.0 );
    mul_lin( p, -3.0 );
    check( p, { -3.0, -2.0, 1.0 }, 1e-5 );
}

TEST( Ferrari, ScaledLeadingCoefficient )
{
    std::vector< real > p{ 2.5 };
    mul_lin( p, -1.5 );
    mul_lin( p, 0.25 );
    mul_lin( p, 0.5 );
    mul_lin( p, 3.0 );
    check( p, { -1.5, 0.25, 0.5, 3.0 }, 1e-8 );
}

TEST( Ferrari, DegeneratesToCubic )
{
    check( { 0.0, 1.0, -6.0, 11.0, -6.0 }, { 1.0, 2.0, 3.0 }, 1e-8 );
}

TEST( Ferrari, RandomSweep )
{
    // random quartics with known factorizations: four real roots,
    // two real + complex pair, or two complex pairs. Near-degenerate
    // root clusters are skipped ( closed-form accuracy limit, documented
    // in the fn_ferrari.hpp caveats ).
    std::mt19937 gen( 7 );
    std::uniform_real_distribution< real > dr( -3.0, 3.0 );
    std::uniform_real_distribution< real > dp( 0.2, 3.0 );

    int cases = 0;
    for ( int k = 0; k < 3000 && cases < 1000; ++k )
    {
        std::vector< real > p{ dr( gen ) };
        if ( std::fabs( p[ 0 ] ) < 0.1 )
        {
            continue;
        }

        std::vector< real > roots;
        const int mode = k % 3;
        if ( mode == 0 )
        {
            for ( int i = 0; i < 4; ++i )
            {
                const real r = dr( gen );
                roots.push_back( r );
                mul_lin( p, r );
            }
        }
        else if ( mode == 1 )
        {
            for ( int i = 0; i < 2; ++i )
            {
                const real r = dr( gen );
                roots.push_back( r );
                mul_lin( p, r );
            }
            mul_quad( p, dr( gen ), dp( gen ) );
        }
        else
        {
            mul_quad( p, dr( gen ), dp( gen ) );
            mul_quad( p, dr( gen ), dp( gen ) );
        }

        bool skip = false;
        for ( size_t i = 0; i < roots.size(); ++i )
        {
            for ( size_t j = i + 1; j < roots.size(); ++j )
            {
                if ( std::fabs( roots[ i ] - roots[ j ] ) < 1e-3 )
                {
                    skip = true;
                }
            }
        }
        if ( skip )
        {
            continue;
        }

        ++cases;
        check( p, roots, 1e-6 );
    }
    EXPECT_GE( cases, 500 );
}
