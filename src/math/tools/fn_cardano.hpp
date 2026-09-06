/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_CARDANO_HPP
#define BELFEM_FN_CARDANO_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "fn_sort.hpp"
#include "fn_sign.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    /**
     * solve a cubic equation
     *     a*x^3 + b*x^2 + c*x + d = 0
     * for its real roots (Cardano / trigonometric method).
     *
     * Returns the distinct real roots in X, sorted ascending: 1 root for
     * positive discriminant, 3 for the casus irreducibilis, 1 or 2 for
     * repeated roots, and 0-2 via the quadratic/linear cascade when the
     * leading coefficients vanish. Complex roots are never returned.
     *
     * Caveats (audited 2026-07-23, Claude Fable + Codex, all live branches
     * verified correct in exact arithmetic):
     * - the epsilon tests on the coefficients and the discriminant are
     *   ABSOLUTE; callers should scale the polynomial to O(1) coefficients
     *   (the EoS and material solvers do).
     * - the exact D == 0.0 repeated-root test only fires for exactly
     *   representable cases; NEAR-multiple roots take the D > 0 or D < 0
     *   branch and lose roughly half the working digits to cancellation —
     *   inherent to closed-form solvers. In the D < 0 branch a near-double
     *   root can push the acos argument marginally outside [-1, 1] through
     *   roundoff; callers needing guaranteed behavior at multiple roots
     *   should prefer a companion-matrix eigensolve.
     * - the quadratic fallback uses the naive formula and loses precision
     *   on the small root when c^2 >> |4*b*d|.
     */

    // forward declaration of the scalar interface, so that the vector
    // wrapper below finds it by ordinary lookup, not just ADL
    template < typename T >
    void
    cardano( const T a, const T b, const T c, const T d, Vector< T > & X );

    template < typename T >
    void
    cardano( const Vector< T > & A, Vector< T > & X )
    {
        BELFEM_ASSERT( A.length() == 4, "coefficient vector needs to have a length of 4");
        cardano( A( 0 ), A( 1 ), A( 2 ), A( 3 ), X );
    }

    template < typename T >
    void
    cardano( const T a, const T b, const T c, const T d, Vector< T > & X )
    {
        if ( std::abs( a ) < BELFEM_EPSILON )
        {
            if( std::abs( b ) < BELFEM_EPSILON )
            {
                if( std::abs( c ) > BELFEM_EPSILON )
                {
                    X.set_size( 1 );

                    X = -d / c;
                }
                else
                {
                    X.set_size( 0 );
                }
            }
            else
            {
                T D = std::pow( c, 2 ) - 4.0 * b * d;

                if( D > 0.0 )
                {
                    // two distinct T roots
                    D = std::sqrt( D );
                    X.set_size( 2 );
                    X( 0 ) = ( -c - D ) / ( 2.0 * b );
                    X( 1 ) = ( -c + D ) / ( 2.0 * b );
                    sort( X );
                }
                else if( D == 0.0 )
                {
                    // one distinct root (double root)
                    X.set_size( 1 );
                    X( 0 ) = -c / ( 2.0 * b );
                }
                else
                {
                    X.set_size( 0 );
                }
            }
        }
        else
        {

            T p = ( 9.0*a*c-3.0*std::pow( b,2) )/(9.0*std::pow( a, 2 ));
            T q = ( 2.0* std::pow( b, 3) - 9.0*a*b*c + 27.0*std::pow(a,2)*d)/(27.0*std::pow( a,3 ));
            T r  = b/(3*a);
            T D = 0.25 * std::pow( q, 2 ) + std::pow( p, 3)/27.0;

            if ( D > 0.0 )
            {
                X.set_size( 1 );
                T u = -0.5*q + std::sqrt( D );
                T v = -0.5*q - std::sqrt( D );
                u = sign( u )*std::pow( std::abs( u ), 1.0/3.0 );
                v = sign( v )*std::pow( std::abs( v ), 1.0/3.0 );

                // note: u and v are real scalars, so x is real by
                // construction and the first branch below is always taken;
                // the complex fallbacks are unreachable defensive code
                // (audit 2026-07-23), kept until a cleanup pass
                std::complex< T > x = u + v;

                if( std::abs( std::imag( x ) ) < BELFEM_EPSILON )
                {
                    X( 0 ) = std::real( x ) - r;
                }
                else
                {
                    const std::complex< T > f1( -0.5,  0.5 * std::sqrt( 3.0 ));
                    const std::complex< T > f2( -0.5, -0.5 * std::sqrt( 3.0 ));

                    x = f1 * u + f2 * v;
                    if( std::abs( std::imag( x ) ) < BELFEM_EPSILON )
                    {
                        X( 0 ) = std::real( x ) - r;
                    }
                    else
                    {
                        x = f2*u + f1 * v;
                        if( std::abs( std::imag( x ) ) < BELFEM_EPSILON )
                        {
                            X( 0 ) = std::real( x ) - r;
                        }
                        else
                        {
                            BELFEM_ERROR( false, "Something went wrong while trying to solve cubic equation" );
                        }
                    }
                }
            }
            else if( D < 0 )
            {
                std::complex< T > u = std::sqrt( -4.0/3.0*p );
                std::complex< T > v = std::acos( -0.5*q*std::sqrt( -27.0/std::pow( p, 3 ) ) )/3.0;
                const std::complex< T > w = 2.0*std::acos( 0.0 )/3.0;

                X.set_size( 3 );
                X( 0 ) = std::real(  u * std::cos( v ) ) -r;
                X( 1 ) = std::real( -u * std::cos( v + w ) ) -r;
                X( 2 ) = std::real( -u * std::cos( v - w ) ) -r;
                sort( X );
            }
            else
            {
                // D == 0: repeated roots
                if( std::abs( p ) < BELFEM_EPSILON )
                {
                    // p = q = 0 → triple root
                    X.set_size( 1 );
                    X( 0 ) = -r;
                }
                else
                {
                    // double root + simple root — return distinct roots only
                    X.set_size( 2 );
                    X( 0 ) = -3.0 * q / ( 2.0 * p ) - r;   // double root
                    X( 1 ) =  3.0 * q / p - r;              // simple root
                    sort( X );
                }
            }
        }

    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_CARDANO_HPP
