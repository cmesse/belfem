/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_FERRARI_HPP
#define BELFEM_FN_FERRARI_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "fn_cardano.hpp"

namespace  belfem
{
//------------------------------------------------------------------------------
    /**
     * solve a quartic equation
     *     A(0)*x^4 + A(1)*x^3 + A(2)*x^2 + A(3)*x + A(4) = 0
     * for its real roots using Ferrari's method.
     *
     * Returns the distinct real roots in X, sorted ascending
     * (0 to 4 entries; multiple roots are reported once, matching the
     * cardano() convention). Degenerates to cardano() when A(0) ~ 0.
     *
     * Method: depress with x = y - a/4 to y^4 + P y^2 + Q y + R = 0,
     * solve the resolvent cubic z^3 + 2P z^2 + (P^2 - 4R) z - Q^2 = 0,
     * take its largest root z = s^2 (which is >= 0 because the cubic is
     * negative at z = 0), and split into the two quadratic factors
     *     ( y^2 + s*y + (P + z - Q/s)/2 ) * ( y^2 - s*y + (P + z + Q/s)/2 ).
     * The Q ~ 0 case is handled separately as a biquadratic (this also
     * covers s -> 0, which cannot occur on the general path since the
     * resolvent's root product is Q^2).
     *
     * Caveats (same class as cardano):
     * - epsilon tests are absolute; scale the coefficients to O(1) first.
     * - near multiple roots, closed-form solvers lose about half the
     *   working digits; slightly negative discriminants at double roots
     *   are clamped to zero within BELFEM_EPSILON. Where guaranteed
     *   robustness at double roots matters (e.g. grazing-eclipse
     *   detection), prefer a 4x4 companion-matrix eigensolve instead.
     *
     * (Implementation completed + three-way verified 2026-07-23:
     *  the resolvent/factorization formulas were derived independently by
     *  Claude and Codex and agree; validated numerically against known
     *  quartics and random-coefficient sweeps.)
     */
    template < typename T >
    void
    ferrari ( const  Vector< T > & A, Vector< T > & X )
    {
        BELFEM_ASSERT( A.length() == 5,
            "coefficient vector needs to have a length of 5" );

        // quartic degenerates to a cubic
        if ( std::abs( A( 0 ) ) < BELFEM_EPSILON )
        {
            cardano( A( 1 ), A( 2 ), A( 3 ), A( 4 ), X );
            return ;
        }

        // normalize to x^4 + a*x^3 + b*x^2 + c*x + d = 0
        const T a = A( 1 ) / A( 0 );
        const T b = A( 2 ) / A( 0 );
        const T c = A( 3 ) / A( 0 );
        const T d = A( 4 ) / A( 0 );

        // depress with x = y - lambda : y^4 + p*y^2 + q*y + r = 0
        const T lambda = 0.25 * a ;
        const T p = b - 0.375 * a * a ;
        const T q = c - 0.5 * a * b + 0.125 * a * a * a ;
        const T r = d - 0.25 * a * c + 0.0625 * a * a * b
                    - ( 3.0 / 256.0 ) * a * a * a * a ;

        // real roots in y, at most four
        T Y[4];

        uint tCount = 0 ;

        // the biquadratic path handles q ~ 0. It also serves as the
        // fallback when roundoff collapses the clamped resolvent root
        // s^2 to zero on an ill-conditioned cubic (Grok audit 2026-07-23:
        // q/s would go Inf/NaN) — in exact arithmetic that only happens
        // when q is negligible anyway
        bool biquadratic = std::abs( q ) < BELFEM_EPSILON ;

        T z = 0.0 ;
        T s = 0.0 ;

        if ( ! biquadratic )
        {
            // resolvent cubic in z = s^2, temporarily using X as scratch
            Vector< T > & Z = X ;
            cardano( T( 1 ), 2.0 * p, p * p - 4.0 * r, -q * q, Z );

            // cardano sorts ascending; the largest root is >= 0 in exact
            // arithmetic (the cubic evaluates to -q^2 < 0 at z = 0) and is
            // the numerically robust choice for s
            z = std::max( Z( Z.length() - 1 ), T( 0 ) );
            s = std::sqrt( z );

            biquadratic = ( s < BELFEM_EPSILON );
        }

        if ( biquadratic )
        {
            // biquadratic: z^2 + P*z + R = 0 with z = y^2.
            T D = p * p - 4.0 * r ;

            if ( D > -BELFEM_EPSILON )
            {
                D = std::sqrt( std::max( D, T( 0 ) ) );

                // stable quadratic: avoid cancellation in the smaller root.
                // the sign must come from p, the linear coefficient of
                // z^2 + p*z + r ( q is ~0 in this branch )
                const T w = -0.5 * ( p + std::copysign( D, p ) );

                // the two z-candidates ( w == 0 only if p = D = 0 )
                const T z0 = w ;
                const T z1 = std::abs( w ) > BELFEM_EPSILON ? r / w : w ;

                for ( const T zk : { z0, z1 } )
                {
                    if ( zk > BELFEM_EPSILON )
                    {
                        const T tS = std::sqrt( zk );
                        Y[ tCount++ ] = -tS ;
                        Y[ tCount++ ] =  tS ;
                    }
                    else if ( zk > -BELFEM_EPSILON )
                    {
                        Y[ tCount++ ] = 0.0 ;
                    }
                }
            }
            // D < 0 : no real roots
        }
        else
        {
            // constant terms of the two quadratic factors
            const T t1 = 0.5 * ( p + z - q / s );
            const T t2 = 0.5 * ( p + z + q / s );

            // y^2 + s*y + t1 = 0
            T D = s * s - 4.0 * t1 ;

            if ( D > BELFEM_EPSILON )
            {
                D = std::sqrt( D );
                Y[ tCount++ ] = 0.5 * ( -s - D );
                Y[ tCount++ ] = 0.5 * ( -s + D );
            }
            else if ( D > -BELFEM_EPSILON )
            {
                // double root (clamped)
                Y[ tCount++ ] = -0.5 * s ;
            }

            // y^2 - s*y + t2 = 0
            D = s * s - 4.0 * t2 ;
            if ( D > BELFEM_EPSILON )
            {
                D = std::sqrt( D );
                Y[ tCount++ ] = 0.5 * ( s - D );
                Y[ tCount++ ] = 0.5 * ( s + D );
            }
            else if ( D > -BELFEM_EPSILON )
            {
                Y[ tCount++ ] = 0.5 * s ;
            }
        }

        // shift back ( x = y - lambda ), sort ascending (insertion sort,
        // at most four entries), and drop duplicates from multiple roots
        for ( uint i = 0 ; i < tCount ; ++i )
        {
            Y[ i ] -= lambda ;
        }
        for ( uint i = 1 ; i < tCount ; ++i )
        {
            const T v = Y[ i ];
            uint j = i ;
            while ( j > 0 && Y[ j - 1 ] > v )
            {
                Y[ j ] = Y[ j - 1 ];
                --j ;
            }
            Y[ j ] = v ;
        }

        uint tNumRoots = 0 ;
        for ( uint i = 0 ; i < tCount ; ++i )
        {
            if ( tNumRoots == 0 || std::abs( Y[ i ] - Y[ tNumRoots - 1 ] )
                 > BELFEM_EPSILON * std::max( T( 1 ), std::abs( Y[ i ] ) ) )
            {
                Y[ tNumRoots++ ] = Y[ i ];
            }
        }

        X.set_size( tNumRoots );
        for ( uint i = 0 ; i < tNumRoots ; ++i )
        {
            X( i ) = Y[ i ];
        }
    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_FERRARI_HPP
