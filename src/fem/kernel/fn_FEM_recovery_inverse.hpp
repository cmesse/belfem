/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_FEM_RECOVERY_INVERSE_HPP
#define BELFEM_FN_FEM_RECOVERY_INVERSE_HPP

#include <algorithm>
#include <cmath>
#include <limits>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_dot.hpp"
#include "fn_gesvd.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * Inverse of the Gram matrix of a patch-recovery fit whose basis is one
         * in slot 0 at the node ( postprocessor_recovery_theory.md §3.2 ).
         *
         * The constant is eliminated exactly through the Schur complement of
         * G( 0, 0 ), the patch weight; only the gradient block
         * S = Grr - gr0 g0r / G00 is inverted by a truncated SVD, dropping
         * singular values at or below sigma_max / aMaxCondition and at or below the
         * round-off of the subtraction that forms S, a multiple of machine
         * epsilon times the largest diagonal entry of Grr. A constant field is
         * therefore always reproduced. The result is the block inverse
         *
         *     [ 1/G00 + a' S+ a   -( S+ a )' ]
         *     [     -S+ a             S+     ]      a = g0r / G00
         *
         * @param[in]  aGram          m x m Gram matrix, m >= 2, G( 0, 0 ) > 0
         * @param[in]  aMaxCondition  truncation limit for the gradient block
         * @param[out] aInverse       m x m, sized by the caller
         * @param[out] aCondition     sigma_max / sigma_min of the computed gradient
         *                            block ( BELFEM_REAL_MAX if sigma_min is zero );
         *                            a block of pure round-off reads as well
         *                            conditioned, the return value is the test
         * @param      aSchur, aCoupling, aSchurCoupling, aSvdU, aSvdVT, aSvdS,
         *             aSvdWork       caller-owned scratch; aSchur is ( m-1 ) x ( m-1 ),
         *                            the vectors have length m-1, the SVD
         *                            containers are sized by gesvd
         * @return true if a direction of the gradient block was dropped
         */
        inline bool
        recovery_block_inverse(
                const Matrix< real > & aGram,
                const real             aMaxCondition,
                Matrix< real >       & aInverse,
                real                 & aCondition,
                Matrix< real >       & aSchur,
                Vector< real >       & aCoupling,
                Vector< real >       & aSchurCoupling,
                Matrix< real >       & aSvdU,
                Matrix< real >       & aSvdVT,
                Vector< real >       & aSvdS,
                Vector< real >       & aSvdWork )
        {
            const uint tNumCoefficients = aGram.n_rows() ;
            const uint tNumGradients = tNumCoefficients - 1 ;
            const real tWeightSum = aGram( 0, 0 );

            BELFEM_ASSERT( tNumCoefficients >= 2 && aGram.n_cols() == tNumCoefficients,
                "recovery Gram matrix must be square with at least two coefficients" );
            BELFEM_ASSERT( tWeightSum > 0.0, "recovery Gram matrix with patch weight %g", tWeightSum );
            BELFEM_ASSERT( aMaxCondition > 0.0, "recovery truncation limit must be positive, is %g", aMaxCondition );
            BELFEM_ASSERT( aSchur.n_rows() == tNumGradients && aSchur.n_cols() == tNumGradients
                && aCoupling.length() == tNumGradients && aSchurCoupling.length() == tNumGradients
                && aInverse.n_rows() == tNumCoefficients && aInverse.n_cols() == tNumCoefficients,
                "recovery_block_inverse: scratch not sized for %u coefficients", tNumCoefficients );

            for ( uint r=0; r<tNumGradients; ++r )
            {
                aCoupling( r ) = aGram( 0, r+1 ) / tWeightSum ;
            }
            for ( uint c=0; c<tNumGradients; ++c )
            {
                for ( uint r=0; r<tNumGradients; ++r )
                {
                    aSchur( r, c ) = aGram( r+1, c+1 ) - aGram( r+1, 0 ) * aCoupling( c );
                }
            }

            // the SVD of a symmetric positive semidefinite block is its eigen
            // decomposition
            gesvd( aSchur, aSvdS, aSvdU, aSvdVT, aSvdWork );

            // an exactly singular S leaves round-off of the size of the
            // subtracted terms, which a cut relative to S alone would invert
            real tScale = 0.0 ;
            for ( uint r=0; r<tNumGradients; ++r )
            {
                tScale = std::max( tScale, std::abs( aGram( r+1, r+1 ) ) );
            }
            const real tCut = std::max( aSvdS( 0 ) / aMaxCondition,
                    64.0 * std::numeric_limits< real >::epsilon() * tScale );

            aCondition = aSvdS( tNumGradients - 1 ) > 0.0 ?
                    aSvdS( 0 ) / aSvdS( tNumGradients - 1 ) : BELFEM_REAL_MAX ;

            bool tTruncated = false ;

            aInverse.fill( 0.0 );
            for ( uint j=0; j<tNumGradients; ++j )
            {
                if ( aSvdS( j ) <= tCut )
                {
                    tTruncated = true ;
                    continue ;
                }
                for ( uint c=0; c<tNumGradients; ++c )
                {
                    for ( uint r=0; r<tNumGradients; ++r )
                    {
                        aInverse( r+1, c+1 ) += aSvdVT( j, r ) * aSvdU( c, j ) / aSvdS( j );
                    }
                }
            }

            for ( uint r=0; r<tNumGradients; ++r )
            {
                aSchurCoupling( r ) = 0.0 ;
                for ( uint c=0; c<tNumGradients; ++c )
                {
                    aSchurCoupling( r ) += aInverse( r+1, c+1 ) * aCoupling( c );
                }
            }
            aInverse( 0, 0 ) = 1.0 / tWeightSum + dot( aCoupling, aSchurCoupling );
            for ( uint r=0; r<tNumGradients; ++r )
            {
                aInverse( 0, r+1 ) = -aSchurCoupling( r );
                aInverse( r+1, 0 ) = -aSchurCoupling( r );
            }

            return tTruncated ;
        }
    }
}

#endif //BELFEM_FN_FEM_RECOVERY_INVERSE_HPP
