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

/**
 * @file
 * @brief Coefficient of determination, R^2.
 * @ingroup grp_linalg
 *
 * @p aExact is the reference data and @p aApproximated the values being judged against
 * it; the two must have the same length. The helper returns exactly 1.0 when either the
 * residual or the total variation falls below BELFEM_EPSILON, rather than dividing two
 * near-zero numbers.
 *
 */

#ifndef BELFEM_FN_R2_HPP
#define BELFEM_FN_R2_HPP


#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "fn_sum.hpp"
#include "fn_norm.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * calculates the R2 coefficient of determination from values
     * of an evaluated function with respect to given samples.
     *
     * See also
     *
     * https://en.wikipedia.org/wiki/Coefficient_of_determination
     *
     * @param[in] aApproximated   values of evaluated function
     * @param[in] aExact          data samples or exact solution
     */
     inline real
     r2( const Vector< real > & aApproximated,
         const Vector< real > & aExact )
         {
         // calculate average of samples
         real tAverage = sum( aExact ) / ( real ) aExact.length();

         // sum of square residuals
         real tRootOfSSres = norm( aExact - aApproximated );

#ifdef BELFEM_BLAZE
         // the minus operator does not exist in blaze. Using a workaround
         index_t tN = aExact.length();
         real tRootOfSStot = 0.0 ;
         for( uint k=0; k<tN; ++k )
         {
             tRootOfSStot += ( aExact( k ) - tAverage ) * ( aExact( k ) - tAverage ) ;
         }
         tRootOfSStot = std::sqrt( tRootOfSStot );
#else
        // total sum of squares
        real tRootOfSStot = norm( aExact - tAverage );
#endif
        if( tRootOfSStot * tRootOfSStot < BELFEM_EPSILON || tRootOfSSres * tRootOfSSres < BELFEM_EPSILON )
        {
            return 1.0;
        }
        else
        {
            return 1.0 - std::pow( tRootOfSSres / tRootOfSStot, 2 );
        }
     }

//------------------------------------------------------------------------------

     inline real
     r2( const Matrix< real > & aApproximated,
         const Matrix< real > & aExact )
     {
         // number of rows
         index_t m = aExact.n_rows();

         // number of columns
         index_t n = aExact.n_cols();

         // calculate average of samples
         real tAverage = 0.0 ;

         // sum of square residuals
         real tSSres = 0.0 ;

         for( index_t i=0; i<m; ++i )
         {
             for( index_t j=0; j<n; ++j )
             {
                 tAverage += aExact( i, j );
                 tSSres +=
                          ( aExact( i, j ) - aApproximated( i, j ) )
                          *  ( aExact( i, j ) - aApproximated( i, j ) );
             }
         }

         tAverage /= m * n ;

         real tSStot = 0.0 ;
         for( index_t i=0; i<m; ++i )
         {
             for( index_t j=0; j<n; ++j )
             {
                 tSStot +=
                         ( aExact( i, j ) - tAverage )
                         *  ( aExact( i, j ) - tAverage );
             }
         }

         if( tSStot < BELFEM_EPSILON || tSSres < BELFEM_EPSILON )
         {
             return 1.0;
         }
         else
         {
             return 1.0 - tSSres / tSStot ;
         }
     }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_R2_HPP
