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
 * @brief Cross product of a normal vector with every column of a matrix.
 * @ingroup grp_linalg
 *
 * Each overload treats every column of @p aA as one input vector and computes
 * @p aN x column(k). The 2D overloads store one scalar z-component per column in a
 * Vector; the 3D overloads store one 3-vector per column in a Matrix.
 *
 * Two conventions differ between the overloads and are easy to get wrong:
 *
 * - the overloads **without** @p aScale **assign** to @p aNxA; those **with** it
 *   **accumulate** into it, so the caller must initialise @p aNxA first;
 * - only the **2D** overloads finish by zeroing result entries whose magnitude is below
 *   BELFEM_EPSILON relative to the norm of the result. That is cosmetic suppression of
 *   rounding dust in what should be exact zeros. The 3D overloads do not do it.
 */

#ifndef BELFEM_FN_CROSSMAT_HPP
#define BELFEM_FN_CROSSMAT_HPP

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "assert.hpp"
namespace belfem
{
    /**
     * @brief Crosses a normal vector with every column of a matrix, 2D.
     * @ingroup grp_linalg
     * @param aN   normal vector, length 2
     * @param aA   matrix whose first two rows hold the vector components; further rows are ignored
     * @param aNxA assigned to, one entry per column of @p aA; must already have that length
     */
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
                  Vector< real > & aNxA )
    {

        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 2,
                      "Length of normal vector must be 2 (is %u)",
                      ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.length() == tN,
                      "Length of solution vector does not match (is %u but expect %u )",
                      ( unsigned int ) aNxA.length(),
                      ( unsigned int ) tN );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( k ) = aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) ;
        }

        // remove dust ( can be removed in the future if too slow, just for beautification )
        real tInvNorm = 1.0/norm( aNxA );
        for( uint k=0; k<tN; ++k )
        {
            if( std::abs(  aNxA( k ) * tInvNorm ) < BELFEM_EPSILON )
            {
                aNxA( k ) = 0.0 ;
            }
        }

    }


    /**
     * @brief Crosses a normal vector with every column of a matrix, 2D, scaled.
     * @ingroup grp_linalg
     * @param aN     normal vector, length 2
     * @param aA     matrix whose first two rows hold the vector components; further rows are ignored
     * @param aScale multiplier applied to each column cross product before accumulation
     * @param aNxA   accumulated into, one entry per column of @p aA; initialise it first
     */
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
            const real             aScale,
            Vector< real > & aNxA )
    {

        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 2,
                       "Length of normal vector must be 2 (is %u)",
                       ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.length() == tN,
                       "Length of solution vector does not match (is %u but expect %u )",
                       ( unsigned int ) aNxA.length(),
                       ( unsigned int ) tN );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( k ) += aScale * ( aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) );
        }

        // remove dust ( can be removed in the future if too slow, just for beautification )
        real tInvNorm = 1.0/norm( aNxA );
        for( uint k=0; k<tN; ++k )
        {
            if( std::abs(  aNxA( k ) * tInvNorm ) < BELFEM_EPSILON )
            {
                aNxA( k ) = 0.0 ;
            }
        }
    }

    /**
     * @brief Crosses a normal vector with every column of a matrix, 3D.
     * @ingroup grp_linalg
     * @param aN     normal vector, length 3
     * @param aA     matrix whose first three rows hold the vector components; further rows are ignored
     * @param aNxA   assigned to; must have 3 rows and one column per column of @p aA
     */
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
                  Matrix< real > & aNxA )
    {
        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 3,
                      "Length of normal vector must be 3 (is %u)",
                      ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.n_cols() == tN,
                      "Number of columns of solution matrix does not match (is %u but expect %u )",
                      ( unsigned int ) aNxA.n_cols(),
                      ( unsigned int ) tN );

        BELFEM_ASSERT( aNxA.n_rows() == 3,
                      "Number of rows of solution matrix does not match (is %u but expect 3)",
                      ( unsigned int ) aNxA.n_rows()  );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( 0, k ) = aN( 1 ) * aA( 2, k ) - aN( 2 ) * aA( 1, k ) ;
            aNxA( 1, k ) = aN( 2 ) * aA( 0, k ) - aN( 0 ) * aA( 2, k ) ;
            aNxA( 2, k ) = aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) ;

        }
    }

    /**
     * @brief Crosses a normal vector with every column of a matrix, 3D, scaled.
     * @ingroup grp_linalg
     * @param aN     normal vector, length 3
     * @param aA     matrix whose first three rows hold the vector components; further rows are ignored
     * @param aScale multiplier applied to each column cross product before accumulation
     * @param aNxA   accumulated into; 3 rows, one column per column of @p aA; initialise it first
     */
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
            const real             aScale,
            Matrix< real > & aNxA )
    {
        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 3,
                      "Length of normal vector must be 3 (is %u)",
                      ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.n_cols() == tN,
                      "Number of columns of solution matrix does not match (is %u but expect %u )",
                      ( unsigned int ) aNxA.n_cols(),
                      ( unsigned int ) tN );

        BELFEM_ASSERT( aNxA.n_rows() == 3,
                      "Number of rows of solution matrix does not match (is %u but expect 3)",
                      ( unsigned int ) aNxA.n_rows()  );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( 0, k ) += aScale * ( aN( 1 ) * aA( 2, k ) - aN( 2 ) * aA( 1, k ) );
            aNxA( 1, k ) += aScale * ( aN( 2 ) * aA( 0, k ) - aN( 0 ) * aA( 2, k ) );
            aNxA( 2, k ) += aScale * ( aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) );
        }
    }
}
#endif //BELFEM_FN_CROSSMAT_HPP
