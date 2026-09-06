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

#ifndef BELFEM_FN_QUATERNION_FROM_ROTATION_MATRIX_HPP
#define BELFEM_FN_QUATERNION_FROM_ROTATION_MATRIX_HPP

#include "cl_Quaternion.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Extract a unit quaternion from a 3x3 rotation matrix
     *
     * Uses the trace-first variant of Shepperd's method: when the trace is
     * positive, w is recovered from the trace (bounded away from zero);
     * otherwise the largest diagonal entry selects which vector component
     * is extracted first.
     *
     * @param aMatrix  A proper rotation matrix (det = +1, orthogonal)
     * @return Unit quaternion representing the same rotation
     */
    template < typename T >
    Quaternion < T >
    quaternion_from_rotation_matrix( const Matrix < T > & aMatrix )
    {
        BELFEM_ASSERT( aMatrix.n_rows() == 3 && aMatrix.n_cols() == 3,
                      "Input matrix must be 3x3" );

        T r00 = aMatrix( 0, 0 );
        T r01 = aMatrix( 0, 1 );
        T r02 = aMatrix( 0, 2 );
        T r10 = aMatrix( 1, 0 );
        T r11 = aMatrix( 1, 1 );
        T r12 = aMatrix( 1, 2 );
        T r20 = aMatrix( 2, 0 );
        T r21 = aMatrix( 2, 1 );
        T r22 = aMatrix( 2, 2 );

        // check that the matrix is a proper rotation (det = +1); the determinant
        // is formed inside the assert so it costs nothing once asserts compile out
        BELFEM_ASSERT( std::abs( r00 * ( r11 * r22 - r12 * r21 )
                               - r01 * ( r10 * r22 - r12 * r20 )
                               + r02 * ( r10 * r21 - r11 * r20 )
                               - T( 1 ) ) < T( 100 ) * BELFEM_EPSILON,
                       "Input must be a proper rotation matrix (det = +1)" );
        T tTrace = r00 + r11 + r22;

        T w, x, y, z;

        if ( tTrace > T( 0 ) )
        {
            // trace > 0: recover w from the trace (w^2 >= 1/4, safe division)
            T s = T( 0.5 ) / std::sqrt( tTrace + T( 1 ) );
            w = T( 0.25 ) / s;
            x = ( r21 - r12 ) * s;
            y = ( r02 - r20 ) * s;
            z = ( r10 - r01 ) * s;
        }
        else if ( r00 > r11 && r00 > r22 )
        {
            // x is the largest component
            T s = T( 2 ) * std::sqrt( T( 1 ) + r00 - r11 - r22 );
            w = ( r21 - r12 ) / s;
            x = T( 0.25 ) * s;
            y = ( r01 + r10 ) / s;
            z = ( r02 + r20 ) / s;
        }
        else if ( r11 > r22 )
        {
            // y is the largest component
            T s = T( 2 ) * std::sqrt( T( 1 ) + r11 - r00 - r22 );
            w = ( r02 - r20 ) / s;
            x = ( r01 + r10 ) / s;
            y = T( 0.25 ) * s;
            z = ( r12 + r21 ) / s;
        }
        else
        {
            // z is the largest component
            T s = T( 2 ) * std::sqrt( T( 1 ) + r22 - r00 - r11 );
            w = ( r10 - r01 ) / s;
            x = ( r02 + r20 ) / s;
            y = ( r12 + r21 ) / s;
            z = T( 0.25 ) * s;
        }

        Quaternion < T > tQ( w, x, y, z );
        tQ.normalize();
        return tQ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_QUATERNION_FROM_ROTATION_MATRIX_HPP
