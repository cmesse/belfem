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

#ifndef BELFEM_FN_QUATERNION_TO_ROTATION_MATRIX_HPP
#define BELFEM_FN_QUATERNION_TO_ROTATION_MATRIX_HPP

#include "cl_Quaternion.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Convert a unit quaternion to a 3x3 rotation matrix
     *
     * Given q = (w, x, y, z), the rotation matrix is:
     *
     *     | 1 - 2(y^2+z^2)    2(xy - wz)      2(xz + wy)   |
     * R = |   2(xy + wz)    1 - 2(x^2+z^2)    2(yz - wx)   |
     *     |   2(xz - wy)      2(yz + wx)    1 - 2(x^2+y^2) |
     *
     * @param aQ      Unit quaternion
     * @param aMatrix Pre-allocated 3x3 output matrix
     */
    template < typename T >
    void
    quaternion_to_rotation_matrix(
            const Quaternion < T > & aQ,
                  Matrix < T >    & aMatrix )
    {
        BELFEM_ASSERT( aMatrix.n_rows() == 3 && aMatrix.n_cols() == 3,
                      "Output matrix must be 3x3" );

        BELFEM_ASSERT( std::abs( aQ.norm() - T( 1 ) ) < T( 100 ) * BELFEM_EPSILON,
                       "Input must be a unit quaternion" );

        T w = aQ.a();
        T x = aQ.b();
        T y = aQ.c();
        T z = aQ.d();

        T xx = x * x;
        T yy = y * y;
        T zz = z * z;
        T xy = x * y;
        T xz = x * z;
        T yz = y * z;
        T wx = w * x;
        T wy = w * y;
        T wz = w * z;

        aMatrix( 0, 0 ) = T( 1 ) - T( 2 ) * ( yy + zz );
        aMatrix( 0, 1 ) = T( 2 ) * ( xy - wz );
        aMatrix( 0, 2 ) = T( 2 ) * ( xz + wy );

        aMatrix( 1, 0 ) = T( 2 ) * ( xy + wz );
        aMatrix( 1, 1 ) = T( 1 ) - T( 2 ) * ( xx + zz );
        aMatrix( 1, 2 ) = T( 2 ) * ( yz - wx );

        aMatrix( 2, 0 ) = T( 2 ) * ( xz - wy );
        aMatrix( 2, 1 ) = T( 2 ) * ( yz + wx );
        aMatrix( 2, 2 ) = T( 1 ) - T( 2 ) * ( xx + yy );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_QUATERNION_TO_ROTATION_MATRIX_HPP
