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

#ifndef BELFEM_FN_QUATERNION_ROTATE_VECTOR_HPP
#define BELFEM_FN_QUATERNION_ROTATE_VECTOR_HPP

#include "cl_Quaternion.hpp"
#include "cl_Vector.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Rotate a 3D vector by a unit quaternion (optimized, in-place)
     *
     * Uses the formula:  v' = v + 2w (u x v) + 2 (u x (u x v))
     *
     * where q = (w, u) with u = (x, y, z).  This avoids constructing
     * intermediate quaternions and is equivalent to q * v * q*.
     *
     * @param aQ       Unit quaternion defining the rotation
     * @param aVector  Input vector (length 3)
     * @param aResult  Output vector (length 3, pre-allocated)
     */
    template < typename T >
    void
    quaternion_rotate_vector(
            const Quaternion < T > & aQ,
            const Vector < T >    & aVector,
                  Vector < T >    & aResult )
    {
        BELFEM_ASSERT( std::abs( aQ.norm() - T( 1 ) ) < T( 100 ) * BELFEM_EPSILON,
                       "Input must be a unit quaternion" );
        BELFEM_ASSERT( aVector.length() == 3, "Input vector must have length 3" );
        BELFEM_ASSERT( aResult.length() == 3, "Output vector must have length 3" );

        T w = aQ.a();
        T ux = aQ.b();
        T uy = aQ.c();
        T uz = aQ.d();

        T vx = aVector( 0 );
        T vy = aVector( 1 );
        T vz = aVector( 2 );

        // t = 2 * (u x v)
        T tx = T( 2 ) * ( uy * vz - uz * vy );
        T ty = T( 2 ) * ( uz * vx - ux * vz );
        T tz = T( 2 ) * ( ux * vy - uy * vx );

        // v' = v + w * t + u x t
        aResult( 0 ) = vx + w * tx + ( uy * tz - uz * ty );
        aResult( 1 ) = vy + w * ty + ( uz * tx - ux * tz );
        aResult( 2 ) = vz + w * tz + ( ux * ty - uy * tx );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_QUATERNION_ROTATE_VECTOR_HPP
