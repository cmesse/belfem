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
 * @brief Closed-form inverse of a 2x2 matrix.
 * @ingroup grp_linalg
 */

#ifndef BELFEM_FN_INV2_HPP
#define BELFEM_FN_INV2_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"


namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Inverse of a 2x2 matrix from the closed-form adjugate.
     * @ingroup grp_linalg
     *
     * Computes the determinant and the adjugate directly, without a backend or LAPACK
     * call. Returning the determinant is deliberate: element Jacobians need it anyway,
     * so it comes for free rather than being recomputed.
     *
     * The singularity check is **relative**, not absolute: in debug builds BELFEM_ASSERT
     * requires det^2 > BELFEM_EPSILON^2 times the product of the squared row norms, so
     * the tolerance scales with the magnitude of @p aA. Release builds compile the check
     * out; a singular or nearly singular matrix then produces invalid or very large
     * floating-point values rather than a diagnostic.
     *
     * @param aA the matrix to invert; not modified
     * @param aB filled with the inverse; must already be 2x2
     * @return the determinant of @p aA
     */
    inline real
    inv2( const Matrix< real > & aA, Matrix< real > & aB )
    {
        // compute determinant
        real aDetA = aA( 0, 0 )*aA( 1,1 )
                   - aA( 1, 0 )*aA( 0,1 );

        aB( 0, 0 ) =   aA( 1, 1 );
        aB( 1, 0 ) = - aA( 1, 0 );
        aB( 0, 1 ) = - aA( 0, 1 );
        aB( 1, 1 ) =   aA( 0, 0 );

        // singularity test on the hadamard ratio |det| / prod( ||row_i|| ),
        // which is dimensionless: an absolute bound on the determinant is a
        // scale test, and rejects well conditioned but small matrices such as
        // any fine mesh jacobian. squared here to avoid the square roots
        BELFEM_ASSERT( aDetA * aDetA > BELFEM_EPSILON * BELFEM_EPSILON
                       * ( aA( 0, 0 ) * aA( 0, 0 ) + aA( 0, 1 ) * aA( 0, 1 ) )
                       * ( aA( 1, 0 ) * aA( 1, 0 ) + aA( 1, 1 ) * aA( 1, 1 ) ),
                       "Can't invert a singular matrix" );

        aB /= aDetA ;

        return aDetA ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_INV2_HPP
