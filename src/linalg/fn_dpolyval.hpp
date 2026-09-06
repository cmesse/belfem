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
 * @brief Evaluates the first derivative of a polynomial.
 * @ingroup grp_linalg
 *
 * **Coefficient order is descending**: `aCoeffs(0)` is the coefficient of the highest
 * power and the last entry is the constant term, so a vector of length n+1 describes a
 * polynomial of degree n. This is the convention Armadillo and MATLAB use, and it is the
 * opposite of the ascending order some other libraries take. Getting it backwards
 * produces a plausible-looking wrong answer rather than an error.
 */

#ifndef BELFEM_FN_DPOLYVAL_HPP
#define BELFEM_FN_DPOLYVAL_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    /**
     * @brief Evaluates the first derivative of a polynomial at one point.
     * @ingroup grp_linalg
     *
     * Differentiates the polynomial described by @p aCoeffs analytically and evaluates
     * the result; nothing is approximated by finite differences.
     *
     * @param aCoeffs coefficients of the polynomial itself, highest power first --
     *                **not** the coefficients of its derivative. Must not be empty.
     * @param aX      the point to evaluate at
     * @return the first derivative at @p aX
     */
    template < typename T >
    T
    dpolyval( const Vector< T > & aCoeffs, const T aX )
    {
        const index_t tN = aCoeffs.length() - 1 ;
        T tPow = ( T ) tN ;
        T aResult = tPow * aCoeffs( 0 );

        for( index_t k=1; k<tN; ++k )
        {
            aResult *= aX;
            tPow -= 1.0 ;
            aResult += tPow * aCoeffs( k );
        }

        return aResult;
    }
}
#endif //BELFEM_FN_DPOLYVAL_HPP
