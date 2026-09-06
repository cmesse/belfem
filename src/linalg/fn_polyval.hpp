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
 * @brief Evaluates a polynomial.
 * @ingroup grp_linalg
 *
 * **Coefficient order is descending**: `aCoeffs(0)` is the coefficient of the highest
 * power and the last entry is the constant term, so a vector of length n+1 describes a
 * polynomial of degree n. This is the convention Armadillo and MATLAB use, and it is the
 * opposite of the ascending order some other libraries take. Getting it backwards
 * produces a plausible-looking wrong answer rather than an error.
 */

#ifndef BELFEM_FN_POLYVAL_HPP
#define BELFEM_FN_POLYVAL_HPP

#ifdef BELFEM_ARMADILLO
#include "armadillo.hpp"
#endif

#include "typedefs.hpp"
#include "assert.h"
#include "cl_Vector.hpp"

namespace belfem
{

//------------------------------------------------------------------------------

    /**
     * @brief Evaluates a polynomial at one point, by Horner's scheme.
     * @ingroup grp_linalg
     * @param aCoeffs coefficients, highest power first Must not be empty.
     * @param aX      the point to evaluate at
     * @return the value of the polynomial at @p aX
     */
    template < typename T >
    T
    polyval( const Vector< T > & aCoeffs, const T aX )
    {
        const index_t tN = aCoeffs.length();
        T aResult = aCoeffs( 0 );

        for( index_t k=1; k<tN; ++k )
        {
            aResult *= aX;
            aResult += aCoeffs( k );
        }

        return aResult;
    }

//------------------------------------------------------------------------------


    template < typename T >
    void
    polyval( const Vector< T > & aCoeffs, const Vector< T > & aX, Vector< T > & aY )
    {
#ifdef BELFEM_ARMADILLO
        aY.vector_data() = arma::polyval( aCoeffs.vector_data(), aX.vector_data() );
#elif BELFEM_BLAZE
        index_t tN = aX.length();
        aY.set_size( tN );

        if( aCoeffs.length() == 0 )
        {
            aY.vector_data() = 0;
            return;
        }

        // Initialize with highest order coefficient
        aY.vector_data() = aCoeffs(0);

        // Horner's method with vectorized operations
        for( index_t i = 1; i < aCoeffs.length(); ++i )
        {
            aY.vector_data() = aY.vector_data() * aX.vector_data() + aCoeffs(i);
        }
#else
        index_t tN = aX.vector_data();
        aY.set_size( tN );
        for( uint k=0; k<tN; ++k )
        {
            aY( k ) = polyval( aCoeffs, aX( k ) );
        }
#endif
    }



//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_POLYVAL_HPP
