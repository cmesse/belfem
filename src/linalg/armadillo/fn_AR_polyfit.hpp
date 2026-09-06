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

#ifndef BELFEM_FN_AR_POLYFIT_HPP
#define BELFEM_FN_AR_POLYFIT_HPP

#include "armadillo.hpp"
#include "cl_AR_Vector.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    void
    polyfit( const Vector< T > & aX, const Vector< T > & aY, const uint & aN, Vector< T > & aCoeffs )
    {
        BELFEM_ASSERT( aX.length() == aY.length(),
                      "Legnths of X and Y vectors do not match ( %lu and %lu )",
                      ( long unsigned int ) aX.length(),
                      ( long unsigned int ) aY.length());

        BELFEM_ASSERT( aX.length() > aN, "not enough samples in vector to create a polynomial of degree %u",
                      ( unsigned int ) aN );

        // get size of vector
        int tN = aX.length();

        // get a reference value to scale the polynomial
        T tXref = ( tN - 1 ) / ( aX( tN - 1 ) - aX( 0 ) );

        // call armadillo
        aCoeffs.vector_data() = arma::polyfit( aX.vector_data() * tXref , aY.vector_data(), aN );

        // scale back vector
        T tScale = 1.0;

        for( int i=aN; i>=0; --i )
        {
            aCoeffs( i ) *= tScale;
            tScale *= tXref;
        }
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_POLYFIT_HPP
