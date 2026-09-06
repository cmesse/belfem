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

#ifndef BELFEM_FN_BZ_POLYFIT_HPP
#define BELFEM_FN_BZ_POLYFIT_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"
#include "fn_gels.hpp"

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
                      ( long unsigned int ) aY.length() );

        BELFEM_ASSERT( aX.length() > aN, "not enough samples in vector to create a polynomial of degree %u",
                      ( unsigned int ) aN );

        // number of entries
        int tN = aN + 1;

        // get the number of samples
        index_t tNumSamples = aX.length();

        // the vandermonde matrix is solved directly rather than through the
        // normal equations, which would square the condition number
        Matrix< T > tVandermonde( tNumSamples, tN );

        // right hand side, also receives the coefficients
        Vector< T > tRHS( tNumSamples );

        // get a reference value to scale the polynomial
        T tXref = ( tNumSamples - 1 ) / ( aX( tNumSamples - 1 ) - aX( 0 ) );

        // loop over all samples
        for( index_t k=0; k<tNumSamples; ++k )
        {
            // scaled value of X to improve condition of matrix
            T tX = aX( k ) * tXref;

            // create polynomial, the highest power sits in the first column
            tVandermonde( k, aN ) = 1.0;

            for( int i=aN-1; i>=0; --i )
            {
                tVandermonde( k, i ) = tVandermonde( k, i+1 ) * tX ;
            }

            tRHS( k ) = aY( k );
        }

        // solve the least squares problem, this destroys the vandermonde
        // matrix and writes the coefficients into the first tN entries of tRHS
        Vector< T > tWork;

        gels( tVandermonde, tRHS, tWork );

        // copy data into coefficient vector
        aCoeffs.set_size( tN );

        T tScale = 1.0;

        for( int i=aN; i>=0; --i )
        {
            aCoeffs( i ) = tRHS( i ) * tScale;
            tScale *= tXref;
        }
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BZ_POLYFIT_HPP
