//
// Created by Christian Messe on 09.12.19.
//

#ifndef BELFEM_FN_BZ_POLYFIT_HPP
#define BELFEM_FN_BZ_POLYFIT_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"
#include "fn_trans.hpp"

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

        // create the vandermonde matrix
        Matrix< real > tVandermonde( tN, tN, 0.0 );

        // help vector
        Matrix< real > tP( tN, 1, 1.0 );

        // right hand side
        Matrix< real > tRHS( tN, 1, 0.0 );

        // get the number of samples
        index_t tNumSamples = aX.length();

        // get a reference value to scale the polynomial
        real tXref = ( aX.length() - 1 ) / ( aX( aX.length() - 1 ) - aX( 0 ) );

        // loop over all samples
        for( index_t k=0; k<tNumSamples; ++k )
        {
            // scaled value of X to improve condition of matrix
            real tX = aX( k ) * tXref;

            // create polynomial
            for( int i=aN-1; i>=0; --i )
            {
                tP( i, 0 ) = tP( i+1, 0 ) * tX ;
            }

            // add polynomial to vandermode matrix
            tVandermonde += tP * trans( tP );

            // add value to right hand side
            tRHS += aY( k ) * tP;
        }

        // solve system using the preconditioned Cholesky solver
        blaze::posv( tVandermonde.matrix_data(), tRHS.matrix_data(), 'U' );

        // copy data into coefficient vector
        aCoeffs.set_size( tN );

        real tScale = 1.0;

        for( int i=aN; i>=0; --i )
        {
            aCoeffs( i ) = tRHS( i, 0 ) * tScale;
            tScale *= tXref;
        }
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BZ_POLYFIT_HPP
