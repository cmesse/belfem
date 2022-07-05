//
// Created by christian on 8/2/21.
//

#ifndef BELFEM_FN_BZ_EIGEN_HPP
#define BELFEM_FN_BZ_EIGEN_HPP
#include <blaze/util/typetraits/IsComplex.h>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_BZ_Vector.hpp"
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    void
    eigen( const Matrix< real > &  aMatrix,
           Vector< real > & aValues )
    {
        BELFEM_ASSERT( aMatrix.n_cols() == aMatrix.n_rows(),
                      "Matrix must be quadratic");

        aValues.set_size( aMatrix.n_cols() );

        blaze::DynamicVector<blaze::complex<real>, blaze::columnVector > tValues ;
        blaze::eigen( aMatrix.matrix_data(), tValues );

        size_t tN = aMatrix.n_cols() ;
        aValues.set_size( aMatrix.n_cols() );
        for( size_t k=0; k<tN; ++k )
        {
            if( std::abs( std::imag( tValues[ k ] ) ) > 1e-15 )
            {
                aValues( k ) = BELFEM_QUIET_NAN ;
            }
            else
            {
                aValues( k ) = std::real( tValues[ k ] );
            }
        }
   }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BZ_EIGEN_HPP
