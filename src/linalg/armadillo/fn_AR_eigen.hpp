//
// Created by christian on 8/2/21.
//

#ifndef BELFEM_FN_AR_EIGEN_HPP
#define BELFEM_FN_AR_EIGEN_HPP
#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_AR_Vector.hpp"
#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    void
    eigen( const Matrix< real > &  aMatrix,
           Vector< real > & aValues )
    {
        BELFEM_ASSERT( aMatrix.n_cols() == aMatrix.n_rows(),
                      "Matrix must be quadratic");

        arma::cx_vec tValues ;


        arma::eig_gen( tValues, aMatrix.matrix_data()  );

        size_t tN = aMatrix.n_cols() ;
        aValues.set_size( aMatrix.n_cols() );
        for( size_t k=0; k<tN; ++k )
        {
            if( std::abs( std::imag( tValues( k ) ) ) > 1e-15 )
            {
                aValues( k ) = BELFEM_QUIET_NAN ;
            }
            else
            {
                aValues( k ) = std::real( tValues( k ) );
            }
        }
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AR_EIGEN_HPP
