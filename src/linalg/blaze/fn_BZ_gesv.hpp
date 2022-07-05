//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_GESV_HPP
#define BELFEM_FN_BZ_GESV_HPP

#include "blaze_config.hpp"

#include <blaze/math/lapack/gesv.h>
#include "assert.hpp"
#include "cl_BZ_Vector.hpp"
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    void
    gesv( Matrix< T > & aA, Vector< T > & aX, Vector< int > & aP )
    {
        BELFEM_ASSERT( aA.n_rows() == aX.length(),
                      "Number of rows of matrix does not match." );
        BELFEM_ASSERT( aA.n_cols() == aX.length(),
                      "Number of cols of matrix does not match." );
        BELFEM_ASSERT( aP.length() >= aX.length(),
                      "Length of pivot vector does not match." );

        blaze::gesv(
                aA.matrix_data(),
                aX.vector_data(),
                aP.data() );
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    gesv( Matrix< T > & aA, Matrix< T > & aX, Vector< int > & aP )
    {
        BELFEM_ASSERT( aA.n_rows() == aX.length(),
                      "Number of rows of matrix does not match." );
        BELFEM_ASSERT( aA.n_cols() == aX.length(),
                      "Number of cols of matrix does not match." );
        BELFEM_ASSERT( aP.length() == aX.length(),
                      "Length of Pivot Vector does not match." );

        blaze::gesv(
                aA.matrix_data(),
                aX.matrix_data(),
                aP.data() );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BZ_GESV_HPP
