//
// Created by Christian Messe on 19.12.20.
//

#ifndef BELFEM_FN_ROTATE_HPP
#define BELFEM_FN_ROTATE_HPP

#include "fn_TR_rotate42.hpp"
#include "cl_Tensor.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    template < typename T >
    inline void
    rotate(  const Tensor< T > & aB, const Matrix< T > & aR, Tensor< T > & aA )
    {
        BELFEM_ASSERT( aB.is_3333(), "Tensor A must be of type 3x3x3x3" );
        BELFEM_ASSERT( aR.n_rows() == 3 && aR.n_cols() == 3,
                      "Matrix B must be of size 3x3" );
        BELFEM_ASSERT( aA.is_3333(), "Tensor C must be of type 3x3x3x3" );

        tensor::rotate42( aB.data(), aR.data(), aA.data() );
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_FN_ROTATE_HPP
