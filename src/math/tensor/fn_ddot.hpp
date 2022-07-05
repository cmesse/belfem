//
// Created by Christian Messe on 19.12.20.
//

#ifndef BELFEM_FN_DDOT_HPP
#define BELFEM_FN_DDOT_HPP

#include "assert.hpp"
#include "fn_TR_contract44.hpp"
#include "cl_Matrix.hpp"
#include "cl_Tensor.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    template < typename T >
    inline void
    ddot(  const Tensor< T > & aA, const Tensor< T > & aB, Tensor< T > & aC )
    {
        BELFEM_ASSERT( aA.is_3333(), "Tensor A must be of type 3x3x3x3" );
        BELFEM_ASSERT( aB.is_3333(), "Tensor B must be of type 3x3x3x3" );
        BELFEM_ASSERT( aC.is_3333(), "Tensor C must be of type 3x3x3x3" );

        tensor::contract44( aA.data(), aB.data(), aC.data() );
    }

//----------------------------------------------------------------------------

    template < typename T >
    inline void
    ddot(  const Tensor< T > & aA, const Matrix< T > & aB, Matrix< T > & aC )
    {
        BELFEM_ASSERT( aA.is_3333(), "Tensor A must be of type 3x3x3x3" );
        BELFEM_ASSERT( aB.n_rows() == 3 && aB.n_cols() == 3,
                      "Matrix B must be of size 3x3" );
        BELFEM_ASSERT( aC.n_rows() == 3 && aC.n_cols() == 3,
                      "Matrix C must be of size 3x3" );

        tensor::contract42( aA.data(), aB.data(), aC.data() );
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_FN_DDOT_HPP
