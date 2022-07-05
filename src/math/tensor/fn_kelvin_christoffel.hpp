//
// Created by Christian Messe on 22.12.20.
//

#ifndef BELFEM_FN_KELVIN_CHRISTOFFEL_HPP
#define BELFEM_FN_KELVIN_CHRISTOFFEL_HPP


#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Tensor.hpp"
#include "fn_TR_kelvin_christoffel.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    template < typename T >
    inline void
    kelvin_christoffel(  const Tensor< T > & aA, const Vector< T > & aB, Matrix< T > & aC )
    {
        BELFEM_ASSERT( aA.is_3333(),
                      "operating tensor must be 3x3x3x3" );

        BELFEM_ASSERT( aB.length() == 3,
                              "argument vector must be allocated as 3x1" );

        BELFEM_ASSERT(    aC.n_rows() == 3
                      && aC.n_cols() == 3,
                      "target matrix must be allocated as 3x3" );

        tensor::kelvin_christoffel( aA.data(), aB.data(), aC.data() );
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_FN_KELVIN_CHRISTOFFEL_HPP
