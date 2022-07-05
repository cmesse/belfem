//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_INV_HPP
#define BELFEM_FN_BZ_INV_HPP

#include "blaze_config.hpp"
#include <blaze/math/functors/Inv.h>
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    auto
    inv( const Matrix< T > & aExpression )
        -> decltype( blaze::inv( aExpression.matrix_data() ) )
    {
        return blaze::inv( aExpression.matrix_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_INV_HPP
