//
// Created by Christian Messe on 08.11.19.
//

#ifndef BELFEM_FN_BZ_DET_HPP
#define BELFEM_FN_BZ_DET_HPP

#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( const T &  A)
     ->decltype( blaze::det( A ) )
    {
        return blaze::det( A );
    }

//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( T &  A )
        ->decltype( blaze::det( A ))
    {
        return blaze::det( A );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_DET_HPP
