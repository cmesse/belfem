//
// Created by Christian Messe on 08.11.19.
//

#ifndef BELFEM_FN_AR_DET_HPP
#define BELFEM_FN_AR_DET_HPP

#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( const T &  A)
        ->decltype( arma::det( A ) )
    {
        return arma::det( A );
    }

//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( T &  A )
        ->decltype( arma::det( A ))
    {
        return arma::det( A );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_DET_HPP
