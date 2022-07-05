//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_AR_INV_HPP
#define BELFEM_FN_AR_INV_HPP

#include "cl_AR_Matrix.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    auto
    inv( const T & aExpression )
        -> decltype( arma::inv( aExpression ) )
    {
        return arma::inv( aExpression );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AR_INV_HPP
