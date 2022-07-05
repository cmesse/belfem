//
// Created by Christian Messe on 2019-01-19.
//

#ifndef BELFEM_FN_INV_HPP
#define BELFEM_FN_INV_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_inv.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_inv.hpp"
#endif

namespace belfem
{

//------------------------------------------------------------------------------
    template< typename T >
    auto
    inv( const Matrix< T > & aA )
    -> decltype( inv( aA.matrix_data())) const
    {
        return inv( aA.matrix_data());
    }

//------------------------------------------------------------------------------

    template< typename T >
    auto
    inv( Matrix< T > & aA )
    -> decltype( inv( aA.matrix_data() ) )
    {
        return inv( aA.matrix_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_INV_HPP
