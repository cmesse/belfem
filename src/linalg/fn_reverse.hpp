//
// Created by christian on 3/12/23.
//

#ifndef BELFEM_FN_REVERSE_HPP
#define BELFEM_FN_REVERSE_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template< typename ET >
    auto
    reverse( ET &  aA )
        ->decltype( arma::reverse( aA ) )
    {
        return arma::reverse( aA );
    }

#elif  BELFEM_BLAZE
    template< typename ET >
    auto
    reverse( ET &  aA )
        ->decltype( blaze::reverse( aA ) )
    {
        return blaze::reverse( aA );
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    reverse( const Vector< T > & aA )
    -> decltype( reverse( aA.vector_data() ) )
    {
        return reverse( aA.vector_data() );
    }

//------------------------------------------------------------------------------

}

#endif //BELFEM_FN_REVERSE_HPP
