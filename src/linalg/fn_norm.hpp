//
// Created by Christian Messe on 2019-01-23.
//

#ifndef BELFEM_FN_NORM_HPP
#define BELFEM_FN_NORM_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template< typename ET >
    auto
    norm( ET &  aA )
        ->decltype( arma::norm( aA, 2 ) )
    {
        return arma::norm( aA, 2 );
    }
#elif  BELFEM_BLAZE
    template< typename ET >
    auto
    norm( ET &  aA )
        ->decltype( blaze::norm( aA ) )
    {
        return blaze::norm( aA );
    }
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    norm( const Vector< T > & aA )
        -> decltype( norm( aA.vector_data() ) )
    {
        return norm( aA.vector_data() );
    }
//------------------------------------------------------------------------------

}
#endif //BELFEM_FN_NORM_HPP
