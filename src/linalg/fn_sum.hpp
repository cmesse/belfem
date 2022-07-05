//
// Created by Christian Messe on 2018-12-26.
//

#ifndef BELFEM_FN_SUM_HPP
#define BELFEM_FN_SUM_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template< typename T>
    T
    sum( const arma::Mat< T > & aA )
    {
        return arma::as_scalar( arma::accu( aA ) );
    }
#elif  BELFEM_BLAZE
    template< typename T>
    auto
    sum( const blaze::DynamicVector< T, blaze::columnVector > & aA )
        -> decltype( blaze::sum( aA ) )
    {
        return blaze::sum( aA );
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T>
    auto
    sum( const Vector< T > & aA )
        -> decltype( sum( aA.vector_data() ) )
    {
        return sum( aA.vector_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_SUM_HPP
