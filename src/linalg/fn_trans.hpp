//
// Created by Christian Messe on 2019-01-05.
//

#ifndef BELFEM_FN_TRANS_HPP
#define BELFEM_FN_TRANS_HPP

#ifdef BELFEM_ARMADILLO
#include "armadillo.hpp"
#elif  BELFEM_BLAZE
#include <blaze/math/functors/Trans.h>
#endif

#include "cl_Matrix.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template < typename ET >
    auto
    trans( ET & aMatrix )
        -> decltype( arma::strans( aMatrix ) )
    {
        return arma::strans( aMatrix );
    }

#elif  BELFEM_BLAZE
    template < typename ET >
    auto
    trans( ET & aMatrix )
        -> decltype( blaze::trans( aMatrix ) )
    {
        return blaze::trans( aMatrix );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename ET >
    auto
    trans( const ET & aMatrix )
        -> decltype( blaze::trans( aMatrix ) )
    {
        return blaze::trans( aMatrix );
    }

#endif
//------------------------------------------------------------------------------

    template < typename T >
    auto
    trans( Matrix < T > & aMatrix )
        -> decltype( trans( aMatrix.matrix_data() ) )
    {
        return  trans( aMatrix.matrix_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    auto
    trans( const Matrix < T > & aMatrix )
        -> decltype( trans( aMatrix.matrix_data() ) )
    {
        return  trans( aMatrix.matrix_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_TRANS_HPP
