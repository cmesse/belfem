//
// Created by Christian Messe on 08.11.19.
//

#ifndef BELFEM_FN_DET_HPP
#define BELFEM_FN_DET_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_det.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_det.hpp"
#endif

namespace belfem
{

//------------------------------------------------------------------------------
    template< typename T >
    auto
    det( const Matrix< T > & aA )
        -> decltype( det( aA.matrix_data())) const
    {
        return det( aA.matrix_data());
    }

//------------------------------------------------------------------------------

    template< typename T >
    auto
    det( Matrix< T > & aA )
        -> decltype( det( aA.matrix_data() ) )
    {
        return det( aA.matrix_data() );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_DET_HPP
