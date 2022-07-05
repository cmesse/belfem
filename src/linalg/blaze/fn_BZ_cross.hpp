//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_CROSS_HPP
#define BELFEM_FN_BZ_CROSS_HPP

#include "cl_BZ_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename ET >
    auto
    cross( const ET & aA, const ET & aB )
        -> decltype( blaze::cross( aA, aB ))
    {
        return blaze::cross( aA, aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T, typename ET >
    auto
    cross( Vector<T> & aA, const ET & aB )
        -> decltype( blaze::cross( aA.vector_data(), aB ))
    {
        return blaze::cross( aA.vector_data(), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename T >
    auto
    cross( const ET & aA, const Vector<T> & aB )
        -> decltype( blaze::cross( aA, aB.vector_data()))
    {
        return blaze::cross( aA.vector_data(), aB );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_CROSS_HPP
