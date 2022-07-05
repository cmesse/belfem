//
// Created by Christian Messe on 26.10.19.
//

#ifndef BELFEM_FN_CROSS_HPP
#define BELFEM_FN_CROSS_HPP

#include "assert.hpp"
#ifdef BELFEM_ARMADILLO
#include "fn_AR_cross.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_cross.hpp"
#endif


namespace belfem
{
    template < typename T >
    auto
    cross( const Vector< T > & aA, const Vector< T > & aB )
        -> decltype( cross( aA.vector_data(), aB.vector_data() ) )
    {
        BELFEM_ASSERT( aA.length() == 3 && aB.length() == 3,
            "Both vectors must have a length of 3");

        return cross( aA.vector_data(), aB.vector_data() );
    }
}
#endif //BELFEM_FN_CROSS_HPP
