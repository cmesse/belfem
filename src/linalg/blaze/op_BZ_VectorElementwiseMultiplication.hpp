//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_OP_BZ_VECTORELEMENTWISEMULTIPLICATION_HPP
#define BELFEM_OP_BZ_VECTORELEMENTWISEMULTIPLICATION_HPP

#include "assert.hpp"
#include "blaze_config.hpp"
#include <blaze/math/blas/gemv.h>

#include "cl_BZ_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    blaze::DynamicVector<T, blaze::columnVector>
    operator%(
            const Vector<T> & aA,
            const Vector<T> & aB )
    {

        BELFEM_ASSERT( aA.length() == aB.length(),
            "Length of vectors does not match ( %lu and %lu )",
                      ( long unsigned int ) aA.length(),
                      ( long unsigned int ) aB.length() );

        blaze::DynamicVector<T, blaze::columnVector> aC( aA.length() );

        const size_t tN = aA.length();

        for( size_t k=0; k<tN; ++k )
        {
            aC[ k ] = aA( k ) * aB( k );
        }
        return aC;
    }

//--------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_VECTORELEMENTWISEMULTIPLICATION_HPP
