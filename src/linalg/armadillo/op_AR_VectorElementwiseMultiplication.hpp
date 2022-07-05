//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_OP_AR_VECTORELEMENTWISEMULTIPLICATION_HPP
#define BELFEM_OP_AR_VECTORELEMENTWISEMULTIPLICATION_HPP

#include "cl_AR_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    template< typename T >
    auto
    operator%(
            const Vector<T> & aA,
            const Vector<T> & aB )
        -> decltype( aA.vector_data() % aB.vector_data())
    {
        return aA.vector_data() % aB.vector_data();
    }

//--------------------------------------------------------------------------------

    template< typename T >
    auto
    operator%( const arma::Mat<T> & aA, const arma::Mat<T> & aB )
        -> decltype( aA.data() % aB.data())
    {
        return aA.data() % aB.data();
    }

//--------------------------------------------------------------------------------

    template< typename T, typename ET >
    auto
    operator%( const Vector<T> & aA, const ET & aB )
    -> decltype( aA.vector_data() % aB )
    {
        return aA.vector_data() % aB;
    }

//--------------------------------------------------------------------------------

    template< typename ET, typename T >
    auto
    operator%( const ET & aA, const Vector<T> & aB )
        -> decltype( aA % aB.vector_data())
    {
        return aA % aB.vector_data();
    }

//--------------------------------------------------------------------------------

}
#endif //BELFEM_OP_AR_VECTORELEMENTWISEMULTIPLICATION_HPP
