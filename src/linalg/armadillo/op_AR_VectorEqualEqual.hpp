//
// Created by Christian Messe on 03.09.19.
//

#ifndef BELFEM_OP_AR_VECTOREQUALEQUAL_HPP
#define BELFEM_OP_AR_VECTOREQUALEQUAL_HPP

#include "cl_AR_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    bool
    operator==( const Vector< T > & aA,
                const Vector< T > & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                arma::accu( aA.vector_data() == aB.vector_data() ) )
               ==  ( std::size_t ) aA.length() ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    bool
    operator==( const T          & aA,
                const Vector< T > & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                arma::accu( aA == aB.vector_data() ) )
               ==  ( std::size_t ) aA.length() ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    bool
    operator==( const Vector< T > & aA,
                const T           & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                arma::accu( aA.vector_data() == aB ) )
               ==  ( std::size_t ) aA.length() ;
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_OP_AR_VECTOREQUALEQUAL_HPP
