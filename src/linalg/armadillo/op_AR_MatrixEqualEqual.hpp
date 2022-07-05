//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_OP_AR_MATRIXEQUALEQUAL_HPP
#define BELFEM_OP_AR_MATRIXEQUALEQUAL_HPP

#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    template< typename T >
    bool
    operator==( const Matrix <T> & aA,
                const Matrix <T> & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                    arma::accu( aA.matrix_data() == aB.matrix_data() ) )
                   ==  ( std::size_t ) aA.n_rows()*aA.n_cols() ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    bool
    operator==( const T & aA,
                const Matrix <T> & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                    arma::accu( aA == aB.matrix_data() ) )
                   ==  ( std::size_t ) aA.n_rows()*aA.n_cols() ;
        return aA == aB.matrix_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    bool
    operator==( const Matrix <T> & aA,
                const T & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                    arma::accu( aA.matrix_data() == aB ) )
                   ==  ( std::size_t ) aA.n_rows()*aA.n_cols() ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_OP_AR_MATRIXEQUALEQUAL_HPP
