//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_OP_MATRIXMINUS_HPP
#define BELFEM_OP_MATRIXMINUS_HPP

#include "cl_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    inline auto
    operator-( const Matrix<T> & aA,
               const Matrix<T> & aB )
    -> decltype( aA.matrix_data() - aB.matrix_data() )
    {
        BELFEM_ASSERT( aA.n_rows() == aB.n_rows(),
                      "Number of rows of matrices does not match ( %lu and %lu ).",
                      ( long unsigned int ) aA.n_rows(),
                      ( long unsigned int ) aB.n_rows() );

        BELFEM_ASSERT( aA.n_cols() == aB.n_cols(),
                      "Number of colums of matrices does not match ( %lu and %lu ).",
                      ( long unsigned int ) aA.n_cols(),
                      ( long unsigned int ) aB.n_cols() );

        return aA.matrix_data() - aB.matrix_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator-( const Matrix< A > & aA,
               const         B  & aB )
        -> decltype( aA.matrix_data() - aB )
    {

        return aA.matrix_data() - aB ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator-( const           A & aA,
               const Matrix< B > & aB )
        -> decltype( aA - aB.matrix_data() )
    {

        return aA - aB.matrix_data() ;
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_OP_MATRIXMINUS_HPP
