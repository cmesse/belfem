//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_OP_MATRIXTIMES_HPP
#define BELFEM_OP_MATRIXTIMES_HPP

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    inline auto
    operator*( const Matrix <T> & aA,
               const Matrix <T> & aB )
        -> decltype( aA.matrix_data() * aB.matrix_data())
    {
        BELFEM_ASSERT( aA.n_cols() == aB.n_rows(),
                      "Number of cols of matrix A must be equal to rows of B ( is %lu and %lu ).",
                      ( long unsigned int ) aA.n_cols(),
                      ( long unsigned int ) aB.n_rows());


        return aA.matrix_data() * aB.matrix_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    inline auto
    operator*( const Matrix <T> & aA,
               const Vector <T> & aB )
        -> decltype( aA.matrix_data() * aB.vector_data())
    {
        BELFEM_ASSERT( aA.n_cols() == aB.length(),
                      "Number of cols of matrix A must be equal to rows of B ( is %lu and %lu ).",
                      ( long unsigned int ) aA.n_cols(),
                      ( long unsigned int ) aB.length());

        return aA.matrix_data() * aB.vector_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator*( const Matrix <A> & aA,
               const B & aB )
        -> decltype( aA.matrix_data() * aB )
    {

        return aA.matrix_data() * aB;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator*( const A & aA,
               const Matrix <B> & aB )
        -> decltype( aA * aB.matrix_data())
    {

        return aA * aB.matrix_data();
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_OP_MATRIXTIMES_HPP
