//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_AR_MAX_HPP
#define BELFEM_FN_AR_MAX_HPP

#include "cl_AR_Vector.hpp"
#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template <typename T >
    inline T
    max( const Vector<T> & aVector )
    {
        return aVector.vector_data().max();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    max( const Matrix<T> & aMatrix )
    {
        return aMatrix.matrix_data().max();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    max( const arma::Mat <T> & aMatrix )
    {
        return aMatrix.max();
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_MAX_HPP
