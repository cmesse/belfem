//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_AR_MIN_HPP
#define BELFEM_FN_AR_MIN_HPP


#include "cl_AR_Vector.hpp"
#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template <typename T >
    inline T
    min( const Vector<T> & aVector )
    {
        return aVector.vector_data().min();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    min( const Matrix<T> & aMatrix )
    {
        return aMatrix.matrix_data().min();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    min( const arma::Mat <T> & aMatrix )
    {
        return aMatrix.min();
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AR_MIN_HPP
