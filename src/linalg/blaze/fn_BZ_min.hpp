//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_MIN_HPP
#define BELFEM_FN_BZ_MIN_HPP

#include "cl_BZ_Vector.hpp"
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    template <typename T >
    inline T
    min( const blaze::Columns< blaze::DynamicMatrix< T, BLAZE_DEFAULT_STORAGE_ORDER >, true, true, false > & aColumn )
    {
        return blaze::min( aColumn );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    min( const blaze::Rows< blaze::DynamicMatrix< T, BLAZE_DEFAULT_STORAGE_ORDER >, false, true, false > & aColumn )
    {
        return blaze::min( aColumn );
    }

//------------------------------------------------------------------------------

    template <typename T >
    inline T
    min( const Vector<T> & aVector )
    {
        return blaze::min( aVector.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    min( const Matrix<T> & aMatrix )
    {
        return blaze::min( aMatrix.matrix_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_MIN_HPP
