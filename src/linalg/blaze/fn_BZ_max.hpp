//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_MAX_HPP
#define BELFEM_FN_BZ_MAX_HPP

#include "cl_BZ_Vector.hpp"
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
    template <typename T >
    inline T
    max( const blaze::Columns< blaze::DynamicMatrix< T, BLAZE_DEFAULT_STORAGE_ORDER >, true, true, false > & aColumn )
    {
        return blaze::max( aColumn );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    max( const blaze::Rows< blaze::DynamicMatrix< T, BLAZE_DEFAULT_STORAGE_ORDER >, false, true, false > & aColumn )
    {
        return blaze::max( aColumn );
    }

//------------------------------------------------------------------------------

    template <typename T >
    inline T
    max( const Vector<T> & aVector )
    {
        return blaze::max( aVector.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    inline T
    max( const Matrix<T> & aMatrix )
    {
        return blaze::max( aMatrix.matrix_data() );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BZ_MAX_HPP
