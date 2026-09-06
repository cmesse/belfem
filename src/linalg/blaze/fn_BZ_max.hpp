/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

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
