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
