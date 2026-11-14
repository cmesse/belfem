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
/**
 * @file fn_is_finite.hpp
 * @brief Finiteness test for a whole matrix or vector, backend-free.
 * @ingroup grp_linalg
 */
#ifndef BELFEM_FN_IS_FINITE_HPP
#define BELFEM_FN_IS_FINITE_HPP

#include <cmath>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief true if every entry of @p aA is finite ( no NaN, no infinity ).
     *
     * Reads through the element accessor, so the result is independent of the
     * backend's storage layout. Meant for setup-grade checks before a decomposition,
     * not for per-integration-point code.
     */
    template< typename T >
    bool
    is_finite( const Matrix< T > & aA )
    {
        const size_t tNumRows = aA.n_rows();
        const size_t tNumCols = aA.n_cols();
        for( size_t j = 0; j < tNumCols; ++j )
        {
            for( size_t i = 0; i < tNumRows; ++i )
            {
                if( ! std::isfinite( aA( i, j ) ) )
                {
                    return false;
                }
            }
        }
        return true;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    bool
    is_finite( const Vector< T > & aA )
    {
        const size_t tN = aA.length();
        for( size_t k = 0; k < tN; ++k )
        {
            if( ! std::isfinite( aA( k ) ) )
            {
                return false;
            }
        }
        return true;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_IS_FINITE_HPP
