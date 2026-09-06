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
 * @file
 * @brief Determinant of a square matrix.
 * @ingroup grp_linalg
 */

/**
 * @fn template<typename T> auto belfem::det( const Matrix<T> & aA )
 * @brief Determinant of a square matrix.
 * @ingroup grp_linalg
 * @param aA square matrix; not modified
 * @return the determinant
 */

#ifndef BELFEM_FN_DET_HPP
#define BELFEM_FN_DET_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_det.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_det.hpp"
#endif

namespace belfem
{

//------------------------------------------------------------------------------
    template< typename T >
    auto
    det( const Matrix< T > & aA )
        -> decltype( det( aA.matrix_data())) const
    {
        return det( aA.matrix_data());
    }

//------------------------------------------------------------------------------

    template< typename T >
    auto
    det( Matrix< T > & aA )
        -> decltype( det( aA.matrix_data() ) )
    {
        return det( aA.matrix_data() );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_DET_HPP
