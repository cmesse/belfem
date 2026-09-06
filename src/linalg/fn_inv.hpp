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
 * @brief Inverse of a square matrix.
 * @ingroup grp_linalg
 *
 * Forwards to the backend's own inverse. For the fixed 2x2 and 3x3 cases prefer
 * belfem::inv2 and belfem::inv3, which use the closed-form adjugate and avoid the
 * backend call entirely.
 */

/**
 * @fn template<typename T> auto belfem::inv( const Matrix<T> & aA )
 * @brief Inverse of a square matrix.
 * @ingroup grp_linalg
 * @param aA square matrix; not modified
 * @return the inverse of @p aA
 */

#ifndef BELFEM_FN_INV_HPP
#define BELFEM_FN_INV_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_inv.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_inv.hpp"
#endif

namespace belfem
{

//------------------------------------------------------------------------------
    template< typename T >
    auto
    inv( const Matrix< T > & aA )
    -> decltype( inv( aA.matrix_data())) const
    {
        return inv( aA.matrix_data());
    }

//------------------------------------------------------------------------------

    template< typename T >
    auto
    inv( Matrix< T > & aA )
    -> decltype( inv( aA.matrix_data() ) )
    {
        return inv( aA.matrix_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_INV_HPP
