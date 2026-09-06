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
 * @brief Euclidean length of a vector.
 * @ingroup grp_linalg
 */

/**
 * @fn template<typename T> auto belfem::norm( const Vector<T> & aA )
 * @brief Euclidean (L2) norm of a vector.
 * @ingroup grp_linalg
 * @param aA vector whose Euclidean norm is computed
 * @return the square root of the sum of the squared entries
 */

#ifndef BELFEM_FN_NORM_HPP
#define BELFEM_FN_NORM_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template< typename ET >
    auto
    norm( ET &  aA )
        ->decltype( arma::norm( aA, 2 ) )
    {
        return arma::norm( aA, 2 );
    }
#elif  BELFEM_BLAZE
    template< typename ET >
    auto
    norm( ET &  aA )
        ->decltype( blaze::norm( aA ) )
    {
        return blaze::norm( aA );
    }
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    norm( const Vector< T > & aA )
        -> decltype( norm( aA.vector_data() ) )
    {
        return norm( aA.vector_data() );
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_FN_NORM_HPP
