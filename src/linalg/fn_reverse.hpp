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
 * @brief Reverses the order of a vector's entries.
 * @ingroup grp_linalg
 *
 * Returns the entries in reverse order. The public wrapper takes its argument by const
 * reference and does not reverse in place.
 *
 */

#ifndef BELFEM_FN_REVERSE_HPP
#define BELFEM_FN_REVERSE_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template< typename T >
    auto
    reverse( arma::Mat <T> &  aA )
        ->decltype( arma::reverse( aA ) )
    {
        return arma::reverse( aA );
    }

#elif  BELFEM_BLAZE
    template< typename T >
    auto
    reverse( blaze::DynamicVector<T, blaze::columnVector> &  aA )
        ->decltype( blaze::reverse( aA ) )
    {
        return blaze::reverse( aA );
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    reverse( const Vector< T > & aA )
    -> decltype( reverse( aA.vector_data() ) )
    {
        return reverse( aA.vector_data() );
    }

//------------------------------------------------------------------------------

}

#endif //BELFEM_FN_REVERSE_HPP
