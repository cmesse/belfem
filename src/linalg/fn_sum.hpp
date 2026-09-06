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
 * @brief Sum of the entries of a vector.
 * @ingroup grp_linalg
 *
 */

#ifndef BELFEM_FN_SUM_HPP
#define BELFEM_FN_SUM_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template< typename T>
    T
    sum( const arma::Mat< T > & aA )
    {
        return arma::as_scalar( arma::accu( aA ) );
    }
#elif  BELFEM_BLAZE
    template< typename T>
    auto
    sum( const blaze::DynamicVector< T, blaze::columnVector > & aA )
        -> decltype( blaze::sum( aA ) )
    {
        return blaze::sum( aA );
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T>
    auto
    sum( const Vector< T > & aA )
        -> decltype( sum( aA.vector_data() ) )
    {
        return sum( aA.vector_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_SUM_HPP
