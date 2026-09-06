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

#ifndef BELFEM_FN_AR_CROSS_HPP
#define BELFEM_FN_AR_CROSS_HPP

#include "cl_AR_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename ET >
    auto
    cross( const ET & aA, const ET & aB )
        -> decltype( arma::cross( aA, aB ))
    {
        return arma::cross( aA, aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T, typename ET >
    auto
    cross( Vector<T> & aA, const ET & aB )
        -> decltype( arma::cross( aA.vector_data(), aB ))
    {
        return arma::cross( aA.vector_data(), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename T >
    auto
    cross( const ET & aA, const Vector<T> & aB )
        -> decltype( arma::cross( aA, aB.vector_data()))
    {
        return arma::cross( aA.vector_data(), aB );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AR_CROSS_HPP
