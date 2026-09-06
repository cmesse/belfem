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

#ifndef BELFEM_FN_AR_DET_HPP
#define BELFEM_FN_AR_DET_HPP

#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( const T &  A)
        ->decltype( arma::det( A ) )
    {
        return arma::det( A );
    }

//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( T &  A )
        ->decltype( arma::det( A ))
    {
        return arma::det( A );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_DET_HPP
