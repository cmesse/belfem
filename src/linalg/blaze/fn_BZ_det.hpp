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

#ifndef BELFEM_FN_BZ_DET_HPP
#define BELFEM_FN_BZ_DET_HPP

#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( const T &  A)
     ->decltype( blaze::det( A ) )
    {
        return blaze::det( A );
    }

//------------------------------------------------------------------------------

    template<typename T >
    auto
    det( T &  A )
        ->decltype( blaze::det( A ))
    {
        return blaze::det( A );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_DET_HPP
