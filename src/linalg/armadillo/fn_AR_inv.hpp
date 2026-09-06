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

#ifndef BELFEM_FN_AR_INV_HPP
#define BELFEM_FN_AR_INV_HPP

#include "cl_AR_Matrix.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    auto
    inv( const T & aExpression )
        -> decltype( arma::inv( aExpression ) )
    {
        return arma::inv( aExpression );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AR_INV_HPP
