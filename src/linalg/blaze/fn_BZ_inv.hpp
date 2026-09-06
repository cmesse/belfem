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

#ifndef BELFEM_FN_BZ_INV_HPP
#define BELFEM_FN_BZ_INV_HPP

#include "blaze_config.hpp"
#include <blaze/math/functors/Inv.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    // blaze::inv() is lazy. Left lazy, blaze restructures B * inv( A ) into
    // trans( solve( trans( A ), trans( B ) ) ) and writes the result through a
    // DMatTransposer, whose resize() forwards ( m, n ) unswapped -- so a
    // non-square B silently transposes the destination. Evaluating here matches
    // the armadillo backend and keeps the operand shapes honest.
    //
    // The IsMatrix_v guard is load bearing: blaze::inv() also has a scalar
    // overload that accepts anything, so an unconstrained template would
    // swallow belfem::Matrix and hide the wrappers in fn_inv.hpp.
    template < typename T, typename = blaze::EnableIf_t< blaze::IsMatrix_v< T > > >
    blaze::DynamicMatrix< blaze::ElementType_t< T >, BLAZE_DEFAULT_STORAGE_ORDER >
    inv( const T & aExpression )
    {
        return blaze::inv( aExpression );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_INV_HPP
