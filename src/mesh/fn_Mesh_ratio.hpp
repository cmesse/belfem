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

#ifndef BELFEM_FN_MESH_RATIO_HPP
#define BELFEM_FN_MESH_RATIO_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    void
    ratio_ar2(
            const real     & aDeltaX0,
            const real     & aLength,
            const index_t  & aNumCells,
            Vector< real > & aX );

//------------------------------------------------------------------------------

    void
    ratio_adx2( const real     & aRatio,
                const real     & aLength,
                const index_t  & aNumCells,
                Vector< real > & aX );

//------------------------------------------------------------------------------

    // test function for ratio_ar2
    real
    _check_ratio(  const real     & aDeltaX0,
                   const real     & aLength,
                   const index_t  & aNumCells,
                   const real     & aRatio );

//------------------------------------------------------------------------------

}

#endif //BELFEM_FN_MESH_RATIO_HPP
