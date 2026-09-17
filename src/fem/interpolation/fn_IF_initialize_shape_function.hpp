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

#ifndef BELFEM_FN_IF_INITIALIZE_SHAPE_FUNCTION_HPP
#define BELFEM_FN_IF_INITIALIZE_SHAPE_FUNCTION_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_IF_InterpolationFunction.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        InterpolationFunction  *
        initialize_shape_function(
                const ElementType      & aElementType,
                const Matrix< real >   & aXi,
                Cell< Matrix< real > > & aN,
                Cell< Matrix< real > > & adNdXi,
                Cell< Matrix< real > > & ad2NdXi2 );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_INITIALIZE_SHAPE_FUNCTION_HPP
