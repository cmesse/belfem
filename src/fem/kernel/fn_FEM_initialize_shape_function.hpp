//
// Created by Christian Messe on 04.11.19.
//

#ifndef BELFEM_FN_FEM_INITIALIZE_SHAPE_FUNCTION_HPP
#define BELFEM_FN_FEM_INITIALIZE_SHAPE_FUNCTION_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "cl_FEM_Element.hpp"
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
