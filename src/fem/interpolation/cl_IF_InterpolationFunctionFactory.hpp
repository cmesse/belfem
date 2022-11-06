//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_INTERPOLATIONFUNCTIONFACTORY_HPP
#define BELFEM_CL_IF_INTERPOLATIONFUNCTIONFACTORY_HPP

#include "Mesh_Enums.hpp"
#include "cl_IF_InterpolationFunction.hpp"

namespace belfem
{
    namespace fem
    {
        class InterpolationFunctionFactory
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            InterpolationFunctionFactory() = default;

//------------------------------------------------------------------------------

            ~InterpolationFunctionFactory() = default;
//------------------------------------------------------------------------------

            InterpolationFunction *
            create_function( const ElementType & aElementType, const InterpolationType aType );

//------------------------------------------------------------------------------

            InterpolationFunction *
            create_lagrange_function( const ElementType & aElementType );

//------------------------------------------------------------------------------

            InterpolationFunction *
            create_hermite_function( const ElementType & aElementType );

//------------------------------------------------------------------------------

            InterpolationFunction *
            create_bernstein_function( const ElementType & aElementType );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_IF_INTERPOLATIONFUNCTIONFACTORY_HPP
