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

#ifndef BELFEM_CL_IF_INTERPOLATIONFUNCTIONFACTORY_HPP
#define BELFEM_CL_IF_INTERPOLATIONFUNCTIONFACTORY_HPP

#include "Mesh_Enums.hpp"
#include "cl_IF_InterpolationFunction.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Creates shape functions by element and interpolation type.
         *
         * @ingroup grp_fem_interpolation
         * @see @ref fem_interpolation_interpolation_usage_guide
         */
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
            create_function( const ElementType aElementType, const InterpolationType aType );

//------------------------------------------------------------------------------

            InterpolationFunction *
            create_lagrange_function( const ElementType aElementType );

//------------------------------------------------------------------------------

            InterpolationFunction *
            create_hermite_function( const ElementType aElementType );

//------------------------------------------------------------------------------

            InterpolationFunction *
            create_bernstein_function( const ElementType aElementType );
            
//------------------------------------------------------------------------------

            InterpolationFunction *
            create_bubble_function( const ElementType aElementType, const uint aFacet );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_IF_INTERPOLATIONFUNCTIONFACTORY_HPP
