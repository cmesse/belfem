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

#ifndef BELFEM_CL_EDGEFUNCTIONFACTORY_HPP
#define BELFEM_CL_EDGEFUNCTIONFACTORY_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"

namespace belfem
{
    namespace fem
    {

        /**
         * @brief Creates edge functions by element type.
         *
         * @ingroup grp_fem_interpolation
         * @see @ref fem_interpolation_nedelec
         */
        class EdgeFunctionFactory
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            EdgeFunctionFactory();

//------------------------------------------------------------------------------

            ~EdgeFunctionFactory() = default;

//------------------------------------------------------------------------------

            EdgeFunction *
            create_edge_function( const ElementType aElementType );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */


#endif //BELFEM_CL_EDGEFUNCTIONFACTORY_HPP
