//
// Created by christian on 12/3/21.
//

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
