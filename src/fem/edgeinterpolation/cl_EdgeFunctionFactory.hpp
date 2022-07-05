//
// Created by christian on 12/3/21.
//

#ifndef BELFEM_CL_EDGEFUNCTIONFACTORY_HPP
#define BELFEM_CL_EDGEFUNCTIONFACTORY_HPP

#include "cl_EF_EdgeFunction.hpp"

namespace belfem
{
    namespace fem
    {
        class EdgeFunctionFactory
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            EdgeFunctionFactory() = default;

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
