//
// Created by christian on 12/3/21.
//

#include "cl_EdgeFunctionFactory.hpp"
#include "cl_EF_TRI3.hpp"
#include "cl_EF_TRI6.hpp"
#include "cl_EF_TET4.hpp"
#include "cl_EF_TET10.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EdgeFunction *
        EdgeFunctionFactory::create_edge_function( const ElementType aElementType )
        {
            switch ( aElementType )
            {
                case ( ElementType::TRI3 ):
                {
                    return new EF_TRI3();
                }
                case ( ElementType::TRI6 ):
                {
                    return new EF_TRI6();
                }
                case ( ElementType::TET4 ):
                {
                    return new EF_TET4();
                }
                case ( ElementType::TET10 ):
                {
                    return new EF_TET10();
                }
                default :
                {
                    BELFEM_ERROR( false,
                                 "no edge function available for selected ElementType" );
                    return nullptr;
                }
            }
        }

//------------------------------------------------------------------------------
    }
}