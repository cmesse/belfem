//
// Created by christian on 12/3/21.
//

#include "cl_EdgeFunctionFactory.hpp"

#include "nedelec/cl_EF_HEX8.hpp"
#include "nedelec/cl_EF_TRI3.hpp"
#include "nedelec/cl_EF_TRI6.hpp"
#include "nedelec/cl_EF_TET4.hpp"
#include "nedelec/cl_EF_TET10.hpp"
#include "nedelec/cl_EF_HEX8TS.hpp"
#include "nedelec/cl_EF_HEX8TB.hpp"
#include "nedelec/cl_EF_QUAD4TS.hpp"
#include "nedelec/cl_EF_PENTA6TS.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EdgeFunctionFactory::EdgeFunctionFactory()
        {

        }

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
                case ( ElementType::PENTA6TS ) :
                {
                    return new EF_PENTA6TS() ;
                }
                case ( ElementType::QUAD4TS ) :
                {
                    return new EF_QUAD4TS() ;
                }
                case ( ElementType::HEX8 ) :
                {
                    return new EF_HEX8() ;
                }
                case ( ElementType::HEX8TS ) :
                {
                    return new EF_HEX8TS() ;
                }
                case ( ElementType::HEX8TB ) :
                {
                    return new EF_HEX8TB() ;
                }
                case ( ElementType::QUAD9TS ) :
                case ( ElementType::PENTA18TS ) :
                {
                    // Quadratic thin-shell Nédélec basis not yet implemented.
                    // Linear shells (QUAD4TS / PENTA6TS / HEX8TS) are
                    // supported; quadratic shells (QUAD9TS / PENTA18TS) fall
                    // through to this error until EF_QUAD9TS / EF_PENTA18TS
                    // are derived from the hierarchical p=2 Nédélec space
                    // and the thin-shell reduction in Messe et al. 2023.
                    // See todo/deferred/edge_function_quadratic_shells.md.
                    BELFEM_ERROR( false,
                        "Edge function for quadratic thin shell %s is not "
                        "implemented. Use linear shells (QUAD4TS/PENTA6TS/HEX8TS) "
                        "or implement EF_%s following Messe et al. 2023.",
                        to_string( aElementType ).c_str(),
                        to_string( aElementType ).c_str() );
                    return nullptr;
                }
                default :
                {
                    BELFEM_ERROR( false,
                        "No edge function available for ElementType %s.",
                        to_string( aElementType ).c_str() );
                    return nullptr;
                }
            }
        }

//------------------------------------------------------------------------------
    }
}