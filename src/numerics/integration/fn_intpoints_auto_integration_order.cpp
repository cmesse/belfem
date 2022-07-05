//
// Created by Christian Messe on 16.06.20.
//

#include "fn_intpoints_auto_integration_order.hpp"
#include "assert.hpp"
#include "meshtools.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    uint
    auto_integration_order( const ElementType aType )
    {
        switch ( mesh::interpolation_order( aType ) )
        {
            case ( InterpolationOrder::CONSTANT ) :
            {
                return 1;
            }
            case ( InterpolationOrder::LINEAR ) :
            {
                return 4;
            }
            case ( InterpolationOrder::SERENDIPITY ) :
            case ( InterpolationOrder::QUADRATIC ) :
            {
                return 7;
            }
            case ( InterpolationOrder::CUBIC ) :
            {
                return 10;
            }
            default:
            {
                BELFEM_ERROR( false,
                             "Don't know what integration order to choose" );
                return 0 ;
            }
        }
    }

//------------------------------------------------------------------------------
}