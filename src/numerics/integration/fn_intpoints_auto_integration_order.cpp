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
            case( InterpolationOrder::QUARTIC ) :
            {
                return 13 ;
            }
            case( InterpolationOrder::QUINTIC ) :
            {
                return 16 ;
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