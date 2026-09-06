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

#include "GT_globals.hpp"
#include "cl_GT_TransportPolyEmpty.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        TransportPolyEmpty::TransportPolyEmpty( const enum TransportPolyType aType ) :
                TransportPoly( aType, 0.0, gTmax, Vector< real >( 1, 0.0 ), TransportPolyKind::EMPTY )
        {

        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::rawpoly( const real T ) const
        {
            BELFEM_ERROR( false, "rawpoly() not supported for empty transport polynomial" );
            return 0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::drawpoly( const real T ) const
        {
            BELFEM_ERROR( false, "drawpoly() not supported for empty transport polynomial" );
            return 0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::ddrawpoly( const real T ) const
        {
            BELFEM_ERROR( false, "ddrawpoly() not supported for empty transport polynomial" );
            return 0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::eval( const real T ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::deval( const real T ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::ddeval( const real T ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */