//
// Created by Christian Messe on 26.08.19.
//

#include "GT_globals.hpp"
#include "cl_GT_TransportPolyEmpty.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        TransportPolyEmpty::TransportPolyEmpty( const enum TransportPolyType aType ) :
                TransportPoly( aType, 0.0, gTmax, Vector< real >( 1, 0.0 ) )
        {

        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::rawpoly( const real & aT ) const
        {
            BELFEM_ERROR( false, "rawpoly() not supported for custom transport polynomial" );
            return 0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::drawpoly( const real & aT ) const
        {
            BELFEM_ERROR( false, "drawpoly() not supported for custom transport polynomial" );
            return 0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::ddrawpoly( const real & aT ) const
        {
            BELFEM_ERROR( false, "drawpoly() not supported for custom transport polynomial" );
            return 0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::eval( const real &  aT ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::deval( const real &  aT ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyEmpty::ddeval( const real &  aT ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */