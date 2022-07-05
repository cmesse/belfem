//
// Created by Christian Messe on 26.08.19.
//

//
// Created by Christian Messe on 26.08.19.
//

#include "GT_globals.hpp"
#include "cl_GT_TransportPolyGlue.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        TransportPolyGlue::TransportPolyGlue(
                const enum TransportPolyType aType,
                const real aTmin,
                const real aTmax,
                const Vector<real> & aCoefficients ) :
                TransportPoly( aType, aTmin, aTmax, aCoefficients )
        {

        }

//------------------------------------------------------------------------------

        real
        TransportPolyGlue::rawpoly( const real & aT ) const
        {
            return    ( (   mCoefficients( 0 )   * aT
                          + mCoefficients( 1 ) ) * aT
                          + mCoefficients( 2 ) ) * aT
                          + mCoefficients( 3 ) ;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyGlue::drawpoly( const real & aT ) const
        {
            return   (   3.0*mCoefficients( 0 )   * aT
                       + 2.0*mCoefficients( 1 ) ) * aT
                           + mCoefficients( 2 );
        }

//------------------------------------------------------------------------------

        real
        TransportPolyGlue::ddrawpoly( const real & aT ) const
        {
            return      6.0*mCoefficients( 1 )  * aT
                    +   2.0*mCoefficients( 2 ) ;
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */