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
                TransportPoly( aType, aTmin, aTmax, aCoefficients, TransportPolyKind::GLUE )
        {

        }

//------------------------------------------------------------------------------

        real
        TransportPolyGlue::rawpoly( const real T ) const
        {
            return    ( (   mCoefficients( 0 )   * T
                          + mCoefficients( 1 ) ) * T
                          + mCoefficients( 2 ) ) * T
                          + mCoefficients( 3 ) ;
        }

//------------------------------------------------------------------------------

        real
        TransportPolyGlue::drawpoly( const real T ) const
        {
            return   (   3.0*mCoefficients( 0 )   * T
                       + 2.0*mCoefficients( 1 ) ) * T
                           + mCoefficients( 2 );
        }

//------------------------------------------------------------------------------

        real
        TransportPolyGlue::ddrawpoly( const real T ) const
        {
            return      6.0*mCoefficients( 0 )  * T
                    +   2.0*mCoefficients( 1 ) ;
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */