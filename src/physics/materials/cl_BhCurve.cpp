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

#include "commtools.hpp"
#include "cl_BhCurve.hpp"


namespace belfem
{
    namespace material
    {
        BhCurve::BhCurve( const string & aPath, const string & aLabel ) :
            mCommRank( comm_rank() ),
            mPath( aPath ),
            mLabel( aLabel )
        {

        }
        real
        BhCurve::nu( const real B ) const
        {
            BELFEM_ERROR( false, "not implemented for baseclass" );
            return BELFEM_QUIET_NAN ;
        }

        real
        BhCurve::mu( const real H ) const
        {
            BELFEM_ERROR( false, "not implemented for baseclass" );
            return BELFEM_QUIET_NAN ;
        }

        void
        BhCurve::dmudH( const real H, real & mu, real & dmudH ) const
        {
            BELFEM_ERROR( false, "not implemented for baseclass" );
        }

    }
}
