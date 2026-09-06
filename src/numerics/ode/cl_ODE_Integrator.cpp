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

#include "cl_ODE_Integrator.hpp"
#include "fn_ODE_RK45.hpp"
#include "fn_ODE_DOP853.hpp"

#include "assert.hpp"

namespace belfem
{
    namespace ode
    {
//------------------------------------------------------------------------------

        Integrator::Integrator(  ODE & aODE, const Type aType ) :
            mODE( aODE ),
            mType( aType )
        {
            // assign integration function
            switch( aType )
            {
                case( Type::RK45 ) :
                {
                    RK45_init( aODE, mWork );
                    mIntegrationFunction = & RK45 ;
                    break ;
                }
                case( Type::DOP853 ) :
                {
                    DOP853_init( aODE, mWork );
                    mIntegrationFunction = & DOP853 ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "unknown ode type" );
                }
            }
        }

//------------------------------------------------------------------------------

        Status
        Integrator::step( real & aT, Vector< real > & aY )
        {
            Status tStatus = ( *mIntegrationFunction )(
                    mODE,
                    aT,
                    aY,
                    mDeltaTime,
                    mWork,
                    mEpsilon,
                    mMaxNumIterations,
                    mMaxTime,
                    mAutoTimestep );

            mTime = aT ;

            return tStatus ;
        }

//------------------------------------------------------------------------------
    }
}