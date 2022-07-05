//
// Created by Christian Messe on 28.12.19.
//

#include "cl_ODE_Integrator.hpp"
#include "fn_ODE_RK45.hpp"
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
            return ( *mIntegrationFunction )(
                    mODE,
                    aT,
                    aY,
                    mDeltaTime,
                    mWork,
                    mEpsilon,
                    mMaxNumIterations,
                    mMaxTime,
                    mAutoTimestep );
        }

//------------------------------------------------------------------------------
    }
}