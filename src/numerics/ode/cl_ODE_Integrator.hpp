//
// Created by Christian Messe on 28.12.19.
//

#ifndef BELFEM_CL_ODE_INTEGRATOR_HPP
#define BELFEM_CL_ODE_INTEGRATOR_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "en_ODE_Type.hpp"
#include "cl_ODE.hpp"

namespace belfem
{
    namespace ode
    {
        class Integrator
        {
            // equation object
            ODE & mODE ;

            // type of this integrator
            Type mType ;

            // memory for integration equation
            Cell< Vector< real > > mWork ;

            // convergence criterion
            real mEpsilon = 1.0e-7 ;

            // maximum number of iterations per step
            uint mMaxNumIterations = 1000 ;

            // maximum T-value
            real mMaxTime = BELFEM_REAL_MAX ;

            // current timestep
            real mTime = 0.0;

            // timestep
            real mDeltaTime = 1.0;

            // flag if timestep is adapted
            bool mAutoTimestep = true ;

            Status
            ( * mIntegrationFunction ) ( ODE         & aODE,
                      real           & aT,
                      Vector< real > & aY,
                      real           & aStep,
                      Cell< Vector< real > > & aWork,
                      const real aEpsilon ,
                      const uint aMaxIterations ,
                      const real aXmax,
                      const bool aAutoTimestep ) ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Integrator( ODE & aObject, const Type aType );

//------------------------------------------------------------------------------

            ~Integrator() = default ;

//------------------------------------------------------------------------------

            /**
             * return the convergence criterion
             */
             inline real &
             epsilon();

//------------------------------------------------------------------------------

            /**
             * set the maximum number of iterations
             */
            inline uint &
            max_num_iterations();

//------------------------------------------------------------------------------

            /**
             * return t-coordinate
             */
            inline real &
            time();

            inline real &
            maxtime();

//------------------------------------------------------------------------------

            /**
             * return t-coordinate
             */
            inline real &
            timestep();

//------------------------------------------------------------------------------

            /**
             * perform one single step
             */
            Status
            step( real & aT, Vector< real > & aY );

//------------------------------------------------------------------------------

            /**
             * return the type
             */
            const Type &
            type() const ;

//------------------------------------------------------------------------------

            /**
             * set the autotimestep flag
             */
             void
             set_auto_timestep( const bool aAutoTimestep=true );

//-----------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------

        inline real &
        Integrator::epsilon()
        {
            return mEpsilon;
        }

//------------------------------------------------------------------------------

        inline uint &
        Integrator::max_num_iterations()
        {
            return mMaxNumIterations;
        }

//------------------------------------------------------------------------------

        inline real &
        Integrator::time()
        {
            return mTime;
        }

//------------------------------------------------------------------------------

        inline real &
        Integrator::maxtime()
        {
            return mMaxTime;
        }

//------------------------------------------------------------------------------

        inline real &
        Integrator::timestep()
        {
            return mDeltaTime;
        }

//------------------------------------------------------------------------------

        inline const Type &
        Integrator::type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------

        inline void
        Integrator::set_auto_timestep( const bool aAutoTimestep )
        {
            mAutoTimestep = aAutoTimestep ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_ODE_INTEGRATOR_HPP
