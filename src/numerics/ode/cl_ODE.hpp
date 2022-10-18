//
// Created by Christian Messe on 27.12.19.
//

#ifndef BELFEM_CL_ODE_HPP
#define BELFEM_CL_ODE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "en_ODE_Status.hpp"
namespace belfem
{
    namespace ode
    {
//------------------------------------------------------------------------------

        /**
         * an object that can create an ordinary differential equation
         */
        class ODE
        {
            const uint mDimension;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ODE( const uint aDimension ) :
                mDimension( aDimension )
            {};

//------------------------------------------------------------------------------

            virtual ~ODE() = default;

//------------------------------------------------------------------------------

            virtual void
            compute(
                    const real           & aT,
                    const Vector< real > & aY,
                          Vector< real > & adYdT ) = 0;

//------------------------------------------------------------------------------

            /**
             * return the number of dimensions of this vector
             */
            inline uint
            dimension() const
            {
                return mDimension;
            }


//------------------------------------------------------------------------------

            virtual real
            check_events(
                    const real           & aT,
                    const Vector< real > & aY )
            {
                return 0.0 ;
            }

        };
    }
}
#endif //BELFEM_CL_ODE_HPP
