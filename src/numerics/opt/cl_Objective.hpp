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

#ifndef BELFEM_CL_OBJECTIVE_HPP
#define BELFEM_CL_OBJECTIVE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace opt
    {
//------------------------------------------------------------------------------

        /**
         * an object that defines a scalar objective function to be optimized.
         *
         * This base class is agnostic of the optimization algorithm: the user
         * derives from it and implements compute_objective(). It plays the same
         * role for the Optimizer that ode::ODE plays for the ode::Integrator.
         */
        class Objective
        {
            // number of design variables
            const uint mDimension;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Objective( const uint aDimension ) :
                mDimension( aDimension )
            {};

//------------------------------------------------------------------------------

            virtual ~Objective() = default;

//------------------------------------------------------------------------------

            /**
             * evaluate the objective at the design point aX.
             *
             * @param aX        design vector, length dimension()
             * @param aGradient gradient output. For derivative-free algorithms
             *                  this vector has length zero and must be ignored.
             *                  For gradient-based algorithms it has length
             *                  dimension() and the implementation must fill it
             *                  with d(objective)/d(aX).
             * @return          the scalar objective value at aX
             */
            virtual real
            compute_objective(
                    const Vector< real > & aX,
                          Vector< real > & aGradient ) = 0;

//------------------------------------------------------------------------------

            /**
             * return the number of design variables
             */
            inline uint
            dimension() const
            {
                return mDimension;
            }

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_OBJECTIVE_HPP
