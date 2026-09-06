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

#include "cl_JcFunction.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace material
    {
        void
        JcFunction::set_dependency( const JcParameter aParameter )
        {
            mDependency.set( static_cast< uint >( aParameter ) );
        }

        real
        JcFunction::eval(
            const real B,    // in T
            const real angle // in rad
            ) const
        {
            BELFEM_ERROR( false, "JcFunction::eval(  const real B, const real angle) is not implemented for this jc-function");
            return BELFEM_QUIET_NAN;
        }

        real
        JcFunction::eval(
            const real B,    // in T
            const real angle, // in rad
            const real T     // in K
            ) const
        {
            BELFEM_ERROR( false, "JcFunction::eval( const real B, const real angle, const real T) is not implemented for this jc-function");
            return BELFEM_QUIET_NAN;
        }

        // the zero defaults are deliberate, not stubs: zero is the exact
        // derivative of a constant function and the conservative fallback
        // for every implementation that has not provided an analytic form
        // ( see the contract on the declarations )

        real
        JcFunction::deval_dB(
            const real B,
            const real angle,
            const real T
            ) const
        {
            return 0.0 ;
        }

        real
        JcFunction::deval_dbeta(
            const real B,
            const real angle,
            const real T
            ) const
        {
            return 0.0 ;
        }

        real
        JcFunction::deval_dT(
            const real B,
            const real angle,
            const real T
            ) const
        {
            return 0.0 ;
        }

    }
}