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

#ifndef BELFEM_CL_GM_EOS_IDGAS_HPP
#define BELFEM_CL_GM_EOS_IDGAS_HPP

#include "typedefs.hpp"
#include "cl_GM_EoS.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * this is the ideal gas equaiton of state
         *
         * \f$ p \, v = R \, T \f$
         *
         * @ingroup grp_physics_gasmodels
         * @see @ref physics_gasmodels_gasmodels_usage_guide
         */
        class EoS_Idgas : public EoS
        {
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            EoS_Idgas( Gas & aParent );

//----------------------------------------------------------------------------

            ~EoS_Idgas() = default;

//----------------------------------------------------------------------------

            void
            remix();

//----------------------------------------------------------------------------
// Thermodynamic State
//----------------------------------------------------------------------------

            real
            p( const real T, const real v ) const;

            real
            v( const real T, const real p ) const;

            real
            T( const real p, const real v ) const;

//----------------------------------------------------------------------------
// State Derivatives
//----------------------------------------------------------------------------

            real
            dpdT( const real T, const real v ) const;

//----------------------------------------------------------------------------

            real
            d2pdT2( const real T, const real v ) const;

//----------------------------------------------------------------------------

            real
            dpdv( const real T, const real v ) const;

//----------------------------------------------------------------------------

            real
            d2pdv2( const real T, const real v ) const;

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

            real
            alpha( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            beta( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            kappa( const real T, const real p ) const;

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

            real
            hdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            cpdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            sdep( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            dhdepdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            real
            dsdepdT( const real T, const real p ) const;

//------------------------------------------------------------------------------

            // pressure derivative of entropy departure
            real
            dsdepdp( const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            hdep0( const real T ) const;

//------------------------------------------------------------------------------

            real
            cpdep0( const real T ) const;

//------------------------------------------------------------------------------

            real
            sdep0( const real T ) const;

//------------------------------------------------------------------------------

            // temperature derivative of entropy departure
            real
            dsdepdT0( const real T ) const;

//------------------------------------------------------------------------------

            void
            eval_critical_point( real & T, real & p, real & v ) const;

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

            real
            v( const uint aIndex, const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            hdep( const uint aIndex, const real T, const real p ) const;

//------------------------------------------------------------------------------

            real
            cpdep( const uint aIndex, const real T, const real p ) const;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */
#endif //BELFEM_CL_GM_EOS_IDGAS_HPP
