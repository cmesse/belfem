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

#ifndef BELFEM_CL_INDUCTOR_HPP
#define BELFEM_CL_INDUCTOR_HPP

#include "cl_TwoTerminals.hpp"
#include "Component_Enums.hpp"
#include "cl_BDF.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        class Inductor : public TwoTerminals
        {
            //! Inductance
            real mValue ;

            //! Time discretized resistance
            real mRL ;

            //! Time discretized current source ( zero = no history )
            real miLh = 0.0 ;

            //! BDF object
            ode::BDF * mBDF  = nullptr;

            //! Time Step Shift Register
            ShiftRegister< real > mH ;

            //! Current Shift Register
            ShiftRegister< real > mI ;

            //! Time step
            real mDeltaTime ;

            //! Order
            uint mOrder = 1 ;

//------------------------------------------------------------------------------

            /**
             * compute mRL and miLh from the current register state; shared by
             * shift() and shift_back() so a revert restores the accepted
             * step's companions exactly
             */
            void
            update_companions() ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Inductor(const real aValue, const uint aOrder, Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

            ~Inductor() override ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            void
            compute_current() override ;

//------------------------------------------------------------------------------

            void
            shift( const real aTime, const real aDeltaTime ) override ;

//------------------------------------------------------------------------------

            void
            shift_back() override ;

//------------------------------------------------------------------------------

            void
            set_timestep( const real aDeltaTime ) override ;

//------------------------------------------------------------------------------

            void
            save_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

            void
            load_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

            real
            get_value() const override ;

//------------------------------------------------------------------------------

            real
            get_discretized_resistance() const override ;

//------------------------------------------------------------------------------

            real
            get_discretized_source() const override ;

//------------------------------------------------------------------------------

        };

        inline ComponentType
        Inductor::component_type() const
        {
            return ComponentType::INDUCTOR ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_INDUCTOR_HPP
