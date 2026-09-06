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

#ifndef BELFEM_CL_SWITCH_HPP
#define BELFEM_CL_SWITCH_HPP

#include "cl_TwoTerminals.hpp"
#include "Component_Enums.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        class Switch : public TwoTerminals
        {
            //! Is the switch closed
            bool mIsClosed ;

            //! Time to switch from on/off states
            real mSwitchTime ;

            //! Has the switch been activated yet
            bool mIsSwitched = false ;

            //! State before the last shift(), for shift_back() on a rejected timestep
            bool mPrevIsClosed ;
            bool mPrevIsSwitched = false ;



//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Switch( const bool aIsClosed , const real aSwitchTime, Cell < ElectricNode* > & aTerminals, string aLabel = "" ) ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            bool
            is_closed() const override ;

//------------------------------------------------------------------------------

            void
            switch_state() override ;

//------------------------------------------------------------------------------

            void
            shift( const real aTime, const real aDeltaTime ) override ;

//------------------------------------------------------------------------------

            void
            shift_back() override ;

//------------------------------------------------------------------------------

            void
            save_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

            void
            load_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

        };

        inline ComponentType
        Switch::component_type() const
        {
            return ComponentType::SWITCH ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SWITCH_HPP
