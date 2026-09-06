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

#ifndef BELFEM_CL_TWOTERMINALS_HPP
#define BELFEM_CL_TWOTERMINALS_HPP

#include "cl_Component.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------
        class TwoTerminals : public Component
        {
        protected:

            //! Current
            real mCurrent = 0.0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TwoTerminals(Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

//------------------------------------------------------------------------------

            real
            get_current() const override ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            void
            set_current( const real aCurrent) override ;

//------------------------------------------------------------------------------

            ElectricNode *
            get_node_plus() const override ;

//------------------------------------------------------------------------------

            ElectricNode *
            get_node_minus() const override ;

//------------------------------------------------------------------------------

        };

        inline ComponentType
        TwoTerminals::component_type() const
        {
            return ComponentType::TERMINALPAIR;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_TWOTERMINALS_HPP
