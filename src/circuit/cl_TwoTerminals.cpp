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

#include "cl_TwoTerminals.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        TwoTerminals::TwoTerminals(Cell < ElectricNode* > & aTerminals, string aLabel) :
                Component(2,aTerminals, aLabel )
        {
        }

//----------------------------------------------------------------------------

        real
        TwoTerminals::get_current() const
        {
            return mCurrent ;
        }

//------------------------------------------------------------------------------

        void
        TwoTerminals::set_current( const belfem::real aCurrent )
        {
            mCurrent = aCurrent ;
        }

//------------------------------------------------------------------------------

        ElectricNode *
        TwoTerminals::get_node_plus() const
        {
            return mTerminals[0];
        }

//------------------------------------------------------------------------------

        ElectricNode *
        TwoTerminals::get_node_minus() const
        {
            return mTerminals[1];
        }

//------------------------------------------------------------------------------
    }
}