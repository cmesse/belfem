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

#include "cl_Resistor.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Resistor::Resistor(const real aValue, Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals,aLabel), mValue(aValue)
        {}

//----------------------------------------------------------------------------

        void
        Resistor::compute_current()
        {
            //Computing the current of the resistor from Ohms law
            this->set_current((mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/mValue);
        }

//------------------------------------------------------------------------------

        real
        Resistor::get_value() const
        {
            return mValue ;
        }

//------------------------------------------------------------------------------

        void
        Resistor::set_value( const real aValue )
        {
            mValue = std::max(aValue, 1e-10) ;
        }

//------------------------------------------------------------------------------
    }
}
