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

#include "cl_Diode.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Diode::Diode(const real aIs, const real aVt, Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals,aLabel), mIs(aIs), mVt(aVt)
        {}

//----------------------------------------------------------------------------

        void
        Diode::compute_current()
        {
            //Shockley diode law
            this->set_current(mIs*(exp((mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/mVt)-1.0));
        }

//------------------------------------------------------------------------------

        real
        Diode::dIdV() const
        {
            return (mIs/mVt)*exp((mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/mVt) ;
        }

//------------------------------------------------------------------------------
    }
}