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

#include "cl_Superconductor.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Superconductor::Superconductor(const real aIc, const real aN, const real aEc, const real aLength, Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals, aLabel), mIc(aIc), mN(aN), mEc(aEc), mLength(aLength)
        {}

//----------------------------------------------------------------------------

        void
        Superconductor::compute_current()
        {
            //E-J power law, inverted for the current
            this->set_current((mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage()>0.0?1.0:-1.0)*
                              mIc*std::pow((std::abs(mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/(mEc*mLength)),(1.0/mN)));
        }

//------------------------------------------------------------------------------

        real
        Superconductor::dIdV() const
        {
            return 1.0*
                   std::min((mIc/(mEc*mLength*mN))*std::pow((std::abs(mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/(mEc*mLength)),(1.0/mN)-1.0),BELFEM_REAL_MAX) ;
        }

//------------------------------------------------------------------------------
    }
}