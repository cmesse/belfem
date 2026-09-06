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

#include "cl_ElectricNode.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        ElectricNode::ElectricNode() :
                graph::Vertex()
        {
            mVoltage = 0.0 ;
        }

//----------------------------------------------------------------------------

        void
        ElectricNode::set_voltage(const real aVoltage)
        {
            mVoltage = aVoltage ;
        }

//-----------------------------------------------------------------------------

        real
        ElectricNode::get_voltage() const
        {
            return mVoltage ;
        }

    }
//----------------------------------------------------------------------------

}
