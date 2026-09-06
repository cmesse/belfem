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

#ifndef BELFEM_CL_DIODE_HPP
#define BELFEM_CL_DIODE_HPP

#include "cl_TwoTerminals.hpp"
#include "Component_Enums.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        class Diode : public TwoTerminals
        {
            //! Saturation current
            real mIs ;

            //! Thermal voltage n·kT/q ( Shockley exponent scale )
            real mVt ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Diode(const real aIs, const real aVt, Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            void
            compute_current() override ;

//------------------------------------------------------------------------------

            real
            dIdV() const override ;

//------------------------------------------------------------------------------

        };

        inline ComponentType
        Diode::component_type() const
        {
            return ComponentType::DIODE ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_DIODE_HPP
