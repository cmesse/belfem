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

#ifndef BELFEM_CL_RESISTOR_HPP
#define BELFEM_CL_RESISTOR_HPP

#include "cl_TwoTerminals.hpp"
#include "Component_Enums.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        class Resistor : public TwoTerminals
        {
            //! Resistance
            real mValue ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Resistor(const real aValue, Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            void
            compute_current() override ;

//------------------------------------------------------------------------------

            real
            get_value() const override ;

//------------------------------------------------------------------------------

            void
            set_value( const real aValue ) ;

//------------------------------------------------------------------------------

        };

        inline ComponentType
        Resistor::component_type() const
        {
            return ComponentType::RESISTOR ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_RESISTOR_HPP
