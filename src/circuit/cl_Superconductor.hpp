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

#ifndef BELFEM_CL_SUPERCONDUCTOR_HPP
#define BELFEM_CL_SUPERCONDUCTOR_HPP

#include "cl_TwoTerminals.hpp"
#include "Component_Enums.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        class Superconductor : public TwoTerminals
        {
            //! Critical current
            real mIc ;

            //! n parameter
            real mN ;

            //! electric field parameter
            real mEc ;

            //! Length of the superconductor
            real mLength ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Superconductor(const real aIc, const real aN, const real aEc, const real aLength, Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

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
        Superconductor::component_type() const
        {
            return ComponentType::SUPERCONDUCTOR ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SUPERCONDUCTOR_HPP
