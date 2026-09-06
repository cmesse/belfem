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

#ifndef BELFEM_EN_FEM_GROUPACTIVATIONMODE_HPP
#define BELFEM_EN_FEM_GROUPACTIVATIONMODE_HPP

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * Defines how a fem::Block or fem::SideSet is activated
         * based on its DomainType
         */
        enum class GroupActivationMode
        {
            GeometryAndDofs = 0,  // Full DOF allocation + calculators (normal operation)
            GeometryOnly    = 1,  // Calculators only, no DOFs (for thin shell source sidesets)
            Inactive        = 2   // Neither DOFs nor calculators (completely dormant)
        };

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_EN_FEM_GROUPACTIVATIONMODE_HPP
