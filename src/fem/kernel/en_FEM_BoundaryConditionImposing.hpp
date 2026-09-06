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

#ifndef BELFEM_EN_FEM_BOUNDARYCONDITIONIMPOSING
#define BELFEM_EN_FEM_BOUNDARYCONDITIONIMPOSING

namespace belfem
{
    namespace fem
    {
        enum class BoundaryConditionImposing
        {
            Free      = 0,
            Dirichlet = 1,     // nodal value such as displacement or temperature
            Neumann   = 2,     // nodal flux, such as heat load or force
            Alpha     = 3,     // special type, for convective heat flux, no RADIATION !
            Lambda    = 4,     // special type for maxwell
            Weak      = 5,     // special type for maxwell
            UNDEFINED = 6
        };
    }
}
#endif // BELFEM_EN_FEM_BOUNDARYCONDITIONIMPOSING
