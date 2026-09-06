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

#ifndef EN_MAXWELL_FORMULATIONS_HPP
#define EN_MAXWELL_FORMULATIONS_HPP


namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            enum class Formulation
            {
                HPhi,
                L2PhiH,
                L2PhiB,
                L2EdgeH,
                UNDEFINED
            };
        }
    }
}
#endif //EN_MAXWELL_FORMULATIONS_HPP
