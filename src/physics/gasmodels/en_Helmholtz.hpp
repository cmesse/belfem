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

#ifndef BELFEM_EN_HELMHOLTZ_HPP
#define BELFEM_EN_HELMHOLTZ_HPP

namespace belfem
{
    enum class HelmholtzModel
    {
        ParaHydrogen   = 0,
        NormalHydrogen = 1,
        OrthoHydrogen  = 2,
        Oxygen         = 3,
        Methane        = 4,
        Nitrogen       = 5,
        UNDEFINED
    };
}
#endif //BELFEM_EN_HELMHOLTZ_HPP
