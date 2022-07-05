//
// Created by Christian Messe on 17.08.20.
//

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
        UNDEFINED
    };
}
#endif //BELFEM_EN_HELMHOLTZ_HPP
