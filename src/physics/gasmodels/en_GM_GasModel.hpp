//
// Created by Christian Messe on 16.09.19.
//

#ifndef BELFEM_EN_GM_GASMODEL_HPP
#define BELFEM_EN_GM_GASMODEL_HPP

namespace belfem
{
    enum class GasModel
    {
        IDGAS,
        SRK,
        PR,
        HELMHOLTZ,
        UNDEFINED
    };
}
#endif //BELFEM_EN_GM_GASMODEL_HPP
