//
// Created by christian on 8/11/21.
//

#ifndef BELFEM_EN_IWG_SIDESETDOFLINKMODE_HPP
#define BELFEM_EN_IWG_SIDESETDOFLINKMODE_HPP

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * there are four ways to link the dofs of a sideset element
         */
        enum class SideSetDofLinkMode
        {
            FacetOnly                = 0,
            FacetAndMaster           = 1,
            FacetAndSlave            = 2,
            MasterAndSlave           = 3,
            Cut                      = 4,
            Shell                    = 5,
            UNDEFINED                = 6
        };

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_EN_IWG_SIDESETDOFLINKMODE_HPP
