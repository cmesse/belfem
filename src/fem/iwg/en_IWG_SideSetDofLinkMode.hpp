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

#ifndef BELFEM_EN_IWG_SIDESETDOFLINKMODE_HPP
#define BELFEM_EN_IWG_SIDESETDOFLINKMODE_HPP

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * how the dofs of a sideset element are linked to its facet, master and slave elements ( six modes plus UNDEFINED )
         */
        enum class SideSetDofLinkMode
        {
            FacetOnly                = 0,
            FacetAndMaster           = 1,
            FacetAndSlave            = 2,
            MasterAndSlave           = 3,
            Cut                      = 4,
            Inactive                 = 5,
            UNDEFINED                = 6
        };

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_EN_IWG_SIDESETDOFLINKMODE_HPP
