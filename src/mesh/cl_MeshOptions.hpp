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

#ifndef BELFEM_CL_MESHOPTIONS_HPP
#define BELFEM_CL_MESHOPTIONS_HPP

#include "typedefs.hpp"

namespace belfem
{
    class MeshOptions
    {
        bool mComputeConnectivities = true ;
        bool mCreateEdges = false ;
        bool mCreateFaces = false ;
        bool mCheckElements = true ;

        string mUnitString = "mm" ;
        real   mUnitScale   = 0.001 ;

    public:

         MeshOptions() = default ;

        ~MeshOptions() = default ;
    };
}

#endif //BELFEM_CL_MESHOPTIONS_HPP
