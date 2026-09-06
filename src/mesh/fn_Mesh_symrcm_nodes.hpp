/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_MESH_SYMRCM_NODES_HPP
#define BELFEM_FN_MESH_SYMRCM_NODES_HPP

#include "cl_Cell.hpp"
#include "cl_Node.hpp"

namespace belfem
{
    namespace mesh
    {
        void
        symrcm( Cell< Node * > & aNodes );
    }
}
#endif //BELFEM_FN_MESH_SYMRCM_NODES_HPP