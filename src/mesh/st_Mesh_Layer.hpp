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

#ifndef BELFEM_ST_LAYER_HPP
#define BELFEM_ST_LAYER_HPP

namespace belfem
{
    namespace mesh
    {
        /**
         *  this is a supporting data type for the tape generation
         *  the mesh entities are passed over to the mesh
         *  and will be deleted by the mesh during destruction
         */
        struct Layer
        {
            // cell with cloned nodes
            Cell< Node * >    Nodes ;
            Cell< Edge * >    Edges ;
            Cell< Face * >    Faces ;
            Cell< Element * > Elements ;
        };
    }

}



#endif //BELFEM_ST_LAYER_HPP
