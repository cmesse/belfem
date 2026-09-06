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

#include "fn_Graph_clear.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
         * tidy up a graph
         */
        void
        clear( Cell< Vertex* > & aGraph )
        {
            for( Vertex* tVertex: aGraph )
            {
                delete tVertex;
            }
            aGraph.clear();
        }

//------------------------------------------------------------------------------
    }
}
