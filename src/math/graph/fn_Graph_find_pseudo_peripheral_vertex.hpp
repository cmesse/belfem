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

#ifndef BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_VERTEX_HPP
#define BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_VERTEX_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {

//------------------------------------------------------------------------------
        /**
        * Finds a pseudo-peripheral vertex in the graph.
        * A pseudo-peripheral vertex is one that has a large eccentricity (distance to farthest vertex).
        * This implementation uses the algorithm that repeatedly performs BFS to find vertices
        * with increasing eccentricity until no further improvement is found.
        *
        * @param aGraph     Cell containing all vertices in the graph
        * @param aStart     Starting vertex for the search (if nullptr, uses first vertex in aGraph)
        * @return           Pointer to a pseudo-peripheral vertex
        */

        Vertex *
        find_pseudo_peripheral_vertex(
                Graph & aGraph,
                      Vertex *     aStart = nullptr );

//------------------------------------------------------------------------------
    }
}
#endif // BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_VERTEX_HPP
