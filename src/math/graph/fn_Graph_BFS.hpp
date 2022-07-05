//
// Created by Christian Messe on 2019-01-20.
//

#ifndef BELFEM_FN_GRAPH_BFS_HPP
#define BELFEM_FN_GRAPH_BFS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        void
        bfs(    Cell< Vertex * > & aGraph,
                uint & aWidth,
                Vertex * aStart = nullptr );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GRAPH_BFS_HPP
