//
// Created by christian on 8/28/26.
//

#ifndef BELFEM_FN_GRAPH_MULTIBFS_HPP
#define BELFEM_FN_GRAPH_MULTIBFS_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
        Vertex *
        multibfs(
            Cell< Vertex * > & aGraph,
            const proc_t aOwnerSubset=gNoOwner ) ;

    }
}
#endif //BELFEM_FN_GRAPH_MULTIBFS_HPP
