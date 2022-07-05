//
// Created by Christian Messe on 2019-01-20.
//

#ifndef BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_NODE_HPP
#define BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_NODE_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        Vertex *
        find_pseudo_peripheral_node(
                Cell< Vertex * > & aGraph,
                      Vertex *     aStart = nullptr );

//------------------------------------------------------------------------------
    }
}
#endif // BELFEM_FN_GRAPH_FIND_PSEUDO_PERIPHERAL_NODE_HPP
