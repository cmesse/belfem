//
// Created by Christian Messe on 2019-01-20.
//

#ifndef BELFEM_OP_GRAPH_NODE_DEGREE_HPP
#define BELFEM_OP_GRAPH_NODE_DEGREE_HPP

#include "cl_Graph_Vertex.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    // comparison object
    struct
    {
        inline bool
        operator()( const graph::Vertex * aA, const graph::Vertex * aB )
        {
            return aA->number_of_vertices() < aB->number_of_vertices();
        }
    } opVertexDegree;

//------------------------------------------------------------------------------
}
#endif //BELFEM_OP_GRAPH_NODE_DEGREE_HPP
