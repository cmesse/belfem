//
// Created by christian on 8/3/22.
//

#ifndef BELFEM_OP_GRAPH_VERTEX_ID_HPP
#define BELFEM_OP_GRAPH_VERTEX_ID_HPP

#include "cl_Graph_Vertex.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    // comparison object
    inline struct
    {
        inline bool
        operator()( const graph::Vertex * aA, const graph::Vertex * aB )
        {
            return aA->id() < aB->id();
        }
    } opVertexID;

//------------------------------------------------------------------------------
}


#endif //BELFEM_OP_GRAPH_VERTEX_ID_HPP
