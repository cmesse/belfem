//
// Created by Christian Messe on 05.11.19.
//

#ifndef BELFEM_OP_GRAPH_VERTEX_INDEX_HPP
#define BELFEM_OP_GRAPH_VERTEX_INDEX_HPP

#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        // comparision object
        struct
        {
            inline bool
            operator()( const Vertex * aA, const Vertex * aB )
            {
                return aA->index() < aB->index();
            }
        } opVertexIndex;

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_OP_GRAPH_VERTEX_INDEX_HPP
