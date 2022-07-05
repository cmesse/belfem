//
// Created by Christian Messe on 2019-01-20.
//

#ifndef BELFEM_FN_GRAPH_SYMCRM_HPP
#define BELFEM_FN_GRAPH_SYMCRM_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        void
        symrcm( Cell< graph::Vertex * > & aGraph,
                      graph::Vertex *     aStart            = nullptr,
                      const size_t        aNumberOfVertices = 0 );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GRAPH_SYMCRM_HPP
