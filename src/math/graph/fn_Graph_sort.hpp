//
// Created by christian on 9/14/21.
//

#ifndef BELFEM_FN_GRAPH_SORT_HPP
#define BELFEM_FN_GRAPH_SORT_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        void
        sort( Cell< Vertex * > & aGraph );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GRAPH_SORT_HPP
