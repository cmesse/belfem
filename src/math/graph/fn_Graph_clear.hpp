//
// Created by Christian Messe on 2019-01-26.
//

#ifndef BELFEM_FN_GRAPH_CLEAR_HPP
#define BELFEM_FN_GRAPH_CLEAR_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
         * tidy up a graph
         */
         void
         clear( Cell< Vertex* > & aGraph );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_GRAPH_CLEAR_HPP
