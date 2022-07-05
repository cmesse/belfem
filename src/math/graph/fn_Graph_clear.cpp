//
// Created by Christian Messe on 2019-01-26.
//

#include "fn_Graph_clear.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
         * tidy up a graph
         */
        void
        clear( Cell< Vertex* > & aGraph )
        {
            for( Vertex* tVertex: aGraph )
            {
                delete tVertex;
            }
            aGraph.clear();
        }

//------------------------------------------------------------------------------
    }
}
