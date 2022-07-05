//
// Created by christian on 9/14/21.
//
#include "fn_Graph_sort.hpp"
#include "op_Graph_Vertex_Index.hpp"
namespace belfem
{
    namespace graph
    {
        void
        sort( Cell< Vertex * > & aGraph )
        {
            belfem::sort( aGraph, graph::opVertexIndex );
        }

    }
}