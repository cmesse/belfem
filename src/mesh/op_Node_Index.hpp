//
// Created by christian on 9/14/21.
//

#ifndef BELFEM_OP_NODE_INDEX_HPP
#define BELFEM_OP_NODE_INDEX_HPP

#include "cl_Node.hpp"

namespace belfem
{
    namespace mesh
    {
        // comparision object
        struct
        {
            inline bool
            operator()( const Node * aA, const Node * aB )
            {
                return aA->index() < aB->index();
            }
        } opNodeIndex;
    }
}


#endif //TBELFEM_OP_NODE_INDEX_HPP
