//
// Created by christian on 6/9/21.
//

#ifndef BELFEM_CL_EDGE_HPP
#define BELFEM_CL_EDGE_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "cl_Vertex.hpp"

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        /**
         * \brief Special Edge class for NEDELEC-Type elements
         */
        class Edge : public Vertex
        {
//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            Edge();

//-----------------------------------------------------------------------------

            ~Edge();

//-----------------------------------------------------------------------------


            EntityType
            entity_type() const ;

//-----------------------------------------------------------------------------
        };
//-----------------------------------------------------------------------------

        inline EntityType
        Edge::entity_type() const
        {
            return EntityType::EDGE ;
        }

//-----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EDGE_HPP
