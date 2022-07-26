//
// Created by christian on 7/26/22.
//

#ifndef BELFEM_CL_MESH_BASIS_HPP
#define BELFEM_CL_MESH_BASIS_HPP

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
    namespace mesh
    {
        class Node ;
        class Edge ;
        class Face ;

        class Basis : public graph::Vertex
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Basis();

            virtual ~Basis() ;

//------------------------------------------------------------------------------

            /**
             * returns the type of this vertex
             */
            virtual EntityType
            entity_type() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_nodes() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_edges() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_faces() const ;

//------------------------------------------------------------------------------

            virtual Node *
            node( const uint aIndex );

//------------------------------------------------------------------------------

            virtual const Node *
            node( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            virtual Edge *
            edge( const uint aIndex );

//------------------------------------------------------------------------------

            virtual const Edge *
            edge( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            virtual Face *
            face( const uint aIndex );

//------------------------------------------------------------------------------

            virtual const Face *
            face( const uint aIndex ) const ;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MESH_BASIS_HPP
