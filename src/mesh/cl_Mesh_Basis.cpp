//
// Created by christian on 7/26/22.
//

#include "cl_Mesh_Basis.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//----------------------------------------------------------------------------

        Basis::Basis() :
                graph::Vertex()
        {
            // mesh data always assume that owner is zero at first
            this->set_owner( 0 );
        }

//----------------------------------------------------------------------------

        Basis::~Basis()
        {

        }

//------------------------------------------------------------------------------

        EntityType
        Basis::entity_type() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : entity_type()");
            return EntityType::UNDEFINED ;
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_nodes() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_nodes()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_edges() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_edges()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_faces() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_faces()" );
            return 0 ;
        }

//------------------------------------------------------------------------------

        Node *
        Basis::node( const uint aIndex )
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : node()" );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        const Node *
        Basis::node( const uint aIndex ) const
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : node() const" );
            return nullptr ;
        }


//------------------------------------------------------------------------------

        Edge *
        Basis::edge( const uint aIndex )
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : edge()" );
            return nullptr ;
        }


//------------------------------------------------------------------------------

        const Edge *
        Basis::edge( const uint aIndex ) const
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : edge() const" );
            return nullptr ;
        }


//------------------------------------------------------------------------------

        Face *
        Basis::face( const uint aIndex )
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : face()" );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        const Face *
        Basis::face( const uint aIndex ) const
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : face() const" );
            return nullptr ;
        }
        
//----------------------------------------------------------------------------
    }
}