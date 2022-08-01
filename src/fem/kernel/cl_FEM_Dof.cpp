//
// Created by Christian Messe on 24.10.19.
//

#include "cl_FEM_Dof.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID, const uint aType, mesh::Node * aNode ) :
            Vertex(),
            mMeshBasis( aNode ),
            mTypeID( aType ),
            mIndexOnField( aNode->index() )
        {
            this->set_id( aID );
            this->set_owner( aNode->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID, const uint aType, mesh::Edge * aEdge, const index_t aDofIndexOnField ) :
                Vertex(),
                mMeshBasis( aEdge ),
                mTypeID( aType ),
                mIndexOnField( aDofIndexOnField )
        {
            this->set_id( aID );
            this->set_owner( aEdge->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID, const uint aType, mesh::Face * aFace, const index_t aDofIndexOnField ) :
                Vertex(),
                mMeshBasis( aFace ),
                mTypeID( aType ),
                mIndexOnField( aDofIndexOnField )
        {
            this->set_id( aID );
            this->set_owner( aFace->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID, const uint aType, mesh::Element * aElement, const index_t aDofIndexOnField ) :
                Vertex(),
                mMeshBasis( aElement ),
                mTypeID( aType ),
                mIndexOnField( aDofIndexOnField )
        {
            this->set_id( aID );
            this->set_owner( aElement->owner() );
        }

//------------------------------------------------------------------------------

        Dof::Dof( const id_t aID, const uint aType, mesh::Facet * aFacet ) :
        Vertex(),
        mMeshBasis( aFacet ),
        mTypeID( aType ),
        mIndexOnField( aFacet->index() )
        {
            this->set_id( aID );
            this->set_owner( aFacet->owner() );
        }

//------------------------------------------------------------------------------
    }
}