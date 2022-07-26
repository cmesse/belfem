//
// Created by christian on 6/15/21.
//

#include "cl_Vertex.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Facet.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
//----------------------------------------------------------------------------

        Vertex::Vertex() :
            Basis()
        {}

//----------------------------------------------------------------------------

        Vertex::~Vertex()
        {
            this->delete_containers();
        }

//----------------------------------------------------------------------------

        void
        Vertex::delete_containers()
        {
            this->reset_vertex_container();
            this->reset_node_container();
            this->reset_edge_container();
            this->reset_facet_container();
            this->reset_element_container();
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_node_container()
        {
            if ( mNodeCounter != 0 )
            {
                delete [] mNodes ;
                mNodeCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_node_container()
        {
            if( mNodeCounter > 0 )
            {
                mNodes = new Node * [ mNodeCounter ];
                mNodeCounter = 0 ;
            }
        }


//------------------------------------------------------------------------------

        void
        Vertex::allocate_node_container( const uint aNumberOfNodes )
        {
            BELFEM_ASSERT( mNodeCounter == 0, "container is already allocated");

            if( aNumberOfNodes > 0 )
            {
                mNodes = new Node * [ aNumberOfNodes ];
                mNodeCounter = aNumberOfNodes ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::add_node( Node * aNode )
        {
            mNodes[ mNodeCounter++ ] = aNode ;
        }

//-----------------------------------------------------------------------------

        void
        Vertex::insert_node( Node * aNode, const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mNodeCounter,
                          "Node index %u is out of bounds, must be < %u.",
                          ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            mNodes[ aIndex ] = aNode ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_edge_container()
        {
            if ( mEdgeCounter != 0 )
            {
                delete [] mEdges ;
                mEdgeCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_edge_container()
        {
            if( mEdgeCounter > 0 )
            {
                mEdges = new Edge * [ mEdgeCounter ];
                mEdgeCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_edge_container( const uint aNumberOfEdges )
        {
            BELFEM_ASSERT( mEdgeCounter == 0, "container is already allocated");

            if( aNumberOfEdges > 0 )
            {
                mEdges = new Edge * [ aNumberOfEdges ];
                mEdgeCounter = aNumberOfEdges ;
            }

        }

//------------------------------------------------------------------------------

        void
        Vertex::add_edge( Edge * aEdge )
        {
            mEdges[ mEdgeCounter++ ] = aEdge ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::insert_edge( Edge * aEdge, const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mEdgeCounter,
                          "Edge index %u is out of bounds, must be < %u.",
                          ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mEdgeCounter );

            mEdges[ aIndex ] = aEdge ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_facet_container()
        {
            if ( mFacetCounter != 0 )
            {
                delete [] mFacets ;
                mFacetCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_facet_container()
        {
            if( mFacetCounter > 0 )
            {
                mFacets = new Facet * [ mFacetCounter ];
                mFacetCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::add_facet( Facet * aFacet )
        {
            mFacets[ mFacetCounter++ ] = aFacet ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_element_container()
        {
            if ( mElementCounter != 0 )
            {
                delete [] mElements ;
                mElementCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_element_container()
        {
            if( mElementCounter > 0 )
            {
                mElements = new Element * [ mElementCounter ];
                mElementCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::add_element( Element * aElement )
        {
            mElements[ mElementCounter++ ] = aElement;
        }

//-----------------------------------------------------------------------------
    }
}