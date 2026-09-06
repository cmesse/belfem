/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

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
            this->reset_face_container();
            this->reset_facet_container();
            this->reset_element_container();
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_node_container()
        {
            if ( mNodeCapacity != 0 )
            {
                free(mNodes) ;
                mNodes = nullptr ;
                mNodeCounter = 0 ;
                mNodeCapacity = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_node_container( const uint aCounter )
        {
            BELFEM_ASSERT( mNodes == nullptr, "container is already allocated");

            uint tCount = aCounter == 0 ? mNodeCounter : aCounter ;

            // the capacity field is narrow; refuse before anything is allocated
            BELFEM_ERROR( tCount <= std::numeric_limits< decltype( mNodeCapacity ) >::max(),
                "node capacity of vertex %lu overflows its counter width ( %u requested )",
                ( long unsigned int ) this->id(), ( unsigned int ) tCount );
            if( tCount > 0 )
            {
                mNodes = ( Node ** ) malloc( tCount * sizeof( Node * ) );
                std::fill( mNodes, mNodes + tCount, nullptr );
                mNodeCounter = 0 ;
                mNodeCapacity = tCount ;
            }
        }


//------------------------------------------------------------------------------

        void
        Vertex::add_node( Node * aNode )
        {
            BELFEM_ASSERT( mNodeCounter < mNodeCapacity, "Node container is full");
            mNodes[ mNodeCounter++ ] = aNode ;
        }

//-----------------------------------------------------------------------------

        void
        Vertex::insert_node( Node * aNode, const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mNodeCapacity,
                          "Error at node %lu, Node index %u is out of bounds, must be < %u.",
                          ( long unsigned int ) this->id(),
                          ( unsigned int ) aIndex,
                          ( unsigned int ) mNodeCapacity );

            if ( mNodes[ aIndex ] == nullptr ) ++mNodeCounter ;

            mNodes[ aIndex ] = aNode ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_edge_container()
        {
            if ( mEdgeCapacity != 0 )
            {
                free( mEdges );
                mEdges = nullptr ;
                mEdgeCounter = 0 ;
                mEdgeCapacity = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_edge_container( const uint aCounter )
        {
            BELFEM_ASSERT( mEdges == nullptr, "container is already allocated");

            uint tCount = aCounter == 0 ? mEdgeCounter : aCounter ;

            // the capacity field is narrow; refuse before anything is allocated
            BELFEM_ERROR( tCount <= std::numeric_limits< decltype( mEdgeCapacity ) >::max(),
                "edge capacity of vertex %lu overflows its counter width ( %u requested )",
                ( long unsigned int ) this->id(), ( unsigned int ) tCount );

            if( tCount > 0 )
            {
                mEdges = ( Edge ** ) malloc( tCount * sizeof( Edge * ) );
                std::fill( mEdges, mEdges + tCount, nullptr );
                mEdgeCounter = 0 ;
                mEdgeCapacity = tCount ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::add_edge( Edge * aEdge )
        {
            BELFEM_ASSERT( mEdgeCounter < mEdgeCapacity, "Edge container is full");
            mEdges[ mEdgeCounter++ ] = aEdge ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::insert_edge( Edge * aEdge, const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mEdgeCapacity,
                          "Edge index %u is out of bounds, must be < %u.",
                          ( unsigned int ) aIndex,
                          ( unsigned int ) mEdgeCapacity );

            if ( mEdges[ aIndex ] == nullptr ) ++mEdgeCounter ;
            mEdges[ aIndex ] = aEdge ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_face_container()
        {
            if ( mFaceCapacity != 0 )
            {
                free( mFaces );
                mFaces = nullptr ;
                mFaceCounter = 0 ;
                mFaceCapacity = 0 ;
            }
        }

        void
        Vertex::allocate_face_container( const uint aCounter )
        {
            BELFEM_ASSERT( mFaces == nullptr, "container is already allocated");

            uint tCount = aCounter == 0 ? mFaceCounter : aCounter ;

            // the capacity field is narrow; refuse before anything is allocated
            BELFEM_ERROR( tCount <= std::numeric_limits< decltype( mFaceCapacity ) >::max(),
                "face capacity of vertex %lu overflows its counter width ( %u requested )",
                ( long unsigned int ) this->id(), ( unsigned int ) tCount );
            if( tCount > 0 )
            {
                mFaces = ( Face ** ) malloc( tCount * sizeof( Face * ) );
                std::fill( mFaces, mFaces + tCount, nullptr );
                mFaceCounter = 0 ;
                mFaceCapacity = tCount ;
            }
        }

        void
        Vertex::add_face( Face * aFace )
        {
            BELFEM_ASSERT( mFaceCounter < mFaceCapacity, "Face container is full");
            mFaces[ mFaceCounter++ ] = aFace ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_facet_container()
        {
            if ( mFacetCapacity != 0 )
            {
                free( mFacets );
                mFacets = nullptr ;
                mFacetCounter = 0 ;
                mFacetCapacity = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_facet_container( uint aNumFacets )
        {
            BELFEM_ASSERT( mFacets == nullptr, "container is already allocated");

            uint tCount = aNumFacets == 0 ? mFacetCounter : aNumFacets ;

            // the capacity field is narrow; refuse before anything is allocated
            BELFEM_ERROR( tCount <= std::numeric_limits< decltype( mFacetCapacity ) >::max(),
                "facet capacity of vertex %lu overflows its counter width ( %u requested )",
                ( long unsigned int ) this->id(), ( unsigned int ) tCount );
            if( tCount > 0 )
            {
                mFacets = ( Facet ** ) malloc( tCount * sizeof( Facet * ) );
                std::fill( mFacets, mFacets + tCount, nullptr );
                mFacetCounter = 0 ;
                mFacetCapacity = tCount ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::add_facet( Facet * aFacet )
        {
            BELFEM_ASSERT( mFacetCounter < mFacetCapacity, "Facet container is full");
            mFacets[ mFacetCounter++ ] = aFacet ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_element_container()
        {
            if ( mElementCapacity != 0 )
            {
                free( mElements );
                mElements = nullptr ;
                mElementCounter = 0 ;
                mElementCapacity = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::allocate_element_container( const uint aCounter )
        {
            BELFEM_ASSERT( mElements == nullptr, "container is already allocated");

            uint tCount = aCounter == 0 ? mElementCounter : aCounter ;

            // the capacity field is narrow; refuse before anything is allocated
            BELFEM_ERROR( tCount <= std::numeric_limits< decltype( mElementCapacity ) >::max(),
                "element capacity of vertex %lu overflows its counter width ( %u requested )",
                ( long unsigned int ) this->id(), ( unsigned int ) tCount );
            if( tCount > 0 )
            {
                mElements = ( Element ** ) malloc( tCount * sizeof( Element * ) );
                std::fill( mElements, mElements + tCount, nullptr );
                mElementCounter = 0 ;
                mElementCapacity = tCount ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::add_element( Element * aElement )
        {
            BELFEM_ASSERT( mElementCounter < mElementCapacity, "Element container is full");
            mElements[ mElementCounter++ ] = aElement;
        }

//------------------------------------------------------------------------------

        void
        Vertex::flag_nodes( const uint8_t aIndex )
        {
            for( uint k=0; k<mNodeCounter; ++k )
            {
                mNodes[ k ]->flag( aIndex ) ;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::unflag_nodes(const uint8_t aIndex )
        {
            for( uint k=0; k<mNodeCounter; ++k )
            {
                mNodes[ k ]->unflag( aIndex ) ;
            }
        }
        
//-----------------------------------------------------------------------------
    }
}