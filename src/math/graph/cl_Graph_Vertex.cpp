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

#include "assert.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Cell.hpp"
#include "fn_Graph_sort.hpp"
#include "op_Graph_Vertex_Index.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        Vertex::~Vertex()
        {
            this->reset_vertex_container();
        }

//------------------------------------------------------------------------------

        void
        Vertex::init_vertex_container()
        {
            if ( mVertices != nullptr ) free( mVertices );
            if( mVertexCounter > 0 )
            {
                mVertices = ( Vertex ** ) malloc( mVertexCounter * sizeof( Vertex * ) );
            }
            else
            {
                mVertices = nullptr;
            }
            mVertexCounter = 0;
        }

//------------------------------------------------------------------------------

        void
        Vertex::init_vertex_container( const uint aSize )
        {
            BELFEM_ASSERT( aSize < std::numeric_limits<uint32_t>::max(), "Too many vertices on Vertex %lu",
                ( long unsigned int ) this->id() );

            if ( mVertices != nullptr ) free( mVertices );
            if ( aSize > 0 )
            {
                mVertices = ( Vertex ** ) malloc( aSize * sizeof( Vertex * ) );
            }
            else
            {
                mVertices = nullptr;
            }
            mVertexCounter = 0;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_vertex_container()
        {
            if ( mVertices != nullptr )
            {
                free( mVertices );
                mVertices = nullptr;
                mVertexCounter = 0;
            }
        }

//------------------------------------------------------------------------------

        void
        Vertex::init_element_container()
        {
            BELFEM_ERROR( false, "init_element_container() not available for graph::Vertex class" );
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_element_container()
        {
            BELFEM_ERROR( false, "reset_element_container() not available for graph::Vertex class" );
        }

//------------------------------------------------------------------------------

        /**
         * sorts the connected vertices according to their index
         */
        void
        Vertex::sort_vertices()
        {
            // copy vertices into temporary Cell
            Graph tVertices( mVertexCounter, nullptr );

            for( uint k=0; k<mVertexCounter; ++k )
            {
                tVertices( k ) = mVertices[ k ];
            }

            // sort cell after index
            sort( tVertices, opVertexIndex );

            // copy vertices back
            for( uint k=0; k<mVertexCounter; ++k )
            {
                mVertices[ k ] = tVertices( k );
            }
        }

        void
        Vertex::reverse_vertices()
        {
            if ( mVertexCounter == 0 ) return;
            Vertex ** tVertices
                = ( Vertex ** ) malloc( mVertexCounter * sizeof( Vertex * ) );
            for ( uint k=0; k<mVertexCounter; ++k )
            {
                tVertices[ k ] = mVertices[ k ];
            }
            uint c = mVertexCounter;
            for ( uint k=0; k<mVertexCounter; ++k )
            {
                mVertices[ k ] = tVertices[ --c ];
            }
            free( tVertices );
        }

//------------------------------------------------------------------------------
    }
}
