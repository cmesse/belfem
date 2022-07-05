//
// Created by Christian Messe on 2019-01-20.
//
#include "assert.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Cell.hpp"
#include "fn_Graph_sort.hpp"
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
            if( mVertexCounter > 0 )
            {
                mVertices = new Vertex * [ mVertexCounter ];
            }
            mVertexCounter = 0;
        }

//------------------------------------------------------------------------------

        void
        Vertex::init_vertex_container( const uint aSize )
        {
            if ( aSize > 0 )
            {
                mVertices = new Vertex * [ aSize ];
            }
            mVertexCounter = 0;
        }

//------------------------------------------------------------------------------

        void
        Vertex::reset_vertex_container()
        {
            if ( mVertexCounter != 0 )
            {
                delete[] mVertices;

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
            Cell< Vertex * > tVertices( mVertexCounter, nullptr );

            for( uint k=0; k<mVertexCounter; ++k )
            {
                tVertices( k ) = mVertices[ k ];
            }

            // sort cell after index
            belfem::graph::sort( tVertices );

            // copy vertices back
            for( uint k=0; k<mVertexCounter; ++k )
            {
                mVertices[ k ] = tVertices( k );
            }
        }

//------------------------------------------------------------------------------
    }
}
