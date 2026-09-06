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

#ifndef BELFEM_CL_GRAPH_NODE_HPP
#define BELFEM_CL_GRAPH_NODE_HPP

#include <limits>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
    namespace graph
    {
        /**
         * @brief Graph node with an adjacency list; the vertex type the graph algorithms operate on.
         *
         * @ingroup grp_math_graph
         * @see @ref math_graph_graph_usage_guide
         */
        class Vertex
        {
//------------------------------------------------------------------------------

            // id of this node
            id_t mID = gNoID;

            // index of this node
            index_t mIndex = gNoIndex;

            // owner of this node
            proc_t mOwner = gNoOwner;

            // level of this node
            index_t mLevel = 0;

            // bitwise flags for this vertex
            uint8_t mFlags = 0;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // counter for connected vertices
            uint32_t mVertexCounter = 0;

            // container for other vertices
            Vertex ** mVertices = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Vertex() = default;

//------------------------------------------------------------------------------

            virtual ~Vertex();

//------------------------------------------------------------------------------

            inline void
            set_id( const id_t aID );

//------------------------------------------------------------------------------

            inline id_t
            id() const;

//------------------------------------------------------------------------------

            virtual inline void
            set_index( const index_t aIndex );

//------------------------------------------------------------------------------

            virtual inline index_t
            index() const;

//------------------------------------------------------------------------------

            virtual inline void
            set_owner( const proc_t aOwner );

//------------------------------------------------------------------------------

            virtual inline proc_t
            owner() const;

//------------------------------------------------------------------------------

            inline index_t
            level() const ;

//------------------------------------------------------------------------------

            inline void
            set_level( const index_t aLevel );

//------------------------------------------------------------------------------

            virtual inline void
            flag( const uint8_t aIndex = 0 );

//------------------------------------------------------------------------------

            virtual inline void
            unflag( const uint8_t aIndex = 0 );

//------------------------------------------------------------------------------

            virtual inline bool
            is_flagged( const uint8_t aIndex = 0 ) const;

//------------------------------------------------------------------------------

            inline void
            increment_vertex_counter();

//------------------------------------------------------------------------------

            void
            init_vertex_container();

//------------------------------------------------------------------------------

            void
            init_vertex_container( const uint aSize );

//------------------------------------------------------------------------------

            void
            reset_vertex_container();

//------------------------------------------------------------------------------

            inline void
            insert_vertex( Vertex * aVertex );

//------------------------------------------------------------------------------

            inline uint
            number_of_vertices() const;

//------------------------------------------------------------------------------

            virtual inline Vertex *
            vertex( const uint aIndex )
            {
                BELFEM_ASSERT( aIndex < mVertexCounter,
                            "Vertex Index %d out of bounds, which must be less than %d",
                            ( int ) aIndex,
                            ( int ) mVertexCounter );

                return mVertices[ aIndex ];
            }

//------------------------------------------------------------------------------

            virtual inline auto
            vertex( const uint aIndex ) const
                -> decltype( this )
            {
                BELFEM_ASSERT( aIndex < mVertexCounter,
                            "Vertex Index %d out of bounds, which must be less than %d",
                            ( int ) aIndex,
                            ( int ) mVertexCounter );

                return mVertices[ aIndex ];
            }

//------------------------------------------------------------------------------

            virtual void
            init_element_container();

//------------------------------------------------------------------------------

            virtual void
            reset_element_container();

//------------------------------------------------------------------------------

            /**
             * sorts the vertices according to their index
             */
            void
            sort_vertices();

//------------------------------------------------------------------------------

            void
            reverse_vertices();

        };

//------------------------------------------------------------------------------

        void
        Vertex::set_id( const id_t aID )
        {
            mID = aID;
        }

//------------------------------------------------------------------------------

        id_t
        Vertex::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

         void
         Vertex::set_index( const index_t aIndex )
         {
            mIndex = aIndex;
         }

//------------------------------------------------------------------------------

        index_t
        Vertex::index() const
        {
            return mIndex;
        }

//------------------------------------------------------------------------------

        void
        Vertex::set_owner( const proc_t aOwner )
        {
            mOwner = aOwner;
        }

//------------------------------------------------------------------------------

        proc_t
        Vertex::owner() const
        {
            return mOwner;
        }

//------------------------------------------------------------------------------

        index_t
        Vertex::level() const
        {
            return mLevel ;
        }

//------------------------------------------------------------------------------

        inline void Vertex::set_level(const index_t aLevel )
        {
            mLevel = aLevel ;
        }

//------------------------------------------------------------------------------

        void
        Vertex::flag( const uint8_t aIndex )
        {
            BELFEM_ASSERT( aIndex < 8, "Flag index %u out of bounds (must be < 8)",
                           ( unsigned int ) aIndex );
            mFlags |= ( 1 << aIndex );
        }

//------------------------------------------------------------------------------

        void
        Vertex::unflag( const uint8_t aIndex )
        {
            BELFEM_ASSERT( aIndex < 8, "Flag index %u out of bounds (must be < 8)",
                           ( unsigned int ) aIndex );
            mFlags &= ~( 1 << aIndex );
        }

//------------------------------------------------------------------------------

        bool
        Vertex::is_flagged( const uint8_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < 8, "Flag index %u out of bounds (must be < 8)",
                           ( unsigned int ) aIndex );
            return mFlags & ( 1 << aIndex );
        }

//------------------------------------------------------------------------------

        void
        Vertex::increment_vertex_counter()
        {
            BELFEM_ASSERT( mVertexCounter < std::numeric_limits<uint32_t>::max(), "Vertex counter overflow." );
            ++mVertexCounter;
        }

//------------------------------------------------------------------------------

        void
        Vertex::insert_vertex( Vertex * aVertex )
        {
            BELFEM_ASSERT( mVertexCounter < std::numeric_limits<uint32_t>::max(), "Vertex counter overflow." );
            mVertices[ mVertexCounter++ ] = aVertex;
        }

//------------------------------------------------------------------------------

        uint
        Vertex::number_of_vertices() const
        {
            return mVertexCounter;
        }

//------------------------------------------------------------------------------
    }

    typedef Cell< graph::Vertex * > Graph;

}
#endif //BELFEM_CL_GRAPH_NODE_HPP
