//
// Created by Christian Messe on 2019-01-20.
//

#ifndef BELFEM_CL_GRAPH_NODE_HPP
#define BELFEM_CL_GRAPH_NODE_HPP

#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace graph
    {
        class Vertex
        {
//------------------------------------------------------------------------------

            // id of this node
            id_t mID = gNoID;

            // index of this node
            index_t mIndex = gNoIndex;

            // owner of this node
            proc_t mOwner = gNoID;

            // level of this node
            uint mLevel = 0;

            // tells if node is flagged or not
            bool mFlag = false;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // counter for connected nodes
            uint mVertexCounter = 0;

            // container for nodes
            Vertex ** mVertices;

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

            inline void
            set_owner( const proc_t aOwner );

//------------------------------------------------------------------------------

            virtual inline proc_t
            owner() const;

//------------------------------------------------------------------------------

            inline uint &
            level();

//------------------------------------------------------------------------------

            inline const uint &
            level() const;

//------------------------------------------------------------------------------

            virtual inline void
            flag();

//------------------------------------------------------------------------------

            virtual inline void
            unflag();

//------------------------------------------------------------------------------

            virtual inline bool
            is_flagged() const;

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
            vertex( const uint & aIndex )
            {
                BELFEM_ASSERT( aIndex < mVertexCounter,
                            "Vertex Index %d out of bounds, which must be less than %d",
                            ( int ) aIndex,
                            ( int ) mVertexCounter );

                return mVertices[ aIndex ];
            }

//------------------------------------------------------------------------------

            virtual inline auto
            vertex( const uint & aIndex ) const
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

        uint &
        Vertex::level()
        {
            return mLevel;
        }

//------------------------------------------------------------------------------

        const uint &
        Vertex::level() const
        {
            return mLevel;
        }

//------------------------------------------------------------------------------

        void
        Vertex::flag()
        {
            mFlag = true;
        }

//------------------------------------------------------------------------------

        void
        Vertex::unflag()
        {
            mFlag = false;
        }

//------------------------------------------------------------------------------

        bool Vertex::is_flagged() const
        {
            return mFlag;
        }

//------------------------------------------------------------------------------

        void
        Vertex::increment_vertex_counter()
        {
            ++mVertexCounter;
        }

//------------------------------------------------------------------------------

        void
        Vertex::insert_vertex( Vertex * aVertex )
        {
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
}
#endif //BELFEM_CL_GRAPH_NODE_HPP
