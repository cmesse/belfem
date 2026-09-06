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

#ifndef BELFEM_CL_ELEMENT_TEMPLATE_HPP
#define BELFEM_CL_ELEMENT_TEMPLATE_HPP

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Bitset.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        /**
          * \brief Lagrange Element templated against
          *
          * N : Number of Nodes
          * C : Number of Corner nodes
          * E : Number of Edges
          * T : Number of Facets
          * F : Number of Faces
          */
        template< uint N, uint C, uint E, uint T, uint F >
        class ElementTemplate : public Element
        {
            //! pointer to nodes
            Node **mNodes;

            //! pointer to edges
            Edge **mEdges;

            //! pointer to faces
            Face **mFaces;

            //! flag telling if edges have been allocated
            bool mHaveEdges = false ;

            //! flag telling if faces have been allocated
            bool mHaveFaces = false ;

            //! bitset stating if the edge orient is plus (true)
            //! or minus (false)
            Bitset<E> mEdgeOrientations ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ElementTemplate( const id_t & aID );

//------------------------------------------------------------------------------

            ~ElementTemplate() override;

//------------------------------------------------------------------------------

            /**
             * return the type of the element
             */
            ElementType
            type() const override;

//------------------------------------------------------------------------------

            uint
            dimension() const override ;

//------------------------------------------------------------------------------

            /**
             * how many nodes does this element have
             */
            uint
            number_of_nodes() const override;

//------------------------------------------------------------------------------

            /**
             * how many corner nodes does this element have
             * ( for linear interpolation )
             */
            uint
            number_of_corner_nodes() const override;

//------------------------------------------------------------------------------

            /**
             * how facets can this element have
             */
            uint
            number_of_facets() const override;

//------------------------------------------------------------------------------

            /**
             * how faces can this element have
             */
            uint
            number_of_faces() const override;

//------------------------------------------------------------------------------

            /**
             * how many edges this element has
             */
            uint
            number_of_edges() const override ;

//------------------------------------------------------------------------------

            /**
             * insert a node to a position in the member array
             */
            void
            insert_node( Node * aNode, const uint aIndex ) override;

//------------------------------------------------------------------------------

            /**
             * return the node pointer
             */
            Node *
            node( const uint aIndex ) override;

            const Node *
            node( const uint aIndex ) const override;

//------------------------------------------------------------------------------

            /**
             * insert an edge to a position in the member array
             */
            void
            insert_edge( Edge * aEdge, const uint aIndex ) override;

//------------------------------------------------------------------------------

            /**
             * insert a face to a position in the member array
             */
            void
            insert_face( Face * aFace, const uint aIndex ) override;

//------------------------------------------------------------------------------

            /**
             * tells if the edge container has been allocated
             */
             bool
             has_edges() const override;

//------------------------------------------------------------------------------

            size_t
            memory() const override;

//------------------------------------------------------------------------------

            /**
             * tells if the face container has been allocated
             */
            bool
            has_faces() const override;

//------------------------------------------------------------------------------

            /**
             * return the edge pointer
             */
            Edge *
            edge( const uint aIndex ) override;

            /**
             * return the edge pointer (const version)
             */
            const Edge *
            edge( const uint aIndex ) const override;

//------------------------------------------------------------------------------

            /**
             * return the face pointer
             */
            Face *
            face( const uint aIndex ) override;

            /**
             * return the face pointer (const version)
             */
            const Face *
            face( const uint aIndex ) const override;

//------------------------------------------------------------------------------

            /**
             * unflag all nodes of this element
             */
             void unflag_nodes( const uint8_t aIndex=0 ) override;

//------------------------------------------------------------------------------

            /**
             * flag all nodes of this element
             */
            void flag_nodes( const uint8_t aIndex=0 ) override;

//------------------------------------------------------------------------------

            /**
             * flag all corner nodes of this element
             */
            void
            flag_corner_nodes(  const uint8_t aIndex=0 ) override;

//------------------------------------------------------------------------------

            /**
             * unflag all corner nodes of this element
             */
            void
            unflag_corner_nodes(  const uint8_t aIndex=0 ) override;

//------------------------------------------------------------------------------

            /**
             * unflag all edges of this element
             */
            void
            unflag_edges() override;

//------------------------------------------------------------------------------

            /**
             * flag all edges of this element
             */
            void
            flag_edges() override;

//------------------------------------------------------------------------------

            /**
             * unflag all faces of this element
             */
            void
            unflag_faces() override;

//------------------------------------------------------------------------------

            /**
             * flag all faces of this element
             */
            void
            flag_faces() override;

//------------------------------------------------------------------------------

            void
            get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes ) override;

//------------------------------------------------------------------------------

            void
            get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes ) override;

//------------------------------------------------------------------------------

            void
            get_edges_of_facet( const uint aFacetIndex, Cell< Edge * > & aEdges ) override;

//------------------------------------------------------------------------------

            void
            get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes ) override;

//------------------------------------------------------------------------------

            /**
              * Initialize the edge container. Called by Edge Factory
              */
            void
            allocate_edge_container() override;

//------------------------------------------------------------------------------

            /**
              * Initialize the face container. Called by Face Factory
              */
            void
            allocate_face_container() override;

//------------------------------------------------------------------------------

            /**
             * display some debug information
             */
             void
             print() const override ;

//------------------------------------------------------------------------------

            void
            set_edge_direction( const uint aEdgeIndex, const bool aIsPlus ) override;

//------------------------------------------------------------------------------


            bool
            edge_direction( const uint aEdgeIndex ) const override ;

//------------------------------------------------------------------------------

            /**
             * Delete the edge container. Called by destructor.
             */
            void
            reset_edge_container() override;

//------------------------------------------------------------------------------

            /**
             * Delete the face container. Called by destructor.
             */
            void
            reset_face_container() override;

//------------------------------------------------------------------------------

            bool
            is_thinshell() const override;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * Initialize the node container. Called by constructor.
             */
            void
            allocate_node_container();

//------------------------------------------------------------------------------

            /**
             * Delete the node container. Called by destructor.
             */
            void
            reset_node_container();

//------------------------------------------------------------------------------

            /**
             * error thrown by the get_*_of_facet() functions if aFacetIndex >= T
             */
            void
            throw_facet_error( const uint aFacetIndex );

//------------------------------------------------------------------------------

            /**
             * error thrown by get_nodes_of_edge if aEdgeIndex >= E
             */
            void
            throw_edge_error( const uint aEdgeIndex );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::allocate_node_container()
        {
            // create the contiainer
            mNodes = new Node * [ N ];

            // populate the members with null pointers
            for ( uint k = 0; k < N; ++k )
            {
                mNodes[ k ] = nullptr;
            }

        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::allocate_edge_container()
        {
            // make sure that container is empty
            BELFEM_ASSERT( ! mHaveEdges, "Edge container for element %lu has already been allocated",
                          ( long unsigned int ) this->id() );

            // create the contiainer
            mEdges = new Edge * [ E ];

            // populate the members with null pointers
            for ( uint k = 0; k < E; ++k )
            {
                mEdges[ k ] = nullptr;
            }

            mHaveEdges = true ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::allocate_face_container()
        {
            // make sure that container is empty
            BELFEM_ASSERT( ! mHaveFaces, "Face container for element %lu has already been allocated",
                          ( long unsigned int ) this->id() );

            // create the contiainer
            mFaces = new Face * [ F ];

            // populate the members with null pointers
            for ( uint k = 0; k < F ; ++k )
            {
                mFaces[ k ] = nullptr;
            }

            mHaveFaces = true ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        ElementTemplate< N, C, E, T, F >::ElementTemplate( const id_t & aID ) :
                Element( aID )
        {
            this->allocate_node_container();

            // default setting, might be overwritten later the curved checker
            mCurvedFlag =
                    mesh::interpolation_order( this->type() )
                    != InterpolationOrder::LINEAR ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        ElementTemplate< N, C, E, T, F >::~ElementTemplate()
        {
            this->reset_face_container();
            this->reset_edge_container();
            this->reset_node_container();
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        ElementType
        ElementTemplate< N, C, E, T, F >::type() const
        {
            BELFEM_ERROR( false, "type() function not implemented for ElementTemplate< %u, %u, %u, %u, %u  > with id %lu",
                         ( unsigned int ) N,
                         ( unsigned int ) C,
                         ( unsigned int ) E,
                         ( unsigned int ) T,
                         ( unsigned int ) F,
                         ( long unsigned int ) this->id() );

            return ElementType::UNDEFINED;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        uint
        ElementTemplate< N, C, E, T, F >::number_of_nodes() const
        {
            return N;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        inline uint
        ElementTemplate< N, C, E, T, F >::number_of_corner_nodes() const
        {
            return C;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        uint
        ElementTemplate< N, C, E, T, F >::number_of_facets() const
        {
            return T;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        uint
        ElementTemplate< N, C, E, T, F >::number_of_faces() const
        {
            return F;
        }


//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        uint
        ElementTemplate< N, C, E, T, F >::number_of_edges() const
        {
            return E;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        size_t
        ElementTemplate< N, C, E, T, F >::memory() const
        {
            // the node array always exists, edges and faces only on demand
            size_t tMem = sizeof( *this ) + this->array_memory()
                        + N * sizeof( Node * );

            if( mHaveEdges )
            {
                tMem += E * sizeof( Edge * );
            }
            if( mHaveFaces )
            {
                tMem += F * sizeof( Face * );
            }

            return tMem ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::reset_node_container()
        {
            // delete the node container
            delete[] mNodes;
        }
//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::reset_face_container()
        {
            // check if container has been allocated
            if ( mHaveFaces )
            {
                // delete the face container
                delete[] mFaces;

                // unset flag
                mHaveFaces = false;
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::reset_edge_container()
        {
            // check if container has been allocated
            if ( mHaveEdges )
            {
                // delete the edge container
                delete[] mEdges;

                // unset flag
                mHaveEdges = false;
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::insert_node( Node * aNode, const uint aIndex )
        {
            // make sure that index is valid
            BELFEM_ASSERT(aIndex < N,
                         "Tried to write node into index %u of %u node element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) N,
                         ( long unsigned int ) this->id() );

            // write node into index
            mNodes[ aIndex ] = aNode;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        Node *
        ElementTemplate< N, C, E, T, F >::node( const uint aIndex )
        {
            // make sure that index is valid
            BELFEM_ASSERT(aIndex < N,
                    "Tried acces node %u of %u node element %lu",
                    ( unsigned int ) aIndex,
                    ( unsigned int ) N,
                    ( long unsigned int ) this->id() );

            // return the node
            return mNodes[ aIndex ];
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        const Node *
        ElementTemplate< N, C, E, T, F >::node( const uint aIndex ) const
        {
            // make sure that index is valid
            BELFEM_ASSERT(aIndex < N,
                         "Tried acces node %u of %u node element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) N,
                         ( long unsigned int ) this->id() );

            // return the node
            return mNodes[ aIndex ];
        }


//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        inline bool
        ElementTemplate< N, C, E, T, F >::has_edges() const
        {
            return mHaveEdges ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        inline bool
        ElementTemplate< N, C, E, T, F >::has_faces() const
        {
            return mHaveFaces ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::insert_edge( Edge * aEdge, const uint aIndex )
        {
            // make sure that index is valid
            BELFEM_ASSERT(aIndex < E,
                         "Tried to write edge into index %u of %u element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) E,
                         ( long unsigned int ) this->id() );

            // write node into index
            mEdges[ aIndex ] = aEdge;
        }

// ------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::insert_face( Face * aFace, const uint aIndex )
        {
            // make sure that index is valid
            BELFEM_ASSERT( aIndex < F,
                         "Tried to write face into index %u of %u element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) F,
                         ( long unsigned int ) this->id() );

            // write node into index
            mFaces[ aIndex ] = aFace;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        Edge *
        ElementTemplate< N, C, E, T, F >::edge( const uint aIndex )
        {

            BELFEM_ASSERT( mHaveEdges, "Edges for element %lu on block %lu have not been allocated",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) this->block_id() );

            // make sure that index is valid
            BELFEM_ASSERT( aIndex < E,
                         "Tried to access edge %u of %u node element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) E,
                         ( long unsigned int ) this->id() );

            // return the node
            return mEdges[ aIndex ];
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        const Edge *
        ElementTemplate< N, C, E, T, F >::edge( const uint aIndex ) const
        {

            BELFEM_ASSERT( mHaveEdges, "Edges for element %lu have not been allocated",
                          ( long unsigned int ) this->id() );

            // make sure that index is valid
            BELFEM_ASSERT(aIndex < E,
                         "Tried to access edge %u of %u node element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) E,
                         ( long unsigned int ) this->id() );

            // return the node
            return mEdges[ aIndex ];
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        Face *
        ElementTemplate< N, C, E, T, F >::face( const uint aIndex )
        {

            BELFEM_ASSERT( mHaveFaces, "Faces for element %lu ( block %lu ) have not been allocated",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) this->block_id() );

            // make sure that index is valid
            BELFEM_ASSERT(aIndex < F,
                         "Tried to access face %u of %u node element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) F,
                         ( long unsigned int ) this->id() );

            // return the node
            return mFaces[ aIndex ];
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        const Face *
        ElementTemplate< N, C, E, T, F >::face( const uint aIndex ) const
        {

            BELFEM_ASSERT( mHaveFaces, "Edges for element %lu have not been allocated",
                          ( long unsigned int ) this->id() );

            // make sure that index is valid
            BELFEM_ASSERT(aIndex < F,
                         "Tried to access face %u of %u node element %lu",
                         ( unsigned int ) aIndex,
                         ( unsigned int ) F,
                         ( long unsigned int ) this->id() );

            // return the node
            return mFaces[ aIndex ];
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::unflag_nodes( const uint8_t aIndex )
        {
            for( uint k=0; k<N; ++k )
            {
                mNodes[ k ]->unflag( aIndex );
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::flag_nodes( const uint8_t aIndex )
        {
            for( uint k=0; k<N; ++k )
            {
                mNodes[ k ]->flag( aIndex );
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::flag_corner_nodes( const uint8_t aIndex )
        {
            for( uint k=0; k<C; ++k )
            {
                mNodes[ k ]->flag( aIndex );
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::unflag_corner_nodes(  const uint8_t aIndex )
        {
            for( uint k=0; k<C; ++k )
            {
                mNodes[ k ]->unflag( aIndex );
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::unflag_edges()
        {
            if( mHaveEdges )
            {
                for ( uint k = 0; k < E; ++k )
                {
                    mEdges[ k ]->unflag();
                }
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::flag_edges()
        {
            if( mHaveEdges )
            {
                for ( uint k = 0; k < E; ++k )
                {
                    mEdges[ k ]->flag();
                }
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::unflag_faces()
        {
            if( mHaveFaces )
            {
                for ( uint k = 0; k < F; ++k )
                {
                    mFaces[ k ]->unflag();
                }
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::flag_faces()
        {
            if( mHaveFaces )
            {
                for ( uint k = 0; k < F; ++k )
                {
                    mFaces[ k ]->flag();
                }
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            BELFEM_ERROR( false,
                    "Function get_nodes_of_facet() not implemented for element %lu",

                          ( long unsigned int ) this->id() );
        }
//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            BELFEM_ERROR( false,
                          "Function get_corner_nodes_of_facet() not implemented for element %lu",
                          ( long unsigned int ) this->id() );
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::get_edges_of_facet( const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            BELFEM_ERROR( false,
                         "invalid call of base class function get_edges_of_facet() from element %lu",
                         ( long unsigned int ) this->id() );
        }


//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // unless this is a 3D element, this function is identical to get_nodes_of_facet
            this->get_nodes_of_facet( aEdgeIndex, aNodes );
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::throw_facet_error( const uint aFacetIndex )
        {
            BELFEM_ERROR( aFacetIndex<T,
                    "invalid facet index %u for element %lu ( must be < %u )",
                         ( unsigned int ) aFacetIndex,
                         ( long unsigned int ) this->id(),
                         ( unsigned int ) T );
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::set_edge_direction( const uint aEdgeIndex, const bool aIsPlus )
        {
            BELFEM_ASSERT( mHaveEdges, "edges have not been allocated for element %lu",
                           ( long unsigned int ) this->id() );

            BELFEM_ASSERT( aEdgeIndex<E,
                          "invalid edge index %u for element %lu ( must be < %u )",
                          ( unsigned int ) aEdgeIndex,
                          ( long unsigned int ) this->id(),
                          ( unsigned int ) E );

            if( aIsPlus )
            {
                mEdgeOrientations.set( aEdgeIndex );
            }
            else
            {
                mEdgeOrientations.reset( aEdgeIndex );
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        uint
        ElementTemplate< N, C, E, T, F >::dimension() const
        {
            BELFEM_ERROR( false, "no dimension implemented for this element" );
            return BELFEM_UINT_MAX ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        bool
        ElementTemplate< N, C, E, T, F >::edge_direction( const uint aEdgeIndex ) const
        {
            BELFEM_ASSERT( mHaveEdges, "edges have not been allocated for element %lu",
                           ( long unsigned int ) this->id() );

            BELFEM_ASSERT( aEdgeIndex<E,
                           "invalid edge index %u for element %lu ( must be < %u )",
                           ( unsigned int ) aEdgeIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) E );
            
            return mEdgeOrientations.test( aEdgeIndex );
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::throw_edge_error( const uint aEdgeIndex )
        {
            BELFEM_ERROR( aEdgeIndex<E,
                         "invalid edge index %u for element %lu ( must be < %u )",
                         ( unsigned int ) aEdgeIndex,
                         ( long unsigned int ) this->id(),
                         ( unsigned int ) E );
        }

//------------------------------------------------------------------------------

        template <uint N, uint C, uint E, uint T, uint F>
        bool ElementTemplate< N, C, E, T, F >::is_thinshell() const
        {
            return false ;
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::print() const
        {
            std::cout << "Element " << this->id() << " of type " <<
                to_string( this->type() ) << std::endl << std::endl

                    << "    Nodes : " << std::endl ;

            for( uint k=0; k<N; ++k )
            {
                std::cout << "     " << k << " " << mNodes[ k ]->id() << std::endl ;
            }
        }


//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_ELEMENT_TEMPLATE_HPP
