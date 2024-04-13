//
// Created by Christian Messe on 2019-07-28.
//

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

            Bitset<E> mEdgeOrientations ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ElementTemplate( const id_t & aID );

//------------------------------------------------------------------------------

            ~ElementTemplate();

//------------------------------------------------------------------------------

            /**
             * return the type of the element
             */
            ElementType
            type() const;

//------------------------------------------------------------------------------

            /**
             * how many nodes does this element have
             */
            uint
            number_of_nodes() const;

//------------------------------------------------------------------------------

            /**
             * how many corner nodes does this element have
             * ( for linear interpolation )
             */
            uint
            number_of_corner_nodes() const;

//------------------------------------------------------------------------------

            /**
             * how facets can this element have
             */
            uint
            number_of_facets() const;

//------------------------------------------------------------------------------

            /**
             * how faces can this element have
             */
            uint
            number_of_faces() const;

//------------------------------------------------------------------------------

            /**
             * how many edges this element has( 3D elements only )
             */
            uint
            number_of_edges() const ;

//------------------------------------------------------------------------------

            /**
             * insert a node to a position in the member array
             */
            void
            insert_node( Node * aNode, const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * return the node pointer
             */
            Node *
            node( const uint aIndex );

            const Node *
            node( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * insert an edge to a position in the member array
             */
            void
            insert_edge( Edge * aEdge, const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * insert a face to a position in the member array
             */
            void
            insert_face( Face * aFace, const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * tells if the edge container has been allocated
             */
             bool
             has_edges() const;

//------------------------------------------------------------------------------

            /**
             * tells if the face container has been allocated
             */
            bool
            has_faces() const;

//------------------------------------------------------------------------------

            /**
             * return the edge pointer
             */
            Edge *
            edge( const uint aIndex );

            /**
             * return the edge pointer (const version)
             */
            const Edge *
            edge( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the face pointer
             */
            Face *
            face( const uint aIndex );

            /**
             * return the face pointer (const version)
             */
            const Face *
            face( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * unflag all nodes of this element
             */
             void unflag_nodes();

//------------------------------------------------------------------------------

            /**
             * flag all nodes of this element
             */
            void flag_nodes();

//------------------------------------------------------------------------------

            /**
             * flag all corner nodes of this element
             */
            virtual void
            flag_corner_nodes();

//------------------------------------------------------------------------------

            /**
             * unflag all corner nodes of this element
             */
            virtual void
            unflag_corner_nodes();

//------------------------------------------------------------------------------

            /**
             * unflag all edges of this element
             */
            void
            unflag_edges();

//------------------------------------------------------------------------------

            /**
             * flag all edges of this element
             */
            void
            flag_edges();

//------------------------------------------------------------------------------

            /**
             * unflag all edges of this element
             */
            void
            unflag_faces();

//------------------------------------------------------------------------------

            /**
             * flag all edges of this element
             */
            void
            flag_faces();

//------------------------------------------------------------------------------

            void
            get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes );

//------------------------------------------------------------------------------

            void
            get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes );

//------------------------------------------------------------------------------

            void
            get_edges_of_facet( const uint aFacetIndex, Cell< Edge * > & aEdges );

//------------------------------------------------------------------------------

            void
            get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes );

//------------------------------------------------------------------------------

            /**
              * Initialize the edge container. Called by Edge Factory
              */
            void
            allocate_edge_container();

//------------------------------------------------------------------------------

            /**
              * Initialize the face container. Called by Face Factory
              */
            void
            allocate_face_container();

//------------------------------------------------------------------------------

            /**
             * display some debug information
             */
             void
             print() const ;

//------------------------------------------------------------------------------

            void
            set_edge_orientation( const uint aEdgeIndex, const bool aIsPlus );

//------------------------------------------------------------------------------


            bool
            edge_orientation( const uint aEdgeIndex ) const ;

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
            delete_node_container();

//------------------------------------------------------------------------------

            /**
             * Delete the node container. Called by destructor.
             */
            void
            delete_edge_container();

//------------------------------------------------------------------------------

            /**
             * Delete the node container. Called by destructor.
             */
            void
            delete_face_container();

//------------------------------------------------------------------------------

            /**
             * error thrown by get_nodes_of_facet if aFacetIndex < F
             */
            void
            throw_facet_error( const uint aFacetIndex );

//------------------------------------------------------------------------------

            /**
             * error thrown by get_nodes_of_edge if aEdgeIndex < E
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
            this->delete_face_container();
            this->delete_edge_container();
            this->delete_node_container();
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
        void
        ElementTemplate< N, C, E, T, F >::delete_node_container()
        {
            // delete the node container
            delete[] mNodes;
        }
//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::delete_face_container()
        {
            // check if container has been allocated
            if ( mHaveFaces )
            {
                // delete the edge container
                delete[] mFaces;

                // unset flag
                mHaveFaces = false;
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::delete_edge_container()
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

            BELFEM_ASSERT( mHaveEdges, "Edges for element %lu have not been allocated",
                          ( long unsigned int ) this->id() );

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

            BELFEM_ASSERT( mHaveFaces, "Faces for element %lu have not been allocated",
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
        ElementTemplate< N, C, E, T, F >::unflag_nodes()
        {
            for( uint k=0; k<N; ++k )
            {
                mNodes[ k ]->unflag();
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::flag_nodes()
        {
            for( uint k=0; k<N; ++k )
            {
                mNodes[ k ]->flag();
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::flag_corner_nodes()
        {
            for( uint k=0; k<C; ++k )
            {
                mNodes[ k ]->flag();
            }
        }

//------------------------------------------------------------------------------

        template< uint N, uint C, uint E, uint T, uint F >
        void
        ElementTemplate< N, C, E, T, F >::unflag_corner_nodes()
        {
            for( uint k=0; k<C; ++k )
            {
                mNodes[ k ]->unflag();
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
        ElementTemplate< N, C, E, T, F >::set_edge_orientation( const uint aEdgeIndex, const bool aIsPlus )
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
        bool
        ElementTemplate< N, C, E, T, F >::edge_orientation( const uint aEdgeIndex ) const
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
