//
// Created by christian on 6/15/21.
//

#ifndef BELFEM_CL_VERTEX_HPP
#define BELFEM_CL_VERTEX_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Graph_Vertex.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
    namespace mesh
    {
        // forward declarations for nodes and edges
        class Node ;
        class Element;
        class Facet;
        class Edge ;

        class Vertex : public graph::Vertex
        {
        protected:
            //! pointer to nodes
            Node **mNodes;

            //! pointer to edges
            Edge **mEdges;

            //! pointer to facets
            Facet **mFacets;

            //! pointer to elements
            Element **mElements;

            //! number of nodes connected to this vertex
            uint mNodeCounter = 0 ;

            //! number of edges connected to this vertex
            uint mEdgeCounter = 0 ;

            //! number of facets connected to this vertex
            uint mFacetCounter = 0 ;

            //! number of elements connected to this vertex
            uint mElementCounter = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Vertex();

            virtual ~Vertex();

//------------------------------------------------------------------------------

            /**
             * returns the type of this mesh vertex
             */
            virtual EntityType
            entity_type() const ;

//------------------------------------------------------------------------------
//  nodes
//------------------------------------------------------------------------------

            void
            reset_node_container();

            void
            increment_node_counter();

            void
            allocate_node_container();

            void
            allocate_node_container( const uint aNumberOfNodes );

            void
            add_node( Node * aNode );

            void
            insert_node( Node * aNode, const uint aIndex );

            virtual uint
            number_of_nodes() const ;

            virtual Node *
            node( const uint aIndex );

            virtual const Node *
            node( const uint aIndex ) const ;

//------------------------------------------------------------------------------
//  edges
//------------------------------------------------------------------------------

            void
            reset_edge_container();

            void
            increment_edge_counter();

            void
            allocate_edge_container();

            void
            allocate_edge_container( const uint aNumberOfEdges );

            void
            add_edge( Edge * aEdge );

            void
            insert_edge( Edge * aEdge, const uint aIndex );

            virtual uint
            number_of_edges() const ;

            virtual Edge *
            edge( const uint aIndex );

            virtual const Edge *
            edge( const uint aIndex ) const ;

//------------------------------------------------------------------------------
//  facets
//------------------------------------------------------------------------------

            void
            reset_facet_container();

            void
            increment_facet_counter();

            void
            allocate_facet_container();

            void
            add_facet( Facet * aFacet );

            virtual uint
            number_of_facets() const ;

            virtual Facet *
            facet( const uint aIndex );

            virtual const Facet *
            facet( const uint aIndex ) const ;

//------------------------------------------------------------------------------
//  elements
//------------------------------------------------------------------------------

            void
            reset_element_container();

            void
            increment_element_counter();

            void
            allocate_element_container();

            void
            add_element( Element * aElement );

            virtual uint
            number_of_elements() const ;

            virtual Element *
            element( const uint aIndex );

            virtual const Element *
            element( const uint aIndex ) const ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * must be called by destructor
             */
            void
            delete_containers();

//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------

        inline void
        Vertex::increment_node_counter()
        {
            ++mNodeCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_edge_counter()
        {
            ++mEdgeCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_facet_counter()
        {
            ++mFacetCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_element_counter()
        {
            ++mElementCounter ;
        }

//------------------------------------------------------------------------------

        inline uint
        Vertex::number_of_nodes() const
        {
            return mNodeCounter ;
        }

//------------------------------------------------------------------------------

        inline uint
        Vertex::number_of_edges() const
        {
            return mEdgeCounter ;
        }

//------------------------------------------------------------------------------

        inline uint
        Vertex::number_of_facets() const
        {
            return mFacetCounter ;
        }

//------------------------------------------------------------------------------

        inline uint
        Vertex::number_of_elements() const
        {
            return mElementCounter ;
        }

//------------------------------------------------------------------------------

        inline Node *
        Vertex::node( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mNodeCounter,
                          "Node index %u is out of bounds for element %u, must be < %u.",
                          ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mNodes[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Node *
        Vertex::node( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mNodeCounter,
                          "Node index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mNodes[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Edge *
        Vertex::edge( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mEdgeCounter,
                          "Edge index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mEdges[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Edge *
        Vertex::edge( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mEdgeCounter,
                          "Edge index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mEdges[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Facet *
        Vertex::facet( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mFacetCounter,
                          "Facet index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mFacets[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Facet *
        Vertex::facet( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mFacetCounter,
                          "Facet index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mFacets[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Element *
        Vertex::element( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mElementCounter,
                          "Element index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Element *
        Vertex::element( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mElementCounter,
                          "Element index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCounter );

            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline EntityType
        Vertex::entity_type() const
        {
            BELFEM_ERROR( false, "invalid call to abstract class");
            return EntityType::UNDEFINED ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_VERTEX_HPP
