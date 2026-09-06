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

#ifndef BELFEM_CL_VERTEX_HPP
#define BELFEM_CL_VERTEX_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Mesh_Basis.hpp"
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

        class Vertex : public Basis
        {
        protected:
            //! pointer to nodes
            Node ** mNodes = nullptr ;

            //! pointer to edges
            Edge ** mEdges = nullptr ;

            //! pointer to faces
            Face ** mFaces = nullptr ;

            //! pointer to facets
            Facet ** mFacets = nullptr ;

            //! pointer to elements
            Element ** mElements = nullptr ;

            //! number of nodes connected to this vertex
            uint8_t mNodeCounter = 0 ;
            uint8_t mNodeCapacity = 0 ;

            //! number of edges connected to this vertex
            uint8_t mEdgeCounter  = 0 ;
            uint8_t mEdgeCapacity = 0 ;

            //! number of faces connected to this vertex
            uint8_t mFaceCounter = 0 ;
            uint8_t mFaceCapacity = 0 ;

            //! number of facets connected to this vertex
            //! this one needs to be larger due to cohomoligies
            uint16_t mFacetCounter = 0 ;
            uint16_t mFacetCapacity = 0 ;

            //! number of elements connected to this vertex
            uint8_t mElementCounter  = 0 ;
            uint8_t mElementCapacity = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Vertex();

            ~Vertex() override;

//------------------------------------------------------------------------------

            /**
             * returns the type of this mesh vertex
             */
            EntityType
            entity_type() const override ;

//------------------------------------------------------------------------------
//  nodes
//------------------------------------------------------------------------------

            void
            reset_node_container();

            void
            increment_node_counter();

            void
            allocate_node_container( const uint aCounter=0 );

            void
            add_node( Node * aNode );

            void
            insert_node( Node * aNode, const uint aIndex );

            uint
            number_of_nodes() const override ;

            Node *
            node( const uint aIndex ) override;

            const Node *
            node( const uint aIndex ) const override ;

//------------------------------------------------------------------------------
//  edges
//------------------------------------------------------------------------------

            void
            reset_edge_container();

            void
            increment_edge_counter();

            void
            allocate_edge_container( const uint aCounter=0 );

            void
            add_edge( Edge * aEdge );

            void
            insert_edge( Edge * aEdge, const uint aIndex );

            uint
            number_of_edges() const override ;

            Edge *
            edge( const uint aIndex ) override;

            const Edge *
            edge( const uint aIndex ) const override ;

//------------------------------------------------------------------------------
//  faces
//------------------------------------------------------------------------------

            void
            reset_face_container();

            void
            increment_face_counter();

            void
            allocate_face_container( const uint aCounter=0 );

            void
            add_face( Face * aFace );

            uint
            number_of_faces() const override ;

            Face *
            face( const uint aIndex ) override;

            const Face *
            face( const uint aIndex ) const override ;

//------------------------------------------------------------------------------
//  facets
//------------------------------------------------------------------------------

            void
            reset_facet_container();

            void
            increment_facet_counter();

            void
            allocate_facet_container( const uint aNumFacets=0 );

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
            reset_element_container() override;

            void
            increment_element_counter();

            void
            allocate_element_container( const uint aCounter=0 );

            void
            add_element( Element * aElement );

            uint
            number_of_elements() const override ;

            Element *
            element( const uint aIndex ) override;

            virtual const Element *
            element( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            virtual void
            flag_nodes( const uint8_t aIndex = 0 );

            virtual void
            unflag_nodes( const uint8_t aIndex = 0 );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * must be called by destructor
             */
            void
            delete_containers();

//------------------------------------------------------------------------------

            //! bytes of the arrays Vertex and its bases own. The five mesh
            //! arrays are charged by their stored capacity, never by
            //! number_of_*(): Facet and Face override those to describe the
            //! wrapped element or the master and slave. Graph vertices, dofs and sources are
            //! fill counts, no capacity is stored for them
            size_t
            array_memory() const ;

//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------

        inline void
        Vertex::increment_node_counter()
        {
            BELFEM_ASSERT(  mNodeCounter < 255, "Node counter overflow." );
            ++mNodeCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_edge_counter()
        {
            BELFEM_ASSERT(  mEdgeCounter < 255, "Edge counter overflow." );
            ++mEdgeCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_face_counter()
        {
            BELFEM_ASSERT(  mFaceCounter < 255, "Face counter overflow." );
            ++mFaceCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_facet_counter()
        {
            BELFEM_ASSERT(  mFacetCounter < 65535, "Facet counter overflow." );
            ++mFacetCounter ;
        }

//------------------------------------------------------------------------------

        inline void
        Vertex::increment_element_counter()
        {
            BELFEM_ASSERT(  mElementCounter < 255, "Element counter overflow." );
            ++mElementCounter ;
        }

//------------------------------------------------------------------------------

        inline uint
        Vertex::number_of_nodes() const
        {
            return mNodeCounter ;
        }

//------------------------------------------------------------------------------

        inline size_t
        Vertex::array_memory() const
        {
            size_t tMem = this->number_of_vertices() * sizeof( graph::Vertex * );

            tMem += this->number_of_dofs() * sizeof( graph::Vertex * );

            // sources and weights are parallel arrays with the same slot count
            tMem += this->number_of_sources() * ( sizeof( Basis * ) + sizeof( real ) );

            tMem += mNodeCapacity    * sizeof( Node * );
            tMem += mEdgeCapacity    * sizeof( Edge * );
            tMem += mFaceCapacity    * sizeof( Face * );
            tMem += mFacetCapacity   * sizeof( Facet * );
            tMem += mElementCapacity * sizeof( Element * );

            return tMem ;
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
            BELFEM_ASSERT( aIndex < ( uint ) mNodeCapacity,
                          "Node index %u is out of bounds for element %u, must be < %u.",
                          ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCapacity );

            return mNodes[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Node *
        Vertex::node( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < ( uint ) mNodeCapacity,
                          "Node index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mNodeCapacity );

            return mNodes[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Edge *
        Vertex::edge( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mEdgeCapacity,
                          "Edge index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mEdgeCapacity );

            return mEdges[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Edge *
        Vertex::edge( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < ( uint ) mEdgeCapacity,
                          "Edge index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mEdgeCapacity );

            return mEdges[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline uint
        Vertex::number_of_faces() const
        {
            return mFaceCounter ;
        }

        inline Face *
        Vertex::face( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mFaceCapacity,
                           "Face index %u is out of bounds for face %u, must be < %u.",
                           ( unsigned int ) aIndex,
                           ( unsigned int ) this->id(),
                           ( unsigned int ) mFaceCapacity );

            return mFaces[ aIndex ];
        }

        inline const Face *
        Vertex::face( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < ( uint ) mFaceCapacity,
                           "Face index %u is out of bounds for face %u, must be < %u.",
                           ( unsigned int ) aIndex,
                           ( unsigned int ) this->id(),
                           ( unsigned int ) mFaceCapacity );

            return mFaces[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Facet *
        Vertex::facet( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mFacetCapacity,
                          "Facet index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mFacetCapacity );

            return mFacets[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Facet *
        Vertex::facet( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < ( uint ) mFacetCapacity,
                          "Facet index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mFacetCapacity );

            return mFacets[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Element *
        Vertex::element( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mElementCapacity,
                          "Element index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mElementCapacity );

            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Element *
        Vertex::element( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < ( uint ) mElementCapacity,
                          "Element index %u is out of bounds for element %u, must be < %u.",
                                  ( unsigned int ) aIndex,
                          ( unsigned int ) this->id(),
                          ( unsigned int ) mElementCapacity );

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
