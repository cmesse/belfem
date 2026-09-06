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

#ifndef BELFEM_MESH_ELEMENT_HPP
#define BELFEM_MESH_ELEMENT_HPP

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"

#include "Mesh_Enums.hpp"
#include "cl_ControlPoint.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Face.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace mesh
    {
        /**
         * \brief Lagrange Element baseclass
         */
        class Element : public Basis
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! flag telling if this element is curved
            bool mCurvedFlag = true ;

            //! geometry group ID as set by gmsh
            //! physical group ID as set by gmsh
            uint16_t mElementTags[ 2 ] = { 0, 0 };

            uint16_t mNumberOfElements = 0;

            //! number of control points connected to this element
            uint8_t mNumberOfControlPoints = 0;

            //! Elements connected to this element
            Element ** mElements = nullptr ;

            //! unlike with the elements, the neighbors are in the order
            //! of the facets and may contain null pointers!
            Element ** mNeighbors = nullptr ;

            // Facets connected to this element (probably not needed anymore)
            Facet ** mFacets = nullptr ;

            //! control points connected to this element, if any
            ControlPoint ** mControlPoints = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Element( const id_t aID );

//------------------------------------------------------------------------------

            ~Element() override ;

//------------------------------------------------------------------------------

            void
            set_geometry_tag( const uint aTag );

//------------------------------------------------------------------------------

            void
            set_physical_tag( const uint aTag );

//------------------------------------------------------------------------------

            void
            set_block_id( const id_t aID );

//------------------------------------------------------------------------------

            id_t
            block_id() const;

//------------------------------------------------------------------------------

            uint
            geometry_tag() const;

//------------------------------------------------------------------------------

            uint
            physical_tag() const ;

//------------------------------------------------------------------------------

            virtual uint
            dimension() const ;

//------------------------------------------------------------------------------

            void
            set_curved_flag();

//------------------------------------------------------------------------------

            void
            unset_curved_flag();

//------------------------------------------------------------------------------

            bool
            is_curved() const;

//------------------------------------------------------------------------------

            uint
            number_of_nodes() const override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual uint
            number_of_corner_nodes() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_corner_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual uint
            number_of_facets() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_facets() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            uint
            number_of_faces() const override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_faces() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            uint
            number_of_edges() const override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_edges() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual bool
            is_thinshell() const { return false; }

//------------------------------------------------------------------------------

            uint
            number_of_control_points() const
            {
                return mNumberOfControlPoints;
            }

//------------------------------------------------------------------------------

            ControlPoint *
            control_point( const uint aIndex )
            {
                BELFEM_ASSERT( aIndex < mNumberOfControlPoints, "invalid control point index" );
                return mControlPoints[ aIndex ];
            }

   //------------------------------------------------------------------------------

            ControlPoint *
            control_point( const uint aIndex ) const
            {
                BELFEM_ASSERT( aIndex < mNumberOfControlPoints, "invalid control point index" );
                return mControlPoints[ aIndex ];
            }

//------------------------------------------------------------------------------

            virtual inline bool
            has_edges() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function has_edges() from Element %lu",
                             ( long unsigned int ) this->id() );

                return false ;
            }

//------------------------------------------------------------------------------

            virtual inline bool
            has_faces() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function has_faces() from Element %lu",
                             ( long unsigned int ) this->id() );

                return false ;
            }

//------------------------------------------------------------------------------

            virtual ElementType
            type() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function type() from Element %lu",
                             ( long unsigned int ) this->id() );
                return ElementType::UNDEFINED;
            }

//------------------------------------------------------------------------------

            Node *
            node( const uint aIndex ) override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function node() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            const Node *
            node( const uint aIndex ) const override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const node() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            Edge *
            edge( const uint aIndex ) override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function edge() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            const Edge *
            edge( const uint aIndex ) const override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const edge() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            Face *
            face( const uint aIndex ) override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function face() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            const Face *
            face( const uint aIndex ) const override
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const face() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }


//------------------------------------------------------------------------------

            virtual void
            insert_node( Node * aNode, const uint aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function insert_node() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            insert_edge( Edge * aEdge, const uint aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function insert_edge() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            insert_face( Face * aFace, const uint aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function insert_face() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            insert_control_point( ControlPoint * aControlPoint, const uint aIndex )
            {
                BELFEM_ASSERT( aIndex < mNumberOfControlPoints, "invalid control point index" );
                mControlPoints[ aIndex ] = aControlPoint;
            }
//------------------------------------------------------------------------------

            virtual void
            unflag_nodes( const uint8_t aIndex=0 )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_nodes( const uint8_t aIndex=0 )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            /**
             * flags nodes for linear representation only
             */
            virtual void
            flag_corner_nodes(  const uint8_t aIndex=0 )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_corner_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            /**
             * flags nodes for linear representation only
             */
            virtual void
            unflag_corner_nodes(  const uint8_t aIndex=0 )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_corner_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            unflag_edges()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_edges() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_edges()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_edges() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            unflag_faces()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_faces() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_faces()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_faces() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            get_nodes_of_facet( const uint aFacetIndex, Cell<Node *> & aNodes )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_nodes_of_facet() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            get_corner_nodes_of_facet( const uint aFacetIndex, Cell<Node *> & aNodes )
            {
                BELFEM_ERROR( false,
                              "invalid call of base class function get_corner_nodes_of_facet() from Element %lu",
                              ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            get_edges_of_facet( const uint aFacetIndex, Cell< Edge * > & aEdges )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_edges_of_facet() from Element %lu",
                             ( long unsigned int ) this->id() );
            }


//------------------------------------------------------------------------------

            virtual void
            get_nodes_of_edge( const uint aEdgeIndex, Cell<Node *> & aNodes )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_nodes_of_edge() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual Node *
            node_of_edge(  const uint aNodeIndex, const uint aEdgeIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_node_of_edge() from Element %lu",
                             ( long unsigned int ) this->id() );
                return nullptr ;
            }

//------------------------------------------------------------------------------

            virtual void
            allocate_edge_container()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function allocate_edge_container() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            allocate_face_container()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function allocate_face_container() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            allocate_control_points_container( const uint aSize )
            {
                BELFEM_ASSERT( mControlPoints == nullptr,
                    "control points already allocated" );

                // the count field is narrow; refuse before anything is allocated
                BELFEM_ERROR( aSize <= std::numeric_limits< decltype( mNumberOfControlPoints ) >::max(),
                    "control point count of element %lu overflows its counter width ( %u requested )",
                    ( long unsigned int ) this->id(), ( unsigned int ) aSize );

                mControlPoints = (ControlPoint **) malloc( aSize * sizeof(ControlPoint *) );

                // we need to initialize with null pointers because sometimes,
                // not all points might be set
                for ( uint i = 0; i < aSize; ++i )
                {
                    mControlPoints[ i ] = nullptr;
                }
                mNumberOfControlPoints = aSize ;
            }

//------------------------------------------------------------------------------

            void
            reset_control_point_container()
            {
                if ( mControlPoints != nullptr )
                {
                    free( mControlPoints );
                    mControlPoints = nullptr ;
                    mNumberOfControlPoints = 0 ;
                }
            }
//------------------------------------------------------------------------------

            virtual void
            set_edge_direction( const uint aEdgeIndex, const bool aIsPlus )
            {
                BELFEM_ERROR( false,
                              "invalid call of base class function set_edge_direction() from Element %lu",
                              ( long unsigned int ) this->id() );
            }


//------------------------------------------------------------------------------

            virtual bool
            edge_direction( const uint aEdgeIndex ) const
            {
                BELFEM_ERROR( false,
                              "invalid call of base class function edge_direction() from Element %lu",
                              ( long unsigned int ) this->id() );

                return false ;
            }

//------------------------------------------------------------------------------

            virtual void
            reset_edge_container()
            {
                BELFEM_ERROR( false,
                              "invalid call of base class function reset_edge_container() from Element %lu",
                              ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            reset_face_container()
            {
                BELFEM_ERROR( false,
                              "invalid call of base class function reset_face_container() from Element %lu",
                              ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            void
            increment_element_counter();

//------------------------------------------------------------------------------

            void
            allocate_element_container( const uint aSize=0 );

//------------------------------------------------------------------------------

            void
            allocate_facet_container( const uint aSize );

//------------------------------------------------------------------------------

            void
            allocate_neighbor_container();

//------------------------------------------------------------------------------

            void
            insert_neighbor( Element * aNeighbor, const uint aIndex );

//------------------------------------------------------------------------------

            void
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            insert_facet( Facet * aFacet, const uint aIndex );

//------------------------------------------------------------------------------

            void
            reset_element_container() override;

//------------------------------------------------------------------------------

            // get the connected element
            Element *
            element( const uint aIndex ) override;

//------------------------------------------------------------------------------

            // get the neighbor, warning: may return nullptr if on domain edge!
            Element *
            neighbor( const uint aIndex );

//------------------------------------------------------------------------------

            // get the facet element
            Facet *
            facet( const uint aIndex );

//------------------------------------------------------------------------------

            uint
            number_of_elements() const override;

//------------------------------------------------------------------------------

            /**
             * returns a cell, needed for DOF handling
             */
            EntityType
            entity_type() const override ;

//------------------------------------------------------------------------------

            virtual void
            print() const ;

//------------------------------------------------------------------------------

            //! bytes of this object and the arrays it owns; the concrete
            //! element overrides this with its own object size
            virtual size_t
            memory() const ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! bytes of the arrays Element and its bases own. Vertices, dofs,
            //! sources and elements are counted by their fill counter ( no
            //! capacity is stored; the two agree once the insert loops that
            //! follow each allocate_*() have run, which is the state memory()
            //! is asked in ). Control points store their capacity as the
            //! count; neighbors and facets are sized by the facet count
            size_t
            array_memory() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        Element::set_block_id( const id_t aID )
        {
            BELFEM_ASSERT( aID < std::numeric_limits<uint16_t>::max(), "Buffer overflow for geometry tag, aka block id" );
            mElementTags[ 0 ] = aID ;
        }

//------------------------------------------------------------------------------

        inline void
        Element::set_geometry_tag( const uint aTag )
        {
            BELFEM_ASSERT( aTag < std::numeric_limits<uint16_t>::max(), "Buffer overflow for geometry tag, aka block id" );
            mElementTags[ 0 ] = aTag;
        }

//------------------------------------------------------------------------------

        inline void
        Element::set_physical_tag( const uint aTag )
        {
            BELFEM_ASSERT( aTag < std::numeric_limits<uint16_t>::max(), "Buffer overflow for physical id" );
            mElementTags[ 1 ] = aTag;
        }

//------------------------------------------------------------------------------

        inline id_t
        Element::block_id() const
        {
            return mElementTags[ 0 ] ;
        }

//------------------------------------------------------------------------------

        inline uint
        Element::geometry_tag() const
        {
            return mElementTags[ 0 ];
        }

//------------------------------------------------------------------------------

        inline uint
        Element::physical_tag() const
        {
            return mElementTags[ 1 ];
        }

//------------------------------------------------------------------------------

        inline uint
        Element::dimension() const
        {
            BELFEM_ERROR( false, "invalid call to Element::dimension()" );
            return BELFEM_UINT_MAX;
        }

//------------------------------------------------------------------------------

        inline void
        Element::set_curved_flag()
        {
            mCurvedFlag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Element::unset_curved_flag()
        {
            mCurvedFlag = false ;
        }

//------------------------------------------------------------------------------

        inline bool
        Element::is_curved() const
        {
            return mCurvedFlag;
        }

//------------------------------------------------------------------------------

        // get the connected element
        inline Element *
        Element::element( const uint aIndex )
        {
            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Facet *
        Element::facet( const uint aIndex )
        {
            return mFacets[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline uint
        Element::number_of_elements() const
        {
            return mNumberOfElements ;
        }

//------------------------------------------------------------------------------

        inline EntityType
        Element::entity_type() const
        {
            return EntityType::CELL ;
        }

//------------------------------------------------------------------------------

        inline Element *
        Element::neighbor( const uint aIndex )
        {
            return mNeighbors[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline size_t
        Element::array_memory() const
        {
            size_t tMem = this->number_of_vertices() * sizeof( graph::Vertex * );

            tMem += this->number_of_dofs() * sizeof( graph::Vertex * );

            // sources and weights are parallel arrays with the same slot count
            tMem += this->number_of_sources() * ( sizeof( Basis * ) + sizeof( real ) );

            tMem += this->number_of_elements() * sizeof( Element * );
            tMem += this->number_of_control_points() * sizeof( ControlPoint * );

            // neighbors and facets exist only once allocated; both carry
            // one slot per facet
            if( mNeighbors != nullptr )
            {
                tMem += this->number_of_facets() * sizeof( Element * );
            }
            if( mFacets != nullptr )
            {
                tMem += this->number_of_facets() * sizeof( Facet * );
            }

            return tMem ;
        }

//------------------------------------------------------------------------------

        inline size_t
        Element::memory() const
        {
            return sizeof( Element ) + this->array_memory() ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_MESH_ELEMENT_HPP
