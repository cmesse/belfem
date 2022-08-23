//
// Created by Christian Messe on 2019-07-22.
//

#ifndef BELFEM_MESH_ELEMENT_HPP
#define BELFEM_MESH_ELEMENT_HPP

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "meshtools.hpp"
#include "Mesh_Enums.hpp"
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

            // flag telling if this element is curved
            bool mCurvedFlag = true ;

            // geometry group ID as set by gmsh
            uint mGeometryTag = 0;

            // physical group ID as set by gmsh
            uint mPhysicalTag = 0;

            // number of elements connected to this element
            uint mNumberOfElements = 0;

            // Elements connected to this element
            Element ** mElements = nullptr ;

            // Facets connected to this element (used for ghost)
            Facet ** mFacets = nullptr ;

            // id of block this element is on
            id_t mBlockID = gNoID ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Element( const id_t aID );

//------------------------------------------------------------------------------

            virtual ~Element() ;

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
            physical_tag() const;

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

            virtual inline uint
            number_of_nodes() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_corner_nodes() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_corner_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_facets() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_facets() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_faces() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_faces() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_edges() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_edges() from Element %lu",
                             ( long unsigned int ) this->id() );

                return 0;
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

            virtual Node *
            node( const uint aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function node() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual const Node *
            node( const uint aIndex ) const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const node() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual Edge *
            edge( const uint aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function edge() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual const Edge *
            edge( const uint aIndex ) const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const edge() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual Face *
            face( const uint aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function face() from Element %lu",
                             ( long unsigned int ) this->id() );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual const Face *
            face( const uint aIndex ) const
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
            unflag_nodes()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_nodes() from Element %lu",
                             ( long unsigned int ) this->id() );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_nodes()
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
            flag_corner_nodes()
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
            unflag_corner_nodes()
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
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            insert_facet( Facet * aFacet, const uint aIndex );

//------------------------------------------------------------------------------

            void
            reset_element_container();

//------------------------------------------------------------------------------

            // get the connected element
            Element *
            element( const uint aIndex );

//------------------------------------------------------------------------------

            // get the facet element
            Facet *
            facet( const uint aIndex );

//------------------------------------------------------------------------------

            uint
            number_of_elements() const;

//------------------------------------------------------------------------------

            /**
             * returns a cell, needed for DOF handling
             */
            EntityType
            entity_type() const ;

//------------------------------------------------------------------------------

            virtual void
            print() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        Element::set_block_id( const id_t aID )
        {
            mBlockID = aID ;
        }

//------------------------------------------------------------------------------

        inline void
        Element::set_geometry_tag( const uint aTag )
        {
            mGeometryTag = aTag;
        }

//------------------------------------------------------------------------------

        inline void
        Element::set_physical_tag( const uint aTag )
        {
            mPhysicalTag = aTag;
        }

//------------------------------------------------------------------------------

        inline id_t
        Element::block_id() const
        {
            return mBlockID ;
        }

//------------------------------------------------------------------------------

        inline uint
        Element::geometry_tag() const
        {
            return mGeometryTag;
        }

//------------------------------------------------------------------------------

        inline uint
        Element::physical_tag() const
        {
            return mPhysicalTag;
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
    }
}
#endif //BELFEM_MESH_ELEMENT_HPP
