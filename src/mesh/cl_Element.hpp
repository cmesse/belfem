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
        class Element
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // id of element
            id_t mID;

            // index in memory
            index_t mIndex = gNoIndex ;

            // proc owner of element
            proc_t mOwner = 0;

            // multi purpose flag
            bool mFlag = false;

            // flag telling if this element is curved
            bool mCurvedFlag = true ;

            // geometry group ID as set by gmsh
            uint mGeometryTag = 0;

            // physical group ID as set by gmsh
            uint mPhysicalTag = 0;

            // number of elements connected to this element
            uint mNumberOfElements = 0;

            // Elements connected to this element
            Element ** mElements ;

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
            set_id( const id_t aID );

//------------------------------------------------------------------------------

            void
            set_index( const index_t aIndex );

//------------------------------------------------------------------------------

            void
            set_geometry_tag( const uint aTag );

//------------------------------------------------------------------------------

            void
            set_physical_tag( const uint aTag );

//------------------------------------------------------------------------------

            void
            set_owner( const proc_t aOwner );

//------------------------------------------------------------------------------

            void
            set_block_id( const id_t aID );

//------------------------------------------------------------------------------

            id_t
            id() const;

//------------------------------------------------------------------------------

            index_t
            index() const;

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

            proc_t
            owner() const;

//------------------------------------------------------------------------------

            void
            flag();

//------------------------------------------------------------------------------

            void
            unflag();

//------------------------------------------------------------------------------

            bool
            is_flagged() const;

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
                             "invalid call of base class function number_of_nodes() from Element %u",
                             ( unsigned int ) mID );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_corner_nodes() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_corner_nodes() from Element %u",
                             ( unsigned int ) mID );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_facets() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_facets() from Element %u",
                             ( unsigned int ) mID );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_faces() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_faces() from Element %u",
                             ( unsigned int ) mID );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline uint
            number_of_edges() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function number_of_edges() from Element %u",
                             ( unsigned int ) mID );

                return 0;
            }

//------------------------------------------------------------------------------

            virtual inline bool
            has_edges() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function has_edges() from Element %u",
                             ( unsigned int ) mID );

                return false ;
            }

//------------------------------------------------------------------------------

            virtual inline bool
            has_faces() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function has_faces() from Element %u",
                             ( unsigned int ) mID );

                return false ;
            }

//------------------------------------------------------------------------------

            virtual ElementType
            type() const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function type() from Element %u",
                             ( unsigned int ) mID );
                return ElementType::UNDEFINED;
            }

//------------------------------------------------------------------------------

            virtual Node *
            node( const uint & aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function node() from Element %u",
                             ( unsigned int ) mID );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual const Node *
            node( const uint & aIndex ) const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const node() from Element %u",
                             ( unsigned int ) mID );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual Edge *
            edge( const uint & aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function edge() from Element %u",
                             ( unsigned int ) mID );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual const Edge *
            edge( const uint & aIndex ) const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const edge() from Element %u",
                             ( unsigned int ) mID );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual Face *
            face( const uint & aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function face() from Element %u",
                             ( unsigned int ) mID );

                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual const Face *
            face( const uint & aIndex ) const
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function const face() from Element %u",
                             ( unsigned int ) mID );

                return nullptr;
            }


//------------------------------------------------------------------------------

            virtual void
            insert_node( Node * aNode, const uint & aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function insert_node() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            insert_edge( Edge * aEdge, const uint & aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function insert_edge() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            insert_face( Face * aFace, const uint & aIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function insert_face() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            unflag_nodes()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_nodes() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_nodes()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_nodes() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            /**
             * flags nodes for linear representation only
             */
            virtual void
            flag_corner_nodes()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_corner_nodes() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            /**
             * flags nodes for linear representation only
             */
            virtual void
            unflag_corner_nodes()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_corner_nodes() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            unflag_edges()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_edges() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_edges()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_edges() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            unflag_faces()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function unflag_faces() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            flag_faces()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function flag_faces() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            get_nodes_of_facet( const uint aFacetIndex, Cell<Node *> & aNodes )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_nodes_of_facet() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            get_edges_of_facet( const uint aFacetIndex, Cell< Edge * > & aEdges )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_edges_of_facet() from Element %u",
                             ( unsigned int ) mID );
            }


//------------------------------------------------------------------------------

            virtual void
            get_nodes_of_edge( const uint aEdgeIndex, Cell<Node *> & aNodes )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_nodes_of_edge() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual Node *
            node_of_edge(  const uint aNodeIndex, const uint aEdgeIndex )
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function get_node_of_edge() from Element %u",
                             ( unsigned int ) mID );
                return nullptr ;
            }

//------------------------------------------------------------------------------

            virtual void
            allocate_edge_container()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function allocate_edge_container() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            virtual void
            allocate_face_container()
            {
                BELFEM_ERROR( false,
                             "invalid call of base class function allocate_face_container() from Element %u",
                             ( unsigned int ) mID );
            }

//------------------------------------------------------------------------------

            void
            increment_element_counter();

//------------------------------------------------------------------------------

            void
            allocate_element_container( const uint aSize=0 );

//------------------------------------------------------------------------------

            void
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            reset_element_container();

//------------------------------------------------------------------------------

            // get the connected element
            Element *
            element( const uint & aIndex );

//------------------------------------------------------------------------------

            uint
            number_of_elements() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        Element::set_id( const id_t aID )
        {
            mID = aID;
        }

//------------------------------------------------------------------------------

        inline void
        Element::set_index( const index_t aIndex )
        {
            mIndex = aIndex;
        }

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

        inline void
        Element::set_owner( const proc_t aOwner )
        {
            mOwner = aOwner;
        }

//------------------------------------------------------------------------------

        inline id_t
        Element::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline index_t
        Element::index() const
        {
            return mIndex;
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

        inline proc_t
        Element::owner() const
        {
            return mOwner;
        }

//------------------------------------------------------------------------------

        inline void
        Element::flag()
        {
            mFlag = true;
        }

//------------------------------------------------------------------------------

        inline void
        Element::unflag()
        {
            mFlag = false;
        }

//------------------------------------------------------------------------------

        inline bool
        Element::is_flagged() const
        {
            return mFlag;
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
        Element::element( const uint & aIndex )
        {
            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline uint
        Element::number_of_elements() const
        {
            return mNumberOfElements ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_MESH_ELEMENT_HPP
