//
// Created by Christian Messe on 2019-07-25.
//

#ifndef BELFEM_CL_FACET_HPP
#define BELFEM_CL_FACET_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
        class Facet : public Vertex
        {
            // wrapped element
            Element * mElement;

            // master element for facet
            Element * mMaster = nullptr;

            // id on master element
            uint mMasterFaceID = BELFEM_UINT_MAX;

            // neighbor element
            Element * mSlave = nullptr;

            // id on slave element
            uint mSlaveFaceID = BELFEM_UINT_MAX;

            // orientation on slave element
            uint mOrientationOnSlave = BELFEM_UINT_MAX ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Facet( Element * aElement );

//------------------------------------------------------------------------------

            ~Facet();

//------------------------------------------------------------------------------

            EntityType
            entity_type() const ;

//------------------------------------------------------------------------------

            /**
             * expose the element pointer
             */
            Element *
            element();

//------------------------------------------------------------------------------

            /**
             * inherited, but forbidden function
             */
            Element *
            element( const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * inherited, but forbidden function
             */
            const Element *
            element( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            id_t
            master_id() const;

//------------------------------------------------------------------------------

            id_t
            slave_id() const;

//------------------------------------------------------------------------------

            auto
            master_index() const -> decltype( mMasterFaceID );

//------------------------------------------------------------------------------

            auto
            slave_index() const -> decltype( mSlaveFaceID );

//------------------------------------------------------------------------------

            /**
             * get a node pointer
             */
            Node *
            node( uint aIndex );

            const Node *
            node( uint aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * get the number of nodes
             */
            uint
            number_of_nodes() const;

//------------------------------------------------------------------------------

            /**
             * get the number of edges
             */
            uint
            number_of_edges() const;

            Edge *
            edge( const uint aIndex );

            const Edge *
            edge( const uint aIndex ) const ;

            Face *
            face( const uint aIndex=0 );

            const Face *
            face( const uint aIndex=0 ) const ;

//------------------------------------------------------------------------------

            /**
             * unflag all nodes that belong to this facet
             */
            void
            unflag_nodes();

//------------------------------------------------------------------------------

            /**
             * flag all nodes that belong to this facet
             */
            void
            flag_nodes();

//------------------------------------------------------------------------------

            /**
             * flag corner nodes nodes that belong to this facet
             */
            void
            flag_corner_nodes();

//------------------------------------------------------------------------------

            void
            set_owner( const proc_t aOwner );

//------------------------------------------------------------------------------

            proc_t
            owner() const;

//------------------------------------------------------------------------------

            void
            set_master( Element * aElement, const uint aIndex, const bool aLinkNodes=true );

//------------------------------------------------------------------------------

            void
            set_slave( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            bool
            has_master() const;

//------------------------------------------------------------------------------

            bool
            has_slave() const;

//------------------------------------------------------------------------------

            Element *
            master();

//------------------------------------------------------------------------------

            Element *
            slave();

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

            virtual void
            set_index( const index_t aIndex );

//------------------------------------------------------------------------------

            virtual index_t
            index() const;

//------------------------------------------------------------------------------

            void
            compute_orientation();

//------------------------------------------------------------------------------

            uint
            orientation_on_slave() const ;

//-----------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline
        EntityType
        Facet::entity_type() const
        {
            return EntityType::FACET ;
        }

//------------------------------------------------------------------------------

        inline
        Element *
        Facet::element()
        {
            return mElement;
        }

//------------------------------------------------------------------------------

        inline Element *
        Facet::element( const uint aIndex )
        {
            BELFEM_ERROR( false, "forbidden call to element( const uint aIndex ) for mesh::Facet");
            return nullptr ;
        }

//------------------------------------------------------------------------------

        inline const Element *
        Facet::element( const uint aIndex ) const
        {
            BELFEM_ERROR( false, "invalid call to element( const uint aIndex ) const for mesh::Facet");
            return nullptr ;
        }

//------------------------------------------------------------------------------

        inline Node *
        Facet::node( uint aIndex )
        {
            return mElement->node( aIndex );
        }

        inline const Node *
        Facet::node( uint aIndex ) const
        {
            return mElement->node( aIndex );
        }

//------------------------------------------------------------------------------

        inline uint
        Facet::number_of_nodes() const
        {
            return mElement->number_of_nodes();
        }

//------------------------------------------------------------------------------

        inline uint
        Facet::number_of_edges() const
        {
            return mElement->number_of_edges() ;
        }

//------------------------------------------------------------------------------

        inline Edge *
        Facet::edge( const uint aIndex )
        {
            return mElement->edge( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Edge *
        Facet::edge( const uint aIndex ) const
        {
            return mElement->edge( aIndex );
        }

//------------------------------------------------------------------------------

        inline Face *
        Facet::face( const uint aIndex )
        {
            return mElement->face( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Face *
        Facet::face( const uint aIndex ) const
        {
            return mElement->face( aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        Facet::unflag_nodes()
        {
            mElement->unflag_nodes();
        }

//------------------------------------------------------------------------------

        inline void
        Facet::flag_nodes()
        {
            mElement->flag_nodes();
        }

//------------------------------------------------------------------------------

        inline void
        Facet::flag_corner_nodes()
        {
            mElement->flag_corner_nodes();
        }

//------------------------------------------------------------------------------

        inline id_t
        Facet::master_id() const
        {
            BELFEM_ASSERT( mMaster != nullptr,
                "Facet %lu has no master element",
                          ( long unsigned int ) this->id() );
            return mMaster->index() + 1;
            //return mMaster->id();
        }

//------------------------------------------------------------------------------

        inline id_t
        Facet::slave_id() const
        {
            BELFEM_ASSERT( mSlave != nullptr,
                          "Facet %lu has no slave element",
                          ( long unsigned int ) this->id() );

            return mSlave->id();
        }

//------------------------------------------------------------------------------
        inline auto
        Facet::master_index() const -> decltype( mMasterFaceID )
        {
            return mMasterFaceID;
        }

//------------------------------------------------------------------------------

        inline auto
        Facet::slave_index() const -> decltype( mSlaveFaceID )
        {
            return mSlaveFaceID;
        }

//------------------------------------------------------------------------------

        inline bool
        Facet::has_master() const
        {
            return mMaster != nullptr;
        }

//------------------------------------------------------------------------------

        inline bool
        Facet::has_slave() const
        {
            return mSlave != nullptr;
        }

//------------------------------------------------------------------------------

        inline Element *
        Facet::master()
        {
            return mMaster;
        }

//------------------------------------------------------------------------------

        inline Element *
        Facet::slave()
        {
            return mSlave;
        }

//------------------------------------------------------------------------------

        inline void
        Facet::set_owner( const proc_t aOwner )
        {
            mElement->set_owner( aOwner );
        }

//------------------------------------------------------------------------------

        inline proc_t
        Facet::owner() const
        {
            return mElement->owner() ;
        }

//------------------------------------------------------------------------------

        inline void
        Facet::flag()
        {
            mElement->flag() ;
        }

//------------------------------------------------------------------------------

        inline void
        Facet::unflag()
        {
            mElement->unflag() ;
        }

//------------------------------------------------------------------------------

        inline bool
        Facet::is_flagged() const
        {
            return mElement->is_flagged() ;
        }

//------------------------------------------------------------------------------

        inline void
        Facet::set_index( const index_t aIndex )
        {
            mElement->set_index( aIndex );
        }

//------------------------------------------------------------------------------

        inline index_t
        Facet::index() const
        {
            return mElement->index() ;
        }
//------------------------------------------------------------------------------

        inline uint
        Facet::orientation_on_slave() const
        {
            return mOrientationOnSlave ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FACET_HPP
