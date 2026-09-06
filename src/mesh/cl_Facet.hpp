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

            // slave element for facet
            Element * mSlave = nullptr;

            // periodic partner
            Facet * mPeriodic = nullptr ;

            // id on master element
            suint mMasterFaceID = BELFEM_SUINT_MAX;

            // id on slave element
            suint mSlaveFaceID  = BELFEM_SUINT_MAX;

            // orientation on slave element
            suint mOrientationOnSlave = BELFEM_SUINT_MAX ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Facet( Element * aElement );

//------------------------------------------------------------------------------

            ~Facet() override;

//------------------------------------------------------------------------------

            EntityType
            entity_type() const override ;

//------------------------------------------------------------------------------

            /**
             * expose the element pointer
             */
            Element *
            element();

            const Element *
            element() const;

//------------------------------------------------------------------------------

            /**
             * returns the master ( aIndex = 0 ) or slave ( aIndex = 1 ) element; see number_of_elements()
             */
            const Element *
            element( const uint aIndex ) const override ;

//------------------------------------------------------------------------------

            id_t
            master_id() const;

//------------------------------------------------------------------------------

            id_t
            slave_id() const;

//------------------------------------------------------------------------------

            auto
            index_on_master() const -> decltype( mMasterFaceID );

//------------------------------------------------------------------------------

            auto
            index_on_slave() const -> decltype( mSlaveFaceID );

//------------------------------------------------------------------------------

            /**
             * get a node pointer
             */
            Node *
            node( uint aIndex ) override ;

            const Node *
            node( uint aIndex ) const override ;

//------------------------------------------------------------------------------

            /**
             * get the number of nodes
             */
            uint
            number_of_nodes() const override ;

//------------------------------------------------------------------------------

            /**
             * get the number of corner nodes
             */
            uint
            number_of_corner_nodes() const ;

//------------------------------------------------------------------------------

            /**
             * get the number of edges
             */
            uint
            number_of_edges() const override ;

            Edge *
            edge( const uint aIndex ) override ;

            const Edge *
            edge( const uint aIndex ) const override ;

            Face *
            face( const uint aIndex=0 ) override ;

            const Face *
            face( const uint aIndex=0 ) const override ;

//------------------------------------------------------------------------------

            /**
             * unflag all nodes that belong to this facet
             */
            void
            unflag_nodes( const uint8_t aIndex=0 ) override ;

//------------------------------------------------------------------------------

            /**
             * flag all nodes that belong to this facet
             */
            void
            flag_nodes( const uint8_t aIndex=0 ) override ;

//------------------------------------------------------------------------------

            /**
             * flag all nodes of the master and the slave element
             */
            void
            flag_master_and_slave_nodes();

//------------------------------------------------------------------------------

            /**
             * flag the corner nodes of the master and the slave element
             */
            void
            flag_master_and_slave_corner_nodes();

//------------------------------------------------------------------------------

            /**
             * flag corner nodes nodes that belong to this facet
             */
            void
            flag_corner_nodes();

//------------------------------------------------------------------------------

            void
            set_owner( const proc_t aOwner ) override ;

//------------------------------------------------------------------------------

            proc_t
            owner() const override ;

//------------------------------------------------------------------------------

            void
            set_master( Element * aElement, const uint aIndex, const bool aLinkNodes=true );

//------------------------------------------------------------------------------

            void
            set_slave( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            void
            set_slave( Element * aElement, const uint aIndex, const uint aOrientation );

//------------------------------------------------------------------------------

            void
            disconnect_slave();

//------------------------------------------------------------------------------

            bool
            has_master() const;

//------------------------------------------------------------------------------

            bool
            has_slave() const;

//------------------------------------------------------------------------------

            Element *
            master();

            const Element *
            master() const ;


//------------------------------------------------------------------------------

            Element *
            slave();

            const Element *
            slave() const ;

//------------------------------------------------------------------------------

            uint
            number_of_elements() const override ;

//------------------------------------------------------------------------------

            Element *
            element( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            void
            flag( const uint8_t aIndex = 0 ) override ;

//------------------------------------------------------------------------------

            void
            unflag( const uint8_t aIndex = 0 ) override ;

//------------------------------------------------------------------------------

            bool
            is_flagged( const uint8_t aIndex = 0 ) const override;

//------------------------------------------------------------------------------

            void
            set_index( const index_t aIndex ) override ;

//------------------------------------------------------------------------------

            index_t
            index() const override ;

//------------------------------------------------------------------------------

            void
            flip();

//------------------------------------------------------------------------------

            void
            compute_orientation();

//------------------------------------------------------------------------------

            uint
            orientation_on_slave() const ;

//----------------------------------------------------------------------------

            bool
            is_curved() const ;

//----------------------------------------------------------------------------

            id_t
            sideset_id() const ;

//----------------------------------------------------------------------------

            void
            set_sideset_id( const id_t aSidesetID ) ;

//-----------------------------------------------------------------------------

            uint
            physical_tag() const ;

//-----------------------------------------------------------------------------

            void
            set_physical_tag( const uint aTag );

//-----------------------------------------------------------------------------

            void
            set_facet_counter( const uint aCount );

//-----------------------------------------------------------------------------

            void
            set_periodic( Facet * aFacet );

//-----------------------------------------------------------------------------

            Facet *
            periodic();

            const Facet *
            periodic() const ;

            bool
            is_periodic() const ;

//-----------------------------------------------------------------------------

            size_t
            memory() const ;

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

        inline
        const Element *
        Facet::element() const
        {
            return mElement;
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
        Facet::number_of_corner_nodes() const
        {
            return mElement->number_of_corner_nodes();
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
        Facet::unflag_nodes( const uint8_t aIndex )
        {
            mElement->unflag_nodes( aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        Facet::flag_nodes( const uint8_t aIndex )
        {
            mElement->flag_nodes( aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        Facet::flag_master_and_slave_nodes()
        {
            if( this->has_master() )
            {
                mMaster->flag_nodes() ;
            }
            if( this->has_slave() )
            {
                mSlave->flag_nodes() ;
            }
        }

//------------------------------------------------------------------------------

        inline void
        Facet::flag_master_and_slave_corner_nodes()
        {
            if( this->has_master() )
            {
                mMaster->flag_corner_nodes() ;
            }
            if( this->has_slave() )
            {
                mSlave->flag_corner_nodes() ;
            }
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
        Facet::index_on_master() const -> decltype( mMasterFaceID )
        {
            return mMasterFaceID;
        }

//------------------------------------------------------------------------------

        inline auto
        Facet::index_on_slave() const -> decltype( mSlaveFaceID )
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

        inline const Element *
        Facet::master() const
        {
            return mMaster;
        }

//------------------------------------------------------------------------------

        inline Element *
        Facet::slave()
        {
            return mSlave;
        }

        inline const Element *
        Facet::slave() const
        {
            return mSlave;
        }


//-----------------------------------------------------------------------------

        inline uint
        Facet::number_of_elements() const
        {
            uint aNumElems = 0 ;
            if ( mMaster != nullptr ) ++ aNumElems ;
            if ( mSlave != nullptr ) ++ aNumElems ;
            return aNumElems ;
        }

//-----------------------------------------------------------------------------

        inline Element *
        Facet::element(const uint aIndex)
        {
            if ( aIndex == 0 )
            {
                if ( mMaster != nullptr ) return mMaster ;
                if ( mSlave != nullptr ) return mSlave ;
                BELFEM_ASSERT( false, "Element index %u out of bounds.", aIndex );
                return nullptr ;
            }
            else
            {
                BELFEM_ASSERT( mMaster != nullptr, "Element index %u out of bounds.", aIndex );
                BELFEM_ASSERT( mSlave != nullptr, "Element index %u out of bounds.", aIndex );
                return mSlave ;
            }
        }

//-----------------------------------------------------------------------------

        inline const Element *
        Facet::element(const uint aIndex) const
        {
            if ( aIndex == 0 )
            {
                if ( mMaster != nullptr ) return mMaster ;
                if ( mSlave != nullptr ) return mSlave ;
                BELFEM_ASSERT( false, "Element index %u out of bounds.", aIndex );
                return nullptr ;
            }
            else
            {
                BELFEM_ASSERT( mMaster != nullptr, "Element index %u out of bounds.", aIndex );
                BELFEM_ASSERT( mSlave != nullptr, "Element index %u out of bounds.", aIndex );
                return mSlave ;
            }
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
        Facet::flag( const uint8_t aIndex )
        {
            mElement->flag( aIndex ) ;
        }

//------------------------------------------------------------------------------

        inline void
        Facet::unflag( const uint8_t aIndex )
        {
            mElement->unflag( aIndex ) ;
        }

//------------------------------------------------------------------------------

        inline bool
        Facet::is_flagged( const uint8_t aIndex ) const
        {
            return mElement->is_flagged( aIndex ) ;
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

//----------------------------------------------------------------------------

        inline bool
        Facet::is_curved() const
        {
            return mElement->is_curved() ;
        }

//----------------------------------------------------------------------------

        inline id_t
        Facet::sideset_id() const
        {
            return mElement->block_id() ;
        }

//----------------------------------------------------------------------------

        inline void
        Facet::set_physical_tag( const uint aTag )
        {
            mElement->set_physical_tag( aTag );
        }

        inline uint
        Facet::physical_tag() const
        {
            return mElement->physical_tag();
        }

//----------------------------------------------------------------------------

        inline void
        Facet::set_sideset_id( const id_t aSidesetID )
        {
            mElement->set_block_id( aSidesetID );
        }

//----------------------------------------------------------------------------

        inline void
        Facet::set_facet_counter( const uint aCount )
        {
            mFacetCounter = aCount ;
        }

//------------------------------------------------------------------------------


        inline void
        Facet::set_periodic( Facet * aFacet )
        {
            mPeriodic = aFacet ;
        }

//------------------------------------------------------------------------------

        inline Facet *
        Facet::periodic()
        {
            return mPeriodic ;
        }

//------------------------------------------------------------------------------

        inline const Facet *
        Facet::periodic() const
        {
            return mPeriodic ;
        }

//------------------------------------------------------------------------------

        inline bool
        Facet::is_periodic() const
        {
            return mPeriodic != nullptr ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FACET_HPP
