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

#include "cl_Facet.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Facet::Facet( Element * aElement ) :
                mElement( aElement )
        {
            this->set_id( aElement->id() );
            this->set_owner( aElement->owner() );
        }

//------------------------------------------------------------------------------

        Facet::~Facet()
        {
            delete mElement;
        }

//------------------------------------------------------------------------------

        void Facet::set_master( Element * aElement, const uint aIndex, const bool aLinkNodes )
        {
            mMaster = aElement;
            mMasterFaceID = aIndex;

            if( aLinkNodes )
            {
                // temporary array containing nodes
                Cell< Node * > tNodes;

                aElement->get_nodes_of_facet( aIndex, tNodes );

                uint tNumNodes = tNodes.size();

                // write nodes in new order onto mesh
                for ( uint k = 0; k < tNumNodes; ++k )
                {
                    mElement->insert_node( tNodes( k ), k );
                }
            }
        }

//------------------------------------------------------------------------------

        void Facet::set_slave( Element * aElement, const uint aIndex )
        {
            mSlave = aElement;
            mSlaveFaceID = aIndex;
        }

//------------------------------------------------------------------------------

        void
        Facet::set_slave( Element * aElement, const uint aIndex, const uint aOrientation )
        {
            mSlave = aElement;
            mSlaveFaceID = aIndex;
            mOrientationOnSlave = aOrientation ;
        }

//------------------------------------------------------------------------------

        void
        Facet::disconnect_slave()
        {
            mSlave = nullptr ;
            mSlaveFaceID = BELFEM_SUINT_MAX ;
            mOrientationOnSlave = BELFEM_SUINT_MAX ;
        }


//------------------------------------------------------------------------------

        void
        Facet::flip()
        {
            if ( ! this->has_master() && ! this->has_slave() ) return;

            std::swap( mMaster, mSlave );
            std::swap( mMasterFaceID, mSlaveFaceID );
            this->set_master( mMaster, mMasterFaceID, true );
            mOrientationOnSlave = BELFEM_SUINT_MAX ;
            this->compute_orientation();
        }

//------------------------------------------------------------------------------

        void
        Facet::compute_orientation()
        {
            // check of orientation has already been computed or provided
            if ( mOrientationOnSlave < BELFEM_SUINT_MAX ) return;

            if( this->has_master() && this->has_slave() )
            {
                Cell< Node * > tNodes ;
                mSlave->get_corner_nodes_of_facet( mSlaveFaceID, tNodes );
                id_t tID = tNodes.first()->original()->id();
                mMaster->get_corner_nodes_of_facet( mMasterFaceID, tNodes );

                uint tNumNodes =  tNodes.size();

                for ( uint t = 0; t < tNumNodes; ++t )
                {
                    if ( tID == tNodes( t )->original()->id() )
                    {
                        mOrientationOnSlave = t + 1 ;
                        return;
                    }
                }

                BELFEM_ERROR( false, "Could not determine slave orientation of facet %lu (master: %lu, slave: %lu)",
                    ( long unsigned int ) this->id(), ( long unsigned int ) mMaster->id(), ( long unsigned int ) mSlave->id() );
            }
        }

//------------------------------------------------------------------------------

        size_t
        Facet::memory() const
        {
            // the wrapped element is a separate object this facet owns
            return sizeof( Facet ) + this->array_memory() + mElement->memory() ;
        }

//------------------------------------------------------------------------------
    }
}
