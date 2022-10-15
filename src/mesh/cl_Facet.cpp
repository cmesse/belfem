//
// Created by Christian Messe on 2019-07-28.
//

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
        Facet::compute_orientation()
        {
            // reset orientation
            mOrientationOnSlave = BELFEM_UINT_MAX ;

            if( this->has_master() && this->has_slave() )
            {
                Cell< Node * > tMasterNodes;
                mMaster->get_nodes_of_facet( this->master_index(), tMasterNodes );

                Cell< Node * > tSlaveNodes;
                mSlave->get_nodes_of_facet( this->slave_index(), tSlaveNodes );

                BELFEM_ASSERT( tMasterNodes.size() == tSlaveNodes.size(), "Facets do not match" );


                // get first id of node on slave side
                id_t tID = tSlaveNodes( 0 )->id();

                // loop over all master nodes
                for ( uint t = 0; t < tMasterNodes.size(); ++t )
                {
                    if ( tID == tMasterNodes( t )->id() )
                    {
                        mOrientationOnSlave = t ;
                        break ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    }
}