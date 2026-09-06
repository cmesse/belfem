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

#include "cl_Face.hpp"
#include "cl_Element.hpp"
#include "cl_Cell.hpp"
#include "fn_to_master_orientation.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Face::Face( Element    * aParent ) :
                Vertex(),
                mMaster( aParent ),
                mIndexOnMaster( 0 ),
                mSlave( nullptr ),
                mIndexOnSlave( gNoIndex ),
                mOrientationOnSlave( BELFEM_UINT_MAX )
        {
            // set the owner to be identical to the one of the primary element
            this->set_owner( aParent->owner() );

            // copy nodes from parent
            uint tNumNodes = aParent->number_of_nodes() ;
            this->allocate_node_container( tNumNodes );
            for( uint k=0; k<tNumNodes; ++k )
            {
                this->insert_node( aParent->node( k ), k );
            }
        }

//------------------------------------------------------------------------------

        Face::Face( Element    * aMaster,
                    const uint   aIndexOnMaster,
                    Element    * aSlave,
                    const uint   aIndexOnSlave,
                    const uint   aOrientationOnSlave ) :
                Vertex(),
                mMaster( aMaster ),
                mIndexOnMaster( aIndexOnMaster ),
                mSlave( aSlave ),
                mIndexOnSlave( aIndexOnSlave ),
                mOrientationOnSlave( ( aOrientationOnSlave == gNoIndex && aSlave != nullptr )?
                        this->compute_orientation( aMaster, aIndexOnMaster, aSlave, aIndexOnSlave ) :
                        aOrientationOnSlave )
        {
            // set the owner to be identical to the one of the primary element

            Cell< Node * > tNodes ;

            if ( aMaster != nullptr )
            {
                aMaster->get_nodes_of_facet( aIndexOnMaster, tNodes );
                this->set_owner( aMaster->owner() );
            }
            else if ( aSlave != nullptr )
            {
                aSlave->get_nodes_of_facet( aIndexOnSlave, tNodes );
                this->set_owner( aSlave->owner() );

                Cell< Node * > tNodesInSlaveOrientation ;

                // first, we grab the nodes from the slave into a temporary container
                aSlave->get_nodes_of_facet( aIndexOnSlave, tNodesInSlaveOrientation );

                // now we align the nodes in the orientation of the master and write them into tNodes
                to_master_orientation( this, tNodesInSlaveOrientation, tNodes );

                // note: this could be wrong but we'll overwrite it later anyways
                this->set_owner( aSlave->owner() );

            }
            else
            {
                return;
            }


            // allocate the nodes
            this->allocate_node_container( tNodes.size() );

            uint tCount = 0 ;

            for( Node * tNode : tNodes )
            {
                this->insert_node( tNode, tCount++ );
            }

        }

//------------------------------------------------------------------------------

        Face::~Face()
        {
            this->delete_containers();
        }

//-----------------------------------------------------------------------------

        uint
        Face::compute_orientation(
                    Element * aMaster,
            const uint aIndexOnMaster,
                    Element * aSlave,
            const uint aIndexOnSlave )
        {
            if( aMaster != nullptr && aSlave != nullptr )
            {
                Cell< Node * > tMasterNodes;
                aMaster->get_corner_nodes_of_facet( aIndexOnMaster, tMasterNodes );

                Cell< Node * > tSlaveNodes;
                aSlave->get_corner_nodes_of_facet( aIndexOnSlave, tSlaveNodes );

                BELFEM_ASSERT( tMasterNodes.size() == tSlaveNodes.size(), "Faces do not match. Master %lu (%u), Slave %lu (%u)",
                    ( long unsigned int ) aMaster->id(), ( unsigned int ) aIndexOnMaster,
                    ( long unsigned int ) aSlave->id(), ( unsigned int ) aIndexOnSlave );

                // get first id of node on slave side
                id_t tID = tSlaveNodes( 0 )->original()->id();

                // loop over all master nodes
                for ( uint aOrientation = 0; aOrientation < tMasterNodes.size(); ++aOrientation )
                {
                    if ( tID == tMasterNodes( aOrientation )->original()->id() )
                    {
                        return aOrientation + 1 ;
                    }
                }

                // catch error
                BELFEM_ERROR( false, "Could not determine orientation of facet %u of element %lu",
                             ( unsigned int ) aIndexOnSlave,
                             ( long unsigned int ) aSlave->id());
            }

            return BELFEM_UINT_MAX ;
        }

//-----------------------------------------------------------------------------

        bool
        Face::edge_direction( const uint aEdgeIndex ) const
        {
            BELFEM_ASSERT( aEdgeIndex < this->number_of_edges(), "Invalid edge index %u at face %lu",
                           ( unsigned int ) aEdgeIndex,
                           ( long unsigned int ) this->id() );

            id_t tA = this->node( aEdgeIndex )->original()->id();
            id_t tB = this->node( aEdgeIndex + 1 < this->number_of_edges() ? aEdgeIndex + 1 : 0 )->original()->id() ;

            const Edge * tEdge = this->edge( aEdgeIndex );

            if( tA == tEdge->node( 0 )->original()->id() && tB == tEdge->node( 1 )->original()->id() )
            {
                return true ;
            }
            else if( tA == tEdge->node( 1 )->original()->id() && tB == tEdge->node( 0 )->original()->id() )
            {
                return false ;
            }
            else
            {
                BELFEM_ERROR( false, "Invalid edge %u at face %lu",
                               ( unsigned int ) aEdgeIndex,
                               ( long unsigned int ) this->id() );
                return false ;
            }
        }

//-----------------------------------------------------------------------------

        uint
        Face::number_of_corner_nodes()
        {
            return mesh::number_of_corner_nodes( mesh::element_type_from_numnodes( 2, this->number_of_nodes() ) ) ;
        }

//-----------------------------------------------------------------------------

        void
        Face::flag_corner_nodes()
        {
            uint n = this->number_of_corner_nodes();

            for( uint k=0; k<n; ++k )
            {
                mNodes[ k ]->flag();
            }
        }

//-----------------------------------------------------------------------------

        void
        Face::set_master( Element *aElement, const uint aIndex )
        {
            mMaster = aElement;
            mIndexOnMaster = aIndex;
        }

//-----------------------------------------------------------------------------

        void Face::set_slave( Element * aElement, const uint aIndex, const uint aOrientation )
        {
            mSlave = aElement;
            mIndexOnSlave = aIndex ;
            mOrientationOnSlave = aOrientation ;
        }

//-----------------------------------------------------------------------------
    }
}
