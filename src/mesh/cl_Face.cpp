//
// Created by christian on 7/12/21.
//
#include "cl_Face.hpp"
#include "cl_Element.hpp"
#include "cl_Cell.hpp"
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
                    const uint   aIndexOnSlave ) :
                Vertex(),
                mMaster( aMaster ),
                mIndexOnMaster( aIndexOnMaster ),
                mSlave( aSlave ),
                mIndexOnSlave( aIndexOnSlave ),
                mOrientationOnSlave(
                        this->compute_orientation( aMaster, aIndexOnMaster, aSlave, aIndexOnSlave ) )
        {
            // set the owner to be identical to the one of the primary element
            this->set_owner( aMaster->owner() );

            // the next few lines copy the node pointers form the master element
            // this can be improved
            Cell< mesh::Node * > tNodes ;

            aMaster->get_nodes_of_facet( aIndexOnMaster, tNodes );

            // allocate the nodes
            this->allocate_node_container( tNodes.size() );

            uint tCount = 0 ;

            for( mesh::Node * tNode : tNodes )
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
            if( aSlave != nullptr )
            {
                Cell< Node * > tMasterNodes;
                aMaster->get_nodes_of_facet( aIndexOnMaster, tMasterNodes );

                Cell< Node * > tSlaveNodes;
                aSlave->get_nodes_of_facet( aIndexOnSlave, tSlaveNodes );

                BELFEM_ASSERT( tMasterNodes.size() == tSlaveNodes.size(), "Facets do not match" );

                // get first id of node on slave side
                id_t tID = tSlaveNodes( 0 )->id();

                // loop over all master nodes
                for ( uint aOrientation = 0; aOrientation < tMasterNodes.size(); ++aOrientation )
                {
                    if ( tID == tMasterNodes( aOrientation )->id() )
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
    }
}