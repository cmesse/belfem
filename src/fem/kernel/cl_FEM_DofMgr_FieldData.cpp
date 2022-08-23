//
// Created by christian on 7/14/21.
//

#include "cl_FEM_DofMgr_FieldData.hpp"

#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"
#include "cl_FEM_DofMgr_FieldData.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "cl_Mesh.hpp"
#include "fn_max.hpp"
#include "fn_entity_type.hpp"
namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//------------------------------------------------------------------------------

            FieldData::FieldData(
                    DofManager * aParent ) :
                    mParent( aParent ),
                    mKernel( aParent->parent() ),
                    mMesh( aParent->parent()->mesh() ),
                    mMyRank( aParent->rank() ),
                    mCommSize( aParent->parent()->number_of_procs() )
            {

            }

//------------------------------------------------------------------------------

            FieldData::~FieldData()
            {
                this->reset();
            }

//------------------------------------------------------------------------------

            void
            FieldData::reset()
            {
                mMyNumberOfOwnedNodes = 0 ;
                mMyNumberOfOwnedGhostElements = 0 ;
                mNodeOwnerList.clear();
                mAllCornerNodeIndices.clear() ;
                mAllNonCornerNodeIndices.clear() ;
            }

//------------------------------------------------------------------------------

            void
            FieldData::update_field_indices( Cell< Dof * > & aDOFs )
            {
                // get the types
                const Vector< index_t > & tTypes = mParent->iwg()->default_dof_types();

                BELFEM_ERROR( tTypes.length() > 0 , "Fields have not been set" );

                const Cell< string > & tLabels = mParent->iwg()->dof_fields() ;

                Vector< index_t > tFieldIndices(  max( tTypes ) + 1, gNoIndex );

                // loop over all fields
                for( uint k=0; k<tLabels.size(); ++k )
                {
                    // get index of field
                    tFieldIndices( tTypes( k ) )
                        = mMesh->field( tLabels( k ) )->index() ;
                }

                for( Dof * tDof : aDOFs )
                {
                    tDof->set_field_index( tFieldIndices( tDof->type_id() ) );
                }
            }

//------------------------------------------------------------------------------

            void
            FieldData::collect_node_owners()
            {
                if( mCommSize > 1 )
                {
                    if ( mMyRank != mKernel->master() )
                    {
                        // count owned nodes
                        mMyNumberOfOwnedNodes = 0;

                        // loop over all nodes on mesh
                        for ( mesh::Node * tNode : mMesh->nodes() )
                        {
                            if ( tNode->owner() == mMyRank )
                            {
                                ++mMyNumberOfOwnedNodes;
                            }
                        }

                        // allocate memory
                        Vector< id_t > tIDs( mMyNumberOfOwnedNodes );
                        index_t tCount = 0;

                        // collect IDs of owned nodes
                        for ( mesh::Node * tNode : mMesh->nodes() )
                        {
                            if ( tNode->owner() == mMyRank )
                            {
                                tIDs( tCount++ ) = tNode->id();
                            }
                        }

                        send( mKernel->master(), tIDs );
                    }
                    else
                    {
                        uint tNumProcs = mKernel->comm_table().length();

                        Cell< Vector< id_t > > tAllIDs( tNumProcs, Vector< id_t >() );

                        receive( mKernel->comm_table(), tAllIDs );

                        // allocate index vector
                        mNodeOwnerList.set_size( tNumProcs, Vector< index_t >() );

                        // loop over all ids
                        for ( uint p = 0; p < tNumProcs; ++p )
                        {
                            Vector< index_t > & tIndices = mNodeOwnerList( p );
                            Vector< id_t > & tIDs = tAllIDs( p );

                            index_t tNumNodes = tIDs.length();
                            tIndices.set_size( tNumNodes );

                            for ( index_t k = 0; k < tNumNodes; ++k )
                            {
                                tIndices( k ) = mMesh->node( tIDs( k ) )->index();
                            }
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            FieldData::collect_ghost_element_owners()
            {
                if( mCommSize > 1 )
                {
                    if ( mMyRank != mKernel->master() )
                    {
                        // count owned nodes
                        mMyNumberOfOwnedGhostElements = 0;

                        // loop over all ghost blocks
                        for( id_t b : mMesh->ghost_block_ids() )
                        {
                            if( mMesh->block_exists( b ) )
                            {
                                // grab elements
                                Cell< mesh::Element * > & tElements = mMesh->block( b )->elements();
                                for ( mesh::Element * tElement: tElements )
                                {
                                    if ( tElement->owner() == mMyRank )
                                    {
                                        ++mMyNumberOfOwnedGhostElements;
                                    }
                                }
                            }
                        }


                        // allocate memory
                        Vector< id_t > tIDs( mMyNumberOfOwnedGhostElements );
                        index_t tCount = 0;

                        // collect IDs of owned nodes
                        for( id_t b : mMesh->ghost_block_ids() )
                        {
                            if( mMesh->block_exists( b ) )
                            {
                                // grab elements
                                Cell< mesh::Element * > & tElements = mMesh->block( b )->elements();
                                for ( mesh::Element * tElement: tElements )
                                {
                                    if ( tElement->owner() == mMyRank )
                                    {
                                       tIDs( tCount++ ) = tElement->id() ;
                                    }
                                }
                            }
                        }

                        send( mKernel->master(), tIDs );
                    }
                    else
                    {
                        uint tNumProcs = mKernel->comm_table().length();

                        Cell< Vector< id_t > > tAllIDs( tNumProcs, Vector< id_t >() );

                        receive( mKernel->comm_table(), tAllIDs );

                        // allocate index vector
                        mGhostElementOwnerList.set_size( tNumProcs, Vector< index_t >() );

                        // loop over all ids
                        for ( uint p = 0; p < tNumProcs; ++p )
                        {
                            Vector< index_t > & tIndices = mGhostElementOwnerList( p );
                            Vector< id_t > & tIDs = tAllIDs( p );

                            index_t tNumNodes = tIDs.length();
                            tIndices.set_size( tNumNodes );

                            for ( index_t k = 0; k < tNumNodes; ++k )
                            {
                                tIndices( k ) = mMesh->element( tIDs( k ) )->index();
                            }
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            FieldData::collect( const string & aLabel )
            {
                EntityType tType = mMesh->field( aLabel )->entity_type() ;

                BELFEM_ASSERT( tType == EntityType::NODE ||
                                       tType == EntityType::ELEMENT,
                    "Field %s is off type %s but must be node to collect",
                    aLabel.c_str(), to_string( mMesh->field( aLabel )->entity_type() ).c_str() );

                if( mCommSize > 1 )
                {
                    // grab field
                    Vector< real > & tField = mMesh->field_data( aLabel );

                    if ( mMyRank != mKernel->master() )
                    {
                        if( tType == EntityType::NODE )
                        {
                            Vector< real > tSubField( mMyNumberOfOwnedNodes );

                            // initialize counter
                            index_t tCount = 0;

                            // collect owned datals
                            for ( mesh::Node * tNode: mMesh->nodes())
                            {
                                if ( tNode->owner() == mMyRank )
                                {
                                    tSubField( tCount++ ) = tField( tNode->index());
                                }
                            }

                            // send data to master
                            send( mKernel->master(), tSubField );
                        }
                        else
                        {
                            Vector< real > tSubField( mMyNumberOfOwnedGhostElements );

                            // initialize counter
                            index_t tCount = 0;

                            for( id_t b : mMesh->ghost_block_ids() )
                            {
                                if( mMesh->block_exists( b ) )
                                {
                                    // grab elements
                                    Cell< mesh::Element * > & tElements = mMesh->block( b )->elements();
                                    for ( mesh::Element * tElement: tElements )
                                    {
                                        if ( tElement->owner() == mMyRank )
                                        {
                                            tSubField( tCount++ ) = tField( tElement->index() );
                                        }
                                    }
                                }
                            }

                            // send data to master
                            send( mKernel->master(), tSubField );

                        }
                    }
                    else
                    {
                        // collect subfields
                        Cell< Vector< real > > tSubFields;
                        receive( mKernel->comm_table(), tSubFields );

                        // assemble
                        uint tNumProcs = mKernel->comm_table().length();

                        if( tType == EntityType::NODE )
                        {
                            for ( uint p = 0; p < tNumProcs; ++p )
                            {
                                // grab index vector
                                const Vector< index_t > & tIndices = mNodeOwnerList( p );
                                const Vector< real > & tSubField = tSubFields( p );

                                index_t tNumNodes = tIndices.length();

                                for ( index_t k = 0; k < tNumNodes; ++k )
                                {
                                    tField( tIndices( k )) = tSubField( k );
                                }
                            }
                        }
                        else
                        {
                            for ( uint p = 0; p < tNumProcs; ++p )
                            {
                                // grab index vector
                                const Vector< index_t > & tIndices = mGhostElementOwnerList( p );
                                const Vector< real > & tSubField = tSubFields( p );

                                index_t tNumNodes = tIndices.length();

                                for ( index_t k = 0; k < tNumNodes; ++k )
                                {
                                    tField( tIndices( k )) = tSubField( k );
                                }
                            }
                        }
                    }
                }
            }

//-----------------------------------------------------------------------------

            void
            FieldData::collect( const Cell< string > & aLabels )
            {
                if( mCommSize > 1 )
                {
                    if ( mMyRank != mKernel->master() )
                    {
                        uint tNumFields = aLabels.size();

                        Matrix< real > tSubFields( mMyNumberOfOwnedNodes, tNumFields );

                        for ( uint f = 0; f < tNumFields; ++f )
                        {
                            BELFEM_ASSERT( mMesh->field( aLabels( f ) )->entity_type() ==  EntityType::NODE,
                                          "Field %s is off type %s but must be node to collect",
                                          aLabels( f ).c_str(), to_string( mMesh->field( aLabels( f ) )->entity_type() ).c_str() );

                            Vector< real > & tField = mMesh->field_data( aLabels( f ) );

                            // initialize counter
                            index_t tCount = 0;

                            // collect owned data
                            for ( mesh::Node * tNode : mMesh->nodes() )
                            {

                                if ( tNode->owner() == mMyRank )
                                {
                                    tSubFields( tCount++, f ) = tField( tNode->index() );
                                }
                            }
                        }
                        // send data to master
                        send( mKernel->master(), tSubFields );
                    }
                    else
                    {
                        // collect subfields
                        Cell< Matrix< real > > tAllSubFields;
                        receive( mKernel->comm_table(), tAllSubFields );

                        uint tNumFields = aLabels.size() ;

                        // assemble
                        uint tNumProcs = mKernel->comm_table().length();

                        for ( uint p = 1; p < tNumProcs; ++p )
                        {
                            // grab index vector
                            const Vector< index_t > & tIndices = mNodeOwnerList( p );
                            const Matrix< real >    & tFields  = tAllSubFields( p );

                            index_t tNumNodes = tIndices.length();

                            for ( uint f = 0; f < tNumFields; ++f )
                            {
                                Vector< real > & tField = mMesh->field_data( aLabels( f ) );

                                for ( index_t k = 0; k < tNumNodes; ++k )
                                {
                                    tField( tIndices( k ) ) = tFields( k, f );
                                }
                            }
                        }
                    }
                }
            }

//-----------------------------------------------------------------------------

            void
            FieldData::initialize_linear_projection_lists()
            {
                // get the block list
                const Vector< id_t > & tSelectedBlocks
                    = mParent->iwg()->selected_blocks();

                Cell< Vector< id_t > > tAllNodeIDs;
                Vector< id_t > tMyNodeIDs;

                // start a timer
                Timer tTimer;

                const bool tIsMaster = mMyRank == mKernel->master() ;

                // get node container of mesh
                Cell< mesh::Node * > & tNodes = mMesh->nodes();

                const uint tNumProcs = mKernel->number_of_procs() ;

                if ( tIsMaster )
                {

                    // make sure that procs are in consecutive order
                    proc_t q = 0;
                    for ( proc_t p : mKernel->comm_table() )
                    {
                        BELFEM_ERROR( p == q++, "proc table must be consecutive" );
                    }

                    // allocate id and index containers
                    tAllNodeIDs.set_size( tNumProcs, {} );
                    mAllCornerNodeIndices.set_size( tNumProcs, {} );

                    // - - - - - - - - - - - - - - - - - - - -
                    // STEP 1 : identify corner nodes per proc
                    // - - - - - - - - - - - - - - - - - - - -

                    for ( proc_t p : mKernel->comm_table() )
                    {
                        mMesh->unflag_all_nodes();
                        Vector< id_t >    & tNodeIDs = tAllNodeIDs( p );
                        Vector< index_t > & tNodeIndices = mAllCornerNodeIndices( p );

                        // find corner nodes
                        for ( index_t tID : tSelectedBlocks )
                        {
                            // grab element list of this block
                            Cell< mesh::Element * > & tElements = mMesh->block( tID )->elements();

                            for ( mesh::Element * tElement : tElements )
                            {
                                if( tElement->owner() == p )
                                {
                                    tElement->flag_corner_nodes();
                                }
                            }
                        }

                        // local counter
                        index_t tCount = 0 ;

                        // count nodes
                        for ( mesh::Node * tNode : tNodes )
                        {
                            // check if node is selected
                            if ( tNode->is_flagged() )
                            {
                                // increment counter
                               ++tCount ;
                            }
                        }

                        if( tCount > 0 )
                        {
                            tNodeIDs.set_size( tCount );
                            tNodeIndices.set_size( tCount );


                            // reset counter
                            tCount = 0;

                            // collect node IDs and indices
                            for ( mesh::Node * tNode: tNodes )
                            {
                                if ( tNode->is_flagged() )
                                {
                                    // add node to list
                                    tNodeIDs( tCount ) = tNode->id() ;
                                    tNodeIndices( tCount++ ) = tNode->index() ;
                                }
                            }
                        }
                    } // end loop over all procs

                    // wait
                    comm_barrier() ;

                    // send IDs to other procs
                    send( mKernel->comm_table(), tAllNodeIDs );

                }
                else
                {
                    // wait
                    comm_barrier() ;

                    // get corner node IDs from master
                    receive( mKernel->master(), tMyNodeIDs );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 2 : identify indices of corner nodes
                //          ( this is the same for all procs )
                // - - - - - - - - - - - - - - - - - - - - - - - -

                Vector< id_t > & tCornerNodeIDs = tIsMaster ? tAllNodeIDs( 0 ) : tMyNodeIDs ;

                // allocate memory
                mMyCornerNodeIndices.set_size( tCornerNodeIDs.length() );

                // initialize counter
                index_t k = 0 ;

                for( id_t tID : tCornerNodeIDs )
                {
                    mMyCornerNodeIndices( k++ ) = mMesh->node( tID )->index() ;
                }

                BELFEM_ASSERT( k == mMyCornerNodeIndices.length(),
                              "Invalid length of mMyCornerNodeIndices. ( is %lu but expect %lu )",
                              ( long unsigned int ) k,
                              ( long unsigned int ) mMyCornerNodeIndices.length() );

                // - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 3 : identify IDs of non corner node indices
                // - - - - - - - - - - - - - - - - - - - - - - - - -

                if( tIsMaster )
                {
                    // make sure that procs are in consecutive order
                    proc_t q = 0;
                    for ( proc_t p : mKernel->comm_table() )
                    {
                        BELFEM_ERROR( p == q++, "proc table must be consecutive" );
                    }

                    // allocate index containers
                    mAllNonCornerNodeIndices.set_size( tNumProcs, {} );

                    // - - - - - - - - - - - - - - - - - - - -
                    // STEP 1 : identify corner nodes per proc
                    // - - - - - - - - - - - - - - - - - - - -

                    for ( proc_t p : mKernel->comm_table() )
                    {
                        mMesh->unflag_all_nodes();

                        Vector< id_t >    & tNodeIDs = tAllNodeIDs( p );
                        Vector< index_t > & tNodeIndices = mAllNonCornerNodeIndices( p ) ;

                        // flag all nodes
                        for ( index_t tID : tSelectedBlocks )
                        {
                            // grab element list of this block
                            Cell< mesh::Element * > & tElements = mMesh->block( tID )->elements();

                            for ( mesh::Element * tElement : tElements )
                            {
                                if( tElement->owner() == p )
                                {
                                    tElement->flag_nodes() ;
                                    tElement->unflag_corner_nodes() ;
                                }
                            }
                        }

                        // local counter
                        index_t tCount = 0 ;

                        // count nodes
                        for ( mesh::Node * tNode : tNodes )
                        {
                            // check if node is selected
                            if ( tNode->is_flagged() )
                            {
                                // increment counter
                                ++tCount ;
                            }
                        }

                        if( tCount > 0 )
                        {
                            tNodeIDs.set_size( tCount );
                            tNodeIndices.set_size( tCount );


                            // reset counter
                            tCount = 0;

                            // collect node IDs and indices
                            for ( mesh::Node * tNode: tNodes )
                            {
                                if ( tNode->is_flagged() )
                                {
                                    // add node to list
                                    tNodeIDs( tCount ) = tNode->id() ;
                                    tNodeIndices( tCount++ ) = tNode->index() ;
                                }
                            }
                        }
                    } // end loop over all procs

                    // wait
                    comm_barrier() ;

                    // send IDs to other procs
                    send( mKernel->comm_table(), tAllNodeIDs );
                }
                else
                {
                    // wait
                    comm_barrier() ;
                    receive( mKernel->master(), tMyNodeIDs );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 4 : correct order of node IDs and indices
                // - - - - - - - - - - - - - - - - - - - - - - - - -

                mMesh->unflag_all_nodes() ;

                Vector< id_t > & tMyNonCornerNodeIDs = tIsMaster ? tAllNodeIDs( 0 ) : tMyNodeIDs ;

                // flag nodes
                for( id_t tID : tMyNonCornerNodeIDs )
                {
                    mMesh->node( tID )->flag() ;
                }

                // reset counter
                k = 0 ;

                // allocate memory
                mMyNonCornerNodeIndices.set_size( tMyNonCornerNodeIDs.length() );

                // allocate container
                mMyNonCornerNodes.set_size( tMyNonCornerNodeIDs.length(), nullptr );

                tMyNodeIDs.set_size( tMyNonCornerNodeIDs.length() );

                for ( index_t tID : tSelectedBlocks )
                {
                    if( mParent->block_exists( tID ) )
                    {
                        Block * tBlock = mParent->block( tID );

                        // get real element type of block ( tBlock->element_type() is the linear one )
                        ElementType tType = tBlock->block()->element_type();

                        // get number of nodes per element
                        uint tNumCornerNodes = mesh::number_of_corner_nodes( tType );
                        uint tNumNodes = mesh::number_of_nodes( tType );

                        // grab elements on block
                        Cell< mesh::Element * > & tElements = tBlock->block()->elements();

                        // loop over all elements on this block
                        for ( mesh::Element * tElement: tElements )
                        {
                            // loop over all non corner nodes of this element
                            for ( uint i = tNumCornerNodes; i < tNumNodes; ++i )
                            {
                                // grab node
                                mesh::Node * tNode = tElement->node( i );

                                // check if node has been processed
                                if ( tNode->is_flagged() )
                                {
                                    mMyNonCornerNodeIndices( k ) = tNode->index();

                                    // store ID for sending
                                    tMyNodeIDs( k ) = tNode->id();

                                    // add node to container
                                    mMyNonCornerNodes( k++ ) = tNode;

                                    // unflag this node, each node is processed only once
                                    tNode->unflag();
                                }
                            }
                        }
                    }
                }

                BELFEM_ASSERT( k == mMyNonCornerNodeIndices.length(),
                              "Invalid length of mMyNonCornerNodeIndices. ( is %lu but expect %lu )",
                              ( long unsigned int ) k,
                              ( long unsigned int ) mMyNonCornerNodeIndices.length() );

                // wait
                comm_barrier() ;

                // - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 5 : send indices in correct order back to master
                // - - - - - - - - - - - - - - - - - - - - - - - - -

                if( tIsMaster )
                {
                    // get IDs from other procs
                    receive( mKernel->comm_table(), tAllNodeIDs );

                    // allocate array
                    mAllNonCornerNodeIndices.set_size( tNumProcs, {} );

                    // create indices
                    for( uint p=1; p<tNumProcs; ++p )
                    {
                        Vector< index_t > & tNodeIDs = tAllNodeIDs( p );

                        // get indices
                        Vector< index_t > & tIndices = mAllNonCornerNodeIndices( p );

                        tIndices.set_size( tNodeIDs.length() );

                        // reset counter
                        k = 0 ;

                        for( id_t tID : tNodeIDs )
                        {
                            tIndices( k++ ) = mMesh->node( tID )->index() ;
                        }
                    }
                }
                else
                {
                    // send IDs to master
                    send( mKernel->master(), tMyNodeIDs );
                }

                // wait
                comm_barrier() ;

                if ( mMyRank == mKernel->master() )
                {
                    message( 4, "    ... time for collecting procjection nodes   : %u ms\n",
                             ( unsigned int ) tTimer.stop() );
                }
            }
//------------------------------------------------------------------------------

            void
            FieldData::communicate_corner_node_data( const Cell< string > & aFieldLabels )
            {
                index_t tNumFields = aFieldLabels.size() ;
                uint tNumProcs = mKernel->number_of_procs() ;

                if ( mMyRank == mKernel->master() )
                {
                    Cell< Vector< real > > tAllData( tNumProcs, {} );

                    // loop over all procs
                    for( uint p=1; p<tNumProcs; ++p )
                    {
                        // grab node indices
                        const Vector< index_t > & tIndices = mAllCornerNodeIndices( p );

                        // grab data
                        Vector< real > & tData = tAllData( p );

                        // number of nodes on this field
                        index_t tNumNodes = tIndices.length() ;

                        // allocate size
                        tData.set_size( tNumFields * tIndices.length() );

                        // initialize counter
                        index_t tCount = 0 ;

                        // loop over all fields
                        for( index_t f=0; f<tNumFields; ++f )
                        {
                            // grab field data
                            Vector< real > & tField = mMesh->field_data( aFieldLabels( f ) );

                            // loop over all nodes
                            for( index_t k=0; k<tNumNodes; ++k )
                            {
                                tData( tCount++ ) = tField( tIndices( k ) );
                            }
                        }
                    }

                    // wait
                    comm_barrier();

                    send( mKernel->comm_table(), tAllData );
                }
                else
                {
                    // Data container
                    Vector< real > tData ;

                    // wait
                    comm_barrier() ;
                    receive( mKernel->master(), tData );

                    // initialize counter
                    index_t tCount = 0 ;

                    index_t tNumNodes = mMyCornerNodeIndices.length() ;

                    // loop over all fields
                    for( index_t f=0; f<tNumFields; ++f )
                    {
                        // grab field data
                        Vector< real > & tField = mMesh->field_data( aFieldLabels( f ) );

                        // loop over all nodes
                        for( index_t k=0; k<tNumNodes; ++k )
                        {
                            tField( mMyCornerNodeIndices( k ) ) = tData( tCount++ );
                        }
                    }

                    // sanity check
                    BELFEM_ASSERT( tCount == tData.length(), "Invalid vector size. Is %lu but expect %lu",
                                  ( long unsigned int ) tCount,
                                  ( long unsigned int ) tData.length() );

                }

                // wait for all procs to be done
                comm_barrier() ;
            }

//------------------------------------------------------------------------------

            void
            FieldData::communicate_noncorner_node_data(
                    const Cell< string > & aFieldLabels,
                    Matrix< real > & aData )
            {
                // get number of fields
                index_t tNumFields = aFieldLabels.size() ;
                uint tNumProcs = mKernel->number_of_procs();

                if( mMyRank == mKernel->master() )
                {

                    // container with other data
                    Cell< Matrix< real > > tAllData( mKernel->number_of_procs(), {} );

                    receive( mKernel->comm_table(), tAllData );

                    // get number of nodes
                    index_t tNumNodes = mMyNonCornerNodeIndices.length() ;

                    // write my own data
                    for( index_t f=0; f<tNumFields; ++f )
                    {
                        // grab field
                        Vector< real > & tData = mMesh->field_data( aFieldLabels( f ) );

                        for( index_t k=0; k<tNumNodes; ++k )
                        {
                            tData( mMyNonCornerNodeIndices( k ) ) = aData( k, f );
                        }
                    }

                    // write data of other procs
                    for( uint p=1; p<tNumProcs; ++p )
                    {
                        // grab data container
                        Matrix< real > & tData = tAllData( p );

                        // grab indices
                        Vector< index_t > & tIndices = mAllNonCornerNodeIndices( p );

                        // get number of nodes
                        index_t tN = tIndices.length() ;

                        // write my own data
                        for( index_t f=0; f<tNumFields; ++f )
                        {
                            // grab field
                            Vector< real > & tField = mMesh->field_data( aFieldLabels( f ) );

                            // write data onto mesh
                            for( index_t k=0; k<tN; ++k )
                            {
                                tField( tIndices( k ) ) = tData( k, f );
                            }
                        }
                    }
                }
                else
                {
                    // send data to master
                    send( mKernel->master(), aData );
                }

                // wait
                comm_barrier() ;
            }

//------------------------------------------------------------------------------

            void
            FieldData::project_linear_field_to_higher_mesh(
                    const Cell< string > & aFieldLabels )
            {
                Timer tTimer ;

                // collect node fields
                Cell< string > tFieldLabels ;
                for( uint f = 0; f < aFieldLabels.size(); ++f )
                {
                    if( mMesh->field( aFieldLabels( f ) )->entity_type() == EntityType::NODE )
                    {
                        tFieldLabels.push( aFieldLabels( f ) );
                    }
                }

                // communicate corner node data
                this->communicate_corner_node_data( tFieldLabels );

                // unflag all nodes on mesh
                mMesh->unflag_all_nodes() ;

                // flag nodes of interest
                for( mesh::Node * tNode : mMyNonCornerNodes )
                {
                    tNode->flag() ;
                }

                // the shape function factory
                InterpolationFunctionFactory tFactory ;

                // number of fields to interpolate
                const index_t tNumFields = tFieldLabels.size() ;

                const index_t tNumNonCornerNodes = mMyNonCornerNodes.size() ;

                Matrix< real > tWork( tNumNonCornerNodes, tNumFields );

                // Node Counter
                index_t tCount = 0 ;

                // loop over all blocks
                for( id_t tBlockID : mParent->iwg()->selected_blocks() )
                {
                    if( mParent->block_exists( tBlockID ) )
                    {
                        // grab block on mesh
                        mesh::Block * tBlock = mParent->block( tBlockID )->block() ;

                        // get element type
                        const ElementType tType = tBlock->element_type();

                        // create the linear shape function
                        InterpolationFunction * tLinShape = tFactory.create_lagrange_function(
                                mesh::linear_element_type( tType ));

                        // create the higher order shape function
                        InterpolationFunction * tHighShape = tFactory.create_lagrange_function(
                                tType );

                        // grab parameter coordinates
                        Matrix< real > tXi;
                        tHighShape->param_coords( tXi );

                        // get number of nodes of element
                        uint tNumNodes = mesh::number_of_nodes( tType );

                        // number of corner nodes
                        uint tNumCornerNodes = mesh::number_of_corner_nodes( tType );

                        // Evaluated function
                        Cell< Matrix< real > > tAllN( tXi.n_cols(), {} );

                        for ( uint k = 0; k < tNumNodes; ++k )
                        {
                            // get vector
                            Matrix< real > & tN = tAllN( k );

                            // allocate size
                            tN.set_size( 1, tNumCornerNodes );

                            // evaluate function
                            tLinShape->N( tXi.col( k ), tN );
                        }

                        // get elements from mesh
                        Cell< mesh::Element * > & tElements = tBlock->elements();

                        Matrix< real > tLinData( tNumCornerNodes, tNumFields );

                        Matrix< real > tNodeData( 1, tNumFields );

                        // loop over all elements
                        for ( mesh::Element * tElement : tElements )
                        {
                            // collect data from fields
                            for ( uint j = 0; j < tNumFields; ++j )
                            {
                                // grab field
                                Vector< real > & tField = mMesh->field_data( tFieldLabels( j ));

                                // grab nodes
                                for ( uint i = 0; i < tNumCornerNodes; ++i )
                                {
                                    tLinData( i, j ) = tField( tElement->node( i )->index() );
                                }
                            }

                            // loop over all non-corner nodes
                            for ( uint k = tNumCornerNodes; k < tNumNodes; ++k )
                            {
                                // grab node
                                mesh::Node * tNode = tElement->node( k );

                                // check if node has been computed
                                if ( tNode->is_flagged() )
                                {
                                    // interpolate data
                                    tNodeData = tAllN( k ) * tLinData;

                                    // write data into work matrix
                                    tWork.set_row( tCount++, tNodeData.row( 0 ) );

                                    // unflag node
                                    tNode->unflag();
                                }
                            }
                        }
                        // delete shape functions
                        delete tLinShape;
                        delete tHighShape;
                    }
                }

                this->communicate_noncorner_node_data( tFieldLabels, tWork );

                if ( mMyRank == mKernel->master() )
                {

                    message( 4, "    ... time for linear data projection         : %u ms\n",
                             ( unsigned int ) tTimer.stop() );

                }
            }

//-----------------------------------------------------------------------------

            void
            FieldData::distribute( const Cell< string > & aFieldLabels )
            {
                // get number of fields
                uint tNumberOfFields = aFieldLabels.size();

                // get master proc
                proc_t tMaster = mKernel->master();

                // check if data must be projected
                if( mParent->enforce_linear_interpolation() )
                {
                    this->project_linear_field_to_higher_mesh( aFieldLabels );
                }

                if ( mMyRank == tMaster )
                {

                    // get number of procs
                    uint tNumberOfProcs = mKernel->number_of_procs();

                    // sanity check
                    if( tNumberOfProcs > 1 )
                    {
                        Vector< uint > tNumFieldsPerProc;
                        receive( mKernel->comm_table(), tNumFieldsPerProc );
                        for ( uint k = 0; k < mKernel->comm_table().length(); ++k )
                        {
                            BELFEM_ERROR( tNumberOfFields == tNumFieldsPerProc( k ) ||
                                                 mKernel->comm_table( k ) == mKernel->master(),
                                         "Number of fields on proc %u does not match. Is %u but expect %u",
                                         ( unsigned int ) mKernel->comm_table( k ),
                                         ( unsigned int ) tNumFieldsPerProc( k ),
                                         ( unsigned int ) tNumberOfFields );

                        }
                    }

                    // container with data to send
                    Cell< Vector< real > > tData( tNumberOfProcs, {} );

                    // loop over all fields
                    for( uint f=0; f<tNumberOfFields; ++f )
                    {
                        mesh::Field * tF = mMesh->field( aFieldLabels( f ) );

                        // grab field data
                        Vector< real > & tField = tF->data() ;

                        index_t tMultiplicity = 0 ;

                        // loop over all procs
                        for( uint p=1; p<tNumberOfProcs; ++p )
                        {
                            const Vector< index_t > & tIndices = this->field_indices( tF->entity_type(), p , tMultiplicity );

                            // get number of nodes
                            index_t tNumberOfEntities = tIndices.length()  ;

                            // get values
                            Vector< real > & tValues = tData( p );

                            // set size for values
                            tValues.set_size( tNumberOfEntities * tMultiplicity );

                            // populate data
                            if( tMultiplicity == 1 )
                            {
                                for( index_t k=0; k<tNumberOfEntities; ++k )
                                {
                                    tValues( k ) = tField( tIndices( k ) );
                                }
                            }
                            else
                            {
                                index_t tCount = 0 ;
                                for( index_t k=0; k<tNumberOfEntities; ++k )
                                {
                                    for( index_t i=0; i<tMultiplicity; ++i )
                                    {
                                        tValues( tCount++ ) = tField( tIndices( k ) * tMultiplicity + i );
                                    }
                                }
                            }

                        }

                        // send data to othere procs
                        send( mKernel->comm_table(), tData );
                    }
                }
                else
                {
                    send( tMaster, tNumberOfFields );

                    // loop over all fields
                    for( uint f=0; f<tNumberOfFields; ++f )
                    {
                        // get data
                        receive( tMaster, mMesh->field_data( aFieldLabels( f ) ) );
                    }
                }
            }

//-----------------------------------------------------------------------------

            const Vector< index_t > &
            FieldData::field_indices( const EntityType  aType,
                                      const proc_t      aTarget,
                                                 uint & aMultiplicity )
            {
                switch ( aType )
                {
                    case( EntityType::NODE ) :
                    {
                        aMultiplicity = 1 ;
                        return mKernel->node_table( aTarget );
                    }
                    case( EntityType::EDGE ) :
                    {
                        aMultiplicity = mParent->iwg()->edge_multiplicity() ;
                        return mKernel->edge_table( aTarget );
                    }
                    case( EntityType::FACE ) :
                    {
                        aMultiplicity = mParent->iwg()->face_multiplicity() ;
                        return mKernel->face_table( aTarget );
                    }
                    case( EntityType::FACET ) :
                    {
                        aMultiplicity = mParent->iwg()->lambda_multiplicity() ;
                        return mKernel->facet_table( aTarget );
                    }
                    case( EntityType::CELL ) :
                    {
                        aMultiplicity = mParent->iwg()->cell_multiplicity() ;
                        return mKernel->element_table( aTarget );
                    }
                    case( EntityType::ELEMENT ) :
                    {
                        aMultiplicity = 1 ;
                        return mKernel->element_table( aTarget );
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "invalid entity type" );
                        aMultiplicity = 0 ;
                        // need to return something here so that the compiler
                        // is happy
                        return mKernel->node_table( 0 );
                    }
                }
            }

//-----------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */
