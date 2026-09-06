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
                    mCommRank( comm_rank() ),
                    mCommSize( comm_size() )
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
                mMyNumberOfOwnedElements = 0 ;
                mNodeOwnerList.clear();
                mElementOwnerList.clear();
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

                    tDof->set_field_index( tFieldIndices( mParent->iwg()->get_field_index( tDof->type_id() ) ) );
                }

            }

//------------------------------------------------------------------------------

            void
            FieldData::collect_node_owners()
            {
                if( mCommSize < 2 ) return ;

                if ( mCommRank != 0 )
                {
                    // count owned nodes
                    mMyNumberOfOwnedNodes = 0;

                    // loop over all nodes on mesh
                    for ( mesh::Node * tNode : mMesh->nodes() )
                    {
                        if ( tNode->owner() == mCommRank )
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
                        if ( tNode->owner() == mCommRank )
                        {
                            tIDs( tCount++ ) = tNode->id();
                        }
                    }

                    comm_barrier() ;
                    send( tIDs );

                }
                else
                {
                    uint tNumProcs = comm_size();

                    Cell< Vector< id_t > > tAllIDs( tNumProcs, {} );

                    comm_barrier() ;
                    belfem::collect( tAllIDs );

                    // allocate index vector
                    mNodeOwnerList.set_size( tNumProcs, {} );

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

//------------------------------------------------------------------------------

            void FieldData::collect_element_owners()
            {
                if( mCommSize < 2 ) return ;

                mMyNumberOfOwnedElements = 0;
                for ( mesh::Element * tElement : mMesh->elements() )
                {
                    if ( tElement->owner() == mCommRank )
                    {
                        ++mMyNumberOfOwnedElements;
                    }
                }

                if ( mCommRank != 0 )
                {
                    Vector< id_t > tIDs( mMyNumberOfOwnedElements );
                    index_t tCount = 0;

                    // collect IDs of owned nodes
                    for ( mesh::Element * tElement : mMesh->elements() )
                    {
                        if ( tElement->owner() == mCommRank )
                        {
                            tIDs( tCount++ ) = tElement->id();
                        }
                    }

                    comm_barrier() ;
                    send( tIDs );
                }
                else
                {
                    uint tNumProcs = comm_size();

                    Cell< Vector< id_t > > tAllIDs( tNumProcs, {} );

                    comm_barrier() ;
                    belfem::collect( tAllIDs );

                    // allocate index vector
                    mElementOwnerList.set_size( tNumProcs, {} );

                    // loop over all ids
                    for ( uint p = 0; p < tNumProcs; ++p )
                    {
                        Vector< index_t > & tIndices = mElementOwnerList( p );
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

//------------------------------------------------------------------------------

            void
            FieldData::collect( const string & aLabel )
            {
                BELFEM_ERROR( mMesh->field( aLabel )->entity_type() == EntityType::NODE ||
                     mMesh->field( aLabel )->entity_type() == EntityType::ELEMENT,
                    "Field %s is off type %s but must be node or element to collect",
                    aLabel.c_str(), to_string( mMesh->field( aLabel )->entity_type() ).c_str() );


                if( mCommSize > 1 )
                {
                    // grab field
                    Vector< real > & tField = mMesh->field_data( aLabel );


                    if ( mCommRank != 0 )
                    {
                        Vector< real > tSubField ;
                        index_t tCount = 0 ;


                        switch ( mMesh->field( aLabel )->entity_type() )
                        {

                            case EntityType::NODE :
                            {
                                tSubField.set_size( mMyNumberOfOwnedNodes );
                                // collect owned data
                                for ( mesh::Node * tNode: mMesh->nodes())
                                {
                                    if ( tNode->owner() == mCommRank )
                                    {
                                        tSubField( tCount++ ) = tField( tNode->index());
                                    }
                                }
                                break ;
                            }
                            case EntityType::ELEMENT :
                            {
                                tSubField.set_size( mMyNumberOfOwnedElements );

                                for ( mesh::Element * tElement: mMesh->elements())
                                {
                                    if ( tElement->owner() == mCommRank )
                                    {
                                        tSubField( tCount++ ) = tField( tElement->index() );
                                    }
                                }
                                break;
                            }
                            default:
                                break;
                        }

                        comm_barrier() ;


                        // send data to master
                        send( tSubField );
                    }
                    else
                    {
                        // collect subfields
                        Cell< Vector< real > > tSubFields;

                        comm_barrier() ;


                        belfem::collect( tSubFields );

                        // assemble
                        uint tNumProcs = comm_size();

                        for ( uint p = 0; p < tNumProcs; ++p )
                        {
                            // grab index vector
                            const Vector< index_t > & tIndices =
                                mMesh->field( aLabel )->entity_type() == EntityType::NODE ? mNodeOwnerList( p ) : mElementOwnerList( p );

                            const Vector< real > & tSubField = tSubFields( p );

                            index_t tNumNodes = tIndices.length();

                            for ( index_t k = 0; k < tNumNodes; ++k )
                            {
                                tField( tIndices( k )) = tSubField( k );
                            }
                        }
                    }
                }

                comm_barrier() ;
            }

//-----------------------------------------------------------------------------

            void
            FieldData::collect( const Cell< string > & aLabels )
            {
                if( mCommSize > 1 )
                {
                    uint tNumFields = aLabels.size();

                    if ( mCommRank != 0 )
                    {
                        // field counters
                        uint tNodeFieldCount = 0 ;
                        uint tElementFieldCount = 0 ;

                        // count fields
                        for ( uint f = 0; f < tNumFields; ++f )
                        {
                            switch( mMesh->field( aLabels( f ) )->entity_type() )
                            {
                                case( EntityType::NODE ) :
                                {
                                    ++tNodeFieldCount ;
                                    break ;
                                }
                                case( EntityType::EDGE ) :
                                {
                                    continue;
                                }
                                case( EntityType::FACE ) :
                                {
                                    continue;
                                }
                                case( EntityType::ELEMENT ) :
                                {
                                    ++tElementFieldCount ;
                                    break ;
                                }
                                default :
                                {
                                    BELFEM_ASSERT( mCommRank != 0,
                                                   "Field %s is off type %s but must be node or element to collect",
                                                   aLabels( f ).c_str(),
                                                   to_string( mMesh->field( aLabels( f ) )->entity_type() ).c_str() );

                                }
                            }
                        }

                        Matrix< real > tNodeData( mMyNumberOfOwnedNodes, tNodeFieldCount );
                        Matrix< real > tElementData( mMyNumberOfOwnedElements, tElementFieldCount );

                        tNodeFieldCount = 0 ;
                        tElementFieldCount = 0 ;

                        for ( uint f = 0; f < tNumFields; ++f )
                        {
                            Vector< real > & tField = mMesh->field_data( aLabels( f ) );

                            index_t tCount = 0;

                            switch ( mMesh->field( aLabels( f ) )->entity_type() )
                            {
                                case EntityType::NODE:
                                {
                                    for ( mesh::Node * tNode: mMesh->nodes())
                                    {
                                        if ( tNode->owner() == mCommRank )
                                        {
                                            tNodeData( tCount++, tNodeFieldCount ) = tField( tNode->index() );
                                        }
                                    }
                                    ++tNodeFieldCount;
                                    break ;
                                }
                                case EntityType::ELEMENT:
                                {
                                    for ( mesh::Element * tElement: mMesh->elements())
                                    {
                                        if ( tElement->owner() == mCommRank )
                                        {
                                            tElementData( tCount++, tElementFieldCount ) = tField( tElement->index() );
                                        }
                                    }
                                    ++tElementFieldCount;
                                    break ;
                                }
                                default:
                                {
                                    break ;
                                }
                            }
                        }

                        comm_barrier() ;

                        // send data to master
                        send( tNodeData );
                        send( tElementData );

                    }
                    else
                    {
                        // collect subfields
                        Cell< Matrix< real > > tAllNodeData ;
                        Cell< Matrix< real > > tAllElementData ;
                        comm_barrier() ;

                        belfem::collect( tAllNodeData );
                        belfem::collect( tAllElementData );

                        // assemble
                        uint tNumProcs = comm_size();

                        for ( uint p = 1; p < tNumProcs; ++p )
                        {
                            // field counters
                            uint tNodeFieldCount = 0 ;
                            uint tElementFieldCount = 0 ;

                            // grab index vector
                            const Vector< index_t > & tNodeIndices    = mNodeOwnerList( p );
                            const Matrix< real >    & tNodeFields     = tAllNodeData( p );

                            const Vector< index_t > & tElementIndices    = mElementOwnerList( p );
                            const Matrix< real >    & tElementFields     = tAllElementData( p );

                            index_t tNumNodes    = tNodeIndices.length();
                            index_t tNumElements = tElementIndices.length();

                            for ( uint f = 0; f < tNumFields; ++f )
                            {
                                switch (  mMesh->field( aLabels( f ) )->entity_type() )
                                {
                                    case EntityType::NODE:
                                    {
                                        Vector< real > & tField = mMesh->field_data( aLabels( f ) );
                                        for ( index_t k = 0; k < tNumNodes; ++k )
                                        {
                                            tField( tNodeIndices( k ) ) = tNodeFields( k, tNodeFieldCount );
                                        }
                                        ++tNodeFieldCount ;
                                        break;
                                    }
                                    case EntityType::ELEMENT:
                                    {
                                        Vector< real > & tField = mMesh->field_data( aLabels( f ) );
                                        for ( index_t k = 0; k < tNumElements; ++k )
                                        {
                                            tField( tElementIndices( k ) ) = tElementFields( k, tElementFieldCount );
                                        }
                                        ++tElementFieldCount ;
                                        break;
                                    }
                                    default:
                                    {
                                        continue;
                                    }
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

                const bool tIsMaster = mCommRank == 0 ;

                // get node container of mesh
                Cell< mesh::Node * > & tNodes = mMesh->nodes();

                const proc_t tNumProcs = comm_size();

                if ( tIsMaster )
                {
                    // allocate id and index containers
                    tAllNodeIDs.set_size( tNumProcs, {} );
                    mAllCornerNodeIndices.set_size( tNumProcs, {} );

                    // - - - - - - - - - - - - - - - - - - - -
                    // STEP 1 : identify corner nodes per proc
                    // - - - - - - - - - - - - - - - - - - - -

                    for ( proc_t p = 1; p < tNumProcs; ++p )
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
                    belfem::distribute( tAllNodeIDs );

                }
                else
                {
                    // wait
                    comm_barrier() ;

                    // get corner node IDs from master
                    receive( tMyNodeIDs );
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
                    // allocate index containers
                    mAllNonCornerNodeIndices.set_size( tNumProcs, {} );

                    // - - - - - - - - - - - - - - - - - - - -
                    // STEP 1 : identify corner nodes per proc
                    // - - - - - - - - - - - - - - - - - - - -

                    for( proc_t p = 1; p < tNumProcs; ++p )
                    {
                        mMesh->unflag_all_nodes();

                        Vector< id_t >    & tNodeIDs     = tAllNodeIDs( p );
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
                    belfem::distribute( tAllNodeIDs );
                }
                else
                {
                    // wait
                    comm_barrier() ;
                    belfem::receive( tMyNodeIDs );
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
                    belfem::collect( tAllNodeIDs );

                    // allocate array
                    mAllNonCornerNodeIndices.set_size( tNumProcs, {} );

                    // create indices
                    for( proc_t p=1; p<tNumProcs; ++p )
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
                    send( tMyNodeIDs );
                }

                // wait
                comm_barrier() ;

                if ( mCommRank == 0 )
                {
                    message( InfoLevel::Verbose, "    ... time for collecting procjection nodes   : %u ms\n",
                             ( unsigned int ) tTimer.stop() );
                }
            }
//------------------------------------------------------------------------------

            void
            FieldData::communicate_corner_node_data( const Cell< string > & aFieldLabels )
            {
                index_t tNumFields = aFieldLabels.size() ;
                uint tNumProcs = mKernel->number_of_procs() ;

                if ( mCommRank == 0 )
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

                    belfem::distribute( tAllData );
                }
                else
                {
                    // Data container
                    Vector< real > tData ;

                    // wait
                    comm_barrier() ;
                    receive( tData );

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

                if( mCommRank == 0 )
                {

                    // container with other data
                    Cell< Matrix< real > > tAllData( mKernel->number_of_procs(), {} );

                    belfem::collect( tAllData );

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
                    send( aData );
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

                if ( mCommRank == 0 )
                {

                    message( InfoLevel::Verbose, "    ... time for linear data projection         : %u ms\n",
                             ( unsigned int ) tTimer.stop() );

                }
            }

//-----------------------------------------------------------------------------

            void
            FieldData::distribute( const Cell< string > & aFieldLabels )
            {
                // get number of fields
                uint tNumberOfFields = aFieldLabels.size();

                if ( mCommRank == 0 )
                {

                    // get number of procs
                    uint tNumberOfProcs = mKernel->number_of_procs();

                    // sanity check
                    if( tNumberOfProcs > 1 )
                    {
                        Vector< uint > tNumFieldsPerProc;
                        belfem::collect( tNumFieldsPerProc );
                        for ( proc_t k = 1; k < comm_size(); ++k )
                        {
                            BELFEM_ERROR( tNumberOfFields == tNumFieldsPerProc( k ),
                                         "Number of fields on proc %u does not match. Is %u but expect %u",
                                         ( unsigned int ) k,
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
                            const Cell< index_t > & tIndices = this->field_indices( tF->entity_type(), p , tMultiplicity );

                            // get number of nodes
                            index_t tNumberOfEntities = tIndices.size()  ;

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

                        // send data to other procs
                        belfem::distribute( tData );
                    }
                }
                else
                {
                    send( tNumberOfFields );

                    // loop over all fields
                    for( uint f=0; f<tNumberOfFields; ++f )
                    {
                        // get data
                        receive( mMesh->field_data( aFieldLabels( f ) ) );
                    }
                }
            }

//-----------------------------------------------------------------------------

            const Cell< index_t > &
            FieldData::field_indices( const EntityType  aType,
                                      const proc_t      aTarget,
                                                 uint & aMultiplicity )
            {
                switch ( aType )
                {
                    case( EntityType::NODE ) :
                    {
                        aMultiplicity = 1 ;
                        return mKernel->comm_table( aTarget )->nodes();
                    }
                    case( EntityType::EDGE ) :
                    {
                        aMultiplicity = mParent->iwg()->edge_multiplicity() ;
                        return mKernel->comm_table( aTarget )->edges();
                    }
                    case( EntityType::FACE ) :
                    {
                        aMultiplicity = mParent->iwg()->face_multiplicity() ;
                        return mKernel->comm_table( aTarget )->faces();
                    }
                    case( EntityType::FACET ) :
                    {
                        aMultiplicity = mParent->iwg()->lambda_multiplicity() ;
                        return mKernel->comm_table( aTarget )->facets();
                    }
                    case( EntityType::CELL ) :
                    {
                        aMultiplicity = mParent->iwg()->cell_multiplicity() ;
                        return mKernel->comm_table( aTarget )->elements();
                    }
                    case( EntityType::ELEMENT ) :
                    {
                        aMultiplicity = 1 ;
                        return mKernel->comm_table( aTarget )->elements();
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "invalid entity type" );
                        aMultiplicity = 0 ;
                        // need to return something here so that the compiler
                        // is happy
                        return mKernel->comm_table( aTarget )->elements();
                    }
                }
            }

//-----------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */
