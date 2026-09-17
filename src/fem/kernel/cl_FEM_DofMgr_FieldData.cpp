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
                const Vector< index_t > & tTypes = mParent->iwg()->default_dof_types();

                BELFEM_ERROR( tTypes.length() > 0 , "Fields have not been set" );

                const Cell< string > & tLabels = mParent->iwg()->dof_fields() ;

                Vector< index_t > tFieldIndices(  max( tTypes ) + 1, gNoIndex );

                for( uint k=0; k<tLabels.size(); ++k )
                {
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
                    mMyNumberOfOwnedNodes = 0;

                    for ( mesh::Node * tNode : mMesh->nodes() )
                    {
                        if ( tNode->owner() == mCommRank )
                        {
                            ++mMyNumberOfOwnedNodes;
                        }
                    }

                    Vector< id_t > tIDs( mMyNumberOfOwnedNodes );
                    index_t tCount = 0;

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

                    mNodeOwnerList.set_size( tNumProcs, {} );

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

                    mElementOwnerList.set_size( tNumProcs, {} );

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

                        send( tSubField );
                    }
                    else
                    {
                        Cell< Vector< real > > tSubFields;

                        comm_barrier() ;

                        belfem::collect( tSubFields );

                        uint tNumProcs = comm_size();

                        for ( uint p = 0; p < tNumProcs; ++p )
                        {
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
                        uint tNodeFieldCount = 0 ;
                        uint tElementFieldCount = 0 ;

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

                        send( tNodeData );
                        send( tElementData );

                    }
                    else
                    {
                        Cell< Matrix< real > > tAllNodeData ;
                        Cell< Matrix< real > > tAllElementData ;
                        comm_barrier() ;

                        belfem::collect( tAllNodeData );
                        belfem::collect( tAllElementData );

                        uint tNumProcs = comm_size();

                        for ( uint p = 1; p < tNumProcs; ++p )
                        {
                            uint tNodeFieldCount = 0 ;
                            uint tElementFieldCount = 0 ;

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
                const Vector< id_t > & tSelectedBlocks
                    = mParent->iwg()->selected_blocks();

                Cell< Vector< id_t > > tAllNodeIDs;
                Vector< id_t > tMyNodeIDs;

                Timer tTimer;

                const bool tIsMaster = mCommRank == 0 ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();

                const proc_t tNumProcs = comm_size();

                if ( tIsMaster )
                {
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

                        for ( index_t tID : tSelectedBlocks )
                        {
                            Cell< mesh::Element * > & tElements = mMesh->block( tID )->elements();

                            for ( mesh::Element * tElement : tElements )
                            {
                                if( tElement->owner() == p )
                                {
                                    tElement->flag_corner_nodes();
                                }
                            }
                        }

                        index_t tCount = 0 ;

                        for ( mesh::Node * tNode : tNodes )
                        {
                            if ( tNode->is_flagged() )
                            {
                               ++tCount ;
                            }
                        }

                        if( tCount > 0 )
                        {
                            tNodeIDs.set_size( tCount );
                            tNodeIndices.set_size( tCount );

                            tCount = 0;

                            for ( mesh::Node * tNode: tNodes )
                            {
                                if ( tNode->is_flagged() )
                                {
                                    tNodeIDs( tCount ) = tNode->id() ;
                                    tNodeIndices( tCount++ ) = tNode->index() ;
                                }
                            }
                        }
                    }
                    comm_barrier() ;

                    belfem::distribute( tAllNodeIDs );

                }
                else
                {
                    comm_barrier() ;

                    receive( tMyNodeIDs );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 2 : identify indices of corner nodes
                //          ( this is the same for all procs )
                // - - - - - - - - - - - - - - - - - - - - - - - -

                Vector< id_t > & tCornerNodeIDs = tIsMaster ? tAllNodeIDs( 0 ) : tMyNodeIDs ;

                mMyCornerNodeIndices.set_size( tCornerNodeIDs.length() );

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
                    mAllNonCornerNodeIndices.set_size( tNumProcs, {} );

                    // - - - - - - - - - - - - - - - - - - - -
                    // STEP 1 : identify corner nodes per proc
                    // - - - - - - - - - - - - - - - - - - - -

                    for( proc_t p = 1; p < tNumProcs; ++p )
                    {
                        mMesh->unflag_all_nodes();

                        Vector< id_t >    & tNodeIDs     = tAllNodeIDs( p );
                        Vector< index_t > & tNodeIndices = mAllNonCornerNodeIndices( p ) ;

                        for ( index_t tID : tSelectedBlocks )
                        {
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

                        index_t tCount = 0 ;

                        for ( mesh::Node * tNode : tNodes )
                        {
                            if ( tNode->is_flagged() )
                            {
                                ++tCount ;
                            }
                        }

                        if( tCount > 0 )
                        {
                            tNodeIDs.set_size( tCount );
                            tNodeIndices.set_size( tCount );

                            tCount = 0;

                            for ( mesh::Node * tNode: tNodes )
                            {
                                if ( tNode->is_flagged() )
                                {
                                    tNodeIDs( tCount ) = tNode->id() ;
                                    tNodeIndices( tCount++ ) = tNode->index() ;
                                }
                            }
                        }
                    }
                    comm_barrier() ;

                    belfem::distribute( tAllNodeIDs );
                }
                else
                {
                    comm_barrier() ;
                    belfem::receive( tMyNodeIDs );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 4 : correct order of node IDs and indices
                // - - - - - - - - - - - - - - - - - - - - - - - - -

                mMesh->unflag_all_nodes() ;

                Vector< id_t > & tMyNonCornerNodeIDs = tIsMaster ? tAllNodeIDs( 0 ) : tMyNodeIDs ;

                for( id_t tID : tMyNonCornerNodeIDs )
                {
                    mMesh->node( tID )->flag() ;
                }

                k = 0 ;

                mMyNonCornerNodeIndices.set_size( tMyNonCornerNodeIDs.length() );

                mMyNonCornerNodes.set_size( tMyNonCornerNodeIDs.length(), nullptr );

                tMyNodeIDs.set_size( tMyNonCornerNodeIDs.length() );

                for ( index_t tID : tSelectedBlocks )
                {
                    if( mParent->block_exists( tID ) )
                    {
                        Block * tBlock = mParent->block( tID );

                        // get real element type of block ( tBlock->element_type() is the linear one )
                        ElementType tType = tBlock->block()->element_type();

                        uint tNumCornerNodes = mesh::number_of_corner_nodes( tType );
                        uint tNumNodes = mesh::number_of_nodes( tType );

                        Cell< mesh::Element * > & tElements = tBlock->block()->elements();

                        for ( mesh::Element * tElement: tElements )
                        {
                            for ( uint i = tNumCornerNodes; i < tNumNodes; ++i )
                            {
                                mesh::Node * tNode = tElement->node( i );

                                if ( tNode->is_flagged() )
                                {
                                    mMyNonCornerNodeIndices( k ) = tNode->index();

                                    tMyNodeIDs( k ) = tNode->id();

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

                comm_barrier() ;

                // - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 5 : send indices in correct order back to master
                // - - - - - - - - - - - - - - - - - - - - - - - - -

                if( tIsMaster )
                {
                    belfem::collect( tAllNodeIDs );

                    mAllNonCornerNodeIndices.set_size( tNumProcs, {} );

                    for( proc_t p=1; p<tNumProcs; ++p )
                    {
                        Vector< index_t > & tNodeIDs = tAllNodeIDs( p );

                        Vector< index_t > & tIndices = mAllNonCornerNodeIndices( p );

                        tIndices.set_size( tNodeIDs.length() );

                        k = 0 ;

                        for( id_t tID : tNodeIDs )
                        {
                            tIndices( k++ ) = mMesh->node( tID )->index() ;
                        }
                    }
                }
                else
                {
                    send( tMyNodeIDs );
                }

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

                    for( uint p=1; p<tNumProcs; ++p )
                    {
                        const Vector< index_t > & tIndices = mAllCornerNodeIndices( p );

                        Vector< real > & tData = tAllData( p );

                        index_t tNumNodes = tIndices.length() ;

                        tData.set_size( tNumFields * tIndices.length() );

                        index_t tCount = 0 ;

                        for( index_t f=0; f<tNumFields; ++f )
                        {
                            Vector< real > & tField = mMesh->field_data( aFieldLabels( f ) );

                            for( index_t k=0; k<tNumNodes; ++k )
                            {
                                tData( tCount++ ) = tField( tIndices( k ) );
                            }
                        }
                    }

                    comm_barrier();

                    belfem::distribute( tAllData );
                }
                else
                {
                    Vector< real > tData ;

                    comm_barrier() ;
                    receive( tData );

                    index_t tCount = 0 ;

                    index_t tNumNodes = mMyCornerNodeIndices.length() ;

                    for( index_t f=0; f<tNumFields; ++f )
                    {
                        Vector< real > & tField = mMesh->field_data( aFieldLabels( f ) );

                        for( index_t k=0; k<tNumNodes; ++k )
                        {
                            tField( mMyCornerNodeIndices( k ) ) = tData( tCount++ );
                        }
                    }

                    BELFEM_ASSERT( tCount == tData.length(), "Invalid vector size. Is %lu but expect %lu",
                                  ( long unsigned int ) tCount,
                                  ( long unsigned int ) tData.length() );

                }

                comm_barrier() ;
            }

//------------------------------------------------------------------------------

            void
            FieldData::communicate_noncorner_node_data(
                    const Cell< string > & aFieldLabels,
                    Matrix< real > & aData )
            {
                index_t tNumFields = aFieldLabels.size() ;
                uint tNumProcs = mKernel->number_of_procs();

                if( mCommRank == 0 )
                {

                    Cell< Matrix< real > > tAllData( mKernel->number_of_procs(), {} );

                    belfem::collect( tAllData );

                    index_t tNumNodes = mMyNonCornerNodeIndices.length() ;

                    for( index_t f=0; f<tNumFields; ++f )
                    {
                        Vector< real > & tData = mMesh->field_data( aFieldLabels( f ) );

                        for( index_t k=0; k<tNumNodes; ++k )
                        {
                            tData( mMyNonCornerNodeIndices( k ) ) = aData( k, f );
                        }
                    }

                    for( uint p=1; p<tNumProcs; ++p )
                    {
                        Matrix< real > & tData = tAllData( p );

                        Vector< index_t > & tIndices = mAllNonCornerNodeIndices( p );

                        index_t tN = tIndices.length() ;

                        for( index_t f=0; f<tNumFields; ++f )
                        {
                            Vector< real > & tField = mMesh->field_data( aFieldLabels( f ) );

                            for( index_t k=0; k<tN; ++k )
                            {
                                tField( tIndices( k ) ) = tData( k, f );
                            }
                        }
                    }
                }
                else
                {
                    send( aData );
                }

                comm_barrier() ;
            }

//------------------------------------------------------------------------------

            void
            FieldData::project_linear_field_to_higher_mesh(
                    const Cell< string > & aFieldLabels )
            {
                Timer tTimer ;

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

                mMesh->unflag_all_nodes() ;

                for( mesh::Node * tNode : mMyNonCornerNodes )
                {
                    tNode->flag() ;
                }

                InterpolationFunctionFactory tFactory ;

                const index_t tNumFields = tFieldLabels.size() ;

                const index_t tNumNonCornerNodes = mMyNonCornerNodes.size() ;

                Matrix< real > tWork( tNumNonCornerNodes, tNumFields );

                index_t tCount = 0 ;

                for( id_t tBlockID : mParent->iwg()->selected_blocks() )
                {
                    if( mParent->block_exists( tBlockID ) )
                    {
                        mesh::Block * tBlock = mParent->block( tBlockID )->block() ;

                        const ElementType tType = tBlock->element_type();

                        InterpolationFunction * tLinShape = tFactory.create_lagrange_function(
                                mesh::linear_element_type( tType ));

                        InterpolationFunction * tHighShape = tFactory.create_lagrange_function(
                                tType );

                        Matrix< real > tXi;
                        tHighShape->param_coords( tXi );

                        uint tNumNodes = mesh::number_of_nodes( tType );

                        uint tNumCornerNodes = mesh::number_of_corner_nodes( tType );

                        Cell< Matrix< real > > tAllN( tXi.n_cols(), {} );

                        for ( uint k = 0; k < tNumNodes; ++k )
                        {
                            Matrix< real > & tN = tAllN( k );

                            tN.set_size( 1, tNumCornerNodes );

                            tLinShape->N( tXi.col( k ), tN );
                        }

                        Cell< mesh::Element * > & tElements = tBlock->elements();

                        Matrix< real > tLinData( tNumCornerNodes, tNumFields );

                        Matrix< real > tNodeData( 1, tNumFields );

                        for ( mesh::Element * tElement : tElements )
                        {
                            for ( uint j = 0; j < tNumFields; ++j )
                            {
                                Vector< real > & tField = mMesh->field_data( tFieldLabels( j ));

                                for ( uint i = 0; i < tNumCornerNodes; ++i )
                                {
                                    tLinData( i, j ) = tField( tElement->node( i )->index() );
                                }
                            }

                            for ( uint k = tNumCornerNodes; k < tNumNodes; ++k )
                            {
                                mesh::Node * tNode = tElement->node( k );

                                if ( tNode->is_flagged() )
                                {
                                    tNodeData = tAllN( k ) * tLinData;

                                    tWork.set_row( tCount++, tNodeData.row( 0 ) );

                                    tNode->unflag();
                                }
                            }
                        }
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
                uint tNumberOfFields = aFieldLabels.size();

                if ( mCommRank == 0 )
                {

                    uint tNumberOfProcs = mKernel->number_of_procs();

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

                    Cell< Vector< real > > tData( tNumberOfProcs, {} );

                    for( uint f=0; f<tNumberOfFields; ++f )
                    {
                        mesh::Field * tF = mMesh->field( aFieldLabels( f ) );

                        Vector< real > & tField = tF->data() ;

                        index_t tMultiplicity = 0 ;

                        for( uint p=1; p<tNumberOfProcs; ++p )
                        {
                            const Cell< index_t > & tIndices = this->field_indices( tF->entity_type(), p , tMultiplicity );

                            index_t tNumberOfEntities = tIndices.size()  ;

                            Vector< real > & tValues = tData( p );

                            tValues.set_size( tNumberOfEntities * tMultiplicity );

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

                        belfem::distribute( tData );
                    }
                }
                else
                {
                    send( tNumberOfFields );

                    for( uint f=0; f<tNumberOfFields; ++f )
                    {
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
