//
// Created by Christian Messe on 02.11.19.
//
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Kernel.hpp"
#include "commtools.hpp"
#include "cl_Timer.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Logger.hpp"
#include "fn_norm.hpp"
#include "fn_max.hpp"
#include "cl_IWG_TATCAD.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_entity_type.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Field::Field(
                Kernel * aParent,
                Mesh * aMesh,
                const Vector< id_t > & aBlockIDs,
                const uint             aNumberOfDOFsPerNode,
                const uint             aNumberOfDOFsPerEdge,
                const bool             aEnforceLinear ) :
                DofManagerBase( DofManagerType::OLD, aParent, aMesh ),
                mBlockIDs( aBlockIDs ),
                mNumberOfDOFsPerNode( aNumberOfDOFsPerNode ),
                mNumberOfDOFsPerEdge( aNumberOfDOFsPerEdge ),
                mEnforceLinear( aEnforceLinear ),
                mBlockIntegrationOrder(
                        aParent->params()->block_integration_order(
                                aMesh->number_of_fields() - aParent->field_offset() ) ),
                mSideSetIntegrationOrder(
                        aParent->params()->sideset_integration_order(
                                aMesh->number_of_fields() - aParent->field_offset() ) ),
                mIntegrationScheme( aParent->params()->integration_scheme() )
        {
            // associate block IDs with indices on mesh
            // distributed meshes only contain blocks
            // with used elements
            this->create_blockmap();

            // create the DOFs of this problem
            // the master proc knows all DOFs, the other procs
            // only have the dofs that are associated with their elements.
            // This routine also populates the DOF map
            this->create_dofs();

            // create the blocks that contain the finite elements
            // that are used for the calculation of the Jacobian
            this->create_blocks();

            // tables tell the master which dofs exist on which proc
            // and in what order they are stored
            this->create_dof_table();

            // create the sidesets that are needed to set boundary conditions.
            // the master proc knows all sidesets, the other procs
            // know only the part of the sideset that is connected to it
            this->create_sidesets();

            // create the bearings that this proc knows. Same logic as sidesets
            this->create_bearings();

            // the empty block is a dummy that is exposed if
            // a set is accessed that doesn't exist on this proc
            mEmptyBlock = new fem::Block( this );

            // the empty sideset is a dummy that is exposed if
            // a set is accesed that doesn't exist on this proc
            Cell< mesh::Facet * > tEmpty;

            mEmptySideset = new fem::SideSet(
                    this,
                    0,
                    tEmpty );

            // same for bearing
            mEmptyBearing = new fem::Bearing( this );

            // in order to collect fields by proc, we need the node owner table
            // this is needed for Neumann boundary conditions
            this->collect_node_owners() ;

            // check if the linear projection flag is on
            if( this->enforce_linear_interpolation() )
            {
                // collect critical nodes
                this->initialize_linear_projection_lists() ;
            }
        }

//------------------------------------------------------------------------------

        Field::~Field()
        {
            if( mSolver != NULL )
            {
                delete mSolver ;
            }

            if ( mJacobian != NULL )
            {
                delete mJacobian;
            }

            if ( mDirichletMatrix != NULL )
            {
                delete mDirichletMatrix;
            }

            for ( Block * tBlock : mBlocks )
            {
                delete tBlock;
            }
            for ( SideSet * tSideSet : mSideSets )
            {
                delete tSideSet;
            }

            for( Bearing * tBearing : mBearings )
            {
                delete tBearing;
            }

            if( mEmptyBlock != NULL )
            {
                delete mEmptyBlock;
            }

            if( mEmptySideset != NULL )
            {
                delete mEmptySideset;
            }

            if( mEmptyBearing != NULL )
            {
                delete mEmptyBearing;
            }

            for ( Dof * tDof : mDOFs )
            {
                delete tDof;
            }

            mMyNonCornerNodes.clear() ;
        }

//------------------------------------------------------------------------------
// User Settings
//------------------------------------------------------------------------------

        void
        Field::set_solver( const SolverType & aSolver )
        {

            // delete solver if it exists already
            if( mSolver != nullptr )
            {
                delete mSolver ;
            }

            // create a new solver
            mSolver = new Solver( aSolver );

            mSolver->set_symmetry_mode( mIWG->symmetry_mode() );
        }

//------------------------------------------------------------------------------

        void
        Field::set_integrated_weak_governing_equation( IWG * aIWG )
        {
            // tidy up memory
            if( mSolver != nullptr )
            {
                mSolver->free();
            }

            // special case if nullpointer is set
            if( aIWG == nullptr )
            {
                mIWG = nullptr ;
                return;
            }

            BELFEM_ERROR( aIWG->number_of_dofs_per_node() == mNumberOfDOFsPerNode,
                         "Number of dofs per node of IWG and field do not match ( %u vs. %u )",
                         ( unsigned int ) aIWG->number_of_dofs_per_node(),
                         ( unsigned int ) mNumberOfDOFsPerNode );


            // associate the IWG with this field
            mIWG = aIWG;

            // associate this field with the IWG
            aIWG->set_field( this );

            // create the fields on the parent mesh
            // and use the names provided by the IWG
            this->create_fields();

            // update indices for the DOFs
            // this->update_field_indices();
        }

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

        void
        Field::init_dofs()
        {
            // set DOF owners
            for( Dof * tDof : mDOFs )
            {
                tDof->set_owner( tDof->node()->owner() );
            }

            // check if user has set the wetted sidesets
            this->detect_wetted_sidesets() ;

            // count how many nodes are wet ( so that the convection table is not needed )
            this->count_wetted_nodes() ;

            // synchronize these fields
            this->distribute_fields( mIWG->dof_fields() );

            this->create_alpha_fields() ;

            // also initialize the boundary conditions
            this->set_boundary_conditions() ;

            this->collect_fields( mIWG->all_fields() );

            this->distribute_fields( mIWG->all_fields() );

            this->init_dof_values() ;

            // wait for other procs
            comm_barrier() ;
        }

//------------------------------------------------------------------------------

        void
        Field::init_dof_values()
        {
            // collect field list from IWG
            const Cell< string > & tLabels = mIWG->dof_fields() ;

            // loop over all number of dofs
            for ( uint k = 0; k < mNumberOfDOFsPerNode; ++k )
            {
                if( mMesh->field( tLabels( k ) )->entity_type() == EntityType::NODE )
                {
                    // grab field
                    Vector< real > & tValues = mMesh->field_data( tLabels( k ) );

                    // loop over all dofs
                    for ( Dof * tDof: mDOFs )
                    {
                        // check of type equals idnex
                        if ( tDof->type_id() == k )
                        {
                            if ( tDof->is_fixed())
                            {
                                tValues( tDof->node()->index()) = tDof->value();
                            }
                            else
                            {
                                tDof->value() = tValues( tDof->node()->index() );
                            }
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::set_boundary_conditions()
        {
            comm_barrier() ;

            // loop over all sidesets and update data
            for ( mesh::SideSet * tSideset : mMesh->sidesets() )
            {
                this->sideset( tSideset->id() )->set_boundary_conditions() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Field::flag_all_dofs()
        {
            for ( Dof * tDOF : mDOFs )
            {
                tDOF->flag();
            }
        }

//------------------------------------------------------------------------------

        void
        Field::flag_used_dofs()
        {
            // unflag all nodes on this mesh
            mMesh->unflag_all_nodes();

            // unflag all edges on this mesh
            mMesh->unflag_all_edges();

            // deselect all dofs on this mesh
            this->unflag_all_dofs();

            if( mNumberOfDOFsPerNode > 0 )
            {
                this->flag_nodes();
            }

            if( mNumberOfDOFsPerEdge > 0 )
            {
                this->flag_edges();
            }

            // loop over all dofs
            for ( Dof * tDOF : mDOFs )
            {
                // test if parent node is flagged
                if ( tDOF->mesh_vertex_is_flagged() )
                {
                    // flag this dof
                    tDOF->flag();
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::flag_nodes()
        {
            // get number of blocks
            uint tNumberOfBlocks = mBlockIDs.length();

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            // loop over all blocks and flag nodes of owned elements
            if( this->enforce_linear_interpolation() )
            {
                // only flag corner nodes on selected blocks
                for ( uint b = 0; b < tNumberOfBlocks; ++b )
                {
                    if ( mBlockMap.key_exists( mBlockIDs( b )))
                    {
                        // get elements from this block
                        Cell< mesh::Element * > & tElements
                                = tBlocks( mBlockMap( mBlockIDs( b )))->elements();

                        // flag all nodes
                        for ( mesh::Element * tElement : tElements )
                        {
                            tElement->flag_corner_nodes();
                        }
                    }
                }
            }
            else
            {
                // flag all nodes on selected blocks
                for ( uint b = 0; b < tNumberOfBlocks; ++b )
                {
                    if ( mBlockMap.key_exists( mBlockIDs( b )))
                    {
                        // get elements from this block
                        Cell< mesh::Element * > & tElements
                                = tBlocks( mBlockMap( mBlockIDs( b )))->elements();

                        // flag all nodes
                        for ( mesh::Element * tElement : tElements )
                        {
                            tElement->flag_nodes();
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::flag_edges()
        {
            // get number of blocks
            uint tNumberOfBlocks = mBlockIDs.length();

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            // loop over all blocks and flag nodes of owned elements
            for ( uint b = 0; b < tNumberOfBlocks; ++b )
            {
                if ( mBlockMap.key_exists( mBlockIDs( b ) ) )
                {
                    // get elements from this block
                    Cell< mesh::Element * > & tElements
                            = tBlocks( mBlockMap( mBlockIDs( b ) ) )->elements();

                    // flag all nodes
                    for ( mesh::Element * tElement : tElements )
                    {
                        tElement->flag_edges();
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::unflag_all_dofs()
        {
            for ( Dof * tDOF : mDOFs )
            {
                tDOF->unflag();
            }
        }

//------------------------------------------------------------------------------

        Block *
        Field::block( const id_t aID )
        {

            if ( mFemBlockMap.key_exists( aID ) )
            {
                return mBlocks( mFemBlockMap( aID ) );
            }
            else
            {
                mEmptyBlock->set_mesh_id( aID );
                return mEmptyBlock;
            }
        }


//------------------------------------------------------------------------------

        bool
        Field::sideset_exists( const id_t & aID ) const
        {
            return mSideSetMap.key_exists( aID ) ;
        }

//------------------------------------------------------------------------------

        SideSet *
        Field::sideset( const id_t aID )
        {
            if ( mSideSetMap.key_exists( aID ) )
            {
                return mSideSets( mSideSetMap( aID ) );
            }
            else
            {
                mEmptySideset->set_mesh_id( aID );
                return mEmptySideset;
            }
        }

//------------------------------------------------------------------------------

        Bearing *
        Field::bearing( const id_t & aID )
        {
            if( mBearingMap.key_exists( aID ) )
            {
                return mBearings( mBearingMap( aID ) );
            }
            else
            {
                return mEmptyBearing;
            }
        }

//------------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------------

        void
        Field::create_fields()
        {
            index_t tNumFields = mIWG->number_of_fields() ;

            for ( index_t f=0; f<tNumFields; ++f )
            {
                const string & tLabel = mIWG->field( f );
                /*BELFEM_ERROR( ! mMesh->field_exists( tLabel ),
                     "Field '%s' already exists on mesh.", tLabel.c_str() ); */

                if( ! mMesh->field_exists( tLabel ) )
                {
                    // check if this could be an edge field
                    EntityType tType = belfem::entity_type( tLabel );

                    // create the field
                    // Vector< real > & tField =
                    mMesh->create_field( tLabel, tType );

                    BELFEM_ERROR( tType == EntityType::NODE, "Todo: Need to adapt field lengths");

                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::create_blockmap()
        {
            index_t tNumBlocks = mMesh->number_of_blocks();

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            for ( index_t b = 0; b < tNumBlocks; ++b )
            {
                mBlockMap[ tBlocks( b )->id() ] = b;
            }

        }

//------------------------------------------------------------------------------

        void
        Field::create_dofs()
        {
            // We assume that each node on the selected block
            // is associated with dofs. However, a field
            // does not have to be present on all blocks of the mesh.
            // This is why we select only relevant dofs

            // grab blocks on mesh
            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            // the counter for the dofs on the used blocks
            index_t tNodeCount = 0;

            // store the maximum node ID here first, then multiply with
            // number of dofs per node to get offset
            mEdgeIdOffset = 0 ;

            if( mNumberOfDOFsPerNode > 0 )
            {
                // unflag all nodes on the mesh
                mMesh->unflag_all_nodes();

                // check for linear flag
                if( this->enforce_linear_interpolation() )
                {
                    // only flag corner nodes on selected blocks
                    // loop over all blocks that are relevant to the field
                    for ( uint b = 0; b < mBlockIDs.length(); ++b )
                    {
                        // check if block exists on current proc
                        if ( mBlockMap.key_exists( mBlockIDs( b )))
                        {
                            // select block
                            mesh::Block * tBlock
                                    = tBlocks( mBlockMap( mBlockIDs( b )));

                            // Flag all nodes on this block
                            tBlock->flag_corner_nodes();
                        }
                    }
                }
                else
                {
                    // flag all nodes on selected blocks
                    // loop over all blocks that are relevant to the field
                    for ( uint b = 0; b < mBlockIDs.length(); ++b )
                    {
                        // check if block exists on current proc
                        if ( mBlockMap.key_exists( mBlockIDs( b )))
                        {
                            // select block
                            mesh::Block * tBlock
                                    = tBlocks( mBlockMap( mBlockIDs( b )));

                            // Flag all nodes on this block
                            tBlock->flag_nodes();
                        }
                    }
                }

                // With the nodes flagged, we now count the number of nodes
                for ( mesh::Node * tNode : mMesh->nodes())
                {
                    if ( tNode->is_flagged() )
                    {
                        // increment counter
                        ++tNodeCount;

                        // remember maximum ID
                        mEdgeIdOffset = tNode->id() > mEdgeIdOffset ? tNode->id() : mEdgeIdOffset ;
                    }
                }
            }

            // finalize offset
            mEdgeIdOffset *= mNumberOfDOFsPerNode ;

            // the counter for the edges on the block
            index_t tEdgeCount = 0 ;

            if( mNumberOfDOFsPerEdge > 0 )
            {
                // unflag all edges on this mesh
                mMesh->unflag_all_edges() ;

                // loop over all blocks that are relevant to the field
                for ( uint b = 0; b < mBlockIDs.length(); ++b )
                {
                    // check if block exists on current proc
                    if ( mBlockMap.key_exists( mBlockIDs( b )))
                    {
                        // select block
                        mesh::Block * tBlock
                                = tBlocks( mBlockMap( mBlockIDs( b )));

                        // Flag all edges on this block
                        tBlock->flag_edges();
                    }
                }

                // With the nodes flagged, we now count the number of nodes
                for ( mesh::Edge * tEdge : mMesh->edges() )
                {
                    if ( tEdge->is_flagged() )
                    {
                        // increment counter
                        ++tEdgeCount;
                    }
                }
            }
            // now we know how big the DOF container needs to be
            mDOFs.set_size( mNumberOfDOFsPerNode * tNodeCount
                               + mNumberOfDOFsPerEdge * tEdgeCount, nullptr );

            // reset counter
            index_t tDofCount = 0 ;

            // create the DOFs in Node-wise order
            if( mNumberOfDOFsPerNode > 0 )
            {
                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    if ( tNode->is_flagged() )
                    {
                        for ( uint k = 0; k < mNumberOfDOFsPerNode; ++k )
                        {
                            // compute DOF id
                            id_t tID = this->calculate_dof_id( tNode, k );

                            // save entry in DOF map
                            mDofMap[ tID ] = tDofCount;

                            // create dof and increment the counter
                            mDOFs( tDofCount++ ) = new Dof( tID, k, tNode );
                        }
                    }
                }
            }

            // create the DOFs in Edge-wise order
            if( mNumberOfDOFsPerEdge > 0 )
            {

                for ( mesh::Edge * tEdge : mMesh->edges() )
                {
                    if ( tEdge->is_flagged())
                    {
                        for ( uint k = 0; k < mNumberOfDOFsPerEdge; ++k )
                        {
                            // compute DOF id
                            id_t tID = this->calculate_dof_id( tEdge, k );

                            // save entry in DOF map
                            mDofMap[ tID ] = tDofCount;

                            // create dof and increment the counter
                            mDOFs( tDofCount++ ) = new Dof( tID, 0, tEdge,
                                                            mNumberOfDOFsPerEdge * tEdge->index() + k );
                        }
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Field::get_field_indices( Vector< index_t > & aFieldIndices, const EntityType aEntityType )
        {
            uint tNumTypes = 0 ;
            switch( aEntityType )
            {
                case( EntityType::NODE ) :
                {
                    tNumTypes = mNumberOfDOFsPerNode ;
                    break ;
                }
                case( EntityType::EDGE ) :
                {
                    tNumTypes = mNumberOfDOFsPerEdge ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown Entity Type");
                }
            }

            aFieldIndices.set_size( tNumTypes );

            const Cell< string > & tLabels = mIWG->dof_fields();

            // initialize counter
            uint tCount = 0 ;

            // sanity check
            for( uint k=0; k<tLabels.size(); ++k )
            {
                if( mMesh->field( tLabels( k ) )->entity_type() == aEntityType )
                {
                    ++tCount ;
                }
            }

            BELFEM_ERROR( tCount == tNumTypes,
                         "Number of DOF Fields does not match FIELD: %u vs IWG: %u",
                         ( unsigned int ) tNumTypes,
                         ( unsigned int ) tCount );

            // reset counter
            tCount = 0 ;

            for( uint k=0; k<tLabels.size(); ++k )
            {
                // get field
                mesh::Field * tField = mMesh->field( tLabels( k ) );

                if( tField->entity_type() == aEntityType )
                {
                    aFieldIndices( tCount++ ) = tField->index() ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::update_field_indices()
        {
            Vector< index_t > tNodeDofFieldIndices ;
            Vector< index_t > tEdgeDofFieldIndices ;

            if( mNumberOfDOFsPerNode > 0 )
            {
                this->get_field_indices( tNodeDofFieldIndices, EntityType::NODE );
            }
            if( mNumberOfDOFsPerEdge > 0 )
            {
                this->get_field_indices( tEdgeDofFieldIndices, EntityType::EDGE );
            }
            for( Dof * tDof : mDOFs )
            {
                if( tDof->is_node() )
                {
                    tDof->set_field_index( tNodeDofFieldIndices( tDof->type_id() ) );
                }
                else if ( tDof->is_edge() )
                {
                    tDof->set_field_index( tEdgeDofFieldIndices( tDof->type_id() ) );
                }
                else
                {
                    BELFEM_ERROR( false, "Unknown DOF type" );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::create_blocks()
        {
            // unflag all nodes on this mesh
            mMesh->unflag_all_nodes();

            // get number of blocks
            uint tNumberOfBlocks = mBlockIDs.length();

            // counter for number of elements per block
            Vector< index_t > tElementsPerBlock( tNumberOfBlocks, 0 );

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            // loop over all blocks and count owned elements
            for ( uint b = 0; b < tNumberOfBlocks; ++b )
            {
                Cell< mesh::Element * > & tElements
                        = tBlocks( mBlockMap( mBlockIDs( b )))->elements();

                for ( mesh::Element * tElement : tElements )
                {
                    if ( tElement->owner() == mMyRank )
                    {
                        ++tElementsPerBlock( b );
                    }
                }
            }

            // count number of blocks that contain at least one element
            uint tCount = 0;
            for ( uint b = 0; b < tNumberOfBlocks; ++b )
            {
                if ( tElementsPerBlock( b ) > 0 )
                {
                    ++tCount;
                }
            }

            // reset block map
            mFemBlockMap.clear() ;

            // allocate block container
            mBlocks.set_size( tCount, nullptr );

            // reset counter
            tCount = 0;

            // create blocks that contain at least one used element
            for ( uint b = 0; b < tNumberOfBlocks; ++b )
            {
                if ( tElementsPerBlock( b ) > 0 )
                {
                    // get mesh block
                    mesh::Block * tBlock = tBlocks( mBlockMap( mBlockIDs( b )));

                    // add entry to map
                    mFemBlockMap[ tBlock->id() ] = tCount ;

                    // add block to container
                    mBlocks( tCount++ ) =
                            new Block( this, tBlock, tElementsPerBlock( b ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::create_dof_table()
        {
            if ( mMyRank != mParent->master() )
            {
                // any other proc collects its DOF IDs and sends them to the master
                index_t tNumberOfDOFs = mDOFs.size();

                // collect IDs of DOFs
                Vector< id_t > tIDs( tNumberOfDOFs );

                for ( index_t k = 0; k < tNumberOfDOFs; ++k )
                {
                    tIDs( k ) = mDOFs( k )->id();
                }

                // send IDs to master
                send( mParent->master(), tIDs );
            }
            else
            {
                // get number of procs
                uint tNumberOfProcs = mParent->comm_table().length() ;

                // initialize container
                mDofIDs.set_size( tNumberOfProcs, {} );

                // collect IDs from all other procs
                receive( mParent->comm_table(), mDofIDs );

                // get ref to my IDs
                this->collect_dof_ids_on_master( mDofIDs( 0 ) );

                // now we can create the dof table that contains the dof indices
                // per proc
                mDofIndices.set_size( tNumberOfProcs, Vector< index_t >() );

                // loop over all procs
                for ( uint p = 0; p < tNumberOfProcs; ++p )
                {
                    // get ID container
                    Vector< id_t > & tIDs = mDofIDs( p );

                    // get index list from dof table
                    Vector< index_t > & tIndices = mDofIndices( p );

                    // get number of DOFs for thos proc
                    index_t tNumberOfDOFs = tIDs.length();

                    // allocate memory for indices
                    tIndices.set_size( tNumberOfDOFs );

                    // populate Indices using map
                    for ( index_t k = 0; k < tNumberOfDOFs; ++k )
                    {
                        tIndices( k ) = mDofMap( tIDs( k ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::collect_dof_ids_on_master( Vector< id_t > & aIDs )
        {
            // the master knows all DOFs that exist
            index_t tNumberOfDOFs = mDOFs.size();
            aIDs.set_size( tNumberOfDOFs );
            for ( index_t k = 0; k < tNumberOfDOFs; ++k )
            {
                aIDs( k ) = mDOFs( k )->id();
            }
        }

//------------------------------------------------------------------------------

        void
        Field::create_sidesets()
        {
            // we have to treat the master proc differently, since it knows
            // the full mesh, but we need only the sidesets that connect to
            // owned elements. For other procs, this work has already be done.
            if ( mMyRank == mParent->master() )
            {
                index_t tSideSetCount = 0;

                Cell< fem::SideSet * > tSideSets(
                        mMesh->number_of_sidesets(),
                        nullptr );

                // loop over all sidesets
                Cell< mesh::SideSet * > & tMeshSideSets = mMesh->sidesets();

                for ( uint s = 0; s < mMesh->number_of_sidesets(); ++s )
                {
                    // get sideset
                    mesh::SideSet * tSideSet = tMeshSideSets( s );

                    // collect facets that connect to master
                    tSideSet->unflag_all_facets();

                    // counter
                    index_t tCount = 0;
                    for ( mesh::Facet * tFacet : tSideSet->facets() )
                    {
                        if ( tFacet->master()->owner() == mMyRank )
                        {
                            tFacet->flag();
                            ++tCount;
                        }
                        else
                        {
                            if ( tFacet->has_slave() )
                            {
                                if ( tFacet->slave()->owner() == mMyRank )
                                {
                                    tFacet->flag();
                                    ++tCount;
                                }
                            }
                        }
                    }

                    if ( tCount > 0 )
                    {
                        // allocate temporary container for facets
                        Cell< mesh::Facet * > tFacets( tCount, nullptr );

                        // reset counter
                        tCount = 0;

                        // collect facets
                        for ( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            if ( tFacet->is_flagged() )
                            {
                                tFacets( tCount++ ) = tFacet;
                            }
                        }

                        // add id to map
                        mSideSetMap[ tSideSet->id() ] = tSideSetCount;

                        // create the sideset
                        tSideSets( tSideSetCount++ ) = new fem::SideSet(
                                this,
                                tSideSet->id(),
                                tFacets );
                    }
                }

                // copy pointers into member cell
                mSideSets.set_size( tSideSetCount, nullptr );

                for ( index_t s = 0; s < tSideSetCount; ++s )
                {
                    mSideSets( s ) = tSideSets( s );
                }
            }
            else
            {
                // for all non-master procs, the facets are already
                // connected to owned elements only
                // this is why the following lines are a bit simpler

                index_t tNumberOfSideSets = mMesh->number_of_sidesets();

                // allocate container
                mSideSets.set_size( tNumberOfSideSets,
                                    nullptr );
                // loop over all sidesets
                Cell< mesh::SideSet * > & tMeshSideSets = mMesh->sidesets();

                // populate container
                for ( index_t s = 0; s < tNumberOfSideSets; ++s )
                {
                    // get sideset
                    mesh::SideSet * tSideSet = tMeshSideSets( s );

                    // create new object
                    mSideSets( s ) = new fem::SideSet(
                            this,
                            tSideSet->id(),
                            tSideSet->facets());

                    // create entry into map
                    mSideSetMap[ tSideSet->id() ] = s;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::create_bearings()
        {
            // get vertices from mesh
            Cell< mesh::Element * > & tVertices = mMesh->vertices();

            index_t tNumberOfVertices = tVertices.size();

            // create bearings
            mBearings.set_size( tNumberOfVertices, nullptr );

            for( index_t k=0; k<tNumberOfVertices; ++k )
            {
                mBearings( k ) = new Bearing(
                        this,
                        tVertices( k )->id(),
                        tVertices( k )->node( 0 ) );

                mBearingMap[ tVertices( k )->id() ] = k;
            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_dof_indices()
        {

            this->flag_used_dofs();

            if (  mMyRank != mParent->master() )
            {
                // get indices
                Vector< index_t > tIndices ;
                receive( mParent->master(), tIndices );

                BELFEM_ASSERT( tIndices.length() == mDOFs.size(),
                              "Length of indices does not match. ( is %lu, expect %lu )",
                              ( long unsigned int ) tIndices.length(),
                              ( long unsigned int ) mDOFs.size() );

                // reset local counters
                mMyNumberOfFreeDofs = 0;
                mMyNumberOfFixedDofs = 0;

                // write indices into dofs
                index_t tCount = 0 ;
                for( Dof * tDof : mDOFs )
                {
                    tDof->set_index( tIndices( tCount++ ) );

                    if( tDof->index() < gNoIndex )
                    {
                        tDof->flag() ;
                        if ( tDof->is_fixed() )
                        {
                            tDof->set_my_index( mMyNumberOfFixedDofs++ );
                        }
                        else
                        {
                            tDof->set_my_index( mMyNumberOfFreeDofs++ ) ;
                        }
                    }
                    else
                    {
                        tDof->unflag() ;
                        tDof->set_my_index( gNoIndex );
                    }
                }

                // send counters to master
                send( mParent->master(), mMyNumberOfFreeDofs );
                send( mParent->master(), mMyNumberOfFixedDofs );

                // receive global counters
                Vector< index_t > tCounters( 2 );
                receive( mParent->master(), tCounters );
                mNumberOfFreeDofs = tCounters( 0 );
                mNumberOfFixedDofs = tCounters( 1 );
            }
            else
            {
                // step 1: count fixed and free dofs
                mNumberOfFreeDofs = 0;
                mNumberOfFixedDofs = 0;

                // loop over all DOFs, note that on master,
                // myindex and index are the same
                for ( Dof * tDOF : mDOFs )
                {
                    if ( tDOF->is_flagged() )
                    {
                        if ( tDOF->is_fixed() )
                        {
                            tDOF->set_my_index( mNumberOfFixedDofs );
                            tDOF->set_index( mNumberOfFixedDofs++ );
                        }
                        else
                        {
                            tDOF->set_my_index( mNumberOfFreeDofs );
                            tDOF->set_index( mNumberOfFreeDofs++ );
                        }
                    }
                    else
                    {
                        tDOF->set_index( gNoIndex );
                        tDOF->set_my_index( gNoIndex );
                    }
                }

                // count number of fixed and free dofs
                mMyNumberOfFreeDofs = mNumberOfFreeDofs;
                mMyNumberOfFixedDofs = mNumberOfFixedDofs;

                // get size of communication table
                uint tNumProcs = mParent->comm_table().length() ;

                // this container holds the global indices for all procs
                Cell< Vector< index_t > > tAllIndices( tNumProcs, {} );

                // populate indices
                for( uint p=1; p<tNumProcs; ++p )
                {
                    // get id vector
                    const Vector< id_t > & tIDs = mDofIDs( p );

                    // get index vector
                    Vector< index_t >    & tIndices = tAllIndices( p );

                    // allocate vector
                    tIndices.set_size( tIDs.length() );

                    // initialize counter
                    index_t tCount = 0 ;

                    // loop over all IDs
                    for( id_t tID : tIDs )
                    {
                        // get index
                        tIndices( tCount++ ) = this->dof( tID )->index() ;
                    }
                }

                // send indices to procs
                send( mParent->comm_table(), tAllIndices );

                // get dof counters from clients
                receive( mParent->comm_table(), mNumberOfFreeDofsPerProc, mMyNumberOfFreeDofs );
                receive( mParent->comm_table(), mNumberOfFixedDofsPerProc, mMyNumberOfFixedDofs );

                // send global counters to other procs
                Vector< index_t > tCounters( 2 );
                tCounters( 0 ) = mNumberOfFreeDofs ;
                tCounters( 1 ) = mNumberOfFixedDofs ;
                send_same( mParent->comm_table(), tCounters );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::create_alpha_fields()
        {
            // check if an alpha BC (special BC for thermal convection) exists on any proc:

            uint tLocalFlag = 0 ;
            uint tGlobalFlag = 0 ;

            for( SideSet * tSideSet : mSideSets )
            {
                for( uint k=0; k<mNumberOfDOFsPerNode; ++k )
                {
                    if( tSideSet->bc_type( k ) == BoundaryConditionImposing::Alpha )
                    {
                        tLocalFlag = 1 ;
                        break ;
                    }
                }
                if( tLocalFlag > 0 )
                {
                    break ;
                }
            }

            if( mMyRank == mParent->master() )
            {
                Vector< uint > tFlags( mParent->comm_table().length(), 0 );
                receive( mParent->comm_table(), tFlags ) ;

                if( max( tFlags ) > 0 || tLocalFlag > 0 )
                {
                    tGlobalFlag = 1 ;
                }

                tFlags.fill( tGlobalFlag );
                send( mParent->comm_table(), tFlags );

            }
            else
            {
                send( mParent->master(), tLocalFlag );
                receive( mParent->master(), tGlobalFlag );
            }


            // fixme: why does broadcast not work here? We use the last send<->receive block as workaround
            // broadcast result to other procs
            // broadcast( mParent->master(), tGlobalFlag );

            if( tGlobalFlag > 0 )
            {
                // add fields to IWG
                mIWG->add_fields( { "alpha", "Tinf" } );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::synchronize_dirichlet_bcs()
        {
            // get comm table
            const Vector< proc_t > & tComm = mParent->comm_table();

            uint tNumProcs = tComm.length();

            if( comm_size() > 1 )
            {
                if ( mMyRank == mParent->master() )
                {
                    // cell with dof indices of fixed dofs
                    Cell< Vector< id_t > > tAllIDs( tNumProcs, {} );

                    // cell with dof values
                    Cell< Vector< real > > tAllValues( tNumProcs, {} );

                    receive( tComm, tAllIDs );
                    receive( tComm, tAllValues );

                    for ( uint p = 1; p < tNumProcs; ++p )
                    {

                        Vector< id_t > & tIDs = tAllIDs( p );
                        Vector< real > & tValues = tAllValues( p );

                        index_t tNumDofs = tIDs.length();

                        for ( index_t k = 0; k < tNumDofs; ++k )
                        {
                            this->dof( tIDs( k ) )->fix( tValues( k ) );
                        }
                    }

                    // now, we make sure that the data is consistent
                    for ( uint p = 1; p < tNumProcs; ++p )
                    {
                        // get ID vector
                        Vector< id_t > & tDofIDs = mDofIDs( p );

                        // count fixed dofs
                        index_t tCount = 0 ;
                        for( id_t tID : tDofIDs )
                        {
                            if( this->dof( tID )->is_fixed() )
                            {
                                ++tCount ;
                            }
                        }

                        // allocate containers
                        Vector< id_t > & tIDs = tAllIDs( p );
                        Vector< real > & tValues = tAllValues( p );

                        tIDs.set_size( tCount );
                        tValues.set_size( tCount );

                        // reset counter
                        tCount = 0 ;

                        // get values and IDs of fixed dofs
                        for( id_t tID : tDofIDs )
                        {
                            // get dof
                            Dof * tDof = this->dof( tID );

                            if( tDof->is_fixed() )
                            {
                                tIDs( tCount ) = tDof->id() ;
                                tValues( tCount++ ) = tDof->value() ;
                            }
                        }


                    }
                    comm_barrier() ;

                    // send containers to other procs
                    send( mParent->comm_table(), tAllIDs ) ;
                    send( mParent->comm_table(), tAllValues );
                }
                else
                {
                    // count fixed dofs
                    index_t tCount = 0;

                    for ( Dof * tDOF : mDOFs )
                    {
                        if ( tDOF->is_fixed() )
                        {
                            ++tCount;
                        }
                    }

                    // container for IDs
                    Vector< id_t > tIDs( tCount );

                    // container for values
                    Vector< real > tValues( tCount );

                    // reset counter
                    tCount = 0;

                    for ( Dof * tDOF : mDOFs )
                    {
                        if ( tDOF->is_fixed() )
                        {
                            tIDs( tCount )    = tDOF->id();
                            tValues( tCount ) = tDOF->value();
                            ++tCount;
                        }
                    }

                    // send fixed IDs and values
                    send( mParent->master(), tIDs );
                    send( mParent->master(), tValues );
                    comm_barrier() ;

                    // receive confirmation from master
                    receive( mParent->master(), tIDs );
                    receive( mParent->master(), tValues );

                    // reset counter
                    tCount = 0 ;

                    // this->free_all_dofs() ;

                    // loop over all IDs
                    for( id_t tID : tIDs )
                    {
                        // get dof
                        Dof * tDof = this->dof( tID );
                        tDof->fix( tValues( tCount ) );
                        ++tCount;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::distribute_alpha_bcs()
        {
            // check that fields exist
            if ( mMesh->field_exists("alpha") && mMesh->field_exists("Tinf") )
            {
                Vector< real > & tAlpha = mMesh->field_data( "alpha" );
                Vector< real > & tTinf  = mMesh->field_data( "Tinf" );

                if( mMyRank == mParent->master() )
                {
                    Cell< Vector< id_t > > tAllNodes( mParent->number_of_procs(), {} );

                    // get the node IDs from the parent
                    receive( mParent->comm_table(), tAllNodes );

                    Cell< Vector< real > > tAllAlpha( mParent->number_of_procs(), {} );
                    Cell< Vector< real > > tAllTinf( mParent->number_of_procs(), {} );

                    // loop over all procs
                    for( proc_t p=0;  mParent->number_of_procs(); ++p )
                    {
                        if( p != mMyRank )
                        {
                            // get the node list
                            Vector< id_t > & tNodeIDs = tAllNodes( p );

                            // get the alpha list
                            Vector< real > & tProcAlpha = tAllAlpha( p ) ;

                            // get the Tinf list
                            Vector< real > & tProcTinf = tAllTinf( p ) ;

                            index_t tNumNodes = tNodeIDs.length() ;

                            // allocate the vectors
                            tProcAlpha.set_size( tNumNodes );
                            tProcTinf.set_size( tNumNodes );

                            // populate data
                            for( index_t k=0; k<tNumNodes; ++k )
                            {
                                // get the node
                                mesh::Node * tNode = mMesh->node( tNodeIDs( k ) ) ;

                                // copy the alpha value
                                tProcAlpha( k ) = tAlpha( tNode->index() );

                                // copy the temperature value
                                tProcTinf( k )  = tTinf( tNode->index() );
                            }
                        }
                    }
                    // send the data to the other procs
                    send( mParent->comm_table(), tAllAlpha );
                    send( mParent->comm_table(), tAllTinf );
                }
                else
                {
                    index_t tCount = 0;

                    // loop over all sidesets
                    for ( fem::SideSet * tSideSet : mSideSets )
                    {
                        if ( tSideSet->bc_type( 0 ) == BoundaryConditionImposing::Alpha )
                        {
                            // add number of nodes to counter
                            tCount += tSideSet->nodes().size();
                        }
                    }

                    // allocate request container
                    Vector< id_t > tNodeIDs( tCount );

                    // collect node indices
                    tCount = 0;
                    for ( fem::SideSet * tSideSet : mSideSets )
                    {
                        if ( tSideSet->bc_type( 0 ) == BoundaryConditionImposing::Alpha )
                        {
                            // get the nodes
                            Cell< mesh::Node * > & tNodes = tSideSet->nodes();

                            for ( mesh::Node * tNode : tNodes )
                            {
                                tNodeIDs( tCount++ ) = tNode->id();
                            }
                        }
                    }

                    // send the request list to the core
                    send( mParent->master(), tNodeIDs );

                    Vector< real > tProcAlpha ;
                    Vector< real > tProcTinf ;

                    // get the field data
                    receive( mParent->master(), tProcAlpha );
                    receive( mParent->master(), tProcTinf );

                    // write data onto the mesh
                    for( index_t k=0; k<tCount; ++k )
                    {
                        // get the node
                        mesh::Node * tNode = mMesh->node( tNodeIDs( tCount ) );

                        tAlpha( tNode->index() ) = tProcAlpha( k ) ;
                        tTinf( tNode->index() )  = tProcTinf( k ) ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_rhs_vector()
        {
            // set size of dof container and reset it
            mRhsVector.set_size( mMyNumberOfFreeDofs, 0.0 );

            // loop over all blocks and elements
            for ( fem::Block * tBlock : mBlocks )
            {
                // check if block has an RHS
                if( tBlock->has_rhs() )
                {
                    // get number of dofs per element
                    uint tN = this->enforce_linear_interpolation() ?
                              mesh::number_of_corner_nodes( tBlock->element_type())
                            : mesh::number_of_nodes( tBlock->element_type() );
                    tN *= mNumberOfDOFsPerNode;
                    tN += mNumberOfDOFsPerEdge * mesh::number_of_edges( tBlock->element_type() );

                    // allocate element RHS
                    Vector< real > tRHS( tN );

                    // link IWG to block
                    mIWG->link_to_group( tBlock );

                    // get elements on block
                    Cell< Element * > & tElements = tBlock->elements();

                    // loop over all elements
                    for ( Element * tElement : tElements )
                    {
                        mIWG->compute_rhs( tElement, tRHS );

                        // add element jacobian to master
                        for ( uint i = 0; i < tN; ++i )
                        {
                            Dof * tRow = tElement->dof( i );

                            if ( !tRow->is_fixed() )
                            {
                                mRhsVector( tRow->my_index() ) += tRHS( i );
                            }
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::detect_wetted_sidesets()
        {
            // get number of wetted sidesets on this proc
            if( mIWG->wetted_sidesets().length() == 0 )
            {
                index_t tCount = 0 ;
                for( SideSet * tSideSet : mSideSets )
                {
                    for ( uint k = 0; k < mNumberOfDOFsPerNode; ++k )
                    {
                        if ( tSideSet->bc_type( k ) == BoundaryConditionImposing::Neumann ||
                             tSideSet->bc_type( k ) == BoundaryConditionImposing::Alpha )
                        {
                            ++tCount;
                            break;
                        }
                    }
                }

                if( tCount > 0 )
                {

                    Vector< id_t > tWettedSidesets( tCount );
                    tCount = 0 ;
                    for( SideSet * tSideSet : mSideSets )
                    {
                        for ( uint k = 0; k < mNumberOfDOFsPerNode; ++k )
                        {
                            if ( tSideSet->bc_type( k ) == BoundaryConditionImposing::Neumann ||
                                 tSideSet->bc_type( k ) == BoundaryConditionImposing::Alpha )
                            {
                                tWettedSidesets( tCount++ ) = tSideSet->id() ;
                                break;
                            }
                        }
                    }

                    mIWG->set_wetted_sidesets( tWettedSidesets );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::count_wetted_nodes()
        {
            // unflag all nodes
            mMesh->unflag_all_nodes() ;

            // loop over all wetted sidesets
            for( id_t tID : mIWG->wetted_sidesets() )
            {
                // select sideset
                if( this->sideset_exists( tID ) )
                {
                    mMesh->sideset( tID )->flag_all_nodes() ;
                }
            }

            // reset node counter
            mNumberOfConvectionNodes = 0 ;

            // count nodes
            for( mesh::Node * tNode : mMesh->nodes() )
            {
                if( tNode->is_flagged() )
                {
                    ++mNumberOfConvectionNodes ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_rhs_matrix()
        {
            // number of cols in rhs matrix
            uint tNumCols = mIWG->num_rhs_cols();

            // set size of dof container and reset it
            mRhsMatrix.set_size( mMyNumberOfFreeDofs, tNumCols, 0.0 );


            // loop over all blocks and elements
            for ( fem::Block * tBlock : mBlocks )
            {
                // get number of dofs per element
                uint tN = this->enforce_linear_interpolation() ?
                        mNumberOfDOFsPerNode * mesh::number_of_corner_nodes( tBlock->element_type() )
                      : mNumberOfDOFsPerNode * mesh::number_of_nodes( tBlock->element_type() );

                // allocate element RHS
                Matrix< real > tRHS( tN, tNumCols );

                // link IWG to block
                mIWG->link_to_group( tBlock );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    mIWG->compute_rhs( tElement, tRHS );

                    // add element jacobian to master
                    for ( uint i = 0; i < tN; ++i )
                    {
                        Dof * tRow = tElement->dof( i );

                        if ( ! tRow->is_fixed()  )
                        {
                            for ( uint j = 0; j < tNumCols; ++j )
                            {
                                mRhsMatrix( tRow->my_index(), j ) += tRHS( i, j );
                            }
                        }
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        real
        Field::residual( const uint aIteration )
        {
            if( mMyRank == mParent->master() )
            {
                BELFEM_ASSERT( mIWG->num_rhs_cols() == 1, "Iwg must have only 1 col for calling Field::residual()" );

                // compute value
                real aResidual = norm( mLhsVector );

                // distribute data
                Vector< real > tResidual( mParent->comm_table().length(), aResidual );

                // catch error
                if( ( aResidual == 0  && aIteration == 0 ) || aResidual > 1E12 || std::isnan( aResidual ) )
                {

                    BELFEM_ERROR( false,
                                 "YOUR SPACECRAFT BLEW UP BECAUSE SOMEONE DEVIDED BY ZERO. \nIN REALITY THIS CAN NOT HAPPEN; BUT NUMERICS ARE DIFFERENT;\nYOU CAN TRY THE FOLLIOWING THINGS: \n    * MAKE A BETTER CONDITIONED MESH\n    * USE A DIFFERENT PRECONDITIONER OR A DIRECT SOLVER\n    * DECREASE THE TIMESTEP\n    * DECREASE THE RELAXATION PARAMETER\n" );
                }

                send( mParent->comm_table(), tResidual );

                return aResidual;
            }
            else
            {
                // get value from master
                real aResidual = BELFEM_REAL_MAX ;

                receive( mParent->master(), aResidual );

                return  aResidual;
            }
        }

//------------------------------------------------------------------------------

        void
        Field::synchronize_rhs_matrix()
        {
            // number of cols in rhs matrix
            uint tNumCols = mIWG->num_rhs_cols();

            if ( mMyRank == mParent->master() )
            {
                Cell< Matrix< real > > tAllRHS;
                receive( mParent->comm_table(), tAllRHS );

                // assemble system
                for ( uint p = 1; p < mParent->comm_table().length(); ++p )
                {
                    // get dof table for this proc
                    Vector< index_t > & tDOFs = mDofIndices( p );
                    Matrix< real > & tRHS = tAllRHS( p );
                    // get number of dofs
                    index_t tNumDOFs = tDOFs.length();

                    // loop over all dofs
                    for ( index_t j = 0; j < tNumCols; ++j )
                    {
                        for ( index_t i = 0; i < tNumDOFs; ++i )
                        {
                            // get dof
                            Dof * tDOF = mDOFs( tDOFs( i ));

                            if ( !tDOF->is_fixed() )
                            {
                                mRhsMatrix( tDOF->index(), j ) += tRHS( i, j );
                            }
                        }
                    }
                }
            }
            else
            {
                // send vector to master
                send( mParent->master(), mRhsMatrix );
            }
        }
//------------------------------------------------------------------------------

        void
        Field::initialize_jacobian()
        {
            // initialize timer
            Timer tTimer;
            if ( mMyRank == mParent->master() )
            {
                message( 4, " Initialize Jacobian ... \n" );
            }

            this->synchronize_dirichlet_bcs();
            this->compute_dof_indices();

            if ( mMyNumberOfFreeDofs > 0 )
            {
                // allocate space for node graph
                Cell< graph::Vertex * > tGraph( mMyNumberOfFreeDofs, nullptr );

                // counter for dofs
                index_t tCount = 0;

                // populate graph
                for ( Dof * tDOF : mDOFs )
                {
                    // test if DOF is free and flagged
                    if ( ! tDOF->is_fixed() && tDOF->is_flagged() )
                    {
                        tGraph( tCount++ ) = tDOF;
                    }
                }

                // first  loop initializes Jacobian matrix
                // second loop initializes Dirichlet matrix, if bearings exist
                // third  loop initializes Bearing Matrix ( in the future )
                uint tNumLoops = ( mMyNumberOfFixedDofs > 0 ) ? 2 : 1;

                // fixed flag is false in first loop and
                // true in second loop
                bool tFlag = false;

                // allocate container for neighbors
                // if you have more than 256 dofs per elemement you have a problem
                Cell< Dof * > tWork( 256, nullptr );

                // start the loop
                for ( uint l = 0; l < tNumLoops; ++l )
                {
                    if ( mNumberOfDOFsPerEdge == 0 )
                    {
                        this->compute_node_based_adjacency( tGraph, tWork, tFlag );
                    }
                    else if ( mNumberOfDOFsPerNode == 0 )
                    {
                        this->compute_edge_based_adjacency( tGraph, tWork, tFlag );
                    }
                    else
                    {
                        this->compute_node_and_edge_based_adjacency( tGraph, tWork, tFlag );
                    }
                    
                    if ( l == 0 )
                    {
                        // with the graph created, we can initialize the Jacobian
                        mJacobian = new SpMatrix( tGraph, mSolver->type() == SolverType::PETSC ?
                        SpMatrixType::CSR : SpMatrixType::CSC,
                                                  mNumberOfFreeDofs, mNumberOfFreeDofs );

                        if ( mMyRank == mParent->master())
                        {
                            message( 4, "    ... time for initializing Jacobian          : %u ms\n",
                                     ( unsigned int ) tTimer.stop());
                        }
                        tTimer.reset();

                        // flip flag for next loop
                        tFlag = true;

                    }
                    else if ( l == 1 )
                    {
                        // for the second loop, we create the Dirichlet Matrix
                        mDirichletMatrix = new SpMatrix( tGraph, SpMatrixType::CSR,
                                                         mNumberOfFreeDofs, mNumberOfFixedDofs );


                        if ( mMyRank == mParent->master() )
                        {
                            message( 4, "    ... time for initializing Dirichlet Matrix  : %u ms\n",
                                     ( unsigned int ) tTimer.stop());
                        }

                        tTimer.reset();
                    }
                    // else if ( l == 2 ) will create the bearing matrix in the future
                } // end matrix loop

                // reset graph containers
                //for ( graph::Vertex * tDof : tGraph )
                //{
                //    tDof->reset_vertex_container();
                //}

            }

            // synchronize indices
            this->create_assembly_table();

            // allocate vectors
            if( mIWG->num_rhs_cols() <= 1 )
            {
                mRhsVector.set_size( mNumberOfFreeDofs, 0.0 );
            }
            else
            {
                mRhsMatrix.set_size( mMyNumberOfFreeDofs,
                        mIWG->num_rhs_cols(), 0.0 );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::reset_jacobian()
        {
            BELFEM_ASSERT( mJacobian != nullptr,
                "Jacobian Matrix was not initialized");

            mJacobian->fill( 0.0 );

            if ( mDirichletMatrix != nullptr )
            {
                mDirichletMatrix->fill( 0.0 );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_jacobian( const bool aReset )
        {
            Timer tTimer;

            if ( aReset )
            {
              this->reset_jacobian() ;
            }

            // loop over all blocks
            for ( fem::Block * tBlock : mBlocks )
            {
                // get number of dofs per element
                uint tN = this->enforce_linear_interpolation() ?
                          mesh::number_of_corner_nodes( tBlock->element_type() )
                        : mesh::number_of_nodes( tBlock->element_type() ) ;
                tN *= mNumberOfDOFsPerNode ;
                tN += mNumberOfDOFsPerEdge * mesh::number_of_edges( tBlock->element_type() );

                // allocate element Jacobian
                Matrix< real > tJ( tN, tN );

                // link IWG to block
                mIWG->link_to_group( tBlock );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    mIWG->compute_jacobian( tElement, tJ );
                    this->assemble_jacobian( tElement, tJ );
                }
            }

            this->synchronize_jacobian();

            if ( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for computing Jacobian             : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_jacobian_and_rhs()
        {
            Timer tTimer;

            BELFEM_ERROR( mIWG->num_rhs_cols() == 1, "Num cols must be 1 for compute_jacobian_and_rhs()" );

            this->reset_jacobian() ;

            SpMatrix & J       =  *mJacobian;
            SpMatrix & D       =  *mDirichletMatrix;

            // residual
            Vector< real > & B = mRhsVector;
            mRhsVector.fill( 0.0 );

            // unflag all nodes on this mesh
            mMesh->unflag_all_nodes();

            // loop over all blocks
            for ( fem::Block * tBlock : mBlocks )
            {
                // get number of dofs per element
                uint tN = this->enforce_linear_interpolation() ?
                          mesh::number_of_corner_nodes( tBlock->element_type() )
                        : mesh::number_of_nodes( tBlock->element_type() ) ;
                tN*= mNumberOfDOFsPerNode ;

                // allocate element Jacobian
                Matrix< real > tJ( tN, tN );

                // allocate element RHS
                Vector< real > tB( tN );

                // link IWG to block
                mIWG->link_to_group( tBlock );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    mIWG->compute_jacobian_and_rhs( tElement, tJ, tB );

                    this->add_to_system_matrices( tElement, tJ, tB, tN, J, D, B );

                }
            } // end loop over all blocks

            if( mIWG->has_alpha() )
            {
                // loop over all sidesets
                for( SideSet * tSideSet : mSideSets )
                {
                    // check if this is an alpha BC
                    if ( tSideSet->bc_type( 0 ) == BoundaryConditionImposing::Alpha )
                    {
                        // link this IWG with the group
                        mIWG->link_to_group( tSideSet ) ;

                        // get number of dofs per element
                        uint tN = this->enforce_linear_interpolation() ?
                                   mesh::number_of_corner_nodes( tSideSet->element_type() )
                                :  mesh::number_of_nodes( tSideSet->element_type() );
                        tN *= mNumberOfDOFsPerNode ;

                        // allocate element Jacobian
                        Matrix< real > tJ( tN, tN );

                        // allocate element RHS
                        Vector< real > tB( tN );

                        // get elements on Sideset
                        Cell< Element * > & tElements = tSideSet->elements();

                        // loop over all elements
                        for ( Element * tElement : tElements )
                        {
                            mIWG->compute_alpha_boundary_condition( tElement, tJ, tB );
                            this->add_to_system_matrices( tElement, tJ, tB, tN, J, D, B );
                        }

                    }
                }
            }

            this->synchronize_jacobian();
            this->collect_rhs_vector();
            if ( mMyRank == mParent->master() )
            {
                // allocate vector
                if( mFieldValues.length() != mNumberOfFreeDofs )
                {
                    mFieldValues.set_size( mNumberOfFreeDofs );
                }

                // collect values for free dofs
                for ( Dof * tDof: mDOFs )
                {
                    if ( ! tDof->is_fixed() )
                    {
                        mFieldValues( tDof->index() ) = tDof->value();
                    }
                }

                message( 4, "    ... time for computing Jacobian and residual : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }
        }

//------------------------------------------------------------------------------

        void
        Field::synchronize_jacobian()
        {
            if ( mMyRank == mParent->master() )
            {

                // get number of procs
                uint tNumberOfProcs = mParent->comm_table().length();

                Cell< Vector< real > > tJacobian;

                receive( mParent->comm_table(), tJacobian );

                // assemble jacobian
                for ( uint p = 1; p < tNumberOfProcs; ++p )
                {
                    // get data
                    Vector< real > & tData = tJacobian( p );

                    // get indices
                    Vector< index_t > & tIndices = mJacobianTable( p );

                    // get number of nonzeros
                    index_t tNNZ = tData.length();

                    BELFEM_ASSERT( tIndices.length() == tData.length(),
                                  "Jacobian Matrix from proc %u has wrong number of nonzeros ( is %lu, expect %lu )",
                                  ( unsigned int ) mParent->comm_table( p ),
                                  ( unsigned int ) tData.length(),
                                  ( unsigned int ) tIndices.length() );

                    // loop over all entries
                    for ( index_t k = 0; k < tNNZ; ++k )
                    {
                        // add entries
                        mJacobian->data( tIndices( k ) ) += tData( k );
                    }
                }

                // tidy up memory a bit
                tJacobian.clear();

                Cell< Vector< real > > tDirichlet( tNumberOfProcs, {} );

                // assemble Dirichlet matrix
                if ( mNumberOfFixedDofs > 0 )
                {
                    // collect data from clients
                    for ( uint p = 1; p < tNumberOfProcs; ++p )
                    {
                        if ( mNumberOfFixedDofsPerProc( p ) > 0 )
                        {
                            // get data
                            Vector< real > & tData = tDirichlet( p );

                            receive( mParent->comm_table( p ), tData );
                        }
                    }
                    // assemble matrix
                    for ( uint p = 1; p < tNumberOfProcs; ++p )
                    {
                        // get data
                        Vector< real > & tData = tDirichlet( p );

                        // get indices
                        Vector< index_t > & tIndices = mDirichletTable( p );

                        // get number of nonzeros
                        index_t tNNZ = tData.length() ;

                        BELFEM_ASSERT( tIndices.length() == tData.length(),
                            "Dirichlet Matrix from proc %u has wrong number of nonzeros ( is %lu, expect %lu )",
                                      ( unsigned int ) mParent->comm_table( p ),
                                      ( unsigned int ) tData.length(),
                                      ( unsigned int ) tIndices.length() );

                        // loop over all entries
                        for ( index_t k = 0; k < tNNZ ; ++k )
                        {
                            // add entries
                            mDirichletMatrix->data( tIndices( k ) ) += tData( k );
                        }
                    }
                }
            }
            else
            {
                // send my data to master
                send( mParent->master(),
                      mJacobian->number_of_nonzeros(),
                      mJacobian->data());

                             if ( mMyNumberOfFixedDofs > 0 )
                {
                    send( mParent->master(),
                         mDirichletMatrix->number_of_nonzeros(),
                         mDirichletMatrix->data() );
                }

            }
        }

//-----------------------------------------------------------------------------
        void
        Field::collect_rhs_vector()
        {
            this->collect_vector( mRhsVector );
        }

//------------------------------------------------------------------------------

        void
        Field::collect_vector( Vector< real > & aVector )
        {

            if ( mMyRank == mParent->master() )
            {
                const Vector< proc_t > & tComm = mParent->comm_table() ;

                Cell< Vector< real > > tAllVectors;
                receive( tComm, tAllVectors );

                // assemble system
                for ( uint p = 1; p < tComm.length(); ++p )
                {
                    // get dof table for this proc
                    Vector< index_t > & tDOFs = mDofIndices( p );

                    Vector< real > & tVector = tAllVectors( p );


                    // sometimes, the vector may be of zero lentgh,
                    // eg, if a proc does not have a wetted surface
                    if ( tVector.length() > 0 )
                    {
                        // get number of dofs fixme: check this
                        index_t tNumDOFs = tDOFs.length();

                        index_t tCount = 0 ;

                        // loop over all dofs
                        for ( index_t i = 0; i < tNumDOFs; ++i )
                        {
                            // get dof
                            Dof * tDOF = mDOFs( tDOFs( i ) );

                            if ( !tDOF->is_fixed() )
                            {
                                aVector( tDOF->index() ) += tVector( tCount++ );
                            }
                        }
                    }
                }
            }
            else
            {
                // send vector to master
                send( mParent->master(), aVector );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::solve()
        {
            // check if a solver exists
            BELFEM_ERROR( mSolver != nullptr,
                "You must set a solver before solving the field" );

            // uint tNumberOfFields = mNumberOfDOFsPerNode + mNumberOfDOFsPerEdge ;

            // populate field data
            const Cell< string > & tFieldLabels = mIWG->all_fields() ;

            uint tNumberOfFields = tFieldLabels.size() ;

            Cell< mesh::Field * > tFields( tNumberOfFields, nullptr );

            for ( uint f = 0; f < tNumberOfFields; ++f )
            {
                tFields( f ) = mMesh->field( tFieldLabels( f ) );
            }

            if ( mMyRank == mParent->master() )
            {
                Timer tTimer;

                // right hand side
                Vector< real > tFixedValues( mNumberOfFixedDofs );

                index_t tCount = 0;

                // loop over all dofs
                for ( Dof * tDof: mDOFs )
                {
                    if ( tDof->is_fixed() )
                    {
                        tFixedValues( tCount++ ) = tDof->value();
                    }
                }

                // add loads over boundary
                if( mConvection.length() > 0 )
                {
                    mRhsVector += mConvection ;
                }

                // add volume loads
                if( mVolumeLoads.length() > 0 )
                {
                    mRhsVector += mVolumeLoads ;
                }

                if ( mNumberOfFixedDofs != 0 )
                {
                    BELFEM_ASSERT( mIWG->num_rhs_cols() == 1,
                        "Can only impose values of RHS is a vector, not a matrix!");

                    mDirichletMatrix->multiply( tFixedValues, mRhsVector, 1.0, 1.0 );
                }

                if(  mIWG->num_rhs_cols() == 1 )
                {
                    if( mIWG->mode() == IwgMode::Iterative )
                    {
                        BELFEM_ASSERT( mFieldValues.length() == mRhsVector.length(),
                            "Length of Field values and RHS vector do not match ( %lu vs. %lu, free dofs: %lu )",
                                      ( long unsigned int ) mFieldValues.length(),
                                      ( long unsigned int ) mRhsVector.length(),
                                      ( long unsigned int ) mNumberOfFreeDofs );

                        mJacobian->multiply( mFieldValues, mRhsVector, 1.0, -1.0 );
                    }

                    // solve the system
                    mSolver->solve( *mJacobian, mLhsVector, mRhsVector ) ;

                    if( mIWG->mode() == IwgMode::Iterative )
                    {
                        for ( Dof * tDof: mDOFs )
                        {
                            // update DOF values
                            if ( !tDof->is_fixed())
                            {
                                tDof->value() -= mIWG->omega() * mLhsVector( tDof->index());
                            }

                            // update value in field
                            tFields( tDof->type_id() )->value( tDof->dof_index_on_field() ) = tDof->value();
                        }
                    }
                    else
                    {
                        // write values into field
                        for ( Dof * tDof: mDOFs )
                        {
                            if ( tDof->is_fixed() )
                            {
                                tFields( tDof->type_id() )->value( tDof->dof_index_on_field() ) = tDof->value();
                            }
                            else
                            {
                                tFields( tDof->type_id() )->value( tDof->dof_index_on_field() )
                                        = mLhsVector( tDof->index() );
                            }
                        }
                    }
                }
                else // right hand side is matrix
                {
                    BELFEM_ERROR( mNumberOfDOFsPerEdge == 0,
                        "solving a right hand side matrix is currently not supported for edge elements");

                    mSolver->solve( *mJacobian, mLhsMatrix, mRhsMatrix );
                    uint k = 0;
                    for( index_t j=0; j<mIWG->num_rhs_cols(); ++j )
                    {
                        for ( index_t i = 0; i < mIWG->number_of_dofs_per_node(); ++i )
                        {
                            // get field
                            Vector< real > & tData = mMesh->field_data( tFieldLabels( k++ ) );

                            // loop over all nodes
                            for( mesh::Node *tNode : mMesh->nodes() )
                            {
                                // get dof
                                Dof * tDOF = this->dof( this->calculate_dof_id( tNode, i ) );

                                // write data of dof into field
                                tData( tDOF->dof_index_on_field() ) = mLhsMatrix( tDOF->index(), j );
                            }
                        }
                    }
                }

                message( 4, "    ... time for solving system of equations    : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }
            else if ( mSolver->type() == SolverType::MUMPS ||
                      mSolver->type() == SolverType::PETSC )
            {
                if ( mIWG->num_rhs_cols() == 1 )
                {
                    mSolver->solve( *mJacobian, mLhsVector, mRhsVector ) ;
                }
                else
                {
                    mSolver->solve( *mJacobian, mLhsMatrix, mRhsMatrix ) ;
                }
            }

            // make field available to other procs
            this->distribute_fields( mIWG->all_fields() );
        }

//------------------------------------------------------------------------------

        void
        Field::create_assembly_table()
        {
            // get master proc
            proc_t tMaster = mParent->master();

            if ( mMyRank == tMaster )
            {
                const Vector< proc_t > & tComm = mParent->comm_table();

                Cell< Vector< int > > tJacobianRows;
                Cell< Vector< int > > tJacobianCols;
                Cell< Vector< int > > tDirichletRows;
                Cell< Vector< int > > tDirichletCols;



                // get data
                receive( tComm, tJacobianRows );
                receive( tComm, tJacobianCols );

                uint tNumberOfProcs = tComm.length();

                // allocate memory
                mJacobianTable.set_size( tNumberOfProcs, {} );

                // create jacobian indices
                for ( uint p = 1; p < tNumberOfProcs; ++p )
                {
                    // get index
                    Vector< index_t > & tIndex = mJacobianTable( p );

                    // get rows
                    Vector< int > & tRows = tJacobianRows( p );

                    // get  cols
                    Vector< int > & tCols = tJacobianCols( p );

                    // get number of nonzeros in this matrix
                    index_t tNNZ = tRows.length();

                    tIndex.set_size( tNNZ );

                    // loop over all entries
                    for ( index_t k = 0; k < tNNZ; ++k )
                    {
                        // compute index
                        tIndex( k ) = mJacobian->index( tRows( k ), tCols( k ) );
                    }
                }

                receive( tComm, tDirichletRows );
                receive( tComm, tDirichletCols );

                if ( mDirichletMatrix != NULL )
                {
                    mDirichletTable.set_size( tNumberOfProcs, {} );

                    // create jacobian indices
                    for ( uint p = 1; p < tNumberOfProcs; ++p )
                    {
                        // get index
                        Vector< index_t > & tIndex = mDirichletTable( p );

                        // get rows
                        Vector< int > & tRows = tDirichletRows( p );

                        // get  cols
                        Vector< int > & tCols = tDirichletCols( p );

                        // get number of nonzeros in this matrix
                        index_t tNNZ = tRows.length();

                        if( tNNZ > 0 )
                        {
                            tIndex.set_size( tNNZ );

                            // loop over all entries
                            for ( index_t k = 0; k < tNNZ; ++k )
                            {
                                // compute index
                                tIndex( k ) = mDirichletMatrix->index( tRows( k ), tCols( k ) );
                            }
                        }
                    }
                }


            }
            else
            {
                // prepare data
                mJacobian->create_coo_indices();

                // send rows
                send( tMaster, mJacobian->number_of_nonzeros(), mJacobian->rows() );
                // send columns
                send( tMaster, mJacobian->number_of_nonzeros(), mJacobian->cols() );

                // tidy up memory
                //mJacobian->free_coo_indices();



                if ( mDirichletMatrix == NULL )
                {

                    index_t tZero = 0;

                    // send zero to master ( for rows )
                    send( tMaster, tZero );
                    // send zero to master ( for cols )
                    send( tMaster, tZero );
                }
                else
                {
                    // prepare data
                    mDirichletMatrix->create_coo_indices();

                    send( tMaster, mDirichletMatrix->number_of_nonzeros(), mDirichletMatrix->rows() );

                    send( tMaster, mDirichletMatrix->number_of_nonzeros(), mDirichletMatrix->cols() );

                    // tidy up memory
                    // mDirichletMatrix->free_coo_indices();
                }
            }


        }

//------------------------------------------------------------------------------

        void
        Field::free_all_dofs()
        {
            for ( Dof * tDOF : mDOFs )
            {
                tDOF->free();
            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_rhs()
        {
            Timer tTimer;

            if( mIWG->num_rhs_cols() == 1 )
            {
                this->compute_rhs_vector();
                this->collect_rhs_vector();
            }
            else if ( mIWG->num_rhs_cols() > 1 )
            {
                this->compute_rhs_matrix();
                this->synchronize_rhs_matrix();
            }


            if( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for computing right hand side      : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }

        }

//------------------------------------------------------------------------------

        void
        Field::distribute_fields( const Cell< string > & aFieldLabels )
        {
            // get number of fields
            uint tNumberOfFields = aFieldLabels.size();

            // get master proc
            proc_t tMaster = mParent->master();

            // check if data must be projected
            if( this->enforce_linear_interpolation() )
            {
                this->project_linear_field_to_higher_mesh( aFieldLabels );
            }

            if ( mMyRank == tMaster )
            {

                // get number of procs
                uint tNumberOfProcs = mParent->number_of_procs();

                // sanity check
                if( mParent->number_of_procs() > 1 )
                {
                    Vector< uint > tNumFieldsPerProc;
                    receive( mParent->comm_table(), tNumFieldsPerProc );
                    for ( uint k = 0; k < mParent->comm_table().length(); ++k )
                    {
                        BELFEM_ERROR( tNumberOfFields == tNumFieldsPerProc( k ) ||
                                     mParent->comm_table( k ) == mParent->master(),
                                     "Number of fields on proc %u does not match. Is %u but expect %u",
                                     ( unsigned int ) mParent->comm_table( k ),
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

                    // loop over all procs
                    for( uint p=1; p<tNumberOfProcs; ++p )
                    {
                        // check if this is a node or an edge field
                        const Vector< index_t > & tIndices
                            = tF->entity_type() == EntityType::EDGE ?
                                    mParent->edge_table( p ) : mParent->node_table( p );

                        // get number of nodes
                        index_t tNumberOfValues = tIndices.length();

                        // get values
                        Vector< real > & tValues = tData( p );

                        // set size for values
                        tValues.set_size( tNumberOfValues );

                        // populate data
                        for( index_t k=0; k<tNumberOfValues; ++k )
                        {
                            tValues( k ) = tField( tIndices( k ) );
                        }

                    }

                    // send data to othere procs
                    send( mParent->comm_table(), tData );
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

//------------------------------------------------------------------------------

        void
        Field::set_all_dofs( const real & aValue )
        {
            for( Dof * tDof : mDOFs )
            {
                tDof->value() = aValue;
            }
        }

//------------------------------------------------------------------------------


        void
        Field::create_convection_table()
        {
            mConvectionTable.clear();

            if( mMyRank == mParent->master() )
            {
                const uint tNumProcs = mParent->comm_table().length();

                // container with node IDs
                Cell< Vector< id_t > > tAllNodeIDs( tNumProcs, {} );

                // get IDs from other procs
                receive( mParent->comm_table(), tAllNodeIDs );

                // allocate convection table
                mConvectionTable.set_size( tNumProcs, {} );

                for( uint p=1; p<tNumProcs; ++p )
                {
                    // get vector with node IDs
                    const Vector< id_t > & tNodeIDs = tAllNodeIDs( p );

                    // get index vector
                    Vector< index_t > & tIndices = mConvectionTable( p );

                    // get number of nodes for this proc
                    index_t tNumNodes = tNodeIDs.length();

                    // set the size of the vector
                    tIndices.set_size( tNumNodes );

                    // loop over all nodes
                    for( index_t k=0; k<tNumNodes; ++k )
                    {
                        // save index from this node
                        tIndices( k ) = mMesh->node( tNodeIDs( k ) )->index();
                    }
                }


                // needed for compatibility
                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    mConvectionMap[ tNode->id() ] = tNode->index();
                }

            }
            else
            {
                // unflag all nodes on mesh
                mMesh->unflag_all_nodes();

                // get wetted sidesets
                const Vector< id_t > & tWettedSidesets = mIWG->wetted_sidesets();

                // check for linear flag
                if ( this->enforce_linear_interpolation() )
                {
                    // loop over all sidesets
                    for ( const id_t & tID : tWettedSidesets )
                    {
                        // get sideset
                        SideSet * tSideSet = this->sideset( tID );

                        // flag all nodes on this sideset
                        for ( Element * tElement : tSideSet->elements())
                        {
                            tElement->element()->flag_corner_nodes();
                        }
                    }
                }
                else
                {
                    // loop over all sidesets
                    for ( const id_t & tID : tWettedSidesets )
                    {
                        // get sideset
                        SideSet * tSideSet = this->sideset( tID );

                        // flag all nodes on this sideset
                        for ( Element * tElement : tSideSet->elements())
                        {
                            tElement->element()->flag_nodes();
                        }
                    }
                }

                // count flagged nodes on mesh
                index_t tCount = 0;
                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        mConvectionMap[ tNode->id() ] = tCount++;
                    }
                }

                // remember the size
                mNumberOfConvectionNodes = tCount;

                // allocate vector with node ids
                Vector< id_t > tIDs;

                if( tCount > 0 )
                {
                    tIDs.set_size( tCount );

                    // reset counter
                    tCount = 0;
                    for ( mesh::Node * tNode : mMesh->nodes())
                    {
                        if ( tNode->is_flagged() )
                        {
                            tIDs( tCount++ ) = tNode->id();
                        }
                    }
                }

                // send IDs to master
                send( mParent->master(), tIDs );

                // tidy up mesh
                mMesh->unflag_all_nodes();
            }
        }

//------------------------------------------------------------------------------

        void
        Field::reset_convection_vector()
        {
            // allocate the convection matrix if it has not been so far
            if( mConvection.length() == 0 )
            {
                if ( mMyRank == mParent->master() )
                {
                    mConvection.set_size(  mNumberOfFreeDofs );
                }
                else if ( mNumberOfConvectionNodes > 0 )
                {
                    mConvection.set_size( mMyNumberOfFreeDofs );
                }
            }

            // reset values
            mConvection.fill( 0.0 );
        }

//------------------------------------------------------------------------------
        // special function to compute convective term
        void
        Field::compute_surface_loads()
        {
            this->reset_convection_vector() ;

            // get wetted sidesets
            const Vector< id_t > & tWettedSidesets = mIWG->wetted_sidesets();

            // loop over all wetted sideset IDs
            for ( const id_t & tID : tWettedSidesets )
            {
                // get the sideset
                SideSet * tSideSet = this->sideset( tID );

                // connect IWG with sideset
                mIWG->link_to_group( tSideSet );

                // get facet container
                Cell< Element * > & tFacets = tSideSet->elements() ;

                // go to next sideset, if this one is empty
                if( tFacets.size() == 0 )
                {
                    continue;
                }

                // number of dofs per element
                uint tNumDofs = this->enforce_linear_interpolation() ?
                          mesh::number_of_corner_nodes(
                                    tFacets( 0 )->element()->type() )
                        : mesh::number_of_nodes(
                                    tFacets( 0 )->element()->type() ) ;
                tNumDofs *= mIWG->number_of_dofs_per_node() ;

                // element contribution
                Vector< real > tConvection( tNumDofs );

                // loop over all facets
                for( Element * tElement : tFacets )
                {
                    // compute contribution of this element
                    mIWG->compute_convection( tElement, tConvection );

                    // add heat load to master
                    for ( uint i = 0; i < tNumDofs; ++i )
                    {
                        Dof * tRow = tElement->dof( i );

                        if ( ! tRow->is_fixed()  )
                        {
                            mConvection( tRow->my_index() ) += tConvection( i ) ;
                        }
                    }
                }
            }

            // scale convection with timestep of iwg
            mConvection *= mIWG->timestep() ;
            this->collect_vector( mConvection );
        }

//------------------------------------------------------------------------------

        void
        Field::collect_node_owners()
        {
            if( comm_size() > 1 )
            {
                if ( mMyRank != mParent->master() )
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

                    send( mParent->master(), tIDs );
                }
                else
                {
                    uint tNumProcs = mParent->comm_table().length();

                    Cell< Vector< id_t > > tAllIDs( tNumProcs, Vector< id_t >() );

                    receive( mParent->comm_table(), tAllIDs );

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
        Field::collect_field( const string & aLabel )
        {
            // grab field
            Vector< real > & tField = mMesh->field_data( aLabel );

            if( mMyRank != mParent->master() )
            {
                Vector< real > tSubField( mMyNumberOfOwnedNodes ) ;

                // initialize counter
                index_t tCount = 0 ;

                // collect owned data
                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    if( tNode->owner() == mMyRank )
                    {
                        tSubField( tCount++ ) = tField( tNode->index() );
                    }
                }

                // send data to master
                send( mParent->master(), tSubField );
            }
            else
            {
                // collect subfields
                Cell< Vector< real > > tSubFields ;
                receive( mParent->comm_table(), tSubFields );

                // assemble
                uint tNumProcs = mParent->comm_table().length() ;

                for( uint p=0; p<tNumProcs; ++p )
                {
                    // grab index vector
                    const Vector< index_t > & tIndices = mNodeOwnerList( p );
                    const Vector< real >    & tSubField = tSubFields( p );

                    index_t tNumNodes = tIndices.length();

                    for ( index_t k = 0; k < tNumNodes; ++k )
                    {
                        tField( tIndices( k ) ) = tSubField( k );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::collect_fields( const Cell< string > & aLabels )
        {
            if( comm_size() > 1 )
            {
                if ( mMyRank != mParent->master() )
                {
                    uint tNumFields = aLabels.size();

                    Matrix< real > tSubFields( mMyNumberOfOwnedNodes, tNumFields );

                    for ( uint f = 0; f < tNumFields; ++f )
                    {
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
                    send( mParent->master(), tSubFields );
                }
                else
                {
                    // collect subfields
                    Cell< Matrix< real > > tAllSubFields;
                    receive( mParent->comm_table(), tAllSubFields );

                    uint tNumFields = tAllSubFields.size() ;

                    // assemble
                    uint tNumProcs = mParent->comm_table().length();

                    for ( uint p = 0; p < tNumProcs; ++p )
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

//------------------------------------------------------------------------------

        void
        Field::add_to_system_matrices(
                      Element         *  aElement,
                const Matrix< real >  & aJel,
                const Vector< real >  & aBel,
                const uint            & aN,
                      SpMatrix        & aJ,
                      SpMatrix        & aD,
                      Vector < real > & aB )
        {
            // add element jacobian to master
            for ( uint i = 0; i < aN; ++i )
            {
                Dof * tRow = aElement->dof( i );

                if ( !tRow->is_fixed() )
                {
                    for ( uint j = 0; j < aN; ++j )
                    {
                        Dof * tCol = aElement->dof( j );
                        if ( !tCol->is_fixed() )
                        {
                            aJ( tRow->index(), tCol->index() ) += aJel( i, j );
                        }
                        else
                        {
                            aD( tRow->index(), tCol->index() ) -= aJel( i, j );
                        }
                    }
                }
            }

            // add residual to vector
            for ( uint i = 0; i < aN; ++i )
            {
                Dof * tRow = aElement->dof( i );

                if ( ! tRow->is_fixed()  )
                {
                    aB( tRow->my_index() ) += aBel( i );
                }
            }
        }
//------------------------------------------------------------------------------

        void
        Field::compute_node_based_adjacency(
                Cell< graph::Vertex * > & aDOFs,
                Cell< Dof * > & aAdjacency,
                const bool aFlag )
        {
            // loop over all DOFs in graph
            for ( graph::Vertex * tDOF : aDOFs )
            {
                // get parent node
                mesh::Node * tNode = reinterpret_cast< Dof * >( tDOF )->node();

                // get number of nodes that are connected to this node
                uint tNumberOfNodes = tNode->number_of_nodes();

                // allocate temporary cell with dofs
                BELFEM_ASSERT( aAdjacency.size() >= mNumberOfDOFsPerNode * ( tNumberOfNodes + 1 ),
                    "Size of work array aAdjacency does not match, is %u but need at least %u",
                              ( unsigned int ) aAdjacency.size(),
                              ( unsigned int ) mNumberOfDOFsPerNode * ( tNumberOfNodes + 1 ) ) ;

                // initialize counter
                uint tCount = 0;

                // add self to list
                for ( uint i = 0; i < mNumberOfDOFsPerNode; ++i )
                {
                    // get neighbor
                    Dof * tNeighbor = this->dof( this->calculate_dof_id( tNode , i ) );

                    // test if neighbor is used and fixed
                    if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                    {
                        aAdjacency( tCount++ ) = tNeighbor;
                    }
                }

                // add neighbors to list
                for ( uint k = 0; k < tNumberOfNodes; ++k )
                {
                    for ( uint i = 0; i < mNumberOfDOFsPerNode; ++i )
                    {
                        // calculate ID of candidate
                        id_t tID = this->calculate_dof_id( tNode->node( k ), i );

                        // test if dof exists
                        if ( mDofMap.key_exists( tID ) )
                        {
                            // get pointer to neighbor
                            Dof * tNeighbor = mDOFs( mDofMap( tID ) );

                            // test if neighbor is used and fixed
                            if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                            {
                                aAdjacency( tCount++ ) = tNeighbor;
                            }
                        }
                    }
                }

                tDOF->reset_vertex_container();

                if ( tCount > 0 )
                {
                    // now we have the relevant nodes in the list
                    // we can now add the to the vertex
                    tDOF->init_vertex_container( tCount );

                    for ( uint k = 0; k < tCount; ++k )
                    {
                        tDOF->insert_vertex( aAdjacency( k ) );
                    }
                }

            }
        }

//------------------------------------------------------------------------------

        void
        Field::compute_edge_based_adjacency(
                Cell< graph::Vertex * > & aDOFs,
                Cell< Dof * > & aAdjacency,
                const bool aFlag )
        {
            // loop over all DOFs in graph
            for ( graph::Vertex * tDOF : aDOFs )
            {
                // get parent node
                mesh::Edge * tEdge = reinterpret_cast< Dof * >( tDOF )->edge();

                // get number of edges that are connected to this node
                uint tNumberOfEdges = tEdge->number_of_edges();

                // allocate temporary cell with dofs
                BELFEM_ASSERT( aAdjacency.size() >= mNumberOfDOFsPerEdge * ( tNumberOfEdges + 1 ),
                              "Size of work array aAdjacency does not match, is %u but need at least %u",
                              ( unsigned int ) aAdjacency.size(),
                              ( unsigned int ) mNumberOfDOFsPerEdge * ( tNumberOfEdges + 1 ) ) ;

                // initialize counter
                uint tCount = 0;

                // add self to list
                for ( uint i = 0; i < mNumberOfDOFsPerEdge; ++i )
                {
                    // get neighbor
                    Dof * tNeighbor = this->dof( this->calculate_dof_id( tEdge , i ) );

                    // test if neighbor is used and fixed
                    if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                    {
                        aAdjacency( tCount++ ) = tNeighbor;
                    }
                }

                // add neighbors to list
                for ( uint k = 0; k < tNumberOfEdges; ++k )
                {
                    for ( uint i = 0; i < mNumberOfDOFsPerEdge; ++i )
                    {
                        // calculate ID of candidate
                        id_t tID = this->calculate_dof_id( tEdge->edge( k ), i );

                        // test if dof exists
                        if ( mDofMap.key_exists( tID ) )
                        {
                            // get pointer to neighbor
                            Dof * tNeighbor = mDOFs( mDofMap( tID ) );

                            // test if neighbor is used and fixed
                            if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                            {
                                aAdjacency( tCount++ ) = tNeighbor;
                            }
                        }
                    }
                }

                tDOF->reset_vertex_container();

                if ( tCount > 0 )
                {
                    // now we have the relevant nodes in the list
                    // we can now add the to the vertex
                    tDOF->init_vertex_container( tCount );

                    for ( uint k = 0; k < tCount; ++k )
                    {
                        tDOF->insert_vertex( aAdjacency( k ) );
                    }
                }

            } // end create graph
        }
//------------------------------------------------------------------------------

        void
        Field::compute_node_and_edge_based_adjacency(
                Cell< graph::Vertex * > & aDOFs,
                Cell< Dof * > & aAdjacency,
                const bool aFlag )
        {
            // loop over all DOFs in graph
            for ( graph::Vertex * tV: aDOFs )
            {
                Dof * tDOF = reinterpret_cast< Dof * >( tV );

                // get parent vertex
                mesh::Vertex * tVertex = tDOF->mesh_vertex();

                // get number of nodes and edges that are connected to this node
                uint tNumberOfNodes = tVertex->number_of_nodes();
                uint tNumberOfEdges = tVertex->number_of_edges();

                // allocate temporary cell with dofs
                BELFEM_ASSERT( aAdjacency.size() >= mNumberOfDOFsPerNode * ( tNumberOfNodes + 1 ) +
                                                   mNumberOfDOFsPerEdge * ( tNumberOfEdges + 1 ),
                              "Size of work array aAdjacency does not match, is %u but need at least %u",
                              ( unsigned int ) aAdjacency.size(),
                              ( unsigned int ) mNumberOfDOFsPerNode * ( tNumberOfNodes + 1 ) +
                              mNumberOfDOFsPerEdge * ( tNumberOfEdges + 1 ) ) ;

                // initialize counter
                uint tCount = 0;

                // add self to list
                if( tDOF->is_node() )
                {
                    // node
                    for ( uint i = 0; i < mNumberOfDOFsPerNode; ++i )
                    {
                        // get neighbor
                        Dof * tNeighbor = this->dof( this->calculate_dof_id( tDOF->node() , i ) );

                        // test if neighbor is used and fixed
                        if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                        {
                            aAdjacency( tCount++ ) = tNeighbor;
                        }
                    }
                }
                else
                {
                    // edge
                    for ( uint i = 0; i < mNumberOfDOFsPerEdge; ++i )
                    {
                        // get neighbor
                        Dof * tNeighbor = this->dof( this->calculate_dof_id( tDOF->edge() , i ) );

                        // test if neighbor is used and fixed
                        if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                        {
                            aAdjacency( tCount++ ) = tNeighbor;
                        }
                    }
                }

                // add node neighbors to list
                for ( uint k = 0; k < tNumberOfNodes; ++k )
                {
                    for ( uint i = 0; i < mNumberOfDOFsPerNode; ++i )
                    {
                        // calculate ID of candidate
                        id_t tID = this->calculate_dof_id( tVertex->node( k ), i );

                        // test if dof exists
                        if ( mDofMap.key_exists( tID ) )
                        {
                            // get pointer to neighbor
                            Dof * tNeighbor = mDOFs( mDofMap( tID ) );

                            // test if neighbor is used and fixed
                            if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                            {
                                aAdjacency( tCount++ ) = tNeighbor;
                            }
                        }
                    }
                }

                // add edge neighbors to list
                for ( uint k = 0; k < tNumberOfEdges; ++k )
                {
                    for ( uint i = 0; i < mNumberOfDOFsPerEdge; ++i )
                    {
                        // calculate ID of candidate
                        id_t tID = this->calculate_dof_id( tVertex->edge( k ), i );

                        // test if dof exists
                        if ( mDofMap.key_exists( tID ) )
                        {
                            // get pointer to neighbor
                            Dof * tNeighbor = mDOFs( mDofMap( tID ) );

                            // test if neighbor is used and fixed
                            if ( tNeighbor->is_flagged() && tNeighbor->is_fixed() == aFlag )
                            {
                                aAdjacency( tCount++ ) = tNeighbor;
                            }
                        }
                    }
                }

                tDOF->reset_vertex_container();

                if ( tCount > 0 )
                {
                    // now we have the relevant nodes in the list
                    // we can now add the to the vertex
                    tDOF->init_vertex_container( tCount );

                    for ( uint k = 0; k < tCount; ++k )
                    {
                        tDOF->insert_vertex( aAdjacency( k ) );
                    }
                }

            }
        }

//------------------------------------------------------------------------------

        void
        Field::initialize_linear_projection_lists()
        {
            Cell< Vector< id_t > > tAllNodeIDs;
            Vector< id_t > tMyNodeIDs;

            // start a timer
            Timer tTimer;

            const bool tIsMaster = mMyRank == mParent->master() ;

            // get node container of mesh
            Cell< mesh::Node * > & tNodes = mMesh->nodes();

            const uint tNumProcs = mParent->number_of_procs() ;

            if ( tIsMaster )
            {

                // make sure that procs are in consecutive order
                proc_t q = 0;
                for ( proc_t p : mParent->comm_table() )
                {
                    BELFEM_ERROR( p == q++, "proc table must be consecutive" );
                }

                // container for nodes per proc
                Vector< index_t > tCount( tNumProcs, 0 );

                // - - - - - - - - - - - - - - - - - - - -
                // STEP 1 : identify corner nodes per proc
                // - - - - - - - - - - - - - - - - - - - -

                // unflag all nodes on mesh
                mMesh->unflag_all_nodes();

                // loop over all blocks
                // find corner nodes
                for ( index_t tID : mBlockIDs )
                {
                    // grab element list of this block
                    Cell< mesh::Element * > & tElements = mMesh->block( tID )->elements();

                    for ( mesh::Element * tElement : tElements )
                    {
                        tElement->flag_corner_nodes();
                    }
                }

                // count flagged nodes
                for ( mesh::Node * tNode : tNodes )
                {
                    // check if node is selected
                    if ( tNode->is_flagged())
                    {
                        // increment counter
                        ++tCount( tNode->owner());
                    }
                }

                // allocate id and index containers
                tAllNodeIDs.set_size( tNumProcs, {} );
                mAllCornerNodeIndices.set_size( tNumProcs, {} );

                for ( uint p = 0; p < tNumProcs; ++p )
                {
                    if ( tCount( p ) > 0 )
                    {
                        tAllNodeIDs( p ).set_size( tCount( p ));
                        mAllCornerNodeIndices( p ).set_size( tCount( p ));
                    }
                }

                // reset counters
                tCount.fill( 0 );

                // collect node indices
                for ( mesh::Node * tNode : tNodes )
                {
                    if ( tNode->is_flagged())
                    {
                        // get owner
                        index_t tOwn = tNode->owner();

                        // add node to list
                        tAllNodeIDs( tOwn )( tCount( tOwn )) = tNode->id();
                        mAllCornerNodeIndices( tOwn )( tCount( tOwn )++ ) = tNode->index();
                    }
                }

                // wait
                comm_barrier();

                // send IDs to other procs
                send( mParent->comm_table(), tAllNodeIDs );
            }
            else
            {

                // wait
                comm_barrier() ;

                // get corner node IDs from master
                receive( mParent->master(), tMyNodeIDs );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2 : identify indices of corner nodes
            //          ( this is the same for all procs )
            // - - - - - - - - - - - - - - - - - - - - - - - -

            Vector< id_t > & tCornerNodeIDs = tIsMaster ? tAllNodeIDs( 0 ) : tMyNodeIDs ;

            // number of nodes
            index_t tNumNodes = tCornerNodeIDs.length() ;

            // allocate memory
            mMyCornerNodeIndices.set_size( tNumNodes );

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

            if( tIsMaster )
            {
                // - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 3 : identify IDs of non corner node indices
                // - - - - - - - - - - - - - - - - - - - - - - - - -

                // unflag
                mMesh->unflag_all_nodes() ;

                // first, we flag all nodes of the elements
                for ( index_t tID : mBlockIDs )
                {
                    // grab element list of this block
                    Cell< mesh::Element * > & tElements = mMesh->block( tID )->elements() ;

                    for ( mesh::Element * tElement : tElements )
                    {
                        tElement->flag_nodes() ;
                    }
                }


                // second, we unflag all corner nodes of the elements
                for ( index_t tID : mBlockIDs )
                {
                    // grab element list of this block
                    Cell< mesh::Element * > & tElements = mMesh->block( tID )->elements() ;

                    for ( mesh::Element * tElement : tElements )
                    {
                        tElement->unflag_corner_nodes() ;
                    }
                }

                // now we introduce a new counter
                Vector< index_t > tCount ( mParent->number_of_procs(), 0 );

                // count nodes per proc
                for( mesh::Node * tNode : tNodes )
                {
                    // check if node is selected
                    if( tNode->is_flagged() )
                    {
                        // increment counter
                        ++tCount( tNode->owner() );
                    }
                }

                // reset container
                tAllNodeIDs.clear() ;
                tAllNodeIDs.set_size( tNumProcs, {} );

                // loop over all procs
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // get node IDs for this proc
                    Vector< id_t > & tNodeIDs = tAllNodeIDs( p );

                    // allocate memory
                    tNodeIDs.set_size( tCount( p ) );
                }

                // reset counter
                tCount.fill( 0 );

                // collect nodes per proc
                for( mesh::Node * tNode : tNodes )
                {
                    // check if node is selected
                    if( tNode->is_flagged() )
                    {
                       tAllNodeIDs( tNode->owner() )( tCount( tNode->owner() ) ++ )
                        = tNode->id() ;
                    }
                }

                // wait
                comm_barrier() ;

                // send data to others
                send( mParent->comm_table(), tAllNodeIDs );
            }
            else
            {
                // wait
                comm_barrier() ;

                receive( mParent->master(), tMyNodeIDs );
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

            for( Block * tBlock : mBlocks )
            {
                // get real element type of block ( tBlock->element_type() is the linear one )
                ElementType tType = tBlock->block()->element_type() ;

                // get number of nodes per element
                uint tNumCornerNodes = mesh::number_of_corner_nodes( tType );
                uint tNumNodes       = mesh::number_of_nodes( tType );

                // grab elements on block
                Cell< mesh::Element * > & tElements = tBlock->block()->elements() ;

                // loop over all elements on this block
                for( mesh::Element * tElement : tElements )
                {
                    // loop over all non corner nodes of this element
                    for( uint i = tNumCornerNodes; i<tNumNodes; ++i )
                    {
                        // grab node
                        mesh::Node * tNode = tElement->node( i );

                        // check if node has been processed
                        if( tNode->is_flagged() )
                        {
                            mMyNonCornerNodeIndices( k ) = tNode->index() ;

                            // store ID for sending
                            tMyNodeIDs( k ) = tNode->id() ;

                            // add node to container
                            mMyNonCornerNodes( k++ ) = tNode ;

                            // unflag this node, each node is processed only once
                            tNode->unflag() ;
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
                receive( mParent->comm_table(), tAllNodeIDs );

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
                send( mParent->master(), tMyNodeIDs );
            }

            // wait
            comm_barrier() ;

            if ( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for collecting procjection nodes   : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }
        }

//------------------------------------------------------------------------------

        void
        Field::communicate_corner_node_data( const Cell< string > & aFieldLabels )
        {
            index_t tNumFields = aFieldLabels.size() ;
            uint tNumProcs = mParent->number_of_procs() ;

            if ( mMyRank == mParent->master() )
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

                send( mParent->comm_table(), tAllData );
            }
            else
            {
                // Data container
                Vector< real > tData ;

                // wait
                comm_barrier() ;
                receive( mParent->master(), tData );

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
        Field::communicate_noncorner_node_data(
                const Cell< string > & aFieldLabels,
                      Matrix< real > & aData )
        {
            // get number of fields
            index_t tNumFields = aFieldLabels.size() ;
            uint tNumProcs = mParent->number_of_procs();

            if( mMyRank == mParent->master() )
            {

                // container with other data
                Cell< Matrix< real > > tAllData( mParent->number_of_procs(), {} );

                receive( mParent->comm_table(), tAllData );

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
                send( mParent->master(), aData );
            }

            // wait
            comm_barrier() ;
        }

//------------------------------------------------------------------------------

        void
        Field::project_linear_field_to_higher_mesh(
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
            for( Block * tBlock : mBlocks )
            {
                // get element type
                const ElementType tType = tBlock->block()->element_type() ;

                if( tType != ElementType::EMPTY )
                {
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
                    Cell< mesh::Element * > & tElements = tBlock->block()->elements();

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

            if ( mMyRank == mParent->master() )
            {

                message( 4, "    ... time for linear data projection         : %u ms\n",
                         ( unsigned int ) tTimer.stop() );

            }
        }

//------------------------------------------------------------------------------
    }
}
