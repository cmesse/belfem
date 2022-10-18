//
// Created by christian on 7/7/21.
//
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "en_FEM_DomainType.hpp"
#include "fn_entity_type.hpp"
#include "cl_IWG.hpp"

namespace belfem
{
    namespace fem
    {
//-----------------------------------------------------------------------------

        DofManager::DofManager(
                Kernel * aParent,
                const index_t aIndex ) :
                DofManagerBase( DofManagerType::NEW, aParent ),
                mIndex( aIndex )
        {

            // create the parameter object
            mParams = new dofmgr::Parameters( aParent,  aIndex );

            // create the dof data object
            mDofData = new dofmgr::DofData( this, mParams ) ;

            // create the block data object
            mBlockData = new dofmgr::BlockData( this );

            // create the sideset data object
            mSideSetData = new dofmgr::SideSetData( this );

            // create the bearing object
            mBearingData = new dofmgr::BearingData( this );

            // create the data object that manages field data
            mFieldData = new dofmgr::FieldData( this );

            // create the solver data object
            mSolverData = new dofmgr::SolverData( this, mDofData, mBlockData, mSideSetData );
        }

//-----------------------------------------------------------------------------

        DofManager::~DofManager()
        {
            // delete the solver object
            delete mSolverData ;

            // delete the field data object
            delete mFieldData ;

            // delete the bearing object
            delete mBearingData ;

            // delete the sideset data object
            delete mSideSetData ;

            // delete the block data object
            delete mBlockData ;

            // delete the dof data object
            delete mDofData ;

            // delete the parameter object
            delete mParams ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::set_equation( IWG * aIWG )
        {
            BELFEM_ERROR( aIWG->is_initialized(), "equation must be initialized before linking to field object" );

            // restore factory settings
            this->reset() ;

            comm_barrier() ;

            // link to equation object
            mIWG = aIWG ;

            // backwards link
            aIWG->set_field( this );

            // create the dofs
            mDofData->create_dofs( aIWG );

            // create the fields
            //this->create_fields( aIWG ) ;

            // create the blocks
            mBlockData->create_blocks() ;

            // create the side sets if they exist
            if( aIWG->selected_sidesets().length() > 0 )
            {
                mSideSetData->create_sidesets() ;

                // create the bearings
                mBearingData->create_bearings() ;
            }

            // create the field information
            mFieldData->collect_node_owners() ;
            mFieldData->collect_ghost_element_owners() ;

            // check if the linear projection flag is on
            if( this->enforce_linear_interpolation() )
            {
                // collect critical nodes
                mFieldData->initialize_linear_projection_lists() ;
            }
        }
//-----------------------------------------------------------------------------

        void
        DofManager::reset()
        {
            mInitializedFlag = false ;

            mDofData->reset() ;
            mBlockData->reset() ;
            mSideSetData->reset() ;
            mBearingData->reset() ;
            mFieldData->reset() ;
            mSolverData->reset() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::print(  const proc_t aRank  )
        {
            if( mMyRank == aRank )
            {
                std::cout << "DofManager" << std::endl;
                mParams->print();
            }
        }

//------------------------------------------------------------------------

        void
        DofManager::create_fields( IWG * aIWG )
        {
            // number of fields on the IWG
            index_t tNumFields = aIWG->number_of_fields() ;

            // loop over all fields on the IWG
            for ( index_t f=0; f<tNumFields; ++f )
            {
                // get the field label
                const string & tLabel = aIWG->field( f );

                // check if this could be an edge field
                EntityType tType = belfem::entity_type( tLabel );

                // create the field if it doesn't exist yet
                Vector< real > & tField = mMesh->field_exists( tLabel ) ?
                        mMesh->field_data( tLabel ) :
                        mMesh->create_field( tLabel, tType );


                // compute size of field
                index_t tFieldSize ;
                switch( tType )
                {
                    case( EntityType::NODE ) :
                    {
                        tFieldSize = mMesh->number_of_nodes() ;
                        break ;
                    }
                    case( EntityType::EDGE ) :
                    {
                        tFieldSize = aIWG->edge_multiplicity() * mMesh->number_of_edges() ;
                        break ;
                    }
                    case( EntityType::FACET ) :
                    {
                        tFieldSize = (
                                mMesh->number_of_facets() * aIWG->lambda_multiplicity()
                                + mMesh->number_of_connectors() ) ;

                        break ;
                    }
                    case( EntityType::ELEMENT ) :
                    {
                        tFieldSize = mMesh->number_of_elements() ;
                        break ;
                    }
                    case( EntityType::FACE ) :
                    {
                        tFieldSize = aIWG->face_multiplicity() * mMesh->number_of_faces() ;
                        break ;
                    }
                    /*case( EntityType::CELL ) :
                    {
                        tFieldSize = aIWG->cell_multiplicity() * mMesh->number_of_elements() ;
                        break ;
                    }*/
                    // entities that sit on shells
                    default:
                    {
                        tFieldSize = 0 ;
                        BELFEM_ERROR( false, "unsupported entity type for creating field");
                    }
                }

                if( tField.length() != tFieldSize )
                {
                    tField.set_size( tFieldSize, 0.0 );
                }
            }

            // with the fields created, we create the map on the dof data
            mDofData->create_field_map( aIWG );

            // tell exodus which fields are not to be written to exodus
            aIWG->hide_fields_from_exodus( mMesh );


        }

//------------------------------------------------------------------------

        void
        DofManager::init_dofs()
        {
            // check if user has set the wetted sidesets
            mSideSetData->detect_wetted_sidesets() ;

            // count how many nodes are wet ( so that the convection table is not needed )
            mSideSetData->count_wetted_nodes() ;

            // synchronize dof-relevant fields
            mFieldData->distribute( mIWG->dof_fields() );

            if( mIWG->selected_sidesets().length() > 0 )
            {
                // check if the fields "alpha" and "Tinf" must be created
                mSideSetData->create_alpha_fields() ;

                // initialize sideset BCs
                mSideSetData->set_boundary_conditions();
            }

            // collect node field data from other procs
            const Cell< string > & tFields = mIWG->all_fields() ;
            uint tNumFields = tFields.size() ;
            for( uint f=0; f<tNumFields; ++f )
            {
                const string & tField = tFields( f );

                if( mMesh->field( tField )->entity_type() == EntityType::NODE ||
                    mMesh->field( tField )->entity_type() == EntityType::ELEMENT )
                {
                    mFieldData->collect( tField );
                }
            }

            // send information back
            mFieldData->distribute( mIWG->all_fields() );

            // set values for dofs
            mDofData->init_dof_values() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::initialize()
        {
            BELFEM_ASSERT( mIWG->is_initialized(), "initialize iwg first");
            BELFEM_ASSERT( mSolverData->solver() != nullptr, "set solver first" );

            // create the field list
            this->create_fields( mIWG );

            // link dofs with fields on mesh
            mFieldData->update_field_indices( mDofData->dofs() );

            this->init_dofs();

            // allocate the matrices
            this->init_matrices();

            // auto set blocks for other fields
            this->auto_set_materials() ;

            mInitializedFlag = true ;

            // inpitialize postprocessors, if they exist
            this->initialize_postprocessors();

        }

//-----------------------------------------------------------------------------

        void
        DofManager::zero()
        {
            mSolverData->reset_matrices() ;

            if( mIWG->num_rhs_cols() == 1 )
            {
                mSolverData->reset_rhs_vector() ;
            }
            else
            {
                mSolverData->reset_rhs_matrix() ;
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::init_matrices()
        {
            // initialize timer
            Timer tTimer;
            if ( mMyRank == mParent->master() )
            {
                message( 4, " Initialize Jacobian ... \n" );
            }

            mDofData->init_dirichlet_bcs() ;

            mDofData->compute_dof_indices();

            mSolverData->allocate_matrices();

            mSolverData->create_assembly_tables();


            if ( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for initializing Jacobian          : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }

        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_jacobian( const bool aReset )
        {
            if( ! mInitializedFlag )
            {
                this->initialize();
            }

            Timer tTimer;

            // in most cases, we want to reset all matrices
            // unless we impose a weak BC first
            if ( aReset )
            {
                mSolverData->reset_matrices();
            }

            if( mIWG->compute_jacobian_on_block() )
            {
                // loop over all blocks
                for ( Block * tBlock : mBlockData->blocks() )
                {
                    // get the number of dofs per element
                    uint tN = mDofData->num_dofs_per_element( tBlock->id() );

                    BELFEM_ASSERT( tN > 0, "No dofs for block %lu", ( long unsigned int ) tBlock->id() );

                    // allocate matrix
                    Matrix< real > tJ( tN, tN );

                    mIWG->link_to_group( tBlock );

                    // get elements on block
                    Cell< Element * > & tElements = tBlock->elements();

                    // loop over all elements
                    for ( Element * tElement : tElements )
                    {

                        // compute element contribution
                        mIWG->compute_jacobian( tElement, tJ );

                        // add contribution to system matrix
                        mSolverData->assemble_jacobian( tElement, tJ );
                    }
                }
            }

            if( mIWG->compute_jacobian_on_sideset() )
            {
                // loop over all sidesets
                for ( SideSet * tSideSet : mSideSetData->sidesets() )
                {
                    // get the number of dofs per element
                    uint tN = mDofData->num_dofs_per_facet( tSideSet->id() );

                    // allocate matrix
                    Matrix< real > tJ( tN, tN );

                    mIWG->link_to_group( tSideSet );

                    // get elements on sideset
                    Cell< Element * > & tElements = tSideSet->elements();

                    // loop over all elements
                    for ( Element * tElement : tElements )
                    {

                        // compute element contribution
                        mIWG->compute_jacobian( tElement, tJ );

                        // add contribution to system matrix
                        mSolverData->assemble_jacobian( tElement, tJ );
                    }
                }
            }

            mSolverData->collect_jacobian();

            if ( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for computing Jacobian             : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_rhs( const bool aReset )
        {
            if( ! mInitializedFlag )
            {
                this->initialize();
            }

            Timer tTimer;

            if( mIWG->num_rhs_cols() == 1 )
            {
                // check if we want to reset the vector
                if( aReset )
                {
                    mSolverData->reset_rhs_vector() ;
                }

                // compute and assemble the vector
                this->compute_rhs_vector() ;

                // wait for other procs to finish
                comm_barrier() ;


                // collect the contributions from the other procs
                mSolverData->collect_rhs_vector() ;
            }
            else
            {
                // check if we want to reset the vector
                if( aReset )
                {
                    mSolverData->reset_rhs_matrix() ;
                }

                // compute and assemble the matrix
                this->compute_rhs_matrix() ;

                // wait for other procs to finish
                comm_barrier() ;

                // collect the contributions from the other procs
                mSolverData->collect_rhs_matrix() ;
            }

            if( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for computing right hand side      : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }

        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_jacobian_and_rhs( const bool aReset )
        {
            if( ! mInitializedFlag )
            {
                this->initialize();
            }

            Timer tTimer;

            // in most cases, we want to reset all matrices
            // unless we impose a weak BC first
            if ( aReset )
            {
                mSolverData->reset_matrices();
                mSolverData->reset_rhs_vector() ;
            }

            // loop over all blocks
            for ( Block * tBlock : mBlockData->blocks() )
            {
                // jump to next block if this block is not used
                if( ! tBlock->is_active() )
                {
                    continue ;
                }

                // get the number of dofs per element
                uint tN = mDofData->num_dofs_per_element( tBlock->id() );

                // allocate matrix
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
                    // compute element contribution
                    mIWG->compute_jacobian_and_rhs( tElement, tJ, tB );

                    // add contribution to system matrix
                    mSolverData->assemble_jacobian_and_rhs( tElement, tJ, tB );
                }
            }

            for ( SideSet * tSideSet : mSideSetData->sidesets() )
            {
                // jump to next sideset if this block is not used
                if( ! tSideSet->is_active() )
                {
                    continue ;
                }

                // link IWG to block
                mIWG->link_to_group( tSideSet );

                // get the number of dofs per element
                uint tN = mIWG->number_of_dofs_per_element( tSideSet );

                // allocate matrix
                Matrix< real > tJ( tN, tN );

                // allocate element RHS
                Vector< real > tB( tN );

                // get elements on block
                Cell< Element * > & tElements = tSideSet->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    // compute element contribution
                    mIWG->compute_jacobian_and_rhs( tElement, tJ, tB );

                    // add contribution to system matrix
                    mSolverData->assemble_jacobian_and_rhs( tElement, tJ, tB );
                }
            }

            if( mIWG->has_alpha() )
            {
                for( SideSet * tSideSet : mSideSetData->sidesets() )
                {
                    if( ! tSideSet->is_active() )
                    {
                        continue ;
                    }

                    // check if this is an alpha BC
                    if ( tSideSet->bc_type( 0 ) == BoundaryConditionImposing::Alpha )
                    {
                        // link this IWG with the group
                        mIWG->link_to_group( tSideSet ) ;

                        // get number of dofs per element
                        uint tN = this->enforce_linear_interpolation() ?
                                  mesh::number_of_corner_nodes( tSideSet->element_type() )
                                   :  mesh::number_of_nodes( tSideSet->element_type() );

                        // allocate element Jacobian
                        Matrix< real > tJ( tN, tN );

                        // allocate element RHS
                        Vector< real > tB( tN );

                        // get elements on Sideset
                        Cell< Element * > & tElements = tSideSet->elements();

                        // loop over all elements on this sideset
                        for ( Element * tElement : tElements )
                        {
                            mIWG->compute_alpha_boundary_condition( tElement, tJ, tB );
                            mSolverData->assemble_jacobian_and_rhs( tElement, tJ, tB );
                        }

                    }
                }
            }

            // unite matrix with values from other procs
            mSolverData->collect_jacobian();

            // unite rhs with values from other procs
            mSolverData->collect_rhs_vector() ;

            // needed for computing the residual later on
            // this field only exists on the master
            mSolverData->update_field_values() ;


            if ( mMyRank == mParent->master() )
            {
                message( 4, "    ... time for computing Jacobian and residual : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_rhs_vector()
        {
            if( ! mInitializedFlag )
            {
                this->initialize();
            }

            for ( Block * tBlock : mBlockData->blocks() )
            {
                // check if block has an RHS
                if ( tBlock->has_rhs() )
                {
                    // get the number of dofs per element
                    uint tN = mDofData->num_dofs_per_element( tBlock->id());

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
                        mSolverData->asseble_rhs( tElement, tRHS );
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_rhs_matrix()
        {
            if( ! mInitializedFlag )
            {
                this->initialize();
            }
            // number of cols in rhs matrix
            uint tNumCols = mIWG->num_rhs_cols();

            for ( Block * tBlock : mBlockData->blocks() )
            {
                // get the number of dofs per element
                uint tN = mDofData->num_dofs_per_element( tBlock->id());

                // allocate element RHS
                Matrix< real > tRHS( tN, tNumCols );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    mIWG->compute_rhs( tElement, tRHS );
                    mSolverData->asseble_rhs( tElement, tRHS );
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::set_solver( const SolverType aSolver )
        {
            mSolverData->set_solver( aSolver );
            comm_barrier() ;
        }


//-----------------------------------------------------------------------------

        void
        DofManager::solve()
        {
            // solve the system
            mSolverData->solve();

            // make result available to other procs
            mFieldData->distribute( mIWG->all_fields() );

            // wait for other procs
            comm_barrier();
        }

//-----------------------------------------------------------------------------

        real
        DofManager::residual( const uint aIteration )
        {
            return mSolverData->residual( aIteration );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::auto_set_materials()
        {
            if( mIndex > 0 )
            {
                // loop over all blocks
                for( Block * tBlock : mBlockData->blocks() )
                {
                    // check if material is empty
                    if( tBlock->material() == nullptr )
                    {
                        // get material of block on first field
                        const Material * tMaterial =  mParent->dofmgr( 0 )->block( tBlock->id() )->material() ;

                        if( tMaterial != nullptr )
                        {
                            // set the material type
                            tBlock->set_material( string_to_material_type( tMaterial->label() ) );
                        }
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::collect_fields( const Cell< string > & aFieldLabels )
        {
            mFieldData->collect( aFieldLabels );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::collect_field( const string & aFieldLabel )
        {
            comm_barrier() ;
            mFieldData->collect( aFieldLabel );
            comm_barrier() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::distribute_fields( const Cell< string > & aFieldLabels )
        {
            mFieldData->distribute( aFieldLabels );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::synchronize_fields( const Cell< string > & aFieldLabels )
        {
            comm_barrier() ;
            mFieldData->collect( aFieldLabels );
            comm_barrier();
            mFieldData->distribute( aFieldLabels );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_volume_loads( const Vector< id_t > & aBlockIDs )
        {

            BELFEM_ASSERT( mIWG->num_rhs_cols() <= 1,
                "can't compute volume loads if RHS is a matrix");

            // get volume vector
            Vector< real > & tVolumeLoads = mSolverData->volume_loads() ;

            // make sure that vector is allocated
            if( tVolumeLoads.length() != mSolverData->my_number_of_free_dofs() )
            {
                tVolumeLoads.set_size( mSolverData->my_number_of_free_dofs() );
            }

            // reset the vector
            tVolumeLoads.fill( 0.0 );

            // loop over all blocks
            for( id_t tBlockID : aBlockIDs )
            {
                // get block
                Block * tBlock = mBlockData->block( tBlockID );

                // check if block has elements
                if( tBlock->number_of_elements() > 0 )
                {
                    // link IWG to block
                    mIWG->link_to_group( tBlock );

                    // rhs vector
                    Vector< real > tRHS( mIWG->number_of_dofs_per_element( tBlock ) );

                    // loop over all elements on block
                    for( Element * tElement : tBlock->elements() )
                    {
                        mIWG->compute_rhs( tElement, tRHS );
                        mSolverData->assemble_volume_loads( tElement, tRHS );
                    }
                }
            }

            // scale loads with timestep of IWG
            tVolumeLoads *= mIWG->timestep() ;

            // wait for other procs
            comm_barrier() ;

            // synchronize data
            mSolverData->collect_vector( tVolumeLoads );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::load_system( const string & aPath )
        {
            mSolverData->load_system( aPath );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::save_system( const string & aPath )
        {
            mSolverData->save_system( aPath );
        }

//-----------------------------------------------------------------------------

#ifdef BELFEM_HDF5
        void
        DofManager::save_system( HDF5 & aFile )
        {
            mSolverData->save_system( aFile );
        }

        void
        DofManager::load_system( HDF5 & aFile )
        {
            mSolverData->load_system( aFile );
        }
#endif

//-----------------------------------------------------------------------------

        void
        DofManager::init_dof_values()
        {
            mDofData->init_dof_values() ;
        }

        void
        DofManager::print_worst_dof()
        {
            mSolverData->print_worst_dof() ;
        }

//-----------------------------------------------------------------------------


        /**
         * compute the matrices for the projections
         */
        void
        DofManager::initialize_postprocessors()
        {
            for ( DofManager * tProjector : mPostprocessors )
            {
                tProjector->initialize() ;
                tProjector->compute_jacobian();
                comm_barrier();
            }
        }

//-----------------------------------------------------------------------------

        /**
         * perform the L2 projections for the secondary fields
         */
        void
        DofManager::postprocess()
        {
            for ( DofManager * tProjector : mPostprocessors )
            {
                //tProjector->compute_jacobian();
                tProjector->compute_rhs() ;
                //tProjector->save_system("projector.hdf5");
                tProjector->solve();
                comm_barrier();
            }
        }


//------------------------------------------------------------------------------

        void
        DofManager::fix_node_dofs_on_ghost_sidesets()
        {
            mMesh->unflag_everything() ;

            for( id_t tID : mIWG->ghost_sideset_ids() )
            {
                if( mMesh->sideset_exists( tID ) )
                {
                    mMesh->sideset( tID )->flag_all_nodes() ;
                }
            }

            // loop over all dofs
            for( Dof * tDof : mDofData->dofs() )
            {
                if( tDof->mesh_basis()->is_flagged() )
                {
                    tDof->fix( 0.0 );
                }
            }
        }

//-----------------------------------------------------------------------------
    }
}