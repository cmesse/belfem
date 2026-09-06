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

#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"
#include "cl_FEM_DofManager.hpp"

#include "cl_FEM_Controller.hpp"
#include "cl_FEM_Kernel.hpp"
#include "en_DomainType.hpp"
#include "fn_entity_type.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_Postprocessor.hpp"
#include "cl_Profiler.hpp"

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
            mParams = new dofmgr::Parameters( aParent );

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

            // crate the eigenvalue object
            mEigenValues = new dofmgr::EigenValues( this );

        }

//-----------------------------------------------------------------------------

        DofManager::~DofManager()
        {
            // delete the postprocessors
            for ( Postprocessor * tProc : mPostprocessors )
            {
                delete tProc;
            }

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

            // delete the eigenvalue object
            delete mEigenValues ;
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

            // create the blocks
            mBlockData->create_blocks() ;

            this->create_element_map() ;

            // create the side sets if they exist
            if( aIWG->selected_sidesets().length() > 0 )
            {
                mSideSetData->create_sidesets() ;
            }

            mBearingData->create_bearings() ;

            // create the field information
            mFieldData->collect_node_owners() ;
            mFieldData->collect_element_owners() ;

            mDofData->collect_hanging_dofs() ;

            // until here, dof indices must not be touched!

            comm_barrier() ;

            mIWG->collect_abstract_node_dofs();

            comm_barrier() ;

            this->consolidate_dofs() ;

            comm_barrier() ;

            this->link_bearings_with_dofs() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::reset()
        {
            mInitializedFlag = false ;
            mJacobianIsUpToDate = false ;

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
            if( mCommRank == aRank )
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
                                mMesh->number_of_facets() * aIWG->lambda_multiplicity() ) ;

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
                    case( EntityType::CELL ) :
                    {
                        tFieldSize = aIWG->cell_multiplicity() * mMesh->number_of_elements() ;
                        break ;
                    }
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
        DofManager::init_dofs( const bool aSeedFreeDofsOnly )
        {
            // check if user has set the wetted sidesets
            mSideSetData->collect_wetted_sidesets() ;

            // count how many nodes are wet ( so that the convection table is not needed )
            mSideSetData->count_wetted_nodes() ;

            // synchronize dof-relevant fields
            mFieldData->distribute( mIWG->dof_fields() );

            mSideSetData->create_alpha_fields() ;

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

            comm_barrier() ;


            // send information back
            mFieldData->distribute( mIWG->all_fields() );

            // set values for dofs
            mDofData->init_dof_values( aSeedFreeDofsOnly ) ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::init_work()
        {
            // loop over all blocks
            for ( Block * tBlock : mBlockData->blocks() )
            {
                if( tBlock->calculator() != nullptr )
                {
                    tBlock->calculator()->allocate() ;
                }
            }

            // loop over all sidesets
            for( SideSet * tSideSet : mSideSetData->sidesets() )
            {
                if( tSideSet->calculator() != nullptr )
                {
                    tSideSet->calculator()->allocate() ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::initialize()
        {
            this->initialize( false );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::initialize( const bool aSeedFreeDofsOnly )
        {
            if( mInitializedFlag )
            {
                return;
            }
            BELFEM_ASSERT( mIWG->is_initialized(), "initialize iwg first");
            BELFEM_ASSERT( mSolverData->solver() != nullptr, "set solver first" );

            // create the field list
            this->create_fields( mIWG );

            // link dofs with fields on mesh
            mFieldData->update_field_indices( mDofData->dofs() );
            mFieldData->update_field_indices( mDofData->hanging_dofs() );

            // initialize degrees of freedom
            this->init_dofs( aSeedFreeDofsOnly );

            // initialize work memory for the groups
            this->init_work() ;

            // allocate the matrices
            this->init_matrices();

            // auto set blocks for other fields
            this->auto_set_materials() ;

            // this information is not needed anymore, tyding up!
            mDofData->disconnect_dofs_from_mesh() ;

            mInitializedFlag = true ;

            // initialize postprocessors, if they exist
            this->initialize_postprocessors();
        }

//-----------------------------------------------------------------------------

        void
        DofManager::extract_abstract_dofs_from_mesh()
        {
            // extract abstract dofs before we disconnect the mesh
            mDofData->extract_abstract_dofs_from_mesh() ;
        }

//------------------------------------------------------------------------------

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
            if ( mCommRank == 0 )
            {
                message( InfoLevel::Verbose, " Initialize Jacobian ... \n" );
            }

            Cell< graph::Vertex * > tFreeDofs ;
            Cell< graph::Vertex * > tFixedDofs ;


            mDofData->synchronize_dirichlet_bcs();
            Vector< id_t > * tGraphData = new Vector< id_t >();
            mSolverData->extract_graph_from_mesh( *tGraphData );
            mDofData->split_dof_container( tFreeDofs, tFixedDofs ) ;

            mDofData->reorder_dofs( *tGraphData, tFreeDofs, tFixedDofs );
            mSolverData->allocate_matrices( *tGraphData, tFreeDofs, tFixedDofs );
            delete tGraphData;
            mDofData->restore_dof_container( tFreeDofs, tFixedDofs );

            mSolverData->create_assembly_tables();
            if ( mCommRank == 0 )
            {
                message( InfoLevel::Verbose, "    ... time for initializing Jacobian          : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }
            comm_barrier();
            mSolverData->compute_memory();
            comm_barrier();
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_jacobian( const bool aReset )
        {
            // initialize the dof manager if it hasn't been done already
            this->initialize();

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
                    Matrix< real > tTJT ;

                    Matrix< real > tTdJT ;

                    mIWG->link_to_group( tBlock );

                    // get elements on block
                    Cell< Element * > & tElements = tBlock->elements();

                    // loop over all elements
                    for ( Element * tElement : tElements )
                    {
                        if( tElement->has_t_matrix() )
                        {
                            mIWG->compute_jacobian( tElement, tJ );
                            tElement->t_matrix()->project( tJ, tTJT );

                            // add contribution to system matrix
                            tElement->t_matrix()->project( mIWG->matrices()->dJdx(), tTdJT );

                            mSolverData->assemble_jacobian( tElement, tTJT);
                            mSolverData->assemble_newton( tElement, tTdJT ) ;
                        }
                        else
                        {
                            // compute element contribution
                            mIWG->compute_jacobian( tElement, tJ );

                            // add contribution to system matrix
                            mSolverData->assemble_jacobian( tElement, tJ);
                            mSolverData->assemble_newton( tElement, mIWG->matrices()->dJdx() ) ;
                        }

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
                    Matrix< real > tTJT ;

                    Matrix< real > tTdJT ;

                    mIWG->link_to_group( tSideSet );

                    // get elements on sideset
                    Cell< Element * > & tElements = tSideSet->elements();

                    // loop over all elements
                    for ( Element * tElement : tElements )
                    {
                        if( tElement->has_t_matrix() )
                        {
                            mIWG->compute_jacobian( tElement, tJ );
                            tElement->t_matrix()->project( tJ, tTJT );
                            // add contribution to system matrix
                            tElement->t_matrix()->project( mIWG->matrices()->dJdx(), tTdJT );

                            mSolverData->assemble_jacobian( tElement, tTJT);
                            mSolverData->assemble_newton( tElement, tTdJT ) ;
                        }
                        else
                        {
                            // compute element contribution
                            mIWG->compute_jacobian( tElement, tJ );

                            // add contribution to system matrix
                            mSolverData->assemble_jacobian( tElement, tJ );
                            mSolverData->assemble_newton( tElement, mIWG->matrices()->dJdx()) ;
                        }
                    }
                }
            }

            mSolverData->collect_matrices();

            if ( mCommRank == 0 )
            {

                message( InfoLevel::Verbose, "    ... time for computing Jacobian             : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }

            mJacobianIsUpToDate = true ;

            // after the Jacobian has been computed,
            // previously stored eigenvalues are void.
            mEigenValues->reset() ;
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

            if( mCommRank == 0 )
            {
                message( InfoLevel::Verbose, "    ... time for computing right hand side      : %u ms\n",
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

            // a fresh assembly voids the two-phase residual handshake:
            // the RHS is rebuilt, so a previously computed residual no
            // longer describes it
            mSolverData->invalidate_residual() ;

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
                Matrix< real > tTJT ;

                // allocate element RHS
                Vector< real > tB( tN );
                Vector< real > tTB ;

                Matrix< real > tTdJT ;

                // link IWG to block
                mIWG->link_to_group( tBlock );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    // compute element contribution
                    mIWG->compute_jacobian_and_rhs( tElement, tJ, tB );

                    if( tElement->has_t_matrix() )
                    {
                        tElement->t_matrix()->project( tJ, tTJT );
                        tElement->t_matrix()->project( tB, tTB );

                        tElement->t_matrix()->project( mIWG->matrices()->dJdx(), tTdJT );
                        mSolverData->assemble_jacobian( tElement, tTJT);
                        mSolverData->assemble_newton( tElement, tTdJT ) ;
                        mSolverData->assemble_rhs( tElement, tTB );
                    }
                    else
                    {
                        // add contribution to system matrix
                        mSolverData->assemble_jacobian( tElement, tJ );
                        mSolverData->assemble_newton( tElement, mIWG->matrices()->dJdx()) ;
                        mSolverData->assemble_rhs( tElement, tB );
                    }
                }
            }

            for ( SideSet * tSideSet : mSideSetData->sidesets() )
            {
                // jump to next sideset if this block is not used
                if( ! tSideSet->is_active() || tSideSet->bc_type( 0 ) == BoundaryConditionImposing::Neumann )
                {
                    continue ;
                }

                // link IWG to block
                mIWG->link_to_group( tSideSet );

                // get the number of dofs per element
                uint tN = mIWG->number_of_dofs_per_element( tSideSet );

                // allocate matrix
                Matrix< real > tJ( tN, tN );
                Matrix< real > tTJT ;

                // allocate element RHS
                Vector< real > tB( tN );
                Vector< real > tTB ;

                Matrix< real > tTdJT ;

                // get elements on block
                Cell< Element * > & tElements = tSideSet->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {


                    if( tElement->has_t_matrix() )
                    {
                        // compute element contribution
                        mIWG->compute_jacobian_and_rhs( tElement, tJ, tB );

                        tElement->t_matrix()->project( tJ, tTJT );
                        tElement->t_matrix()->project( tB, tTB );

                        tElement->t_matrix()->project( mIWG->matrices()->dJdx(), tTdJT );

                        mSolverData->assemble_jacobian( tElement, tTJT);
                        mSolverData->assemble_newton( tElement, tTdJT ) ;
                        mSolverData->assemble_rhs( tElement, tTB );
                    }
                    else
                    {
                        // compute element contribution
                        mIWG->compute_jacobian_and_rhs( tElement, tJ, tB );

                        // add contribution to system matrix
                        mSolverData->assemble_jacobian( tElement, tJ );
                        mSolverData->assemble_newton( tElement, mIWG->matrices()->dJdx() ) ;
                        mSolverData->assemble_rhs( tElement, tB );
                    }
                }
            }

            if( mIWG->has_convection() )
            {
                for( SideSet * tSideSet : mSideSetData->sidesets() )
                {
                    if( ! tSideSet->is_active() )
                    {
                        continue ;
                    }

                    // check if this is an alpha BC
                    if ( tSideSet->bc_type( 0 ) == BoundaryConditionImposing::Neumann )
                    {
                        // link this IWG with the group
                        mIWG->link_to_group( tSideSet ) ;

                        // get number of dofs per element
                        uint tN = mesh::number_of_nodes( tSideSet->element_type() );

                        // allocate element RHS
                        Vector< real > tB( tN );
                        Vector< real > tTB ;

                        // get elements on Sideset
                        Cell< Element * > & tElements = tSideSet->elements();

                        // loop over all elements on this sideset
                        for ( Element * tElement : tElements )
                        {
                            if( tElement->has_t_matrix() )
                            {
                                mIWG->compute_convection( tElement, tB );
                                tElement->t_matrix()->project( tB, tTB );
                                mSolverData->asseble_surface_loads( tElement, tTB );
                            }
                            else
                            {
                                mIWG->compute_convection( tElement, tB );
                                mSolverData->asseble_surface_loads( tElement, tB );
                            }
                        }

                    }
                }
            }

            // unite matrix with values from other procs
            mSolverData->collect_matrices();

            // unite rhs with values from other procs
            mSolverData->collect_rhs_vector() ;

            // needed for computing the residual later on
            // this field only exists on the master
            mSolverData->update_field_values() ;


            if ( mCommRank == 0 )
            {
                message( InfoLevel::Verbose, "    ... time for computing Jacobian and residual : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }

            mJacobianIsUpToDate = true ;

            // after the Jacobian has been computed,
            // previously stored eigenvalues are void.
            mEigenValues->reset() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_full_matrices()
        {
            mSolverData->reset_matrices( true );

            index_t n = mSolverData->my_number_of_free_dofs() ;
            index_t m = n + mSolverData->my_number_of_fixed_dofs() ;
            index_t off = mSolverData->number_of_free_dofs() ;

            Cell< Dof * > & tDofs = mDofData->dofs() ;
            for ( index_t k=n; k<m; ++k )
            {
                Dof * tDof = tDofs( k ) ;
                tDof->set_index( tDof->index() + off );
                tDof->set_my_index( tDof->my_index() + n );
            }

            Timer tTimer;
            if ( mCommRank == 0 )
            {
                message( InfoLevel::Default, "    computing full matrices...\n" );
            }
            for ( Block * tBlock : mBlockData->blocks() )
            {
                // jump to next block if this block is not used
                if( ! tBlock->is_active() )
                {
                    continue ;
                }

                Matrix< real > tTKT ;
                Matrix< real > tTMT ;


                // link IWG to block
                mIWG->link_to_group( tBlock );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    // compute element contribution
                    mIWG->compute_mkf( tElement );

                    if( tElement->has_t_matrix() )
                    {
                        tElement->t_matrix()->project( mIWG->matrices()->M(), tTMT );
                        tElement->t_matrix()->project( mIWG->matrices()->K(), tTKT );

                        mSolverData->assemble_full_matrices( tElement, tTMT, tTKT );
                    }
                    else
                    {
                        // add contribution to system matrix
                        mSolverData->assemble_full_matrices( tElement, mIWG->matrices()->M(), mIWG->matrices()->K() );
                    }
                }
            }

            for ( index_t k=n; k<m; ++k )
            {
                Dof * tDof = tDofs( k ) ;
                tDof->set_index( tDof->index() - off );
                tDof->set_my_index( tDof->my_index() - n );
            }

            comm_barrier();

            mSolverData->collect_matrices( true );

            if ( mCommRank == 0 )
            {
                message( InfoLevel::Default, "    ... time for computing full matrices        : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::reset_convection()
        {
            if( ! mInitializedFlag )
            {
                this->initialize();
            }
            mSolverData->reset_convection() ;
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
                    Vector< real > tTRHS ;

                    // link IWG to block
                    mIWG->link_to_group( tBlock );

                    // get elements on block
                    Cell< Element * > & tElements = tBlock->elements();

                    // loop over all elements
                    for ( Element * tElement : tElements )
                    {
                        if( tElement->has_t_matrix() )
                        {
                            mIWG->compute_rhs( tElement, tRHS );
                            tElement->t_matrix()->project( tRHS, tTRHS );
                            mSolverData->assemble_rhs( tElement, tTRHS );
                        }
                        else
                        {
                            mIWG->compute_rhs( tElement, tRHS );
                            mSolverData->assemble_rhs( tElement, tRHS );
                        }

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
                Matrix< real > tTRHS ;

                // link IWG to block
                mIWG->link_to_group( tBlock );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements();

                // loop over all elements
                for ( Element * tElement : tElements )
                {
                    if( tElement->has_t_matrix() )
                    {
                        mIWG->compute_rhs( tElement, tRHS );

                        Vector< real > tG ;
                        tTRHS.set_size( tElement->number_of_dofs(), tNumCols );

                        for( uint k=0; k<tNumCols; ++k )
                        {
                            tElement->t_matrix()->project( tRHS.col( k ), tG );
                            tTRHS.set_col( k, tG );
                        }
                        mSolverData->assemble_rhs( tElement, tTRHS );
                    }
                    else
                    {
                        mIWG->compute_rhs( tElement, tRHS );
                        mSolverData->assemble_rhs( tElement, tRHS );
                    }
                }
            }

            // after the Jacobian has been computed,
            // previously stored eigenvalues are void.
            mEigenValues->reset() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::set_solver( const SolverParameters & aParams )
        {
            mSolverData->set_solver( aParams );
            comm_barrier() ;
        }


//-----------------------------------------------------------------------------

        void
        DofManager::solve()
        {
            BELFEM_ERROR( mJacobianIsUpToDate,
                "You need to compute or update the Jacobian matrix before calling solve()!" );

            // solve the system
            mSolverData->solve();

            // make result available to other procs
            mFieldData->distribute( mIWG->all_fields() );

            // wait for other procs
            comm_barrier();

            // the solve/distribute exchange ends here
            comm_drain_check( "DofManager::solve" );

            mJacobianIsUpToDate = false ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::compute_residual()
        {
            BELFEM_ERROR( mJacobianIsUpToDate,
                "You need to compute or update the Jacobian matrix before calling compute_residual()!" );

            // phase 1 of the two-phase iterative solve ( certified
            // exit ): residual of the committed state, no solve, no update.
            // Deliberately does NOT distribute fields and does NOT consume
            // mJacobianIsUpToDate — the assembly stays valid for the solve
            mSolverData->compute_residual();
        }

//-----------------------------------------------------------------------------

        void
        DofManager::solve_from_residual()
        {
            BELFEM_ERROR( mJacobianIsUpToDate,
                "You need to compute or update the Jacobian matrix before calling solve_from_residual()!" );

            // phase 2: linear solve + update on the residual phase 1 left
            mSolverData->solve_from_residual();

            // make result available to other procs
            mFieldData->distribute( mIWG->all_fields() );

            // wait for other procs
            comm_barrier();

            // the solve/distribute exchange ends here
            comm_drain_check( "DofManager::solve_from_residual" );

            mJacobianIsUpToDate = false ;
        }

//-----------------------------------------------------------------------------

        real
        DofManager::residual( const uint aIteration )
        {
            return mSolverData->residual( aIteration );
        }

//-----------------------------------------------------------------------------

        real
        DofManager::absolute_residual() const
        {
            return mSolverData->absolute_residual();
        }

//-----------------------------------------------------------------------------

        real
        DofManager::pre_update_residual() const
        {
            return mSolverData->pre_update_residual();
        }

//-----------------------------------------------------------------------------

        real
        DofManager::fixed_point_residual() const
        {
            return mSolverData->fixed_point_residual();
        }

//-----------------------------------------------------------------------------

        bool
        DofManager::solve_failed()
        {
            return mSolverData->solver()->wrapper()->failed();
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
                            tBlock->set_material( tMaterial->label() );
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
                    Vector< real > tTRHS;

                    // loop over all elements on block
                    for( Element * tElement : tBlock->elements() )
                    {
                        if( tElement->has_t_matrix() )
                        {
                            mIWG->compute_rhs( tElement, tRHS );
                            tElement->t_matrix()->project( tRHS, tTRHS );
                            mSolverData->assemble_volume_loads( tElement, tTRHS );
                        }
                        else
                        {
                            mIWG->compute_rhs( tElement, tRHS );
                            mSolverData->assemble_volume_loads( tElement, tRHS );
                        }
                    }
                }
            }

            // scale loads with timestep of IWG
            tVolumeLoads *= mIWG->delta_time() ;

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

//-----------------------------------------------------------------------------

        void
        DofManager::seed_dof_values()
        {
            // seeding reads the fields the initialize() pass has already
            // linked; calling it earlier would copy unpopulated data
            BELFEM_ASSERT( mInitializedFlag,
                "seed_dof_values() called before initialize()" );

            mDofData->init_dof_values( true ) ;
        }

//-----------------------------------------------------------------------------

        /**
         * compute the matrices for the projections
         */
        void
        DofManager::initialize_postprocessors()
        {
            mPostprocessorSourceFields.clear() ;
            mPostprocessorTargetFields.clear() ;

            for ( Postprocessor * tPostproc : mPostprocessors )
            {
                if ( tPostproc->projector() != nullptr )
                {
                    tPostproc->projector()->initialize() ;
                    tPostproc->projector()->disconnect_dofs_from_mesh();
                    //tPostproc->projector()->initialize() ;
                }
                for ( const string & tField : tPostproc->source_fields() )
                {
                    mPostprocessorSourceFields.push( tField );
                }
                for ( const string & tField : tPostproc->target_fields() )
                {
                    mPostprocessorTargetFields.push( tField );
                }
                comm_barrier();
            }
            unique( mPostprocessorSourceFields );
            unique( mPostprocessorTargetFields );

        }

//-----------------------------------------------------------------------------

        /**
         * perform the L2 projections for the secondary fields
         */
        void
        DofManager::postprocess()
        {
            this->distribute_fields( mPostprocessorSourceFields );
            mIWG->custom_postprocess();
            for ( Postprocessor * tProjector : mPostprocessors )
            {
                tProjector->run();
                comm_barrier();
            }
            //this->collect_fields( mPostprocessorTargetFields );

            // where the historical stray message was born -- the postprocess
            // exchanges must leave the fabric empty
            comm_drain_check( "DofManager::postprocess" );
        }

//-----------------------------------------------------------------------------

        void
        DofManager::consolidate_dofs()
        {
            for ( Block * tBlock: mBlockData->blocks())
            {
                for ( Element * tElement: tBlock->elements())
                {
                    this->consolidate_dofs( tElement );
                }
            }
            for ( SideSet * tSideSet: mSideSetData->sidesets())
            {
                for ( Element * tElement: tSideSet->elements())
                {
                    this->consolidate_dofs( tElement );
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::consolidate_dofs( Element * aElement )
        {
            // collect local dofs
            Cell< Dof * > tLocalDofs( aElement->number_of_dofs(), nullptr );
            for( uint k=0; k<aElement->number_of_dofs(); ++k )
            {
                tLocalDofs( k ) = aElement->dof( k ) ;
            }

            // count sources and hanging dofs
            index_t tNumFreeDofs = 0 ;
            index_t tNumHangingDofs = 0 ;
            index_t tNumSources = 0 ;
            for( Dof * tDof : tLocalDofs )
            {
                if( tDof->is_hanging() )
                {
                    ++tNumHangingDofs ;
                    tNumSources += tDof->number_of_sources() ;
                }
                else
                {
                    ++tNumFreeDofs ;
                }
            }
            if( tNumHangingDofs == 0 )
            {
                return;
            }

            // collect global dofs
            Cell< Dof * > tGlobalDofs( tNumFreeDofs + tNumSources, nullptr );
            index_t tCount = 0 ;
            for( Dof * tDof : tLocalDofs )
            {
                if( tDof->is_hanging() )
                {
                    for( uint s=0; s<tDof->number_of_sources(); ++s )
                    {
                        tGlobalDofs( tCount++ ) = tDof->source( s );
                    }
                }
                else
                {
                    tGlobalDofs( tCount++ ) = tDof ;
                }
            }
            unique( tGlobalDofs );
            Map< id_t, uint > tIndex ;
            tCount = 0 ;
            for( Dof * tDof : tGlobalDofs )
            {
                tIndex[ tDof->id()] = tCount++ ;
            }


            Matrix< real > tT( tLocalDofs.size(), tGlobalDofs.size(), 0. );

            Cell< Dof * > tSources ;
            Vector< real > tCoefficients ;

            uint tRow = 0 ;
            for( Dof * tDof : tLocalDofs )
            {
                if( ! tDof->is_hanging() )
                {
                    tT( tRow, tIndex[ tDof->id() ] ) = 1.0 ;
                }
                else
                {
                    for( uint k=0; k<tDof->number_of_sources(); ++k )
                    {
                        tT( tRow, tIndex[ tDof->source( k )->id() ] ) = tDof->weight( k );
                    }
                } // end is hanging
                ++tRow;
            }

            aElement->relink_dofs( tGlobalDofs, tLocalDofs, tT );

        }

//-----------------------------------------------------------------------------

        void
        DofManager::link_bearings_with_dofs()
        {

            Vector< id_t > tTable ;

            index_t tCount ;

            if( mCommRank == 0 )
            {
                // only on the master proc, the nodes know which dofs
                // they are connected to. We need to collect and share
                // this information

                // count the memory size for the dof table
                tCount = 3 * mMesh->vertices().size() ;
                for( mesh::Element * tVertex : mMesh->vertices() )
                {
                    tCount += tVertex->node( 0 )->number_of_dofs() ;
                }

                // allocate lookup table
                tTable.set_size( tCount, 0 );
                tCount = 0 ;

                // populate message
                for( mesh::Element * tVertex : mMesh->vertices() )
                {
                    // get the node

                    mesh::Node * tNode = tVertex->node( 0 );

                    tTable( tCount++ ) = tVertex->id() ;
                    tTable( tCount++ ) = tNode->id() ;
                    tTable( tCount++ ) = tNode->number_of_dofs() ;
                    for( index_t k = 0; k<tNode->number_of_dofs(); ++k )
                    {
                        tTable( tCount++ ) = tNode->dof( k )->id() ;
                    }
                }
            }

            broadcast( tTable );
            tCount = 0 ;

            // loop over all entries
            while( tCount < tTable.length() )
            {
                // get the bearing id
                id_t tVertexID = tTable( tCount++ );

                id_t tNodeID   = tTable( tCount++ );

                // get the number of bearings
                index_t tNumDofs = tTable( tCount++ );

                // check if we know this node
                if( mMesh->node_exists( tNodeID ) && tNumDofs > 0 )
                {
                    Bearing * tBearing = mBearingData->bearing( tVertexID );

                    if( tBearing == nullptr )
                    {
                        tCount += tNumDofs ;
                        continue;
                    }
                    if( mMesh->node( tNodeID )->number_of_dofs() > 0 )
                    {
                        // link dofs with bearing
                        tBearing->allocate_dof_container( tNumDofs );
                        for( index_t k=0; k<tNumDofs; ++k )
                        {
                            tBearing->insert_dof( mDofData->dof( tTable( tCount++ )), k );
                        }
                    }
                    else
                    {
                        // jump to next dof
                        tCount += tNumDofs ;
                    }
                }
                else
                {
                    // jump to next dof
                    tCount += tNumDofs ;
                }


            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::create_element_map()
        {
            for ( Block * tBlock : mBlockData->blocks() )
            {
                for ( Element * tElement : tBlock->elements() )
                {
                    mElementMap[ tElement->id() ] = tElement ;
                }
                for ( Element * tElement : tBlock->aura_elements() )
                {
                    mElementMap[ tElement->id() ] = tElement ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        DofManager::disconnect_dofs_from_mesh()
        {
            mDofData->disconnect_dofs_from_mesh() ;
        }

//-----------------------------------------------------------------------------

        void
        DofManager::full_lhs( Vector< real > & aLHS )
        {
            Cell< Dof * > & tDofs = mDofData->dofs();

            aLHS.set_size( tDofs.size() );

            index_t tCount = 0 ;
            for ( Dof * tDof : tDofs )
            {
                aLHS( tCount++ ) = tDof->value() ;
            }
        }

    }

//-----------------------------------------------------------------------------
}
