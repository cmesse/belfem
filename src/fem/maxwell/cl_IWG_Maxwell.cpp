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
#include "matrices/mt_maxwell_phi.hpp"
#include "matrices/mt_maxwell_h.hpp"
#include "matrices/mt_maxwell_background.hpp"
#include "matrices/mt_maxwell_symmetry.hpp"
//#include "matrices/mt_maxwell_interface.hpp"

#include "cl_IWG_Maxwell.hpp"

#include "cl_FEM_Controller.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Kernel.hpp"
#include "assert.hpp"
namespace belfem
{
    namespace fem
    {
        IWG_Maxwell::IWG_Maxwell(
            const maxwell::Formulation aFormulation,
            const ModelDimensionality aDimensionality,
                  const bool aHigherOrder,
                  const bool aUseEnrichment ) :
            IWG_Timestep(
                belfem::IwgType::Maxwell,
                aDimensionality,
                IwgMode::Iterative,
                SymmetryMode::Unsymmetric,
                DofMode::BlockSpecific,
                SideSetDofLinkMode::MasterAndSlave ),
            mFormulation( aFormulation ),
            mNumberOfDimensions( aDimensionality == ModelDimensionality::ThreeD ? 3 : 2 ),
            mHigherOrder( aHigherOrder ),
            mUseEnrichment( aUseEnrichment ),
            mFields( mDofMap, mDofFields, mOtherFields, mHiddenFields, mAllFields )
        {
            // first two used by ghost, third used by gauging. These are
            // FALLBACKS for paths that never reach Controller::set_params
            // ( the hand-built test stacks ), not the production defaults:
            // since 2026-09-01 the controller writes every slot from the
            // deck, and an absent block means 0 for eta and for chi
            // ( opt-in; see fn_FEM_ghost_switch.hpp and the input reference )
            mPenalty = { 4.0 , 1e-3, 1.0e-4 };

            mEnrichSideSets = mUseEnrichment ;

            if( mFormulation == maxwell::Formulation::HPhi )
            {
                mUseEdges = true ;

                // main dofs for h-phi
                mFields.Air.push( "phi");

                mFields.Conductor.push( "edge_h") ;

                if( mHigherOrder )
                {
                    mFields.Conductor.push( "face_h");
                    mEdgeDofMultiplicity = 2 ;
                    mFaceDofMultiplicity = 2 ;
                }
                else
                {
                    mEdgeDofMultiplicity = 1 ;
                    mFaceDofMultiplicity = 0 ;
                }

                mFields.Ferro.push("phi");

                mFields.NonDof.push( "element_rho" );

                mTimeStepMatrices->set_flag( MatrixFlag::M ) ;
                mTimeStepMatrices->set_flag( MatrixFlag::K ) ;
                mTimeStepMatrices->set_flag( MatrixFlag::F ) ;
                mTimeStepMatrices->set_flag( MatrixFlag::dKdX_times_x ) ;
                mTimeStepMatrices->set_flag( MatrixFlag::dMdX_times_x ) ;
                mTimeStepMatrices->set_flag( MatrixFlag::dMdX_times_h ) ;
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::set_currents( Vector< real > & aI )
        {
            if (aI.length() == 0) return ;

            // aI( k ) goes onto abstract node dof k; that ordering is the
            // factory's contract and cannot be rechecked here.
            //
            // The loop is deliberately tolerant of a SHORT container: abstract
            // nodes are collected on the master only ( MaxwellFactory calls
            // set_abstract_nodes inside a rank guard ), while every rank is
            // handed the full current vector. A worker therefore has no dofs to
            // fix and must fix none. An earlier version of this function raised
            // a hard error on that mismatch and aborted every parallel run.
            index_t tCount = 0 ;

            for( Dof * tDof : mAbstractNodeDofs )
            {
                if (tCount < aI.length())
                {
                    BELFEM_ASSERT( tDof != nullptr,
                                   "set_currents: abstract node dof %lu is null",
                                   ( long unsigned int ) tCount );

                    tDof->fix( aI( tCount++ ));
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::init_activation_maps()
        {
            // Call base implementation first
            IWG::init_activation_maps();

            // Thin shell sidesets need GeometryOnly mode
            // (calculators for normal H-field computation, but no DOFs)
            mSideSetActivationModes[ DomainType::ThinShell ] = GroupActivationMode::GeometryOnly;

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::initialize()
        {
            // wait for other procs
            if( mDofMode == DofMode::BlockSpecific )
            {
                proc_t tCommSize = comm_size() ;

                if( tCommSize > 1 )
                {
                    if( mCommRank == 0 )
                    {
                        Vector< id_t > tData( 2 * ( mBlockTypes.size() + mSideSetTypes.size() ) + 2 );
                        uint tCount = 0 ;
                        tData( tCount++ ) = mBlockTypes.size() ;
                        tData( tCount++ ) = mSideSetTypes.size() ;

                        for( auto tPair : mBlockTypes )
                        {
                            tData( tCount++ ) = tPair.first ;
                            tData( tCount++ ) = static_cast< id_t >( tPair.second );
                        }
                        for( auto tPair : mSideSetTypes )
                        {
                            tData( tCount++ ) = tPair.first ;
                            tData( tCount++ ) = static_cast< id_t >( tPair.second );
                        }


                        comm_barrier() ;
                        broadcast( tData );
                        broadcast( mBlockIDs );
                        broadcast( mSideSetIDs );

                    }
                    else
                    {
                        mBlockTypes.clear() ;
                        mSideSetTypes.clear() ;

                        Vector< id_t > tData ;

                        comm_barrier() ;
                        broadcast( tData );
                        broadcast( mBlockIDs );
                        broadcast( mSideSetIDs );

                        index_t tCount = 0 ;
                        index_t tNumBlocks   = tData( tCount++ );
                        index_t tNumSideSets = tData( tCount++ );

                        for( index_t b=0; b<tNumBlocks; ++b )
                        {
                            id_t tID = tData( tCount++ ) ;
                            DomainType tType = static_cast< DomainType >( tData( tCount++ ) ) ;
                            mBlockTypes[ tID ] = tType ;
                        }
                        for( index_t s=0; s<tNumSideSets; ++s )
                        {
                            id_t tID = tData( tCount++ ) ;
                            DomainType tType = static_cast< DomainType >( tData( tCount++ ) ) ;
                            mSideSetTypes[ tID ] = tType ;
                        }
                    }
                }

                // initialize the Dof list
                mFields.initialize( this );

                if( mFormulation == maxwell::Formulation::HPhi )
                {
                    mAbstractDofType = mFields.doftype("phi");
                }

                mFields.collect_block_dofs( mBlockIDs, mBlockTypes, mDofsPerBlock );
                mFields.collect_sideset_dofs( mSideSetIDs, mSideSetTypes, mDofsPerSideSet );
            }
            else
            {
                BELFEM_ERROR( mDofFields.size() > 0, "No dofs assigned for this IWG" );
            }

            // call parent function to set
            this->select_blocks( mBlockIDs );
            this->select_sidesets( mSideSetIDs );

            // call init function from parent
            IWG::initialize();

            // some dandy fix
            mNumberOfRhsDofsPerEdge =
                    mNumberOfRhsDofsPerEdge < mEdgeDofMultiplicity ?
                    mEdgeDofMultiplicity : mNumberOfRhsDofsPerEdge ;

            mNumberOfRhsDofsPerFace =
                    mNumberOfRhsDofsPerFace < mFaceDofMultiplicity ?
                    mFaceDofMultiplicity : mNumberOfRhsDofsPerFace ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_mkf(
            Element * aElement)
        {
            // link calculator with element
            mGroup->calculator()->link( aElement );

            mTimeStepMatrices->reset() ;

            // call the function
            ( *mFunMKF )( mGroup->calculator(), mTimeStepMatrices );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::link_to_group( Group  * aGroup )
        {
            IWG::link_to_group(aGroup);

            mGroup = aGroup ;
            mCalc = aGroup->calculator() ;
            if( aGroup->number_of_elements() == 0 )
            {
                return;

            }

            // a kernel without a Controller has no thermal kernel by
            // construction; assign on both guard arms so a
            // relink of a non-empty group can never leave stale state
            // ( the empty-group return above still skips the write )
            Kernel * tKernel = aGroup->parent()->parent() ;
            mHaveThermal = tKernel->has_controller()
                    && tKernel->controller()->thermal_kernel() != nullptr ;

            // size the timestep matrices ( M, K, f, dMdx/dKdx blocks ) for
            // this group's element size — dropped during the dispatch rewrite,
            // bdf1 dies on the unsized f() otherwise
            uint tNumDofs = mGroup->parent() != nullptr ?
                    mGroup->elements()( 0 )->number_of_local_dofs() : 0 ;
            mTimeStepMatrices->initialize( tNumDofs ) ;

            switch ( aGroup->domain_type() )
            {
                case DomainType::Air :
                case DomainType::Buffer :
                {
                    switch ( aGroup->element_type() )
                    {
                        case ElementType::TRI3 :
                        {
                            mFunMKF = & maxwell::phi_tri3 ;
                            break ;
                        }
                        case ElementType::TET4 :
                        {
                            mFunMKF = & maxwell::phi_tet4 ;
                            break ;
                        }
                        case ElementType::TRI6 :
                        case ElementType::TET10 :
                        {
                            mFunMKF = & maxwell::phi_tri6_tet10 ;
                            break ;
                        }
                        default:
                        {
                            mFunMKF = & maxwell::phi ;
                            break ;
                        }
                    }
                    break ;
                }
                case DomainType::Ferro :
                {
                    mFunMKF = this->algorithm() == SolverAlgorithm::NewtonRaphson ?
                          & maxwell::phi_ferro_newton
                        : & maxwell::phi_ferro_picard ;
                    break ;
                }
                case DomainType::InterfaceFerroAir :
                {
                    // there is no ferro-air interface weak form, in 2d or in 3d:
                    // the nodes are duplicated for visualization only. Under h-phi
                    // the sideset is set Inactive after the hanging-edge pass
                    // ( MaxwellFactory::create_hanging_edges_and_facets ), so an
                    // active group here means the deactivation was bypassed
                    BELFEM_ERROR( false, "Ferro-Air Interfaces must be disabled!");
                    break ;
                }
                case DomainType::Conductor :
                case DomainType::ThinShell :
                {
                    Material * tMat = aGroup->material();

                    // material specifics ( rho law, T source, mu source ) are
                    // resolved inside MaxwellData; only the mu tangent needs a
                    // kernel distinction: constant mu skips the dMdx blocks
                    if ( this->algorithm() == SolverAlgorithm::NewtonRaphson )
                    {
                        mFunMKF = tMat->is_constant( MaterialProperty::mu ) ?
                              & maxwell::h_newton_mu0
                            : & maxwell::h_newton_mu ;
                    }
                    else
                    {
                        mFunMKF = & maxwell::h_picard ;
                    }
                    break ;
                }
                case DomainType::InterfaceCondAir :
                {
                    BELFEM_ERROR( false, "Conductor-Air Interfaces must be disabled!");
                    break ;
                }
                case DomainType::InterfaceCondFerro :
                {
                    // the h-side is coupled to the phi-side by the hanging-edge
                    // condensation, not by a weak form. Under h-phi the sideset
                    // is set Inactive once that pass is done
                    BELFEM_ERROR( false, "Conductor-Ferro Interfaces must be disabled!");
                    break ;
                }
                case DomainType::InterfaceAirCoil :
                {
                    BELFEM_ERROR( false, "Air-Coil Interfaces must be disabled!");
                    break ;
                }
                case DomainType::InterfaceFerroCoil :
                {
                    BELFEM_ERROR( false, "Ferro-Coil Interfaces must be disabled!");
                    break ;
                }
                case DomainType::AirAntiSymmetry :
                case DomainType::BufferAntiSymmetry :
                {
                    BELFEM_ERROR( false, "Anti Symmetry fields must be disabled!");
                    break ;
                }
                case DomainType::AirSymmetry :
                case DomainType::FerroSymmetry :
                case DomainType::BufferSymmetry :
                {
                    if( mesh::dimension( aGroup->element_type() ) == 1 )
                    {
                        mFunMKF = & maxwell::symmetry_phi_2d ;
                    }
                    else if( mesh::dimension( aGroup->element_type() ) == 2 )
                    {
                        mFunMKF = & maxwell::symmetry_phi_3d ;
                    }
                    else
                    {
                        BELFEM_ERROR( false, "Sideset dimension error");
                    }
                    break ;
                }
                case DomainType::FerroAntiSymmetry :
                {
                    BELFEM_ERROR( false, "Anti Symmetry fields must be disabled!");
                    break ;
                }
                case DomainType::ConductorAntiSymmetry :
                {
                    BELFEM_ERROR( false, "Anti Symmetry fields must be disabled!");
                    break ;
                }
                case DomainType::ConductorSymmetry :
                {
                    if( mesh::dimension( aGroup->element_type() ) == 2 )
                    {
                        mFunMKF = & maxwell::h_symmetry_3d ;
                    }
                    else
                    {
                        mFunMKF = & maxwell::h_symmetry_2d ;
                    }
                    break ;
                }
                case DomainType::BackgroundField :
                {
                    mFunMKF = & maxwell::background_phi ;
                    break ;
                }
                case DomainType::Ghost :
                {
                    mFunMKF = & maxwell::h_ghost ;
                    break ;
                }
                case DomainType::LeftCoating :
                case DomainType::RightCoating :
                {
                    // constant mu is enforced in the MaxwellData ctor, so no
                    // mu-kernel distinction is needed here
                    mFunMKF = this->algorithm() == SolverAlgorithm::NewtonRaphson ?
                          & maxwell::h_side_connector_newton
                        : & maxwell::h_side_connector ;
                    break ;
                }
                default:
                {
                    std::cout << to_string(aGroup->domain_type()) << std::endl ;
                    BELFEM_ERROR( false, "Not implemented");
                }
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::create_custom_vectors_and_matrices( Calculator * aCalc )
        {
            ElementType tType = aCalc->group()->type() == GroupType::BLOCK ?
                    aCalc->group()->element_type() : aCalc->group()->master_type() ;

            uint e = mesh::number_of_nedelec_dofs( tType ) ;
            aCalc->create_vector( "nedelec_h", e );

            // number of nodes per element
            uint n = mesh::number_of_nodes( tType );

            // number of dimensions
            uint d = mesh::dimension(  tType );

            aCalc->create_vector( "phi0", n );

            uint tNumDimensions = mesh::dimension( tType );

            aCalc->create_vector("j", tNumDimensions );
            aCalc->create_vector("b", tNumDimensions );
            aCalc->create_vector("h", tNumDimensions );

            aCalc->create_vector( "phi_m", n );
            aCalc->create_vector( "phi_s", n );
            aCalc->create_vector( "theta", n );

            aCalc->create_vector( "normal",  d );

            // the next three are needed for side connector wall elements
            aCalc->create_vector( "tangent",  d );
            aCalc->create_vector( "binomial", d );
            aCalc->create_vector( "Tseam", 4 );

            aCalc->create_vector( "ht", d );
            aCalc->create_vector( "hn", d );
            aCalc->create_vector("hm", tNumDimensions );
            aCalc->create_vector("hs", tNumDimensions );

            // per-point scratch of compute_h_trace ( thin-shell normal recovery )
            aCalc->create_vector( "hk", d );

            // slave twin of "nedelec_h" ( sized for the master above ): the
            // edge dofs of a facet's slave volume, read by
            // Calculator::nedelec_data_slave_h when that volume is an
            // h-conductor
            if (    aCalc->group()->type() != GroupType::BLOCK
                 && aCalc->group()->slave_type() != ElementType::EMPTY
                 && aCalc->group()->slave_type() != ElementType::UNDEFINED )
            {
                aCalc->create_vector( "nedelec_h_s",
                    mesh::number_of_nedelec_dofs( aCalc->group()->slave_type() ) );
            }

            aCalc->create_vector( "bt", d );
            aCalc->create_vector( "bn", d );
            
            aCalc->create_matrix( "Ctj", e, 1 );

            // workspaces for the gauge-tangent drho channel
            // ( h_newton_mu0 / h_newton_mu ). Gq holds G*q, whose row count
            // is set by the element type's G operator ( TET4: 6, TET10: 12,
            // TRI3/PENTA6TS: 3 ) — sized to the largest and resized by the
            // first assignment
            aCalc->create_vector( "Gq", 12 );
            aCalc->create_matrix( "GtGq", e, 1 );

            // workspaces for the rho(|B|,beta) field-derivative tangent
            // ( add_rho_field_tangent in mt_maxwell_h.cpp )
            aCalc->create_vector( "vE", d );
            aCalc->create_vector( "vC", d );
            aCalc->create_matrix( "Rw", e, 1 );

            // workspaces for the signed history contraction of the
            // nonlinear-mu mass tangent ( phi_ferro_newton, h_newton_mu );
            // the row count of B/E on every nonlinear-mu path is the mesh
            // spatial dimension ( equals mesh::dimension( tType ) for bulk
            // and thin-shell layer blocks alike; groups where it differs
            // use constant-mu kernels and never touch these )
            uint tNumSpaceDim = aCalc->mesh()->number_of_dimensions() ;
            aCalc->create_vector( "hcur",  tNumSpaceDim );
            aCalc->create_vector( "hhist", tNumSpaceDim );

            if ( aCalc->group()->domain_type() == DomainType::Ghost )
            {
                aCalc->create_matrix( "K++", e, e );
                aCalc->create_matrix( "K+-", e, e );
                aCalc->create_matrix( "K-+", e, e );
                aCalc->create_matrix( "K--", e, e );

                aCalc->create_matrix( "D+", d, e );
                aCalc->create_matrix( "D-", d, e );
            }

            // todo: move the lines below into calculator class
            if ( ! mUseEnrichment || ! aCalc->group()->has_enrichment() ) return ;

            // number of enriched functions
            uint m = aCalc->group()->enrichment_data( 0 )->function()->number_of_bases() ;

            // node interpolation matrix for bubbles
            aCalc->create_matrix("BubbleN1", 1, n*m );

            // gradient operator matrix for bubbles
            aCalc->create_matrix("BubbleB1", d, n*m );

            // container for n'*B
            aCalc->create_matrix("BubbleNtB1", 1, n*m );

            // node interpolation matrix for bubbles
            aCalc->create_matrix("BubbleN2", 1, n*m );

            // gradient operator matrix for bubbles
            aCalc->create_matrix("BubbleB2", d, n*m );

            // container for n'*B
            aCalc->create_matrix("BubbleNtB2", 1, n*m );

            // help matrix to compute derivatives
            aCalc->create_matrix("BubbleH", d, n*m );

            // other help matrices
            aCalc->create_matrix("BubbleU", n, n*m );

            aCalc->create_matrix("BubbleV", n*m, n );
            aCalc->create_matrix("BubbleW", n*m, n*m );
            aCalc->create_matrix("BubbleInvW", n*m, n*m );


        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::custom_postprocess()
        {
            if ( mCommRank == 0 )
            {
                mMesh->unflag_all_nodes() ;
                for ( mesh::Block * tBlock : mMesh->blocks() )
                {
                    if ( mField->block_exists( tBlock->id() ) )
                    {
                        if ( mField->block( tBlock->id() )->domain_type() == DomainType::Air || mField->block( tBlock->id() )->domain_type() == DomainType::Buffer )
                        {
                            tBlock->flag_nodes() ;
                        }
                    }
                }

                Vector< real > & Bx = mMesh->field_exists( "Bx") ? mMesh->field_data("Bx") : mMesh->create_field( "Bx");
                Vector< real > & Hx = mMesh->field_data("Hx");



                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    if ( tNode->is_flagged() )
                    {
                        Bx( tNode->index() ) = constant::mu0 * Hx( tNode->index() ) ;
                    }
                }

                Vector< real > & By = mMesh->field_exists( "By") ? mMesh->field_data("By") : mMesh->create_field( "By");
                Vector< real > & Hy = mMesh->field_data("Hy");

                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    if ( tNode->is_flagged() )
                    {
                        By( tNode->index() ) = constant::mu0 * Hy( tNode->index() ) ;
                    }
                }

                Vector< real > & Bz = mMesh->field_exists( "Bz") ? mMesh->field_data("Bz") : mMesh->create_field( "Bz");

                if ( mMesh->number_of_dimensions() == 3 )
                {
                    Vector< real > & Hz = mMesh->field_data("Hz");
                    for ( mesh::Node * tNode : mMesh->nodes() )
                    {
                        if ( tNode->is_flagged() )
                        {
                            Bz( tNode->index() ) = constant::mu0 * Hz( tNode->index() ) ;
                        }
                    }
                }
            }

            comm_barrier() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::collect_abstract_node_dofs()
        {
            mAbstractNodeDofs.set_size( mAbstractNodes.size(), nullptr );

            index_t tCount = 0 ;

            for( mesh::Node * tNode : mAbstractNodes )
            {
                // grab dof
                Dof * tDof = reinterpret_cast< Dof * > ( tNode->dof( 0 ) );

                mAbstractNodeDofs( tCount++ ) = tDof ;
            }
        }

//------------------------------------------------------------------------------

    }
}
