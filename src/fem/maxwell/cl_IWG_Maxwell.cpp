//
// Created by christian on 12/16/21.
//

#include "cl_IWG_Maxwell.hpp"
#include "commtools.hpp"
#include "meshtools.hpp"
#include "assert.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

#include "cl_EdgeFunctionFactory.hpp"
#include "cl_MaxwellMaterial.hpp"
#include "cl_FEM_DofManager.hpp"
#include "fn_crossmat.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell::IWG_Maxwell(
                const ElementType aElementType,
                const IwgType aType,
                const IwgMode aMode,
                const SymmetryMode aSymmetryMode,
                const SideSetDofLinkMode aSideSetDofLinkMode,
                const bool         aUseEdges ) :
                IWG_Timestep(
                        aType,
                        aMode,
                        aSymmetryMode,
                        DofMode::BlockSpecific,
                        aSideSetDofLinkMode ),
                mElementType( aElementType ),
                mNumberOfDimensions( mesh::dimension( aElementType ) ),
                mUseEdges( aUseEdges ),
                mFields( mDofMap, mDofFields, mOtherFields, mHiddenFields, mAllFields )
        {
            mNumberOfSpatialDimensions = mNumberOfDimensions ;
            mNumberOfDerivativeDimensions =  mNumberOfDimensions ;

            // create the edge function if this IWG has edge elements
            if( aUseEdges )
            {
                // create the factory
                EdgeFunctionFactory tFactory ;

                // create the function
                mEdgeFunction = tFactory.create_edge_function( mElementType );
            }

            // make sure that this is a backwards implicit scheme
            BELFEM_ERROR( this->theta() == 1.0,
                         "Maxwell IWGs must be Backward Implicit (theta=1)" );

            // set the computation flags for the dof manager
            mComputeJacobianOnBlock   = true ;
            mComputeJacobianOnSideset = true ;
        }

//------------------------------------------------------------------------------

        IWG_Maxwell::~IWG_Maxwell()
        {
            // delete edge function if it was set
            if( mEdgeFunction != nullptr )
            {
                delete mEdgeFunction ;
            }
            if( mEdgeFunctionTS != nullptr )
            {
                delete mEdgeFunctionTS ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::initialize()
        {
            // wait for other procs
            comm_barrier() ;
            if( mDofMode == DofMode::BlockSpecific )
            {
                // create commlist
                proc_t tNumProcs = comm_size();

                Vector< id_t > tAllBlockIDs( mBlockIDs );
                Vector< id_t > tAllSideSetIDs( mSideSetIDs );

                if( tNumProcs > 1 )
                {
                    if( comm_rank() == 0 )
                    {
                        // create commlist
                        Vector< proc_t > tCommList( tNumProcs - 1 );
                        for( proc_t p=1; p<tNumProcs; ++p )
                        {
                            tCommList( p-1 ) = p ;
                        }

                        // send lists to others
                        send_same( tCommList, tAllBlockIDs );
                        send_same( tCommList, tAllSideSetIDs );
                    }
                    else
                    {
                        receive( 0, tAllBlockIDs );
                        receive( 0, tAllSideSetIDs );
                    }
                }

                // initialize the Dof list
                mFields.initialize( this );
                mFields.collect_block_dofs( tAllBlockIDs, mBlockTypes, mDofsPerBlock );
                mFields.collect_sideset_dofs( tAllSideSetIDs, mSideSetTypes, mMagfieldTypeMap, mDofsPerSideSet );
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
        IWG_Maxwell::link_to_group( Group * aGroup )
        {
            // skip empty blocks and sidesets
           if( aGroup->number_of_elements() == 0 ) return;

           // call function form parentcd
           IWG::link_to_group( aGroup );

           // link the edge collector
           this->link_collect_nedelec_data_function( aGroup );

           // link the material law, if this is a superconducting block
           this->link_sc_matrix_function( aGroup ) ;

           // link the material law, if this is a ferro block
           this->link_ferro_matrix_function( aGroup ) ;

           // now link the matrix functions that are IWG specific
           this->link_jacobian_function( aGroup );

           // precompute some data for the integration
           if( aGroup->type() == GroupType::BLOCK )
           {
               mEdgeFunction->precompute( aGroup->integration_points() );
           }
           else if ( aGroup->type() == GroupType::SIDESET || aGroup->type() == GroupType::SHELL )
           {
               // special function for second order thin shell
               if( mEdgeFunctionTS != nullptr )
               {
                   mEdgeFunctionTS->precompute( aGroup->integration_points() );
               }

               // set integration function
               switch( mesh::geometry_type( aGroup->slave_type() ) )
               {
                   case( GeometryType::TRI ) :
                   {
                       mGetInterfaceData = & IWG_Maxwell::integration_data_tri ;
                       break ;
                   }
                   case( GeometryType::TET ) :
                   {
                       mGetInterfaceData = & IWG_Maxwell::integration_data_tet ;
                       break ;
                   }
                   default :
                   {
                       mGetInterfaceData = nullptr ;
                   }
               }
           }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::allocate_work_matrices( Group * aGroup )
        {
            IWG::allocate_work_matrices( aGroup );

            // check matrix function for superconductors
            switch( aGroup->domain_type() )
            {
                case(  DomainType::SuperConductor ) :
                {
                    // allocate container for node data, used for temperature
                    aGroup->work_tau().set_size( mNumberOfNodesPerElement );

                    // allocate container for edge data
                    aGroup->work_nedelec().set_size( mNumberOfRhsEdgeDofsPerElement );

                    // allocate container for all dofs
                    aGroup->work_phi().set_size( mNumberOfDofsPerElement );

                    // stiffness matrix
                    aGroup->work_K().set_size( mNumberOfDofsPerElement,
                                               mNumberOfDofsPerElement, 0.0 );

                    break ;
                }
                case( DomainType::Air ) :
                case( DomainType::FerroMagnetic ) :
                case( DomainType::Coil ) :
                {
                    // allocate container for all dofs, phi or a field
                    aGroup->work_phi().set_size( mNumberOfNodesPerElement );

                    // for A-Vector
                    if( aGroup->parent()->mesh()->number_of_dimensions() == 3 )
                    {
                        aGroup->work_Phi().set_size(  mNumberOfNodesPerElement, 3 );
                    }

                    // stiffness matrix
                    aGroup->work_K().set_size( mNumberOfDofsPerElement,
                                               mNumberOfDofsPerElement, 0.0 );

                    // help matrix for vector
                    aGroup->work_chi().set_size( mNumberOfSpatialDimensions );

                    break ;
                }
                case( DomainType::Cut ) :
                {
                    // stiffness matrix
                    aGroup->work_K().set_size( mNumberOfDofsPerElement,
                                               mNumberOfDofsPerElement, 0.0 );
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                case( DomainType::InterfaceScFm )  :
                case( DomainType::InterfaceFmAir ) :
                {
                    mNumberOfNodesPerElement = mesh::number_of_nodes( aGroup->slave_type() );
                    mNumberOfEdgesPerElement = mesh::number_of_edges( aGroup->master_type() );
                    mNumberOfFacesPerElement = mesh::number_of_faces( aGroup->master_type() );

                    mNumberOfEdgeDofsPerElement =   mNumberOfEdgesPerElement * mEdgeDofMultiplicity
                                                  + mNumberOfFacesPerElement * mFaceDofMultiplicity ;

                    // work vector with dofs
                    aGroup->work_phi().set_size( mNumberOfDofsPerElement );

                    aGroup->work_K().set_size( mNumberOfDofsPerElement,
                                               mNumberOfDofsPerElement, 0.0 );

                    aGroup->work_J().set_size( mNumberOfDimensions, mNumberOfDimensions );
                    aGroup->work_invJ().set_size( mNumberOfDimensions, mNumberOfDimensions );

                    aGroup->work_C().set_size( mNumberOfDimensions, mNumberOfNodesPerElement );


                    // data for nxB or nxN in 2D
                    aGroup->work_sigma().set_size( mNumberOfNodesPerElement );

                    // data for nxE or nxC in 2D
                    aGroup->work_tau().set_size( mNumberOfEdgeDofsPerElement  );

                    // data for nxB or nxN in 3D
                    aGroup->work_Sigma().set_size( mNumberOfDimensions, mNumberOfNodesPerElement );

                    // data for nxE or nxC in 3D
                    aGroup->work_Tau().set_size( mNumberOfDimensions, mNumberOfEdgeDofsPerElement );

                    // data for node coordinates of surface element ( 3D only )
                    aGroup->work_X().set_size( mesh::number_of_nodes( aGroup->element_type() ), mNumberOfDimensions  );

                    // contains the expression N_xi * X for surface
                    aGroup->work_L().set_size( 2, 3 );

                    // for perpendicular coupling term
                    aGroup->work_Chi().set_size( mNumberOfNodesPerElement, mNumberOfEdgeDofsPerElement );
                    break ;
                }
                case ( DomainType::Boundary ) :
                {
                    // we don't need this data for the wave, but for the farfield (magnetic wall)
                    mNumberOfNodesPerElement = mesh::number_of_nodes( aGroup->master_type() );

                    aGroup->work_X().set_size( mNumberOfNodesPerElement, mNumberOfDimensions );
                    aGroup->work_J().set_size( mNumberOfDimensions, mNumberOfDimensions );

                    break ;
                }
                case( DomainType::ThinShell ) :
                {
                    mNumberOfNodesPerElement = mesh::number_of_nodes( aGroup->slave_type() );

                    // container for node coordinates
                    aGroup->work_Xi().set_size( mNumberOfNodesPerElement, mNumberOfDimensions );
                    aGroup->work_Eta().set_size( mNumberOfNodesPerElement, mNumberOfDimensions );

                    // data for node coordinates of surface element
                    aGroup->work_X().set_size( mesh::number_of_nodes( aGroup->element_type() ), mNumberOfDimensions  );

                    // stiffness matrix
                    aGroup->work_K().set_size( mNumberOfDofsPerElement,
                                               mNumberOfDofsPerElement, 0.0 );

                    // work vector for B-Matrix
                    aGroup->work_D().set_size( 1, 2*mNumberOfNodesPerElement );

                    aGroup->work_B().set_size( 1, 2 * mNumberOfNodesPerElement );
                    aGroup->work_C().set_size( 1, mNumberOfDofsPerEdge + 1 );

                    // data for nxB or nxN in 2D
                    aGroup->work_Sigma().set_size( 1,mNumberOfDofsPerElement, 0.0 );

                    // transformation matrix n*B or n*N in 2D
                    aGroup->work_Tau().set_size( 1,mNumberOfDofsPerElement, 0.0 );

                    aGroup->work_phi().set_size( mNumberOfNodesPerElement );
                    aGroup->work_psi().set_size( mNumberOfNodesPerElement );

                    // the geometry jacobian for the master and slave elements
                    aGroup->work_J().set_size( mNumberOfDimensions, mNumberOfDimensions );

                    // help matrix for mass
                    uint tOrder =  mesh::interpolation_order_numeric( aGroup->element_type() ) ;
                    uint tN = tOrder == 2 ? 6 : 2 ;

                    aGroup->work_M().set_size( tN, tN, 0.0 );
                    aGroup->work_L().set_size( tN, tN, 0.0 );

                    // data for edge dofs
                    aGroup->work_sigma().set_size( tN , 0.0 ) ;

                    // todo: add face dofs for 3D
                    if( mNumberOfDofsPerFace > 0 && mMesh->number_of_dimensions() == 3 )
                    {
                        aGroup->work_tau().set_size( 2  * mFaceDofMultiplicity * mNumberOfDofsPerFace, 0.0 );
                    }

                    // data for temperatures
                    //aGroup->work_tau().set_size( this->number_of_ghost_sidesets(), 0.0 );

                    aGroup->work_nedelec().set_size( mNumberOfDofsPerElement );

                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false,
                                 "Unknown Domain type for %s %lu",
                                 aGroup->type() == GroupType::BLOCK ? "block" : "sideset",
                                 ( long unsigned int ) aGroup->id() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false, "Invalid call to IWG_Maxwell::compute_jacobian_and_rhs");
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::save( const string & aPath )
        {
#ifdef BELFEM_HDF5

            if( mField->is_master() )
            {

                // create a new HDF5 file
                HDF5 tFile( aPath, FileMode::NEW );

                tFile.create_group("MeshInfo");
                tFile.save_data( "numNodes", mMesh->number_of_nodes() );
                tFile.save_data( "numEdges", mMesh->number_of_edges() );
                tFile.save_data( "numElements", mMesh->number_of_elements() );
                tFile.close_active_group();

                tFile.create_group("TimeInfo");
                tFile.save_data( "timestamp", mMesh->time_stamp() );
                tFile.save_data( "timestep", mMesh->time_step() );
                tFile.close_active_group();

                tFile.create_group("Fields");

                for( string tDofLabel : mAllFields )
                {
                    tFile.save_data( tDofLabel, mMesh->field_data( tDofLabel ) );
                }

                // check if temperature exists
                if( mMesh->field_exists("T") )
                {
                    tFile.save_data( "T", mMesh->field_data( "T" ) );
                }

                tFile.close_active_group() ;

                // also save matrix
                reinterpret_cast< DofManager * > ( mField )->save_system( tFile );

                tFile.close() ;
            }
#endif
        }

//------------------------------------------------------------------------------

        int
        IWG_Maxwell::load( const string & aPath, Mesh * aMesh )
        {
#ifdef BELFEM_HDF5

            uint tFlag = 1 ;

            if( mField->is_master() )
            {
                // load the HDF5 file
                HDF5 tFile( aPath, FileMode::OPEN_RDONLY );

                tFile.select_group( "MeshInfo" );
                index_t tNumNodes ;
                index_t tNumEdges ;
                index_t tNumElements ;
                tFile.load_data( "numNodes", tNumNodes );
                tFile.load_data( "numEdges", tNumEdges );
                tFile.load_data( "numElements", tNumElements );
                tFile.close_active_group() ;

                // check if this is the correct file
                bool tFileIsOK = tNumNodes == mMesh->number_of_nodes() ;
                tFileIsOK = tFileIsOK && tNumEdges == mMesh->number_of_edges() ;
                tFileIsOK = tFileIsOK && tNumElements == mMesh->number_of_elements() ;

                // check if all datasets exist
                hid_t tGroup = tFile.select_group( "Fields" );
                for ( string tDofLabel: mDofFields )
                {
                    tFileIsOK = tFileIsOK && hdf5::dataset_exists( tGroup, tDofLabel );
                }
                tFile.close_active_group() ;

                // send flag to other procs
                tFlag = tFileIsOK ? 1 : 0 ;
                Vector< uint > tData( mField->parent()->comm_table().length(),
                                      tFlag );
                comm_barrier() ;
                send( mField->parent()->comm_table(), tData );
                comm_barrier() ;

                if( tFileIsOK )
                {
                    tFile.select_group( "TimeInfo" );
                    tFile.load_data( "timestamp", aMesh->time_stamp());
                    tFile.load_data( "timestep", aMesh->time_step());

                    tFile.close_active_group();

                    tFile.select_group( "Fields" );

                    for ( string tDofLabel: mAllFields )
                    {
                        tFile.load_data( tDofLabel, aMesh->field_data( tDofLabel ) );
                    }

                    // check if temperature exists
                    if ( aMesh->field_exists( "T" ))
                    {
                        tFile.load_data( "T", aMesh->field_data( "T" ));
                    }

                    tFile.close_active_group();
                    tFile.close();

                    comm_barrier();
                    const Vector< proc_t > & tCommTable = mField->parent()->comm_table();

                    Vector< real > tTime( tCommTable.length(), aMesh->time_stamp());
                    Vector< uint > tTimeCount( tCommTable.length(), aMesh->time_step());
                    send( tCommTable, tTime );
                    send( tCommTable, tTimeCount );
                }
            }
            else
            {
                comm_barrier() ;
                receive( mField->parent()->master(), tFlag );
                comm_barrier() ;

                // check if file is OK
                if( tFlag == 1 )
                {
                    comm_barrier();
                    receive( mField->parent()->master(), aMesh->time_stamp());
                    receive( mField->parent()->master(), aMesh->time_step());
                }
            }

            // make sure that fields are synchronized
            if( gComm.size() > 1 && tFlag == 1)
            {
                if ( aMesh->field_exists( "T" ) )
                {
                    Cell< string > tDofFields = mAllFields ;
                    tDofFields.push( "T" );
                    mField->distribute_fields( tDofFields );
                }
                else
                {
                    mField->distribute_fields( mAllFields );
                }
            }



            // return error code
            if( tFlag == 1 )
            {
                return  0 ;
            }
            else
            {
                return 1 ;
            }
#else
            BELFEM_ERROR( false, "Trying to load an hdf5 file, but we are not linked against HDF5");
#endif
        }

//------------------------------------------------------------------------------
// protected
//------------------------------------------------------------------------------

        void
        IWG_Maxwell::link_jacobian_function( Group * aGroup )
        {
            BELFEM_ERROR( false, "Invalid call to IWG_Maxwell::link_jacobian_function");
        }

//------------------------------------------------------------------------------
// private
//------------------------------------------------------------------------------

        void
        IWG_Maxwell::link_collect_nedelec_data_function( Group * aGroup )
        {
            switch( mesh::interpolation_order( aGroup->element_type() ) )
            {
                case( InterpolationOrder::LINEAR ) :
                {
                    mCollectNedelecData = & IWG_Maxwell::collect_nedelec_data_linear ;
                    break ;
                }
                case( InterpolationOrder::QUADRATIC ) :
                {
                    mCollectNedelecData = & IWG_Maxwell::collect_nedelec_data_qaudratic ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false,
                                 "Invalid interpolation order for group %lu",
                                 ( long unsigned int ) aGroup->id() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::link_sc_matrix_function( Group * aGroup )
        {
            if (  ( aGroup->domain_type() == DomainType::SuperConductor && aGroup->type() == GroupType::BLOCK ) )
            //    || ( aGroup->domain_type() == DomainType::ThinShell && aGroup->type() == GroupType::SHELL ) )
            {
                // get material of block
                const Material * tMat = aGroup->material() ;

                BELFEM_ASSERT(
                        tMat != nullptr,
                        "no material associated with block %lu", aGroup->id() );

                // make sure that permeability law is constant
                BELFEM_ERROR(
                        tMat->permeability_law() == belfem::PermeabilityLaw::Constant,
                        "superconductor material %s requires constant permeability",
                        tMat->label().c_str() );

                if( aGroup->element_type() == ElementType::TRI3 || aGroup->element_type() == ElementType::TET4 )
                {
                    switch ( tMat->resistivity_law() )
                    {
                        case( belfem::ResistivityLaw::Constant ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_const_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_t_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependB ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_b_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependBT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_bt_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJ ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_j_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_jt_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJB ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_jb_first_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJBT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_jbt_first_order ;
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid material law for block %lu.",
                                         ( long unsigned int ) aGroup->id() );
                        }
                    }
                }
                else if( aGroup->element_type() == ElementType::TRI6 || aGroup->element_type() == ElementType::TET10 )
                {
                    switch ( tMat->resistivity_law() )
                    {
                        case( belfem::ResistivityLaw::Constant ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_const_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_t_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependB ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_b_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependBT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_bt_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJ ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_j_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_jt_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJB ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_jb_higher_order ;
                            break ;
                        }
                        case( belfem::ResistivityLaw::DependJBT ) :
                        {
                            mComputeScMatrices = & IWG_Maxwell::sc_matrices_jbt_higher_order ;
                            break ;
                        }

                        default :
                        {
                            BELFEM_ERROR( false, "Invalid material law for block %lu.",
                                         ( long unsigned int ) aGroup->id() );
                        }
                    }
                }
                else
                {
                    mComputeScMatrices = nullptr ;
                }

            }
            else
            {
                mComputeScMatrices = nullptr ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::link_ferro_matrix_function( Group * aGroup )
        {
            if( aGroup->domain_type() == DomainType::FerroMagnetic && aGroup->type() == GroupType::BLOCK )
            {
                switch( aGroup->element_type() )
                {
                    case( ElementType::TRI3 ) :
                    {
                        mComputeFerroMatrix = & IWG_Maxwell::ferro_stiffness_tri3 ;
                        break ;
                    }
                    case( ElementType::TRI6 ) :
                    {
                        mComputeFerroMatrix = & IWG_Maxwell::ferro_stiffness_tri6 ;
                        break ;
                    }
                    case( ElementType::TET4 ) :
                    {
                        mComputeFerroMatrix = & IWG_Maxwell::ferro_stiffness_tet4 ;
                        break ;
                    }
                    case( ElementType::TET10 ) :
                    {
                        mComputeFerroMatrix = & IWG_Maxwell::ferro_stiffness_tet10 ;
                        break ;
                    }
                    default :
                    {
                        mComputeFerroMatrix = nullptr ;
                    }
                }
            }
            else
            {
                mComputeFerroMatrix = nullptr ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_const_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE  ;
            }

            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C();

            // compute contribution for stiffness matrix
            aK = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J() * tMat->rho_el() )
                  * trans( tC ) * tC ;

            // scale mass matrix
            aM *= constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_j_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE  ;
            }

            // scale mass matrix
            aM *= constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;

            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C( );

            // compute stiffness matrix
            aK = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J()
                    * tMat->rho_el( norm( tC * mGroup->work_nedelec() ) ) )
                    * trans( tC ) * tC ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_jt_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );


            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C();

            real tK = 0.0 ;
            real tJ = norm( tC * mGroup->work_nedelec() ) ; // current density

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE ;

                // compute contribution for stiffness matrix
                tK += tW( k ) * tMat->rho_el(
                        tJ,             // current density
                        dot( mGroup->n( k ), mGroup->work_tau() ) ); // temperature
            }

            // compute stiffness matrix
            aK = ( tK * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

            // scale mass matrix
            aM *=  constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_jb_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // grab material
            const Material * tMat = mGroup->material() ;

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // permeability constant
            const real tMu = constant::mu0 * tMat->mu_r();

            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C( 0 );

            // help parameter for stiffness
            real tK = 0.0 ;

            // compute current
            real tJ = norm( tC * mGroup->work_nedelec() ) ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE ;

                // compute contribution for stiffness matrix
                tK += tW( k ) * tMat->rho_el(
                        tJ,   //current density
                        0, // temperature
                        tMu * norm( tE * mGroup->work_nedelec() ) ) ;
            }

            // compute stiffness matrix
            aK = ( tK * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

            // scale mass matrix
            aM *= tMu * mEdgeFunction->abs_det_J() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_jbt_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            // permeability constant
            const real tMu = constant::mu0 * tMat->mu_r();


            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C();

            // current density
            real tJ = norm( tC * mGroup->work_nedelec() );

            // help parameter for stiffness
            real tK = 0.0 ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE  ;

                // compute contribution for stiffness matrix
                tK += tW( k ) * tMat->rho_el(
                                tJ,             // current density
                                dot( mGroup->n( k ), mGroup->work_tau() ),   // temperature
                                tMu * norm( tE * mGroup->work_nedelec() ) )   ;   // magnetic flux density
            }

            // compute stiffness matrix
            aK = ( tK * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

            // scale stiffness matrix
            aM *=  tMu * mEdgeFunction->abs_det_J() ;
        }
//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_t_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C();

            // help constant for stiffness
            real tK = 0.0 ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE ;

                // compute contribution for stiffness matrix
                tK += tW( k ) * tMat->rho_el(
                        BELFEM_QUIET_NAN,             // current density
                        dot( mGroup->n( k ), mGroup->work_tau() ) ); // temperature
            }

            // compute stiffness matrix
            aK = ( tK * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

            // scale mass matrix
            aM *=  constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
        }

//-----------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_b_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C();

            // help constant for stiffness
            real tK = 0.0 ;

            real tMu = constant::mu0 * tMat->mu_r() ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE ;

                // compute contribution for stiffness matrix
                tK += tW( k ) * tMat->rho_el(
                        BELFEM_QUIET_NAN,             // current density
                        BELFEM_QUIET_NAN, // temperature
                        tMu * norm( tE * mGroup->work_nedelec() )
                ); // temperature
            }

            // compute stiffness matrix
            aK = ( tK * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

            // scale mass matrix
            aM *=  tMu * mEdgeFunction->abs_det_J() ;
        }

//-----------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_bt_first_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // get operator for curl function
            const Matrix< real > & tC = mEdgeFunction->C();

            // help constant for stiffness
            real tK = 0.0 ;

            real tMu = constant::mu0 * tMat->mu_r() ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get operator for edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // compute contribution for mass matrix
                aM += tW( k ) * trans( tE ) * tE ;

                // compute contribution for stiffness matrix
                tK += tW( k ) * tMat->rho_el(
                        BELFEM_QUIET_NAN,             // current density
                        dot( mGroup->n( k ), mGroup->work_tau() ),
                        tMu * norm( tE * mGroup->work_nedelec() )
                        ); // temperature
            }

            // compute stiffness matrix
            aK = ( tK * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

            // scale mass matrix
            aM *=  tMu * mEdgeFunction->abs_det_J() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_const_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {

                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tC ) * tC ;

                }

                // scale mass matrix
                aK *= tMat->rho_el() ;

                // scale stiffness matrix
                aM *= constant::mu0 * tMat->mu_r();
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += tW( k ) * trans( tC ) * tC ;
                }

                // scale mass matrix
                aK *= tMat->rho_el() * mEdgeFunction->abs_det_J() ;

                // scale stiffness matrix
                aM *= constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_j_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE  ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el( norm( tC * mGroup->work_nedelec() ) ) )
                            * trans( tC ) * tC ;
                }

                // scale mass matrix
                aM *= constant::mu0 * tMat->mu_r();
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE  ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * tMat->rho_el( norm( tC * mGroup->work_nedelec() ) ) )
                            * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *= constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_jt_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el(
                                    norm( tC * mGroup->work_nedelec() ),             // current density
                                    dot( mGroup->n( k ), mGroup->work_tau() ) ) ) // temperature
                                    * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  constant::mu0 * tMat->mu_r();
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k )
                            * tMat->rho_el(
                                    norm( tC * mGroup->work_nedelec() ),             // current density
                                    dot( mGroup->n( k ), mGroup->work_tau() ) ) ) // temperature
                                            * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_jb_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // permeability constant
            const real tMu = constant::mu0 * tMat->mu_r();

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE  ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el(
                                    norm( tC * mGroup->work_nedelec() ),   //current density
                                    0, // temperature
                                    tMu * norm( tE * mGroup->work_nedelec() ) ) ) // magnetic flux density
                                    * trans( tC )  * tC  ;
                }

                // scale mass matrix
                aM *= tMu ;
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE  ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * tMat->rho_el(
                                    norm( tC * mGroup->work_nedelec() ),   //current density
                                    0, // temperature
                                    tMu * norm( tE * mGroup->work_nedelec() ) ) ) // magnetic flux density
                                            * trans( tC )  * tC  ;
                }

                aM *= tMu * mEdgeFunction->abs_det_J() ;
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_jbt_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            // permeability constant
            const real tMu = constant::mu0 * tMat->mu_r();

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el(
                                    norm( tC * mGroup->work_nedelec() ),             // current density
                                    dot( mGroup->n( k ), mGroup->work_tau() ),   // temperature
                                    tMu * norm( tE * mGroup->work_nedelec() ) )  )     // magnetic flux density
                                    *  trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  tMu ;
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k )
                            * tMat->rho_el(
                                    norm( tC * mGroup->work_nedelec() ),             // current density
                                    dot( mGroup->n( k ), mGroup->work_tau() ),   // temperature
                                    tMu * norm( tE * mGroup->work_nedelec() ) )  )     // magnetic flux density
                           *  trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  tMu * mEdgeFunction->abs_det_J() ;

                // scale stiffness matrix
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_t_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el( BELFEM_QUIET_NAN,             // current density
                                    dot( mGroup->n( k ), mGroup->work_tau() ) ) ) // temperature
                                            * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *= constant::mu0 * tMat->mu_r();
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k )
                            * tMat->rho_el(
                                    BELFEM_QUIET_NAN,             // current density
                                    dot( mGroup->n( k ), mGroup->work_tau() ) ) ) // temperature
                                            * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *= constant::mu0 * tMat->mu_r() * mEdgeFunction->abs_det_J() ;
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }
//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_b_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            real tMu = constant::mu0 * tMat->mu_r();

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el( BELFEM_QUIET_NAN,             // current density
                                            BELFEM_QUIET_NAN, // temperature
                                            tMu * norm( tE * mGroup->work_nedelec()  ) // b-field
                    ) )
                          * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  constant::mu0 * tMat->mu_r();
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k )
                            * tMat->rho_el(
                            BELFEM_QUIET_NAN,             // current density
                            BELFEM_QUIET_NAN,
                            tMu * norm( tE * mGroup->work_nedelec() ) ) ) // temperature
                          * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  tMu * mEdgeFunction->abs_det_J() ;
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::sc_matrices_bt_higher_order(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // grab temperature data
            this->collect_node_data( aElement, "T", mGroup->work_tau() );

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            real tMu = constant::mu0 * tMat->mu_r();

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * tMat->rho_el( BELFEM_QUIET_NAN,             // current density
                                            dot( mGroup->n( k ), mGroup->work_tau() ), // temperature
                                            tMu * norm( tE * mGroup->work_nedelec()  ) // b-field
                                            ) )
                          * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  constant::mu0 * tMat->mu_r();
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get operator for edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get operator for curl function
                    const Matrix< real > & tC = mEdgeFunction->C( k );

                    // compute contribution for mass matrix
                    aM += tW( k ) * trans( tE ) * tE ;

                    // compute contribution for stiffness matrix
                    aK += ( tW( k )
                            * tMat->rho_el(
                            BELFEM_QUIET_NAN,             // current density
                            dot( mGroup->n( k ), mGroup->work_tau() ),
                            tMu * norm( tE * mGroup->work_nedelec() ) ) ) // temperature
                          * trans( tC ) * tC  ;
                }

                // scale mass matrix
                aM *=  tMu * mEdgeFunction->abs_det_J() ;
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::ferro_stiffness_tri3(
                Element * aElement,
                Matrix< real > & aK )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // get curl matrix ( constant for this element )
            const Matrix< real > & tC = mEdgeFunction->CA(  );

            // grab node data
            this->collect_node_data( aElement,
                                        "az",
                                        mGroup->work_phi() );

            aK = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J()
                    * mMaterial->nu_s( norm( tC * mGroup->work_phi() ) ) )
                 * trans( tC ) * tC ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::ferro_stiffness_tri6(
                Element * aElement,
                Matrix< real > & aK )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // grab node data
            Vector< real > & tAz = mGroup->work_phi() ;

            this->collect_node_data( aElement, "az", tAz );

            // reset matrices
            aK.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get curl matrix
                    const Matrix< real > & tC = mEdgeFunction->CA( k );

                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * mMaterial->nu_s( norm( tC * tAz ) ) )
                            * trans( tC ) * tC ;
                }
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get curl matrix
                    const Matrix< real > & tC = mEdgeFunction->CA( k );

                    aK += ( tW( k ) * mMaterial->nu_s( norm( tC * tAz ) ) )
                                    * trans( tC ) * tC ;
                }
                aK *= mEdgeFunction->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::ferro_stiffness_tet4(
                Element * aElement, Matrix< real > & aK )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // get curl matrix ( constant for this elemet )
            const Matrix< real > & tC = mEdgeFunction->CA(  );

            // grab node data
            Matrix< real > & tA = mGroup->work_Phi() ;

            // populate with data from mesh
            this->collect_node_data( aElement, mFields.Ferro, tA );

            // compute the B-Field
            const Vector< real > & tB = this->curl_a( 0, tA );

            aK = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J()
                    * mMaterial->nu_s( norm( tB ) ) )  * trans( tC ) * tC ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::ferro_stiffness_tet10(
                Element * aElement, Matrix< real > & aK )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // reset matrix
            aK.fill( 0.0 );

            // grab node data
            Matrix< real > & tA = mGroup->work_Phi() ;

            // populate with data from mesh
            this->collect_node_data( aElement, mFields.Ferro, tA );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            if( aElement->element()->is_curved() )
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get curl matrix
                    const Matrix< real > & tC = mEdgeFunction->CA( k );

                    // compute the B-Field
                    const Vector< real > & tB = this->curl_a( k, tA );

                    aK += ( tW( k ) * mEdgeFunction->abs_det_J()
                            * mMaterial->nu_s( norm( tB ) ) )
                            * trans( tC ) * tC ;
                }
            }
            else
            {
                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get curl matrix
                    const Matrix< real > & tC = mEdgeFunction->CA( k );

                    // compute the B-Field
                    const Vector< real > & tB = this->curl_a( k, tA );

                    aK += ( tW( k ) * mMaterial->nu_s( norm( tB ) ) )
                                    * trans( tC ) * tC ;
                }
                aK *= mEdgeFunction->abs_det_J();
            }
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_jacobian_and_rhs_air_phi_linear(
            Element        * aElement,
            Matrix< real > & aJacobian,
            Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, true, false );

            // get gradient operator matrix ( constant for this element )
            const Matrix< real > & tB = mEdgeFunction->B(  );

            // grab node data from last timestep
            this->collect_node_data( aElement,
                                     "phi0",
                                     mGroup->work_phi() );

            aJacobian = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J()
                    * constant::mu0 ) * trans( tB ) * tB ;

            aRHS = aJacobian * mGroup->work_phi() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_jacobian_and_rhs_air_phi_higher_order(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, true, false );


            // grab node data from last timestep
            this->collect_node_data( aElement,
                                     "phi0",
                                     mGroup->work_phi() );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            aJacobian.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get gradient operator matrix
                const Matrix< real > & tB = mEdgeFunction->B( k );

                // integrate mass matrix
                aJacobian += ( tW( k ) * mEdgeFunction->abs_det_J() ) * trans( tB ) * tB ;
            }

            aJacobian *= constant::mu0 ;
#
            // compute right hand side
            aRHS = aJacobian * mGroup->work_phi() ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell::collect_q0_hphi_2d( Element * aElement )
        {
            // old dof vector
            Vector< real > & aQ0      = mGroup->work_phi() ;

            this->collect_nedelec_data( aElement->master(), NedelecField::H0, aQ0 );
            uint tCount = mNumberOfEdgeDofsPerElement ;
            this->collect_node_data( aElement->slave(), "phi0", aQ0, tCount );
            this->collect_lambda_data( aElement, "lambda0", aQ0( tCount++ ) );

            BELFEM_ASSERT( tCount == mNumberOfDofsPerElement, "number of dofs does not match" );

            return aQ0 ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell::collect_q0_hphi_3d( Element * aElement )
        {
            // old dof vector
            Vector< real > & aQ0      = mGroup->work_phi() ;


            this->collect_nedelec_data( aElement->master(), NedelecField::H0, aQ0 );
            uint tCount = mNumberOfEdgeDofsPerElement ;
            this->collect_node_data( aElement->slave(), "phi0", aQ0, tCount );

            this->collect_lambda_data( aElement, "lambda_x0", aQ0( tCount++ ) );
            this->collect_lambda_data( aElement, "lambda_y0", aQ0( tCount++ ) );
            this->collect_lambda_data( aElement, "lambda_z0", aQ0( tCount++ ) );

            BELFEM_ASSERT( tCount == mNumberOfDofsPerElement, "number of dofs does not match" );

            return aQ0 ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell::collect_q0_ha_2d( Element * aElement )
        {
            // old dof vector
            Vector< real > & aQ0      = mGroup->work_phi() ;


            this->collect_nedelec_data( aElement->master(), NedelecField::H0, aQ0 );
            uint tCount = mNumberOfEdgeDofsPerElement ;
            this->collect_node_data( aElement->slave(), "az0", aQ0, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofsPerElement, "number of dofs does not match" );

            return aQ0 ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell::collect_q0_ha_3d( Element * aElement )
        {
            // old dof vector
            Vector< real > & aQ0      = mGroup->work_phi() ;

            this->collect_nedelec_data( aElement->master(), NedelecField::H0, aQ0 );
            uint tCount = mNumberOfEdgeDofsPerElement ;
            this->collect_node_data( aElement->slave(), "ax0", aQ0, tCount );
            this->collect_node_data( aElement->slave(), "ay0", aQ0, tCount );
            this->collect_node_data( aElement->slave(), "az0", aQ0, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofsPerElement, "number of dofs does not match" );

            return aQ0 ;
        }

//------------------------------------------------------------------------------


        const Vector< real > &
        IWG_Maxwell::collect_q0_aphi_2d( Element * aElement )
        {
            Vector< real > & aQ0    = mGroup->work_phi() ;

            uint tCount = 0 ;

            this->collect_node_data( aElement->master(), "az0", aQ0, tCount );
            this->collect_node_data( aElement->slave(), "phi0", aQ0, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofsPerElement, "number of dofs does not match" );

            return aQ0 ;
        }
//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell::collect_q0_aphi_3d( Element * aElement )
        {
            Vector< real > & aQ0    = mGroup->work_phi() ;

            uint tCount = 0 ;
            this->collect_node_data( aElement->master(), "ax0", aQ0, tCount );
            this->collect_node_data( aElement->master(), "ay0", aQ0, tCount );
            this->collect_node_data( aElement->master(), "az0", aQ0, tCount );
            this->collect_node_data( aElement->slave(), "phi0", aQ0, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofsPerElement, "number of dofs does not match" );

            return aQ0 ;
        }


//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_interface_ha_tri3(
                Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK )
        {
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tExn    = mGroup->work_tau() ;

            // compute the normal
            const Vector< real > & tn = this->normal_straight_2d( aElement );

            // reset matrix
            aK.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tE = mEdgeFunction->E( k );
                tExn( 0 ) = tn( 1 ) * tE( 0, 0 ) - tn( 0 ) * tE( 1, 0 );
                tExn( 1 ) = tn( 1 ) * tE( 0, 1 ) - tn( 0 ) * tE( 1, 1 );
                tExn( 2 ) = tn( 1 ) * tE( 0, 2 ) - tn( 0 ) * tE( 1, 2 );

                const Vector< real > & tN = tNodeFunction->phi( k );

                aK( 3, 0 ) += tN( 0 ) * tExn( 0 ) * tW( k );
                aK( 4, 0 ) += tN( 1 ) * tExn( 0 ) * tW( k );
                aK( 5, 0 ) += tN( 2 ) * tExn( 0 ) * tW( k );
                aK( 3, 1 ) += tN( 0 ) * tExn( 1 ) * tW( k );
                aK( 4, 1 ) += tN( 1 ) * tExn( 1 ) * tW( k );
                aK( 5, 1 ) += tN( 2 ) * tExn( 1 ) * tW( k );
                aK( 3, 2 ) += tN( 0 ) * tExn( 2 ) * tW( k );
                aK( 4, 2 ) += tN( 1 ) * tExn( 2 ) * tW( k );
                aK( 5, 2 ) += tN( 2 ) * tExn( 2 ) * tW( k );
            }

            aK *= mGroup->work_det_J() ;
            aM = -trans( aK );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_interface_ha_tri6(
                Element        * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {

            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            aK.fill( 0.0 );

            if( aElement->element()->is_curved() )
            {
                // grab node coords for normal vector
                aElement->get_node_coors( mGroup->work_X() );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // compute the normal
                    const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                    real tOmega = tW( k ) * mGroup->work_det_J() ;

                    // get the edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get the node function
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // assemble stiffness matrix
                    for ( uint j = 0; j<8; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            aK( i+8, j ) += tOmega * tN( i ) *
                                    ( tE( 0, j ) * tn( 1 )
                                    - tE( 1, j ) * tn( 0 ) );
                        }
                    }
                }
            }
            else
            {
                // compute the normal
                const Vector< real > & tn = this->normal_straight_2d( aElement );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    real tOmega = tW( k ) * mGroup->work_det_J() ;

                    // get the edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // get the node function
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // assemble stiffness matrix
                    for ( uint j = 0; j<8; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            aK( i+8, j ) += tOmega * tN( i ) *
                                    ( tE( 0, j ) * tn( 1 )
                                    - tE( 1, j ) * tn( 0 ) );
                        }
                    }
                }
            }

            aM = -trans( aK );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_interface_ha_tet4(
                Element        * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Matrix< real > & tnxE    = mGroup->work_Tau() ;

            // compute the normal
            const Vector< real > & tn = this->normal_straight_3d( aElement );

            // reset matrix
            aK.fill( 0.0 );

            uint I ;
            real tOmega ;

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // edge function
                const Matrix< real > & tE = mEdgeFunction->E( k );

                // node function
                const Vector< real > & tN = tNodeFunction->phi( k );

                // compute n x E (note: must multiply with -1 later, since nxE+Exn=0)
                crossmat( tn, tE, tnxE );

                // integration increment
                tOmega = tW( k ) * mGroup->work_det_J() ;

                for( uint j=0; j<6; ++j )
                {
                    I = 6;
                    for ( uint l = 0; l < 3; ++l )
                    {
                        for( uint i=0; i<4; ++i )
                        {
                            aK( I++, j ) += tOmega * tN( i ) * tnxE( l, j );
                        }
                    }
                }
            }

            aM = trans( aK );

            // nxE -> Exn
            aK *= -1. ;
        }
//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_interface_ha_tet10(
                Element        * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // collect the node coords for the surface element
            aElement->get_node_coors( mGroup->work_X() );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Matrix< real > & tnxE    = mGroup->work_Tau() ;

            // reset matrix
            aK.fill( 0.0 );

            real tOmega ;
            uint I ;

            if( aElement->element()->is_curved() )
            {
                // grab node coords for normal vector
                aElement->get_node_coors( mGroup->work_X() );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get the normal
                    const Vector < real > & tn = this->normal_curved_3d( aElement, k );

                    // get the edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // compute tnxE
                    crossmat( tn, tE,  tnxE );

                    // compute the scaling parameter
                    tOmega = tW( k ) * mGroup->work_det_J() ;

                    // get the node function
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // assemble matrix
                    for( uint j = 0; j<20; ++j )
                    {
                        I = 20 ;
                        for( uint d=0; d<3; ++d )
                        {
                            for( uint i=0; i<10; ++i )
                            {
                                aK( I++, j ) += tOmega * tnxE( d, j ) * tN( i );
                            }
                        }
                    }
                }
            }
            else
            {
                // get the normal
                const Vector < real > & tn = this->normal_straight_3d( aElement );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get the edge function
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // compute tnxE
                    crossmat( tn, tE,  tnxE );

                    // compute the scaling parameter ( negative so that nxE->Exn )
                    tOmega = tW( k ) * mGroup->work_det_J() ;

                    // get the node function
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // assemble matrix
                    for( uint j = 0; j<20; ++j )
                    {
                        I = 20 ;
                        for( uint d=0; d<3; ++d )
                        {
                            for( uint i=0; i<10; ++i )
                            {
                                aK( I++, j ) += tOmega * tnxE( d, j ) * tN( i );
                            }
                        }
                    }
                }
            }

            aM = trans( aK );

            // nxE -> Exn
            aK *= -1. ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::shift_fields()
        {
            // loop over all defined dofs
            for( string tDof : mDofFields )
            {
               // create name for old dof
                string tOldDof = tDof + "0" ;

                // check if field for old dofs exists
                if( mMesh->field_exists( tOldDof ) )
                {
                    // shift fields
                    mMesh->field_data( tOldDof ) = mMesh->field_data(tDof );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::update_group_types()
        {
            // loop over all blocks
            for( id_t tID : mBlockIDs )
            {
                if( mField->block_exists( tID ) )
                {
                    mField->block( tID )->set_domain_type( mBlockTypes( tID ) ) ;
                }
            }

            for( id_t tID : mMesh->ghost_block_ids() )
            {
                if( mField->block_exists( tID ) )
                {
                    mField->block( tID )->set_domain_type( DomainType::Ghost ) ;
                }
            }


            // loop over all sidesets
            for( id_t tID : mBlockIDs )
            {
                if( mField->sideset_exists( tID ) )
                {
                    mField->sideset( tID )->set_domain_type( mSideSetTypes( tID ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::compute_jacobian_and_rhs_cut(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {


            aJacobian( 0, 0 ) =  0.0 ;
            aJacobian( 1, 0 ) =  0.0 ;
            aJacobian( 2, 0 ) =  1.0 ;
            aJacobian( 0, 1 ) =  0.0 ;
            aJacobian( 1, 1 ) =  0.0 ;
            aJacobian( 2, 1 ) = -1.0 ;
            aJacobian( 0, 2 ) =  1.0 ;
            aJacobian( 1, 2 ) = -1.0 ;
            aJacobian( 2, 2 ) =  0.0 ;


            aRHS( 0 ) = 0.0 ;
            aRHS( 1 ) = 0.0 ;
            this->collect_lambda_data( aElement, "lambda_I", aRHS( 2 ) );

            // compute residual vector
            aRHS *= mDeltaTime ;

            //  finalize Jacobian ( M + delta_t * theta * K )
            aJacobian *= mDeltaTime ;

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::set_nan_values()
        {
            uint tDim = mesh::dimension( mElementType );

            if( mField->is_master() )
            {
                // flag current nodes
                mMesh->unflag_all_nodes();
                for ( id_t tID: mBlockIDs )
                {
                    if ( mMesh->block_exists( tID ) &&
                            (   mBlockTypes( tID ) == DomainType::SuperConductor
                             || mBlockTypes( tID ) == DomainType::Coil ) )
                    {
                        mMesh->block( tID )->flag_nodes();
                    }
                }

                // reset the current field
                if ( tDim == 3 )
                {
                    Vector< real > & tJx = mMesh->field_data( "jx" );
                    tJx.fill( BELFEM_QUIET_NAN );
                    Vector< real > & tJy = mMesh->field_data( "jy" );
                    tJy.fill( BELFEM_QUIET_NAN );
                    Vector< real > & tJz = mMesh->field_data( "jy" );
                    tJz.fill( BELFEM_QUIET_NAN );
                    for ( mesh::Node * tNode: mMesh->nodes())
                    {
                        if ( tNode->is_flagged())
                        {
                            tJx( tNode->index()) = 0.0;
                            tJy( tNode->index()) = 0.0;
                            tJz( tNode->index()) = 0.0;
                        }
                    }

                }
                else
                {
                    Vector< real > & tJz = mMesh->field_data( "jz" );
                    tJz.fill( BELFEM_QUIET_NAN );
                    uint tCount = 0 ;
                    for ( mesh::Node * tNode: mMesh->nodes() )
                    {
                        if ( tNode->is_flagged() )
                        {
                            tJz( tNode->index() ) = 0.0;
                            ++tCount ;
                        }
                    }
                }

                if ( is_maxwell_hphi( mType ))
                {
                    // flag air nodes
                    mMesh->unflag_all_nodes();
                    for ( id_t tID: mBlockIDs )
                    {
                        if ( mMesh->block_exists( tID ) && mBlockTypes( tID ) == DomainType::Air )
                        {
                            mMesh->block( tID )->flag_nodes();
                        }
                    }

                    // reset phi field
                    Vector< real > & tPhi = mMesh->field_data( "phi" );
                    tPhi.fill( BELFEM_QUIET_NAN );
                    for ( mesh::Node * tNode: mMesh->nodes())
                    {
                        if ( tNode->is_flagged())
                        {
                            tPhi( tNode->index()) = 0.0;
                        }
                    }

                    // flag ferro nodes
                    // flag air nodes
                    mMesh->unflag_all_nodes();
                    for ( id_t tID: mBlockIDs )
                    {
                        if ( mMesh->block_exists( tID ) && mBlockTypes( tID ) == DomainType::FerroMagnetic )
                        {
                            mMesh->block( tID )->flag_nodes();
                        }
                    }
                }
                else if ( is_maxwell_ha( mType ))
                {
                    // flag air nodes
                    mMesh->unflag_all_nodes();
                    for ( id_t tID: mBlockIDs )
                    {
                        if ( mMesh->block_exists( tID ) &&
                             ( mBlockTypes( tID ) == DomainType::Air
                               || mBlockTypes( tID ) == DomainType::FerroMagnetic
                               || mBlockTypes( tID ) == DomainType::Coil ))
                        {
                            mMesh->block( tID )->flag_nodes();
                        }
                    }
                }

                if ( tDim == 3 )
                {
                    Vector< real > & tAx = mMesh->field_data( "ax" );
                    tAx.fill( BELFEM_QUIET_NAN );
                    Vector< real > & tAy = mMesh->field_data( "ay" );
                    tAy.fill( BELFEM_QUIET_NAN );
                    Vector< real > & tAz = mMesh->field_data( "az" );
                    tAz.fill( BELFEM_QUIET_NAN );
                    for ( mesh::Node * tNode: mMesh->nodes())
                    {
                        if ( tNode->is_flagged())
                        {
                            tAx( tNode->index()) = 0.0;
                            tAy( tNode->index()) = 0.0;
                            tAz( tNode->index()) = 0.0;
                        }
                    }
                }
                else
                {
                    Vector< real > & tAz = mMesh->field_data( "az" );
                    tAz.fill( BELFEM_QUIET_NAN );
                    for ( mesh::Node * tNode: mMesh->nodes())
                    {
                        if ( tNode->is_flagged())
                        {
                            tAz( tNode->index()) = 0.0;
                        }
                    }
                }
            }

            comm_barrier() ;

        }

//------------------------------------------------------------------------------

        // specifically for thin shells
        const Vector< real > &
        IWG_Maxwell::collect_q0_hphi_thinshell( Element * aElement )
        {
            // grab the output vector
            Vector< real > & aQ0 = mGroup->work_nedelec() ;

            BELFEM_ASSERT( mGroup->domain_type() == DomainType::ThinShell,
                "function IWG_Maxwell::collect_thin_shell_hphi0 can only be applied to a thin shell" );

            // grab field data from mesh
            const Vector< real > & tPhi  = mMesh->field_data( "phi0" );
            const Vector< real > & tHe   = mMesh->field_data("edge_h0" );
            const Vector< real > & tHf   = mNumberOfDofsPerFace > 0 ? mMesh->field_data("face_h0" ) : tHe ;


            uint tLambdaCount = 0 ;

            // loop over all dofs
            for( uint k=0; k<mNumberOfDofsPerElement; ++k )
            {
                // grab dof
                Dof * tDof = aElement->dof( k );

                switch( tDof->mesh_basis()->entity_type() )
                {
                    case( EntityType::NODE ) :
                    {
                        aQ0( k ) = tPhi( tDof->mesh_basis()->index() );
                        break ;
                    }
                    case( EntityType::EDGE ) :
                    {
                        aQ0( k ) = tHe( tDof->mesh_basis()->index() );
                        break ;
                    }
                    case( EntityType::FACE ) :
                    {
                        aQ0( k ) = tHf( tDof->mesh_basis()->index() );
                        break ;
                    }
                    case( EntityType::FACET ) :
                    {
                        // get lambda field
                        const Vector< real > & tL = mMesh->field_data( mFields.ThinShellLast( tLambdaCount++ ) );
                        aQ0( k ) = tL( tDof->mesh_basis()->index() );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid dof type" );
                    }
                }
            }

            return aQ0 ;
        }

//------------------------------------------------------------------------------

        // called by factory to create the index map for the thin shells
        void
        IWG_Maxwell::create_ghost_map(  Mesh * aMesh,
                                        const Vector< id_t >   & aThinShellSideSets,
                                        const Vector< proc_t > & aCommTable )
        {
            // reset the map
            mGhostElementMap.clear() ;

            // initialize the counter
            index_t tCount = 0 ;

            proc_t tCommSize = comm_size() ;

            if( comm_rank() == 0 )
            {
                // collect the sideset ids
                for( id_t tS : aThinShellSideSets )
                {
                    tCount += aMesh->sideset( tS )->number_of_facets() ;
                }

                // wait for other procs
                comm_barrier() ;

                // send the number of elements
                send( aCommTable, tCount );

                // exit if there is nothing to do
                if( tCount == 0 )
                {
                    return ;
                }

                // collect the facets
                Cell< mesh::Facet * > tFacets( tCount, nullptr );
                Vector< id_t > tElementIDs( tCount );

                tCount = 0 ;

                // collect the sideset ids
                for( id_t tS : aThinShellSideSets )
                {
                    if ( aMesh->sideset_exists( tS ) )
                    {
                        Cell< mesh::Facet * > & tSFacets = aMesh->sideset( tS )->facets() ;
                        for( mesh::Facet * tFacet : tSFacets )
                        {
                            tFacets( tCount++ ) = tFacet ;
                        }
                    }
                }

                // counter for layers
                uint l = 0 ;
                uint tNumLayers = mGhostBlockIDs.length() ;

                // count the elements per proc
                Vector< index_t > tOwnerCount( tCommSize, 0 );


                // loop over all ghost blocks
                for( id_t b : mGhostBlockIDs )
                {
                    // grab the elements on the block
                    Cell< mesh::Element * > & tElements = aMesh->block( b )->elements() ;

                    // count owned elements

                    // sanity check
                    BELFEM_ASSERT( tElements.size() == tCount, "number of elements does not match ( is %lu expect %lu)",
                                   ( long unsigned int ) tElements.size(), tCount );

                    // reset the counter
                    tCount = 0 ;

                    // loop over all facets
                    for( mesh::Facet * tFacet : tFacets )
                    {
                        // create a key
                        luint tKey = tFacet->id() * tNumLayers + l ;

                        // get element
                        mesh::Element * tElement = tElements( tCount++ ) ;

                        // make sure that ownership is correct
                        BELFEM_ASSERT( tFacet->owner() == tElement->owner(),
                        "Owner of facet %lu is %u but the one of ghost element %lu is %u. They must be the same",
                                       ( long unsigned int ) tFacet->id(),
                                       ( unsigned int ) tFacet->owner(),
                                       ( long unsigned int ) tElement->id(),
                                       ( unsigned int ) tElement->owner() );

                        // increment the owner
                        ++tOwnerCount( tFacet->owner() );

                        // link key to element on mesh
                        mGhostElementMap[ tKey ] = tElement ;
                    }

                    // increment layer counter
                    ++l ;
                }

                // we are done if we run in serial
                if( tCommSize == 1 )
                {
                    return ;
                }

                // for all other procs we must now create the subtables
                Cell< Vector< long unsigned int > > tAllKeys( tCommSize, {} );
                Cell< Vector< id_t > > tAllIDs( tCommSize, {} );

                // allocate memory
                for( proc_t p=0; p<tCommSize; ++p )
                {
                    if( tOwnerCount( p ) > 0 )
                    {
                        tAllKeys( p ).set_size( tOwnerCount( p ));
                        tAllIDs( p ).set_size( tOwnerCount( p ));
                    }
                }

                // reset the counters
                tOwnerCount.fill( 0 ) ;
                l = 0 ;

                for( id_t b : mGhostBlockIDs )
                {
                    tCount = 0 ;

                    // grab the elements on the block
                    Cell< mesh::Element * > & tElements = aMesh->block( b )->elements() ;

                    // loop over all facets
                    for ( mesh::Facet * tFacet: tFacets )
                    {

                        // get the id container
                        Vector< luint > & tKeys = tAllKeys( tFacet->owner() );
                        Vector< id_t >  & tIDs  = tAllIDs( tFacet->owner() );

                        // store values
                        tKeys( tOwnerCount( tFacet->owner() ) )   = tFacet->id() * tNumLayers + l;
                        tIDs( tOwnerCount( tFacet->owner() )++ ) = tElements( tCount++ )->id() ;
                    }

                    // increment layer counter
                    ++l ;
                }

                // wait
                comm_barrier() ;

                // distribute values
                send( aCommTable, tAllKeys );
                send( aCommTable, tAllIDs );

            }
            else
            {
                comm_barrier() ;

                receive( 0, tCount );

                // exit if there is nothing to do
                if( tCount == 0 )
                {
                    return ;
                }

                // wait
                comm_barrier() ;

                // get the keys from the master
                Vector< luint > tMyKeys ;
                receive( 0, tMyKeys );

                // get the ids from the master
                Vector< id_t > tMyIDs ;
                receive( 0, tMyIDs );

                // create the map
                tCount = 0 ;
                for( luint tKey : tMyKeys )
                {
                    mGhostElementMap[ tKey ] = aMesh->element( tMyIDs( tCount++ ) );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell::set_magfield_bc_type( const id_t aSideSetID, const MagfieldBcType aType )
        {
            mMagfieldTypeMap[ aSideSetID ] = aType ;
        }
        
//------------------------------------------------------------------------------
    }
}