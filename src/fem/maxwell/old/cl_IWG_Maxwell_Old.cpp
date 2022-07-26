//
// Created by christian on 10/4/21.
//
#include "assert.hpp"
#include "cl_Timer.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_IWG_Maxwell_Old.hpp"
#include "cl_MaxwellMaterial.hpp"
#include "meshtools.hpp"
#include "geometrytools.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "cl_HDF5.hpp"
#include "cl_FEM_DofManager.hpp"
#include "fn_FEM_compute_element_current.hpp"
#include "en_FEM_MagfieldBcType.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"
#include "fn_inv.hpp"
#include "fn_det.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

    IWG_Maxwell_Old::IWG_Maxwell_Old(
            const uint aNumberOfDimensions,
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
            mNumberOfDimensions( aNumberOfDimensions ),
            mUseEdges( aUseEdges )
    {
        mNumberOfSpatialDimensions = aNumberOfDimensions ;
        mNumberOfDerivativeDimensions = aNumberOfDimensions ;
        mWorkN.set_size( aNumberOfDimensions, BELFEM_QUIET_NAN );
        mWorkT.set_size( aNumberOfDimensions, BELFEM_QUIET_NAN );

        mComputeJacobianOnBlock   = true ;
        mComputeJacobianOnSideset = true ;
    }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::set_blocks(
                const Vector< id_t >      & aBlockIDs,
                const Cell< DomainType >  & aBlockTypes )
        {
            mBlockIDs = aBlockIDs ;

            mBlockTypes.clear() ;

            index_t tCount = 0 ;
            for( id_t tID : aBlockIDs )
            {
                mBlockTypes[ tID ] = aBlockTypes( tCount++ );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::set_sidesets(
                const Vector< id_t > & aSideSetIDs,
                const Cell< DomainType > & aSideSetTypes )
        {
            mSideSetIDs = aSideSetIDs;
            index_t tCount = 0 ;
            for ( id_t tID: aSideSetIDs )
            {
                mSideSetTypes[ tID ] = aSideSetTypes( tCount++ );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::collect_dof_fields()
        {
            if( mDofFields.size() == 0 )
            {
                // append dofs
                mDofFields.clear() ;

                for ( string tDof : mScDofs )
                {
                    mDofFields.push( tDof );
                    mInterfaceDofsScFm.push( tDof );
                    mInterfaceDofsScAir.push( tDof );
                }

                for ( string tDof : mFerroDofs )
                {
                    mDofFields.push( tDof );
                    mInterfaceDofsScFm.push( tDof );
                    mInterfaceDofsFmAir.push( tDof );
                }

                for ( string tDof : mCoilDofs )
                {
                    mDofFields.push( tDof );
                }

                for ( string tDof : mCutDofs )
                {
                    mDofFields.push( tDof );
                    mInterfaceDofsCut.push( tDof );
                    mLambdaDofs.push( tDof );
                }

                for ( string tDof : mAirDofs )
                {
                    mDofFields.push( tDof );
                    mInterfaceDofsFmAir.push( tDof );
                    mInterfaceDofsScAir.push( tDof );
                    mInterfaceDofsCut.push( tDof );
                }

                for ( string tDof : mInterfaceDofs )
                {
                    mDofFields.push( tDof );
                    mInterfaceDofsScAir.push( tDof );
                    mLambdaDofs.push( tDof );
                }

                for ( string tDof : mFarfieldDofs )
                {
                    mDofFields.push( tDof );
                    mLambdaDofs.push( tDof );
                }
                for ( string tDof : mBoundaryDofs )
                {
                    mDofFields.push( tDof );
                    mLambdaDofs.push( tDof );
                }
                for( string tDof:  mSymmetryDofs )
                {
                    mDofFields.push( tDof );
                    mLambdaDofs.push( tDof );
                }

                // make lists unique
                unique( mDofFields );
                unique( mInterfaceDofsScAir );
                unique( mInterfaceDofsScFm );
                unique( mInterfaceDofsFmAir );
                unique( mInterfaceDofsCut );
                unique( mLambdaDofs );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::save( const string & aPath )
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

        /**
         * load data from HDF5 field
         */
        int
        IWG_Maxwell_Old::load( const string & aPath, Mesh * aMesh  )
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

                    for ( string tDofLabel: mDofFields )
                    {
                        tFile.load_data( tDofLabel, aMesh->field_data( tDofLabel ));
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
                    Cell< string > tDofFields = mDofFields ;
                    tDofFields.push( "T" );
                    mField->distribute_fields( tDofFields );
                }
                else
                {
                    mField->distribute_fields( mDofFields );
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

        void
        IWG_Maxwell_Old::collect_doftypes_per_block_and_sideset()
        {
            // create a temporary map for the dofs
            Map< string, index_t > tDofMap ;

            index_t tCount = 0 ;

            for ( string tDof : mDofFields )
            {
                tDofMap[ tDof ] = tCount++ ;
            }

            // add interface dofs
            for ( string tDof : mInterfaceDofsScAir )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }
            for ( string tDof : mInterfaceDofsScFm )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }
            for ( string tDof : mInterfaceDofsFmAir )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }

            // add cut dofs
            for ( string tDof : mInterfaceDofsCut )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }

            // add farfield dofs
            for ( string tDof : mFarfieldDofs )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }

            // add stmmetry dofs
            for ( string tDof : mSymmetryDofs )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }

            // add prescribed boundary dofs
            for ( string tDof : mBoundaryDofs )
            {
                // check if key already exists
                // should always be the case except for lambda dofs
                if( ! tDofMap.key_exists( tDof ) )
                {
                    tDofMap[ tDof ] = tCount++ ;
                }
            }

            // determine the number of blocks that are used
            uint tNumBlocks = mBlockIDs.length();
            mDofsPerBlock.set_size( tNumBlocks, Vector< id_t >() );

            // loop over all blocks
            for ( uint b=0; b<tNumBlocks; ++b )
            {

                // get dof vector
                Vector< id_t > & tDofsPerBlock = mDofsPerBlock( b );

                // check type of block
                switch ( mBlockTypes( mBlockIDs( b ) ) )
                {
                    case( DomainType::SuperConductor ) :
                    {
                        tDofsPerBlock.set_size( mScDofs.size() ) ;
                        for( index_t k=0; k<mScDofs.size(); ++k )
                        {
                            tDofsPerBlock( k ) = tDofMap( mScDofs( k ) );
                        }
                        break ;
                    }
                    case( DomainType::Coil ) :
                    {
                        tDofsPerBlock.set_size( mCoilDofs.size() ) ;
                        for( index_t k=0; k<mCoilDofs.size(); ++k )
                        {
                            tDofsPerBlock( k ) = tDofMap( mCoilDofs( k ) );
                        }
                        break ;
                    }
                    case( DomainType::FerroMagnetic ) :
                    {
                        tDofsPerBlock.set_size( mFerroDofs.size() ) ;
                        for( index_t k=0; k<mFerroDofs.size(); ++k )
                        {
                            tDofsPerBlock( k ) = tDofMap( mFerroDofs( k ) );
                        }
                        break ;
                    }
                    case( DomainType::Air ) :
                    {

                        tDofsPerBlock.set_size( mAirDofs.size() ) ;
                        for( index_t k=0; k<mAirDofs.size(); ++k )
                        {
                            tDofsPerBlock( k ) = tDofMap( mAirDofs( k ) );
                        }
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid block type" );
                    }
                }
            }

            this->select_blocks( mBlockIDs );

            // set dofs per sideset
            uint tNumSideSets = mSideSetIDs.length() ;
            mDofsPerSideSet.set_size( tNumSideSets, Vector< id_t >() );
            mSideSetIndices.clear() ;

            for ( uint s=0; s<tNumSideSets; ++s )
            {
                // write id into list
                mSideSetIndices[ mSideSetIDs( s ) ] = s ;

                // get dof vector
                Vector< index_t > & tDofsPerSideSet = mDofsPerSideSet( s );

                switch ( mSideSetTypes( mSideSetIDs( s ) ) )
                {
                    case( DomainType::InterfaceScAir ) :
                    {
                        tDofsPerSideSet.set_size( mInterfaceDofsScAir.size() ) ;
                        for( index_t k=0; k<mInterfaceDofsScAir.size(); ++k )
                        {
                            tDofsPerSideSet( k ) = tDofMap( mInterfaceDofsScAir( k ) );
                        }
                        break ;
                    }
                    case( DomainType::InterfaceScFm ) :
                    {
                        tDofsPerSideSet.set_size( mInterfaceDofsScFm.size() ) ;
                        for( index_t k=0; k<mInterfaceDofsScFm.size(); ++k )
                        {
                            tDofsPerSideSet( k ) = tDofMap( mInterfaceDofsScFm( k ) );
                        }
                        break ;
                    }
                    case( DomainType::InterfaceFmAir ) :
                    {
                        tDofsPerSideSet.set_size( mInterfaceDofsFmAir.size() ) ;
                        for( index_t k=0; k<mInterfaceDofsFmAir.size(); ++k )
                        {
                            tDofsPerSideSet( k ) = tDofMap( mInterfaceDofsFmAir( k ) );
                        }
                        break ;
                    }
                    case( DomainType::Cut ) :
                    {
                        tDofsPerSideSet.set_size( mInterfaceDofsCut.size() ) ;
                        for( index_t k=0; k<mInterfaceDofsCut.size(); ++k )
                        {
                            tDofsPerSideSet( k ) = tDofMap( mInterfaceDofsCut( k ) );
                        }

                        break ;
                    }
                    case( DomainType::Boundary ) :
                    {
                        switch( mMagfieldTypeMap( mSideSetIDs( s ) ) )
                        {
                            case( MagfieldBcType::Farfied ) :
                            {
                                tDofsPerSideSet.set_size( mAirDofs.size() + mFarfieldDofs.size() );
                                for( index_t k=0; k<mAirDofs.size(); ++k )
                                {
                                    tDofsPerSideSet( k ) = tDofMap( mAirDofs( k ) );
                                }
                                index_t c = mAirDofs.size() ;
                                for( index_t k=0; k<mFarfieldDofs.size(); ++k )
                                {
                                    tDofsPerSideSet( c++ ) = tDofMap( mFarfieldDofs( k ) );
                                }
                                break ;
                            }
                            case( MagfieldBcType::Symmetry ) :
                            {
                                tDofsPerSideSet.set_size( mAirDofs.size() + mSymmetryDofs.size() );
                                for( index_t k=0; k<mAirDofs.size(); ++k )
                                {
                                    tDofsPerSideSet( k ) = tDofMap( mAirDofs( k ) );
                                }
                                index_t c = mAirDofs.size() ;
                                for( index_t k=0; k<mSymmetryDofs.size(); ++k )
                                {
                                    tDofsPerSideSet( c++ ) = tDofMap( mSymmetryDofs( k ) );
                                }
                                break ;
                            }
                            case( MagfieldBcType::Wave ) :
                            {
                                tDofsPerSideSet.set_size( mAirDofs.size() + mBoundaryDofs.size() );
                                for( index_t k=0; k<mAirDofs.size(); ++k )
                                {
                                    tDofsPerSideSet( k ) = tDofMap( mAirDofs( k ) );
                                }
                                index_t c = mAirDofs.size() ;
                                for( index_t k=0; k<mBoundaryDofs.size(); ++k )
                                {
                                    tDofsPerSideSet( c++ ) = tDofMap( mBoundaryDofs( k ) );
                                }
                                break ;
                            }
                            default :
                            {
                                BELFEM_ERROR( false, "Invalid BC type for sideset %lu", (long unsigned int ) mSideSetIDs( s )) ;
                            }
                        }
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid sideset type" );
                    }
                }
            }

            this->select_sidesets( mSideSetIDs );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::initialize()
        {
            BELFEM_ASSERT( ! this->is_initialized(), "iwg has already been initialized" );

            this->collect_dof_fields();
            this->collect_doftypes_per_block_and_sideset() ;

            // call init function from parent
            IWG::initialize();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::set_nan_values()
        {
            // reset phi values
            if( mMesh->field_exists( "phi" ) )
            {
                Vector< real > & tPhi = mMesh->field_data( "phi" );
                tPhi.fill( BELFEM_QUIET_NAN );

                mMesh->unflag_all_nodes();
                for( id_t tID : mBlockIDs )
                {
                    if( mMesh->block_exists( tID ) && mBlockTypes( tID ) == DomainType::Air )
                    {
                        mMesh->block( tID )->flag_nodes() ;
                    }
                }

                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        tPhi( tNode->index() ) = 0.0 ;
                    }
                }

                // flag nodes for B
                for( id_t tID : mBlockIDs )
                {
                    if( mMesh->block_exists( tID ) &&
                        (   mBlockTypes( tID ) == DomainType::FerroMagnetic
                         || mBlockTypes( tID ) == DomainType::SuperConductor ) )
                    {
                        mMesh->block( tID )->flag_nodes() ;
                    }
                }

                Cell< string > tFields = { "bx", "by", "bz" };

                for( string tField: tFields )
                {
                    if( mMesh->field_exists( tField ) )
                    {
                        Vector< real > & tB = mMesh->field_data( tField );
                        tB.fill( BELFEM_QUIET_NAN );
                        for( mesh::Node * tNode : mMesh->nodes() )
                        {
                            if( tNode->is_flagged() )
                            {
                                tB( tNode->index() ) = 0.0 ;
                            }
                        }
                    }
                }

                // flag nodes for a
                mMesh->unflag_all_nodes() ;

                for( id_t tID : mBlockIDs )
                {
                    if( mMesh->block_exists( tID ) && mBlockTypes( tID ) == DomainType::FerroMagnetic )
                    {
                        mMesh->block( tID )->flag_nodes() ;
                    }
                }

            }
            else
            {
                for( id_t tID : mBlockIDs )
                {
                    if( mMesh->block_exists( tID ) &&
                        (      mBlockTypes( tID ) == DomainType::Air
                            || mBlockTypes( tID ) == DomainType::Coil
                            || mBlockTypes( tID ) == DomainType::FerroMagnetic) )
                    {
                        mMesh->block( tID )->flag_nodes() ;
                    }
                }

            }

            Cell< string > tFields = { "ax", "ay", "az" };

            for( string tField: tFields )
            {
                if( mMesh->field_exists( tField ) )
                {
                    Vector< real > & tA = mMesh->field_data( tField );
                    tA.fill( BELFEM_QUIET_NAN );
                    for( mesh::Node * tNode : mMesh->nodes() )
                    {
                        if( tNode->is_flagged() )
                        {
                            tA( tNode->index() ) = 0.0 ;
                        }
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::add_boundary_condition( BoundaryCondition * aBoundaryCondition )
        {
            IWG::add_boundary_condition( aBoundaryCondition );

            if( aBoundaryCondition->physics() == BoundaryConditionPhysics::Magfield )
            {
                // get subtype of bc
                MagfieldBcType tType
                    = reinterpret_cast< MaxwellBoundaryConditionMagfield * >( aBoundaryCondition )->subtype() ;

                for( id_t tID: aBoundaryCondition->sidesets() )
                {
                    mMagfieldTypeMap[ tID ] = tType ;
                }
            }
        }

//------------------------------------------------------------------------------

        int
        IWG_Maxwell_Old::check_mesh ( Mesh * aMesh, const proc_t aMasterRank )
        {
            // we flag all elements that sit on a superconducting block
            if( aMesh == nullptr )
            {
                return 1;
            }
            else
            {
                Timer tTimer ;

                proc_t tMyRank = comm_rank();

                if( tMyRank == aMasterRank )
                {
                    message( 4, "Checking mesh ...");
                }

                // first, we unflag all elements on the mesh
                aMesh->unflag_all_elements();

                // now we loop over all blocks ...
                for( id_t tBlockID : mBlockIDs )
                {
                    DomainType tType = mBlockTypes( tBlockID );
                    if( aMesh->block_exists( tBlockID ) )
                    {
                        // ... and check if this block is a solid
                        if(       tType == DomainType::SuperConductor
                               || tType == DomainType::Coil
                               || tType == DomainType::FerroMagnetic )
                        {
                            // if so, we flag all elements on the block
                            aMesh->block( tBlockID )->flag_elements() ;
                        }
                    }
                }

                // now, we loop over all relevant sidesets ...
                for( id_t tSideSetID : mSideSetIDs )
                {
                    DomainType tType = mSideSetTypes( tSideSetID );

                    if( aMesh->sideset_exists( tSideSetID ) )
                    {
                        // ... and check if this sideset is an interface
                        if( tType == DomainType::InterfaceScAir
                         || tType == DomainType::InterfaceFmAir
                         || tType == DomainType::InterfaceScFm )
                        {
                            // if so, we grab the facets of this sideset ...
                            Cell< mesh::Facet * > & tFacets = aMesh->sideset( tSideSetID )->facets();


                            // now we make sure that the orientation of the face is correct
                            Cell< mesh::Node * > tNodes ;
                            uint tCount ;
                            for( mesh::Facet * tFacet : tFacets )
                            {

                                tFacet->master()->get_nodes_of_facet( tFacet->master_index(), tNodes );

                                mesh::Element * tElement = tFacet->element();
                                tCount = 0 ;
                                for( mesh::Node * tNode : tNodes )
                                {
                                    tElement->insert_node( tNode, tCount++ );
                                }
                            }


                            // now we do the same thing for the edges
                            if( aMesh->edges_exist() )
                            {
                                Cell< mesh::Edge * > tEdges ;
                                for( mesh::Facet * tFacet : tFacets )
                                {
                                    tFacet->master()->get_edges_of_facet( tFacet->master_index(), tEdges );

                                    mesh::Element * tElement = tFacet->element();

                                    // make sure that edge container has been allocated
                                    if( ! tElement->has_edges() )
                                    {
                                        tElement->allocate_edge_container();
                                    }

                                    tCount = 0 ;
                                    for( mesh::Edge * tEdge : tEdges )
                                    {
                                        tElement->insert_edge( tEdge, tCount++ );
                                    }
                                }
                            }
                        }
                    }
                }


                if ( tMyRank == aMasterRank )
                {
                    message( 4, "    ... time for checking mesh                  : %u ms\n",
                             ( unsigned int ) tTimer.stop() );
                }

                return 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::allocate_work_matrices( Group    * aGroup )
        {
            // Jacobian
            aGroup->work_invJ().set_size( mNumberOfDimensions, mNumberOfDimensions );
            aGroup->work_J().set_size( mNumberOfDimensions, mNumberOfDimensions );

            if ( aGroup->type() == GroupType::SIDESET )
            {
                SideSet *  tSideSet = reinterpret_cast< SideSet * >( aGroup );

                switch( mSideSetTypes( tSideSet->id() ) )
                {
                    case( DomainType::InterfaceScAir ) :
                    case( DomainType::InterfaceFmAir ) :
                    case( DomainType::InterfaceScFm ) :
                    {
                        // get number of nodes per element
                        mNumberOfNodesPerElement = aGroup->parent()->enforce_linear_interpolation() ?
                                mesh::number_of_corner_nodes( tSideSet->slave_type() ) :
                                mesh::number_of_nodes( tSideSet->slave_type() );

                        // get number of edges per element
                        mNumberOfEdgesPerElement = mesh::number_of_edges( tSideSet->master_type() );

                        // integration points on master
                        aGroup->work_Tau().set_size( mNumberOfDimensions, aGroup->integration_weights().length() );

                        // integration points on slave
                        aGroup->work_Sigma().set_size( mNumberOfDimensions, aGroup->integration_weights().length() );

                        // all dofs from last timestep
                        aGroup->work_psi().set_size( mNumberOfDofsPerElement );

                        // all dofs from current timestep
                        aGroup->work_phi().set_size( mNumberOfDofsPerElement );

                        // help value for edge data
                        aGroup->work_chi().set_size( mNumberOfEdgesPerElement );

                        // help matrix for node data
                        aGroup->work_sigma().set_size( mNumberOfNodesPerElement * mAirDofs.size() );

                        aGroup->work_N().set_size( 1, mNumberOfDofsPerElement, 0.0 );

                        // help matrix for shape functions
                        aGroup->work_L().set_size( mNumberOfDimensions == 2 ? 1 : 3, mNumberOfDofsPerElement, 0.0 );

                        // help matrix for derivative
                        aGroup->work_H().set_size( mNumberOfDimensions, mNumberOfDofsPerElement );

                        break ;
                    }
                    case( DomainType::Cut ) :
                    {
                        mNumberOfNodesPerElement = 2 ;
                        mNumberOfEdgesPerElement = 0 ;

                        break ;
                    }
                    case( DomainType::Boundary ) :
                    {
                        // get number of nodes per element
                        mNumberOfNodesPerElement = aGroup->parent()->enforce_linear_interpolation() ?
                                                   mesh::number_of_corner_nodes( tSideSet->master_type() ) :
                                                   mesh::number_of_nodes( tSideSet->master_type() );

                        // get number of edges per element
                        mNumberOfEdgesPerElement = mesh::number_of_edges( tSideSet->master_type() );

                        // integration points on master
                        aGroup->work_Tau().set_size( mNumberOfDimensions, aGroup->integration_weights().length() );

                        // integration points on slave ( acutually not needed )
                        aGroup->work_Sigma().set_size( mNumberOfDimensions, aGroup->integration_weights().length() );


                        // help vector for node data
                        aGroup->work_chi().set_size( mNumberOfNodesPerElement );

                        // help matrix for node data and B-Field
                        aGroup->work_Chi().set_size(  mNumberOfNodesPerElement, mNumberOfDimensions);

                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "unsupported sideset type");
                    }
                }


            }
            else if ( aGroup->type() == GroupType::BLOCK )
            {
                // get number of nodes per element
                mNumberOfNodesPerElement = aGroup->parent()->enforce_linear_interpolation() ?
                                           mesh::number_of_corner_nodes( aGroup->element_type() ) :
                                           mesh::number_of_nodes( aGroup->element_type() );

                // get number of edges per element
                mNumberOfEdgesPerElement = mesh::number_of_edges( aGroup->element_type() );

                switch( mBlockTypes( aGroup->id() ) )
                {
                    case( DomainType::SuperConductor ) :
                    {
                        // edge values for last timestep
                        aGroup->work_psi().set_size( mNumberOfDofsPerElement );

                        // node temperature values for last timestep
                        aGroup->work_sigma().set_size( mNumberOfNodesPerElement );

                        // node temperature for current timestep
                        aGroup->work_tau().set_size( mNumberOfNodesPerElement );

                        break ;
                    }
                    case( DomainType::Air ) :
                    case( DomainType::Coil ) :
                    case( DomainType::FerroMagnetic ) :
                    {
                        // node values for last timestep
                        aGroup->work_psi().set_size( mNumberOfDofsPerElement );

                        // node values for current timestep
                        aGroup->work_phi().set_size( mNumberOfDofsPerElement );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "unsupported block type");
                    }
                }
            }

            // for grad phi
            aGroup->work_B().set_size( mNumberOfDimensions,
                                       mNumberOfNodesPerElement,
                                       0.0 );

            // for curl a
            aGroup->work_D().set_size( mNumberOfDimensions,
                                       mNumberOfDimensions == 2
                                        ? mNumberOfNodesPerElement :
                                          mNumberOfDimensions * mNumberOfNodesPerElement,
                                       0.0 );

            // gradient matrix for edges
            aGroup->work_dEdX().set_size(
                    mNumberOfDimensions,
                    mNumberOfEdgesPerElement );

            // gradient matrix for edges
            aGroup->work_dEdXi().set_size(
                    mNumberOfDimensions,
                    mNumberOfEdgesPerElement );

            mNablaXi.set_size( mNumberOfDimensions );
            mNablaEta.set_size( mNumberOfDimensions );
            mNablaZeta.set_size( mNumberOfDimensions );
            mNablaTau.set_size( mNumberOfDimensions );


        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::link_to_group( Group * aGroup )
        {
            IWG::link_to_group( aGroup );

            // skip empty blocks and sidesets
            if( aGroup->number_of_elements() == 0 ) return;


              ElementType tType = aGroup->type() == GroupType::CUT ? ElementType::EMPTY : aGroup->type() == GroupType::SIDESET ?
                    reinterpret_cast< SideSet * >( aGroup )->master_type() : aGroup->element_type() ;

            switch ( tType )
            {
                case( ElementType::EMPTY ) :
                {
                    break ;
                }
                case( ElementType::TRI3 ) :
                case( ElementType::TRI6 ) :
                {
                    mFunNabla   = & IWG_Maxwell_Old::nabla_tri3 ;
                    mFunE       = & IWG_Maxwell_Old::E_tri3 ;
                    mFunE2D     = & IWG_Maxwell_Old::E_2d_tri3 ;
                    mFunCurlH   = & IWG_Maxwell_Old::curl_h_tri3 ;
                    mFunCurlA   = & IWG_Maxwell_Old::curl_a_tri3 ;
                    mFunGradPhi = & IWG_Maxwell_Old::grad_phi_tri3 ;
                    mFunNormal  = & IWG_Maxwell_Old::normal_linear ;
                    mFunProjectIntpoints = & IWG_Maxwell_Old::project_intpoints_tri ;
                    break ;
                }
                case( ElementType::TET4 ) :
                case( ElementType::TET10 ) :
                {
                    mFunNabla   = & IWG_Maxwell_Old::nabla_tet4 ;
                    mFunE       = & IWG_Maxwell_Old::E_tet4 ;
                    mFunE3D     = & IWG_Maxwell_Old::E_3d_tet4 ;
                    mFunCurlH   = & IWG_Maxwell_Old::curl_h_tet4 ;
                    mFunCurlA   = & IWG_Maxwell_Old::curl_a_tet4 ;
                    mFunGradPhi = & IWG_Maxwell_Old::grad_phi_tet4 ;
                    mFunNormal  = & IWG_Maxwell_Old::normal_linear ;
                    mFunProjectIntpoints = & IWG_Maxwell_Old::project_intpoints_tet ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "unsupported element type: %s");
                }
            }

            // check matrix function for superconductors
            if( aGroup->type() == GroupType::BLOCK )
            {
                if( mBlockTypes( aGroup->id() ) == DomainType::SuperConductor )
                {
                    // get material of block
                    const Material * tMat = aGroup->material();

                    // make sure that permeability law is constant
                    BELFEM_ERROR(
                        tMat->permeability_law() == PermeabilityLaw::Constant,
                        "superconductor material %s requires constant permeability",
                        tMat->label().c_str() );

                    switch( tMat->resistivity_law() )
                    {
                        case( ResistivityLaw::Constant ) :
                        {
                            mFunScMatrices = & IWG_Maxwell_Old::sc_matrices_const ;
                            break ;
                        }
                        case( ResistivityLaw::DependJ ) :
                        {
                            mFunScMatrices = & IWG_Maxwell_Old::sc_matrices_powerlaw_ej ;
                            break ;
                        }
                        case( ResistivityLaw::DependJT ) :
                        {
                            mFunScMatrices = & IWG_Maxwell_Old::sc_matrices_powerlaw_ejt ;
                            break ;
                        }
                        case( ResistivityLaw::DependJB ) :
                        {
                            mFunScMatrices = & IWG_Maxwell_Old::sc_matrices_powerlaw_ejb ;
                            break ;
                        }
                        case( ResistivityLaw::DependJBT ) :
                        {
                            mFunScMatrices = & IWG_Maxwell_Old::sc_matrices_powerlaw_ejbt ;
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "unsupported material law for material %s for block %lu",
                                         tMat->label().c_str(),
                                         ( long unsigned int ) aGroup->id() );
                        }
                    }
                }
            }

            // now link the matrix functions that are IWG specific
            this->link_jacobian_function( aGroup );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::nabla_tri3( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {

                Matrix< real > & tJ = mGroup->work_J();


                tJ( 0, 0 ) = aElement->element()->node( 0 )->x() - aElement->element()->node( 2 )->x();
                tJ( 1, 0 ) = aElement->element()->node( 1 )->x() - aElement->element()->node( 2 )->x();
                tJ( 0, 0 ) = aElement->element()->node( 0 )->y() - aElement->element()->node( 2 )->y();
                tJ( 1, 0 ) = aElement->element()->node( 1 )->y() - aElement->element()->node( 2 )->y();

                mGroup->work_det_J() =
                        tJ( 0, 0 ) * tJ( 1, 1 )
                        - tJ( 0, 1 ) * tJ( 1, 0 );

                mDomainIncrement = std::abs( mGroup->work_det_J() );

                mGroup->work_det_invJ() = 1.0 / mGroup->work_det_J();

                Matrix< real > & tNabla = mGroup->work_invJ();

                tNabla( 0, 0 ) =  tJ( 1, 1 );
                tNabla( 1, 0 ) = -tJ( 1, 0 );
                tNabla( 0, 1 ) = -tJ( 0, 1 );
                tNabla( 1, 1 ) =  tJ( 0, 0 );

                tNabla *= mGroup->work_det_invJ();

                mNablaXi = tNabla.col( 0 );
                mNablaEta = tNabla.col( 1 );
                mNablaZeta( 0 ) = -mNablaXi( 0 ) - mNablaEta( 0 );
                mNablaZeta( 1 ) = -mNablaXi( 1 ) - mNablaEta( 1 );

            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::nabla_tet4( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & tX = mGroup->node_coords() ;

                for( uint k=0; k<4; ++k )
                {
                    tX( 0, k ) = aElement->element()->node( k )->x() ;
                    tX( 1, k ) = aElement->element()->node( k )->y() ;
                    tX( 2, k ) = aElement->element()->node( k )->z() ;
                }

                Matrix< real > & tJ = mGroup->work_J();

                tJ( 0, 0 ) = tX( 0, 0 ) - tX( 3, 0 );
                tJ( 1, 0 ) = tX( 1, 0 ) - tX( 3, 0 );
                tJ( 2, 0 ) = tX( 2, 0 ) - tX( 3, 0 );

                tJ( 0, 1 ) = tX( 0, 1 ) - tX( 3, 1 );
                tJ( 1, 1 ) = tX( 1, 1 ) - tX( 3, 1 );
                tJ( 2, 1 ) = tX( 2, 1 ) - tX( 3, 1 );

                tJ( 0, 2 ) = tX( 0, 2 ) - tX( 3, 2 );
                tJ( 1, 2 ) = tX( 1, 2 ) - tX( 3, 2 );
                tJ( 2, 2 ) = tX( 2, 2 ) - tX( 3, 2 );

                mGroup->work_det_J() =
                        tJ( 0, 0 ) *
                        ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ))
                        + tJ( 0, 1 ) *
                          ( tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 ))
                        + tJ( 0, 2 ) *
                          ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ));

                mDomainIncrement = std::abs( mGroup->work_det_J());

                mGroup->work_det_invJ() = 1.0 / mGroup->work_det_J();

                Matrix< real > & tNabla = mGroup->work_invJ();

                tNabla( 0, 0 ) = tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 );
                tNabla( 1, 0 ) = tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 );
                tNabla( 2, 0 ) = tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 );

                tNabla( 0, 1 ) = tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 );
                tNabla( 1, 1 ) = tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 );
                tNabla( 2, 1 ) = tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 );

                tNabla( 0, 2 ) = tJ( 0, 1 ) * tJ( 1, 2 ) - tJ( 0, 2 ) * tJ( 1, 1 );
                tNabla( 1, 2 ) = tJ( 0, 2 ) * tJ( 1, 0 ) - tJ( 0, 0 ) * tJ( 1, 2 );
                tNabla( 2, 2 ) = tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 );

                tNabla *= mGroup->work_det_invJ();

                mNablaXi = tNabla.col( 0 );
                mNablaEta = tNabla.col( 1 );
                mNablaZeta = tNabla.col( 2 );
                mNablaTau( 0 ) = -mNablaXi( 0 ) - mNablaEta( 0 ) - mNablaZeta( 0 );
                mNablaTau( 1 ) = -mNablaXi( 1 ) - mNablaEta( 1 ) - mNablaZeta( 1 );
                mNablaTau( 2 ) = -mNablaXi( 2 ) - mNablaEta( 2 ) - mNablaZeta( 2 );
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_Old::curl_h_tri3( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & aC = mGroup->work_C();

                aC( 0, 0 ) = aElement->edge_direction( 0 ) ? 1.0 : -1.0 ;
                aC( 0, 1 ) = aElement->edge_direction( 1 ) ? 1.0 : -1.0 ;
                aC( 0, 2 ) = aElement->edge_direction( 2 ) ? 1.0 : -1.0 ;

                aC *= mGroup->work_det_invJ() + mGroup->work_det_invJ() ;

                return aC ;
            }
            else
            {
                return mGroup->work_C() ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_Old::curl_h_tet4( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & aC = mGroup->work_C();

                Matrix< real > & tX = mGroup->node_coords() ;

                for( uint k=0; k<4; ++k )
                {
                    tX( 0, k ) = aElement->element()->node( k )->x() ;
                    tX( 1, k ) = aElement->element()->node( k )->y() ;
                    tX( 2, k ) = aElement->element()->node( k )->z() ;
                }

                if( aElement->edge_direction( 0 ) )
                {
                    aC( 0, 0 ) = tX( 3, 0 ) - tX( 2, 0 ) ;
                    aC( 1, 0 ) = tX( 3, 1 ) - tX( 2, 1 ) ;
                    aC( 2, 0 ) = tX( 3, 2 ) - tX( 2, 2 ) ;
                }
                else
                {
                    aC( 0, 0 ) = tX( 2, 0 ) - tX( 3, 0 ) ;
                    aC( 1, 0 ) = tX( 2, 1 ) - tX( 3, 1 ) ;
                    aC( 2, 0 ) = tX( 2, 2 ) - tX( 3, 2 ) ;
                }
                if( aElement->edge_direction( 1 ) )
                {
                    aC( 0, 1 ) = tX( 3, 0 ) - tX( 0, 0 ) ;
                    aC( 1, 1 ) = tX( 3, 1 ) - tX( 0, 1 ) ;
                    aC( 2, 1 ) = tX( 3, 2 ) - tX( 0, 2 ) ;
                }
                else
                {
                    aC( 0, 1 ) = tX( 0, 0 ) - tX( 3, 0 ) ;
                    aC( 1, 1 ) = tX( 0, 1 ) - tX( 3, 1 ) ;
                    aC( 2, 1 ) = tX( 0, 2 ) - tX( 3, 2 ) ;
                }
                if( aElement->edge_direction( 2 ) )
                {
                    aC( 0, 2 ) = tX( 3, 0 ) - tX( 1, 0 ) ;
                    aC( 1, 2 ) = tX( 3, 1 ) - tX( 1, 1 ) ;
                    aC( 2, 2 ) = tX( 3, 2 ) - tX( 1, 2 ) ;
                }
                else
                {
                    aC( 0, 2 ) = tX( 1, 0 ) - tX( 3, 0 ) ;
                    aC( 1, 2 ) = tX( 1, 1 ) - tX( 3, 1 ) ;
                    aC( 2, 2 ) = tX( 1, 2 ) - tX( 3, 2 ) ;
                }
                if( aElement->edge_direction( 3 ) )
                {
                    aC( 0, 3 ) = tX( 2, 0 ) - tX( 1, 0 ) ;
                    aC( 1, 3 ) = tX( 2, 1 ) - tX( 1, 1 ) ;
                    aC( 2, 3 ) = tX( 2, 2 ) - tX( 1, 2 ) ;
                }
                else
                {
                    aC( 0, 3 ) = tX( 1, 0 ) - tX( 2, 0 ) ;
                    aC( 1, 3 ) = tX( 1, 1 ) - tX( 2, 1 ) ;
                    aC( 2, 3 ) = tX( 1, 2 ) - tX( 2, 2 ) ;
                }
                if( aElement->edge_direction( 4 ) )
                {
                    aC( 0, 4 ) = tX( 0, 0 ) - tX( 2, 0 ) ;
                    aC( 1, 4 ) = tX( 0, 1 ) - tX( 2, 1 ) ;
                    aC( 2, 4 ) = tX( 0, 2 ) - tX( 2, 2 ) ;
                }
                else
                {
                    aC( 0, 4 ) = tX( 2, 0 ) - tX( 0, 0 ) ;
                    aC( 1, 4 ) = tX( 2, 1 ) - tX( 0, 1 ) ;
                    aC( 2, 4 ) = tX( 2, 2 ) - tX( 0, 2 ) ;
                }

                if( aElement->edge_direction( 5 ) )
                {
                    aC( 0, 5 ) = tX( 1, 0 ) - tX( 0, 0 ) ;
                    aC( 1, 5 ) = tX( 1, 1 ) - tX( 0, 1 ) ;
                    aC( 2, 5 ) = tX( 1, 2 ) - tX( 0, 2 ) ;
                }
                else
                {
                    aC( 0, 5 ) = tX( 0, 0 ) - tX( 1, 0 ) ;
                    aC( 1, 5 ) = tX( 0, 1 ) - tX( 1, 1 ) ;
                    aC( 2, 5 ) = tX( 0, 2 ) - tX( 1, 2 ) ;
                }

                aC *= mGroup->work_det_invJ() + mGroup->work_det_invJ() ;

                return aC ;
            }
            else
            {
                return mGroup->work_C() ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_Old::curl_a_tri3( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & aC = mGroup->work_D();
                const Matrix< real > & tInvJ = mGroup->work_invJ() ;

                aC( 0, 0 ) =  tInvJ( 1, 0 ) ;
                aC( 1, 0 ) = -tInvJ( 0, 0 );
                aC( 0, 1 ) =  tInvJ( 1, 1 );
                aC( 1, 1 ) = -tInvJ( 0, 1 ) ;
                aC( 0, 2 ) = -tInvJ( 1, 0 ) - tInvJ( 1, 1 ) ;
                aC( 1, 2 ) =  tInvJ( 0, 0 )  +  tInvJ( 0, 1 ) ;

                return aC ;
            }
            else
            {
                return mGroup->work_D() ;
            }
        }
//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_Old::curl_a_tet4( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & aC = mGroup->work_D();
                const Matrix< real > & tInvJ = mGroup->work_invJ() ;

                aC( 1, 0 ) =  tInvJ( 2, 0 );
                aC( 2, 0 ) = -tInvJ( 1, 0 );
                aC( 0, 1 ) = -tInvJ( 2, 0 );
                aC( 2, 1 ) =  tInvJ( 0, 0 );
                aC( 0, 2 ) =  tInvJ( 1, 0 );
                aC( 1, 2 ) = -tInvJ( 0, 0 );
                aC( 1, 3 ) =  tInvJ( 2, 1 );
                aC( 2, 3 ) = -tInvJ( 1, 1 );
                aC( 0, 4 ) = -tInvJ( 2, 1 );
                aC( 2, 4 ) =  tInvJ( 0, 1 );
                aC( 0, 5 ) =  tInvJ( 1, 1 );
                aC( 1, 5 ) = -tInvJ( 0, 1 );
                aC( 1, 6 ) =  tInvJ( 2, 2 );
                aC( 2, 6 ) = -tInvJ( 1, 2 );
                aC( 0, 7 ) = -tInvJ( 2, 2 );
                aC( 2, 7 ) =  tInvJ( 0, 2 );
                aC( 0, 8 ) =  tInvJ( 1, 2 );
                aC( 1, 8 ) = -tInvJ( 0, 2 );
                aC( 1, 9 ) = -tInvJ( 2, 0 ) - tInvJ( 2, 1 ) - tInvJ( 2, 2 );
                aC( 2, 9 ) =  tInvJ( 1, 0 ) + tInvJ( 1, 1 ) + tInvJ( 1, 2 );
                aC( 0, 10 ) =  tInvJ( 2, 0 ) + tInvJ( 2, 1 ) + tInvJ( 2, 2 );
                aC( 2, 10 ) = -tInvJ( 0, 0 ) - tInvJ( 0, 1 ) - tInvJ( 0, 2 );
                aC( 0, 11 ) = -tInvJ( 1, 0 ) - tInvJ( 1, 1 ) - tInvJ( 1, 2 );
                aC( 1, 11 ) =  tInvJ( 0, 0 ) + tInvJ( 0, 1 ) + tInvJ( 0, 2 );

                return aC ;
            }
            else
            {
                return mGroup->work_D() ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_Old::grad_phi_tri3( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & aB = mGroup->work_B() ;

                aB( 0, 0 ) = -mNablaXi( 0 );
                aB( 1, 0 ) = -mNablaXi( 1 );

                aB( 0, 1 ) = -mNablaEta( 0 );
                aB( 1, 1 ) = -mNablaEta( 1 );

                aB( 0, 1 ) = -mNablaZeta( 0 );
                aB( 1, 1 ) = -mNablaZeta( 1 );

                return aB ;
            }
            else
            {
                return mGroup->work_B() ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_Old::grad_phi_tet4( Element * aElement, const uint aIndex )
        {
            if( aIndex == 0 )
            {
                Matrix< real > & aB = mGroup->work_B() ;

                aB( 0, 0 ) = -mNablaXi( 0 );
                aB( 1, 0 ) = -mNablaXi( 1 );
                aB( 2, 0 ) = -mNablaXi( 2 );

                aB( 0, 1 ) = -mNablaEta( 0 );
                aB( 1, 1 ) = -mNablaEta( 1 );
                aB( 2, 1 ) = -mNablaEta( 2 );

                aB( 0, 2 ) = -mNablaZeta( 0 );
                aB( 1, 2 ) = -mNablaZeta( 1 );
                aB( 2, 2 ) = -mNablaZeta( 2 );

                aB( 0, 3 ) = -mNablaTau( 0 );
                aB( 1, 3 ) = -mNablaTau( 1 );
                aB( 2, 3 ) = -mNablaTau( 2 );

                return aB ;
            }
            else
            {
                return mGroup->work_B() ;
            }
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell_Old::normal_linear( Element * aElement, const uint aIndex )
        {
            // assume that the normal has already been computed by
            // project_intpoints
            return mWorkN ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::print_dofs( Element * aElement )
        {
            std::cout << "Element : " << aElement->id() << std::endl ;

            unsigned int tCount = 0 ;
            uint tNumDofs = aElement->number_of_dofs() ;

            for( uint d=0; d<tNumDofs; ++d )
            {
                Dof * tDof = aElement->dof( d );

                // get label
                const string & tLabel = mDofFields( tDof->type_id() ) ;

                int tN = 6 - tLabel.length() ;
                tN = tN < 0 ? 0 : tN ;
                string tSpace = "";
                for ( int i=0; i<tN; ++i )
                {
                    tSpace += " ";
                }
                string tType ;
                if( tDof->is_node() )
                {
                    tType = "node" ;
                }
                else if( tDof->is_edge() )
                {
                    tType = "edge" ;
                }
                else if( tDof->is_face() )
                {
                    tType = "face" ;
                }
                else if( tDof->is_lambda() )
                {
                    tType = "lambda" ;
                }
                else
                {
                    tType = "????" ;
                }

                fprintf( stdout, "    %2u  %s%s id %8u     index %8u     %s %8u\n",
                         tCount++,
                         tLabel.c_str(),
                         tSpace.c_str(),
                         ( unsigned int ) tDof->id(),
                         ( unsigned int ) tDof->index(),
                         tType.c_str(),
                         ( unsigned int ) tDof->mesh_basis()->id()
                );
            }
        }


//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::link_jacobian_function( Group * aGroup )
        {
            BELFEM_ERROR( false,
                         "compute_jacobian_and_rhs() not implemented for abstract base class");
        }


//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::sc_matrices_const(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;


            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // update nabla functions
                this->update_nabla( aElement, k );

                // edge function
                const Matrix< real > & tE = this->E( aElement, k );

                // curl operator
                const Matrix< real > & tC = this->curl_h( aElement, k );

                // compute mass matrix
                aM += tW( k ) * trans( tE ) * tE * mDomainIncrement ;

                // compute stiffness matrix
                aK += tW( k ) * trans( tC ) * tMat->rho_el() * tC * mDomainIncrement ;
            }

            aM *= constant::mu0 * tMat->mu_r();
        }
//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::sc_matrices_powerlaw_ej(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // vector with field at current time
            Vector< real > & tH = mGroup->work_phi() ;

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // update nabla functions
                this->update_nabla( aElement, k );

                // edge function
                const Matrix< real > & tE = this->E( aElement, k );

                // curl operator
                const Matrix< real > & tC = this->curl_h( aElement, k );

                // compute mass matrix
                aM += tW( k ) * trans( tE ) * tE * mDomainIncrement ;

                // compute stiffness matrix
                aK += tW( k ) * trans( tC ) * tMat->rho_el( norm( tC * tH ) )
                        * tC * mDomainIncrement ;
            }

            aM *= constant::mu0 * tMat->mu_r();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::sc_matrices_powerlaw_ejt(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // vector with field at current time
            Vector< real > & tH = mGroup->work_phi() ;

            // temperature from current timestep
            Vector< real > & tT = mGroup->work_tau() ;
            this->collect_node_data( aElement, "T", tT );

            // compute new balance if this is not a backwards Euler
            if( this->theta() != 1.0 )
            {
                tT *= this->theta() ;

                // collect data from last timestep
                this->collect_node_data( aElement, "T0", mGroup->work_sigma() );
                tT += ( 1.0 - this->theta() ) * mGroup->work_sigma() ;
            }

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // update nabla functions
                this->update_nabla( aElement, k );

                // edge function
                const Matrix< real > & tE = this->E( aElement, k );

                // curl operator
                const Matrix< real > & tC = this->curl_h( aElement, k );

                // compute mass matrix
                aM += tW( k ) * trans( tE ) * tE * mDomainIncrement ;

                // compute stiffness matrix
                aK += tW( k ) * trans( tC ) * tMat->rho_el(
                        norm( tC * tH ),
                        dot( mGroup->n( k ), tT ) ) * tC * mDomainIncrement ;
            }

            aM *= constant::mu0 * tMat->mu_r();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::sc_matrices_powerlaw_ejb(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // vector with field at current time
            Vector< real > & tH = mGroup->work_phi() ;

            // permeability constant
            const real tMu = constant::mu0 * tMat->mu_r();

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // update nabla functions
                this->update_nabla( aElement, k );

                // edge function
                const Matrix< real > & tE = this->E( aElement, k );

                // curl operator
                const Matrix< real > & tC = this->curl_h( aElement, k );

                // compute mass matrix
                aM += tW( k ) * trans( tE ) * tE * mDomainIncrement ;

                // compute stiffness matrix
                aK += tW( k ) * trans( tC )
                        * tMat->rho_el(
                            norm( tC * tH ),
                            0.0,
                            tMu * norm( tE * tH  ) ) * tC * mDomainIncrement ;
            }

            aM *= tMu ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::sc_matrices_powerlaw_ejbt(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK )
        {
            aM.fill( 0.0 );
            aK.fill( 0.0 );

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights();

            // grab material
            const Material * tMat = mGroup->material() ;

            // vector with field at current time
            Vector< real > & tH = mGroup->work_phi() ;

            // temperature from current timestep
            Vector< real > & tT = mGroup->work_tau() ;
            this->collect_node_data( aElement, "T", tT );

            // permeability constant
            const real tMu = constant::mu0 * tMat->mu_r();

            // compute new balance if this is not a backwards Euler
            if( this->theta() != 1.0 )
            {
                tT *= this->theta() ;

                // collect data from last timestep
                this->collect_node_data( aElement, "T0", mGroup->work_sigma() );
                tT += ( 1.0 - this->theta() ) * mGroup->work_sigma() ;
            }

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // update nabla functions
                this->update_nabla( aElement, k );

                // edge function
                const Matrix< real > & tE = this->E( aElement, k );

                // curl operator
                const Matrix< real > & tC = this->curl_h( aElement, k );

                // compute mass matrix
                aM += tW( k ) * trans( tE ) * tE * mDomainIncrement ;

                // compute stiffness matrix
                aK += tW( k ) * trans( tC ) * tMat->rho_el(
                        norm( tC * tH ),
                        dot( mGroup->n( k ), tT ),
                        tMu * norm( tE * tH ) ) * tC * mDomainIncrement ;
            }

            aM *= tMu ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::project_intpoints_tri(
                Element * aElement,
                      Matrix< real > & aIntpointsOnMaster,
                      Matrix< real > & aIntPointsOnSlave )
        {
            const Matrix< real > & tIntpoints = mGroup->integration_points() ;

            // check if a slave side exists
            uint n = aElement->slave() == nullptr ? 1 : 2 ;

            // loop over both sides
            for( uint s=0; s<n; ++s )
            {
                // projected coords
                Matrix< real > & tXi = s==0 ? aIntpointsOnMaster : aIntPointsOnSlave ;

                // get indes
                uint tSideIndex = s==0 ? aElement->facet()->master_index() : aElement->facet()->slave_index();

                // get the number of points
                uint tNumPoints = tIntpoints.n_cols() ;

                switch( tSideIndex )
                {
                    case( 0 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) = 0.5 - 0.5 * tIntpoints( 0, k ) ;

                            // eta
                            tXi( 1, k ) = 0.5 * tIntpoints( 0, k ) + 0.5 ;
                        }

                        break;
                    }
                    case( 1 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) = 0.0 ;

                            // eta
                            tXi( 1, k ) =  0.5 - 0.5 * tIntpoints( 0, k ) ;
                        }

                        break;
                    }
                    case( 2 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) =  0.5 * tIntpoints( 0, k ) + 0.5 ;

                            // eta
                            tXi( 1, k ) = 0.0 ;
                        }

                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid side of triangle");
                    }
                }

                if( s == 1 )
                {
                    tXi *= -1 ;
                    tXi += 1.0 ;
                }
            }

            // compute tangent, normal and domain increment
            const mesh::Node * tA = nullptr ;
            const mesh::Node * tB = nullptr ;

            // pick node indices for edge
            switch( aElement->facet()->master_index() )
            {
                case( 0 ) :
                {
                    tA = aElement->master()->element()->node( 0 );
                    tB = aElement->master()->element()->node( 1 );
                    break ;
                }
                case( 1 ) :
                {
                    tA = aElement->master()->element()->node( 1 );
                    tB = aElement->master()->element()->node( 2 );
                    break ;
                }
                case( 2 ) :
                {
                    tA = aElement->master()->element()->node( 2 );
                    tB = aElement->master()->element()->node( 0 );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid side of triangle");
                }
            }

            // compute tangent
            mWorkT( 0 ) = tB->x() - tA->x() ;
            mWorkT( 1 ) = tB->y() - tA->y() ;

            // compute length of element, tangent and normal
            mDomainIncrement = norm( mWorkT );

            // scale tangent
            mWorkT /= mDomainIncrement ;

            // since the sum of the weights is 2, we must
            // scale the increment
            mDomainIncrement *= 0.5 ;

            // finally, we can compute the normal
            mWorkN( 0 ) =  mWorkT( 1 );
            mWorkN( 1 ) = -mWorkT( 0 );

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::project_intpoints_tet(
                Element * aElement,
                Matrix< real > & aIntpointsOnMaster,
                Matrix< real > & aIntPointsOnSlave )
        {
            const Matrix< real > & tIntpoints = mGroup->integration_points() ;
            // check if a slave side exists
            uint n = aElement->slave() == nullptr ? 1 : 2 ;

            // loop over both sides
            for( uint s=0; s<n; ++s )
            {
                // projected coords
                Matrix< real > & tXi = s==0 ? aIntpointsOnMaster : aIntPointsOnSlave ;

                // get indes
                uint tSideIndex = s==0 ? aElement->facet()->master_index() : aElement->facet()->slave_index() ;

                // get the number of points
                uint tNumPoints = tIntpoints.n_cols() ;

                switch( tSideIndex )
                {
                    case( 0 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) = tIntpoints( 0, k ) ;

                            // eta
                            tXi( 1, k ) = tIntpoints( 1, k ) ;

                            // zeta
                            tXi( 2, k ) = 0.0 ;
                        }

                        break;
                    }
                    case( 1 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) = 0.0 ;

                            // eta
                            tXi( 1, k ) = tIntpoints( 0, k ) ;

                            // zeta
                            tXi( 2, k ) = tIntpoints( 1, k ) ;
                        }

                        break;
                    }
                    case( 2 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) = tIntpoints( 1, k ) ;

                            // eta
                            tXi( 1, k ) = 0.0 ;

                            // zeta
                            tXi( 2, k ) = tIntpoints( 0, k ) ;
                        }

                        break;
                    }
                    case( 3 ) :
                    {
                        for( uint k=0; k<tNumPoints; ++k )
                        {
                            // xi
                            tXi( 0, k ) = tIntpoints( 0, k ) ;

                            // eta
                            tXi( 1, k ) = tIntpoints( 1, k ) ;

                            // zeta
                            tXi( 2, k ) = 1.0 - tIntpoints( 0, k ) - tIntpoints( 1, k ) ;
                        }

                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid side of tetrahedron");
                    }
                }

                if( s == 1 )
                {
                    tXi *= -1 ;
                    tXi += 1.0 ;
                }
            }

            // compute tangent, normal and domain increment
            const mesh::Node * tA = nullptr ;
            const mesh::Node * tB = nullptr ;
            const mesh::Node * tC = nullptr ;

            // pick node indices for edge
            switch( aElement->facet()->master_index() )
            {
                case( 0 ) :
                {
                    tA = aElement->master()->element()->node( 0 );
                    tB = aElement->master()->element()->node( 1 );
                    tC = aElement->master()->element()->node( 3 );
                    break ;
                }
                case( 1 ) :
                {
                    tA = aElement->master()->element()->node( 1 );
                    tB = aElement->master()->element()->node( 2 );
                    tC = aElement->master()->element()->node( 3 );
                    break ;
                }
                case( 2 ) :
                {
                    tA = aElement->master()->element()->node( 2 );
                    tB = aElement->master()->element()->node( 0 );
                    tC = aElement->master()->element()->node( 3 );
                    break ;
                }
                case( 3 ) :
                {
                    tA = aElement->master()->element()->node( 0 );
                    tB = aElement->master()->element()->node( 1 );
                    tC = aElement->master()->element()->node( 2 );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid side of triangle");
                }
            }

            // direction of first edge
            mWorkS( 0 ) = tB->x() - tA->x() ;
            mWorkS( 1 ) = tB->y() - tA->y() ;
            mWorkS( 2 ) = tB->z() - tA->z() ;

            // direction of second edge
            mWorkT( 0 ) = tC->x() - tB->x() ;
            mWorkT( 1 ) = tC->y() - tB->y() ;
            mWorkT( 2 ) = tC->z() - tB->z() ;

            // normal vector
            mWorkN = cross( mWorkS, mWorkT );

            // = 2 * triangle surface, but the sum of all integration points
            //   is 0.5, so this is correct
            mDomainIncrement = norm( mWorkN );

            // scale normal vector
            mWorkN /= mDomainIncrement ;


        }
//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::compute_element_current_2d( Element * aElement, real & aJz )
        {
            BELFEM_ASSERT(  mField->mesh()->number_of_dimensions() == 2,
                "this function must be called for 2D-meshes only" );

            Vector< real > & tPsi = mGroup->work_psi() ;

            // get edge dofs from master element
            this->collect_edge_data( aElement ,
                                     "edge_h", tPsi );


            this->update_nabla( aElement, 0 );

            aJz = dot( this->curl_h( aElement, 0 ).row( 0 ), tPsi.vector_data() );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::compute_element_current_3d( Element * aElement, real & aJx, real & aJy, real & aJz )
        {
            BELFEM_ASSERT(  mField->mesh()->number_of_dimensions() == 3,
                             "this function must be called for 3D-meshes only" );

            Vector< real > & tPsi = mGroup->work_psi() ;

            // get edge dofs from master element
            this->collect_edge_data( aElement ,"edge_h", tPsi );

            this->update_nabla( aElement, 0 );

            // interpolation function
            const Matrix< real > & tC = this->curl_h( aElement, 0 );

            // any 3d vector will do
            mWorkS = tC * tPsi ;

            aJx = mWorkS( 0 );
            aJy = mWorkS( 1 );
            aJz = mWorkS( 2 );
        }

//----------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::ferro_stiffness_2d(  Element * aElement,
                                          Matrix< real > & aK )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // current potential
            Vector< real > & tAz       = mGroup->work_phi() ;
            Vector< real > & tAz0      = mGroup->work_psi() ;

            // reset result vectors
            aK.fill( 0.0 );

            // collect field
            this->collect_node_data( aElement, "a0z", tAz0 );
            this->collect_node_data( aElement, "az", tAz );
            if( this->theta() != 1.0 )
            {
                tAz *= this->theta() ;
                tAz += ( 1.0 - this->theta() ) * tAz0 ;
            }

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute curl operator
                const Matrix< real > & tC = this->curl_a( aElement, k ) ;


                aK += tW( k ) * trans( tC ) * mMaterial->nu_s( norm( tC * tAz ) ) * tC * this->domain_increment() ;
            }
        }

//----------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::ferro_stiffness_3d(  Element * aElement,
                                          Matrix< real > & aK )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // current potential
            Vector< real > & tQ       = mGroup->work_phi() ;
            Vector< real > & tQ0      = mGroup->work_psi() ;

            // reset result vectors
            aK.fill( 0.0 );

            Matrix< real > & tA       = mGroup->work_Phi();
            Matrix< real > & tA0      = mGroup->work_Psi();

            // collect field
            this->collect_node_data( aElement, mAfields0, tA0 );
            this->collect_node_data( aElement, mAfields, tA );

            uint tCount = 0 ;

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tQ( tCount++ ) = tA( i, 0 );
                tQ( tCount++ ) = tA( i, 1 );
                tQ( tCount++ ) = tA( i, 2 );
            }
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tQ0( tCount++ ) = tA0( i, 0 );
                tQ0( tCount++ ) = tA0( i, 1 );
                tQ0( tCount++ ) = tA0( i, 2 );
            }

            if( this->theta() != 1.0 )
            {
                tQ *= this->theta() ;
                tQ += ( 1.0 - this->theta() ) * tQ0 ;
            }

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute curl operator
                const Matrix< real > & tC = this->curl_a( aElement, k ) ;


                aK += tW( k ) * trans( tC ) * mMaterial->nu_s( norm( tC * tQ ) ) * tC * this->domain_increment() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::compute_Nxn_2d(
                const Matrix< real > & aXi,
                const uint             aIntpoint,
                const uint             aOffset,
                const Vector< real > & aNormal,
                Matrix< real > & aN,
                Matrix< real > & aNxn )
        {
            // get node shape function
            Matrix< real > & tN = mGroup->work_N() ;

            // populate node shapes
            // may need to change this when introducing #facedof
            uint tCount = aOffset ;
            tN( 0, tCount++ ) = aXi( 0, aIntpoint );
            tN( 0, tCount++ ) = aXi( 1, aIntpoint );
            tN( 0, tCount++ ) = 1. - aXi( 0, aIntpoint ) - aXi( 1, aIntpoint );

            // expression cross(n, N )
            tCount = aOffset ;
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                aNxn( 0, tCount ) = -aNormal( 1 ) * tN( 0, tCount );
                aNxn( 1, tCount ) =  aNormal( 0 ) * tN( 0, tCount );
                ++tCount;
            }
        }

//--------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::compute_B_2d(
                const uint       aOffset,
                Matrix< real > & aB )
        {
            uint tCount = aOffset ;
            const Matrix< real > & tiJ = mGroup->work_invJ() ;

            aB( 0, tCount   ) =  tiJ( 0, 0 );
            aB( 1, tCount++ ) =  tiJ( 1, 0 );
            aB( 0, tCount   ) =  tiJ( 0, 1 );
            aB( 1, tCount++ ) =  tiJ( 1, 1 );
            aB( 0, tCount   ) = -tiJ( 0, 0 )-tiJ( 0, 1 );
            aB( 1, tCount++ ) = -tiJ( 1, 0 )-tiJ( 1, 1 );
        }

//----------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::interface_ha_2d( Element * aElement,
                      Matrix< real > & aM, Matrix< real > & aK )

        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 1: determine integration points for both elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // integration points on master element ( superconductor )
            Matrix< real > & tXiH = mGroup->work_Tau() ;

            // integration points on slave element ( air )
            Matrix< real > & tXiA = mGroup->work_Sigma() ;

            // project the integration points on the edge
            this->project_intpoints( aElement, tXiH, tXiA );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 2: get containers
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // reset mass matrix
            aM.fill( 0.0 );

            Matrix< real > & tN = mGroup->work_N() ;

            // get edge shape function
            Matrix< real > & tE = mGroup->work_E() ;

            // help matrix for Nxn
            Matrix< real > & tNxn = mGroup->work_Psi() ;


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 3: compute matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // populate node shape function
                this->compute_Nxn_2d( tXiA, k, mNumberOfEdgesPerElement, this->normal( aElement, k ), tN, tNxn );

                // populate edge shape function
                tE = this->E( aElement, tXiH( 0, k ), tXiH( 1, k ) );

                // stiffness
                aM += tW( k ) * trans( tE ) * tNxn * this->domain_increment() ;
            }

            aK = -trans( aM );
        }

//----------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::interface_phia_2d( Element * aElement,
                        Matrix< real > & aM, Matrix< real > & aK )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 1: determine integration points for both elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // integration points on master element ( superconductor )
            Matrix< real > & tXiH = mGroup->work_Tau() ;

            // integration points on slave element ( air )
            Matrix< real > & tXiA = mGroup->work_Sigma() ;

            // project the integration points on the edge
            this->project_intpoints( aElement, tXiH, tXiA );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 2: grab matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // reset mass matrix
            aM.fill( 0.0 );

            // get node shape function
            Matrix< real > & tN = mGroup->work_N() ;

            // get function
            Matrix< real > & tB = mGroup->work_B() ;

            // help matrix for Nxn
            Matrix< real > & tNxn = mGroup->work_Psi() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 3: compute matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // change this if introducing #facedof
            this->compute_B_2d( mNumberOfNodesPerElement, tB );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // populate node shape function
                this->compute_Nxn_2d( tXiA, k, 0, this->normal( aElement, k ), tN, tNxn );

                // stiffness
                aM += tW( k ) * trans( tB ) * tNxn * this->domain_increment() ;


            }
            aM *= constant::mu0 ;

            aK = -trans( aM );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::interface_hphi_2d( Element * aElement, Matrix< real > & aK, Vector< real > & aRHS )
        {

            BELFEM_ASSERT(
                    aElement->master()->element()->type() == ElementType::TRI3,
                    "Element must be a linear triangle" );

            // dot product
            aRHS.fill( 0.0 );
            aK.fill( 0.0 );

            // integration points on master element ( air )
            Matrix< real > & tXiH = mGroup->work_Tau() ;
            Matrix< real > & tXiPhi = mGroup->work_Sigma() ;

            this->update_nabla( aElement->slave(), 0.0 );

            Matrix< real > & tB = mGroup->work_B() ;
            this->compute_B_2d( mNumberOfEdgesPerElement, tB );

            this->update_nabla( aElement->master(), 0.0 );
            this->project_intpoints( aElement, tXiH, tXiPhi );

            const Vector< real > & tn = this->normal( aElement, 0 );
            const real nx = tn( 0 );
            const real ny = tn( 1 );

            Vector< real > & tN = mGroup->work_phi() ;

            const Vector< real > & tW = mGroup->integration_weights() ;

            uint n = mNumberOfNodesPerElement + mNumberOfEdgesPerElement ;

            /*Vector< real > H( 3 );
            this->collect_edge_data( aElement->master(), "edge_h", H );

            Vector< real > Phi( 3 );
            this->collect_node_data( aElement->slave(), "phi", Phi ); */


            for( uint k=0; k<tW.length(); ++k )
            {
                tN( 3 ) = tXiPhi( 0, k );
                tN( 4 ) = tXiPhi( 1, k );
                tN( 5 ) = 1-tXiPhi( 0, k )-tXiPhi( 1, k );

                const Matrix< real > & tE = this->E( aElement->master(), tXiH( 0, k ), tXiH( 1, k ) );

                //real Bx = tE( 0, 0 ) * H( 0 ) + tE( 0, 1 ) * H( 1 ) + tE( 0, 2 ) * H( 2 );
                //real By = tE( 1, 0 ) * H( 0 ) + tE( 1, 1 ) * H( 1 ) + tE( 1, 2 ) * H( 2 );
                //Bx *= constant::mu0;
                //By *= constant::mu0;

                /*
                real bx = tB( 0, 3 ) * Phi( 0 ) + tB( 0, 4 ) * Phi( 1 ) + tB( 0, 5 ) * Phi( 2 );
                real by = tB( 1, 3 ) * Phi( 0 ) + tB( 1, 4 ) * Phi( 1 ) + tB( 1, 5 ) * Phi( 2 );
                bx *= -constant::mu0;
                by *= -constant::mu0;

                if( aElement->id() == 25 )
                {


                    H.print("H");
                    Phi.print("Phi");

                    std::cout << "b " << bx << " " << by << std::endl;
                    std::cout << "B " << Bx << " " << By << std::endl << std::endl;
                }*/

                for( uint i=0; i<mNumberOfEdgesPerElement; ++i )
                {
                    aK( i, n ) += tW( k ) * ( nx*tE( 1, i ) -ny*tE( 0, i ) ) * this->domain_increment();
                }
                for( uint i=mNumberOfEdgesPerElement; i<n; ++i )
                {
                    aK( i, n ) += tW( k ) * ( nx*tB( 1, i ) -ny*tB( 0, i ) ) * this->domain_increment() ;
                    //aRHS( i ) -= tW( k ) * tN( i ) *  ( nx*Bx +ny*By ) * this->domain_increment() ;
                    for( uint j=0; j<mNumberOfEdgesPerElement; ++j )
                    {
                      aK( i, j ) -= tW( k ) * tN( i ) * ( nx * tE( 0, j ) + ny * tE( 1, j ) ) * this->domain_increment();
                    }
                }
            }
            aK.set_row( n, aK.col( n )) ;
            aK *= constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Old::projection_jacobian_b2d( Element * aElement, Matrix< real> & aJacobian )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // reset result vectors
            aJacobian.fill( 0.0 );

            uint tCount;

            // loop over all integration points
            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute shape
                const Vector< real > & tn = mGroup->n( k );
                Matrix< real > & tN = mGroup->work_N();

                tCount = 0;
                for ( uint i = 0; i < mNumberOfNodesPerElement; ++i )
                {
                    tN( 0, tCount++ ) = tn( i );
                    tN( 1, tCount++ ) = tn( i );
                }

                aJacobian += tW( k ) * trans( tN ) * tN * this->domain_increment();
            }
        }

//------------------------------------------------------------------------------

    } /* end namespace fem */
}  /* end namespace belfem */