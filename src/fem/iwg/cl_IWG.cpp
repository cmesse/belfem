//
// Created by Christian Messe on 08.11.19.
//

#include "cl_IWG.hpp"
#include "en_FEM_DomainType.hpp"
#include "assert.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_BoundaryCondition.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Kernel.hpp"
#include "meshtools.hpp"
#include "commtools.hpp"
#include "cl_Material.hpp"
#include "fn_entity_type.hpp"


namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG::IWG(  const IwgType aType,
                   const IwgMode aMode,
                   const SymmetryMode aSymmetryMode,
                   const DofMode aDofMode,
                   const SideSetDofLinkMode aSideSetDofLinkMode ) :
            mRank( comm_rank() ),
            mType( aType ),
            mMode( aMode ),
            mSymmetryMode( aSymmetryMode ),
            mDofMode( aDofMode ),
            mSideSetDofLinkMode( aSideSetDofLinkMode )
        {
            // set default value for solver algorithm
            switch( aMode )
            {
                case( IwgMode::Direct ) :
                {
                    mSolverAlgorithm = SolverAlgorithm::Direct ;
                    break ;
                }
                case( IwgMode::Iterative ) :
                {
                    mSolverAlgorithm = SolverAlgorithm::NewtonRaphson ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unsupported IWG Mode");
                    mSolverAlgorithm = SolverAlgorithm::UNDEFINED ;
                    break ;
                }
            }
        }

//------------------------------------------------------------------------------

        IWG::~IWG()
        {
            // delete block dof tables if they have been created
            this->delete_block_dof_tables();

            // delete sideset dof tables if they have been created
            this->delete_sideset_dof_tables();

            // delete boundary conditions if they have been set
            this->delete_boundary_conditions();
        }

//------------------------------------------------------------------------------

        void
        IWG::set_algorithm( const SolverAlgorithm aAlgorithm  )
        {
            BELFEM_ERROR( mMode == IwgMode::Iterative,
                         "algorithm can only be set for iterative solvers" );



            BELFEM_ERROR( aAlgorithm !=  SolverAlgorithm::Direct,
                         "forbidden algorithm for iterative solvers: SolverAlgorithm::Direct" );

            mSolverAlgorithm = aAlgorithm ;
        }

//------------------------------------------------------------------------------

        SymmetryMode
        IWG::symmetry_mode() const
        {
            return mSymmetryMode ;
        }

//------------------------------------------------------------------------------

        void
        IWG::set_num_rhs_cols( const uint & aNumRhsCols )
        {
            mNumberOfRhsCols = aNumRhsCols ;
        }

//------------------------------------------------------------------------------

        void
        IWG::set_wetted_sidesets( const Vector< id_t > & aSideSets )
        {
            // count sidesets that exist on thos proc
            uint tCount = 0;

            for ( id_t tID : aSideSets )
            {

                // check of sideset exists on this proc
                if ( mField->sideset_exists( tID ) && mField->mesh()->sideset_exists( tID ) )
                {
                    // increment counter
                    ++tCount;
                }
            }

            // allocate container
            mWettedSidesets.set_size( tCount );

            // reset counter
            tCount = 0;

            for ( id_t tID : aSideSets )
            {
                // check of sideset exists on this proc
                if ( mField->sideset_exists( tID ) && mField->mesh()->sideset_exists( tID ) )
                {
                    // increment counter
                    mWettedSidesets( tCount++ ) = tID;
                }
            }

            this->collect_nodes_on_wetted_sitdesets( mField->mesh(), aSideSets );
        }

//------------------------------------------------------------------------------

        void
        IWG::collect_nodes_on_wetted_sitdesets( Mesh * aMesh, const Vector< id_t > & aSideSets )
        {
            // unflag all nodes on the mesh
            aMesh->unflag_all_nodes() ;

            // loop over all wedded sidesets
            for( id_t tID : aSideSets )
            {
                if( aMesh->sideset_exists( tID ) )
                {
                    aMesh->sideset( tID )->flag_all_nodes() ;
                }
            }

            // get node contailer
            Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

            // initialize counter
            index_t tCount = 0 ;

            // count flagged nodes
            for( mesh::Node * tNode : tNodes )
            {
                // check if node is flagged
                if( tNode->is_flagged() )
                {
                    // increment counter
                    ++tCount ;
                }
            }

            // allocate memory
            mNodesOnWettedSidesets.set_size( tCount );

            // reset counter
            tCount = 0 ;

            // now loop again and add the indices to the container
            for( mesh::Node * tNode : tNodes )
            {
                // check if node is flagged
                if( tNode->is_flagged() )
                {
                   mNodesOnWettedSidesets( tCount++ ) = tNode->index() ;
                }
            }

        }

//------------------------------------------------------------------------------

        const Vector< id_t >  &
        IWG::wetted_sidesets() const
        {
            return mWettedSidesets;
        }

//------------------------------------------------------------------------------

        IwgType
        IWG::type() const
        {
            return mType;
        }

//------------------------------------------------------------------------------

        IwgMode
        IWG::mode() const
        {
            return mMode;
        }

//------------------------------------------------------------------------------

        void
        IWG::select_block( const id_t aBlockID )
        {
            // remember block IDs
            mBlockIDs.set_size( 1, aBlockID );

            // populate map
            mBlockIndices.clear() ;

            index_t tCount = 0 ;
            for( id_t tID : mBlockIDs )
            {
                mBlockIndices[ tID ] = tCount++;
            }

            if( mDofMode == DofMode::AllBlocksEqual )
            {
                this->assign_dofs_per_block( mBlockIDs );
            }

            // allocate memory
            this->create_block_dof_tables( tCount );
        }

//------------------------------------------------------------------------------

        void
        IWG::select_blocks( const Vector< id_t > & aBlockIDs )
        {
            // remember block IDs
            mBlockIDs = aBlockIDs ;

            // populate map
            mBlockIndices.clear() ;

            index_t tCount = 0 ;
            for( id_t tID : aBlockIDs )
            {
                mBlockIndices[ tID ] = tCount++;
            }

            if( mDofMode == DofMode::AllBlocksEqual )
            {
                this->assign_dofs_per_block( aBlockIDs );
            }

            // allocate memory
            this->create_block_dof_tables( tCount );
        }

//------------------------------------------------------------------------------

        void
        IWG::select_sidesets( const Vector< id_t > & aSideSetIDs )
        {
            // remember sideset IDs
            mSideSetIDs = aSideSetIDs ;

            // populate map
            mSideSetIndices.clear() ;

            index_t tCount = 0 ;
            for( id_t tID : aSideSetIDs )
            {
                mSideSetIndices[ tID ] = tCount++;
            }

            // also set ghost sidesets if selected
            for( id_t tID : mGhostSideSetIDs )
            {
                std::cout << "ghost set " << tID << " " << tCount << std::endl ;

                mSideSetIndices[ tID ] = tCount++ ;
            }

            if( mDofMode == DofMode::AllBlocksEqual )
            {
                this->assign_dofs_per_sideset( aSideSetIDs );
            }

            // allocate memory
            this->create_sideset_dof_tables( aSideSetIDs.length(),
                                             mGhostSideSetIDs.length() );
        }

//------------------------------------------------------------------------------

        void
        IWG::set_thin_shell_link_mode(
                const SideSetDofLinkMode aLinkMode )
        {
            mThinShellDofLinkMode = aLinkMode ;
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            BELFEM_ERROR( false, "compute_jacobian() not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false, "compute_rhs() [ Vector ] not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_rhs(
                Element        * aElement,
                Matrix< real > & aRHS )
        {
            BELFEM_ERROR( false, "compute_rhs() [ Matrix ] not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false, "compute_jacobian_and_rhs() not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_convection(
                Element        * aElement,
                Vector< real > & aConvection )
        {
            BELFEM_ERROR( false, "compute_convection() not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        void
        IWG::link_to_group( Group * aGroup )
        {
            mGroup    = aGroup;
            mField    = aGroup->parent();
            mMaterial = aGroup->material();

            if ( mGroup->element_type() == ElementType::EMPTY )
            {
                mNumberOfDofsPerElement = 0 ;
                mNumberOfIntegrationPoints = 0;
                mNumberOfEdgesPerElement = 0 ;
                mNumberOfNodesPerElement = 0 ;
            }
            else
            {
                mNumberOfNodesPerElement = aGroup->number_of_nodes_per_element() ;
                mNumberOfEdgesPerElement = aGroup->number_of_edges_per_element() ;
                mNumberOfFacesPerElement = aGroup->number_of_faces_per_element() ;

                switch ( mDofMode )
                {
                    case ( DofMode::AllBlocksEqual ) :
                    {


                        mNumberOfNodeDofsPerElement =  mNumberOfNodesPerElement *
                                                       mNumberOfDofsPerNode;

                        mNumberOfEdgeDofsPerElement =
                                  mNumberOfEdgesPerElement * mNumberOfDofsPerEdge
                                + mNumberOfFacesPerElement * mNumberOfDofsPerFace ;

                        mNumberOfDofsPerElement = mNumberOfNodeDofsPerElement
                                + mNumberOfEdgeDofsPerElement ;

                        break;
                    }
                    case ( DofMode::BlockSpecific ) :
                    {
                        switch ( aGroup->type() )
                        {
                            case ( GroupType::BLOCK ) :
                            {
                                mNumberOfDofsPerNode = this->number_of_dofs_per_node( aGroup->id() );

                                mNumberOfDofsPerEdge = this->number_of_dofs_per_edge( aGroup->id() );
                                mNumberOfDofsPerFace = this->number_of_dofs_per_face( aGroup->id() );

                                mNumberOfNodeDofsPerElement =  mNumberOfNodesPerElement *
                                                               mNumberOfDofsPerNode;

                                mNumberOfEdgeDofsPerElement =
                                        mNumberOfEdgesPerElement * mNumberOfDofsPerEdge
                                        + mNumberOfFacesPerElement * mNumberOfDofsPerFace ;

                                // the element may contain lambda-dofs,
                                // so we look it up in our table
                                mNumberOfDofsPerElement = this->number_of_dofs_per_element(
                                        reinterpret_cast< Block * > ( aGroup ) );
                                break ;
                            }
                            case( GroupType::CUT ) :
                            {
                                mNumberOfDofsPerNode = 1 ;
                                mNumberOfDofsPerEdge = 0 ;
                                mNumberOfDofsPerFace = 0 ;
                                mNumberOfDofsPerElement = 3 ;

                                // the following lines are not needed,
                                // but just for the sake of my sanity
                                mNumberOfNodesPerElement = 2 ;
                                mNumberOfNodeDofsPerElement = 2 ;
                                mNumberOfEdgeDofsPerElement = 0 ;

                                break ;
                            }
                            case ( GroupType::SIDESET ) :
                            case( GroupType::SHELL ) :
                            {
                                index_t tIndex = mSideSetIndices( aGroup->id() );

                                mNumberOfDofsPerNode = mSideSetDofs( tIndex )->Node.length() ;

                                mNumberOfDofsPerEdge = mSideSetDofs( tIndex )->Edge.length() * mEdgeDofMultiplicity ;
                                mNumberOfDofsPerFace = mSideSetDofs( tIndex )->Face.length() * mFaceDofMultiplicity ;

                                mNumberOfNodeDofsPerElement = mNumberOfNodesPerElement *
                                                              mNumberOfDofsPerNode;

                                if( aGroup->domain_type() == DomainType::ThinShell )
                                {
                                    mNumberOfEdgeDofsPerElement =
                                            ( mNumberOfEdgesPerElement * mNumberOfDofsPerEdge
                                            + mNumberOfFacesPerElement * mNumberOfDofsPerFace ) * aGroup->number_of_thin_shell_layers() ;
                                }
                                else
                                {
                                    mNumberOfEdgeDofsPerElement =
                                            mNumberOfEdgesPerElement * mNumberOfDofsPerEdge
                                            + mNumberOfFacesPerElement * mNumberOfDofsPerFace;
                                }


                                // the element may contain lambda-dofs,
                                // so we look it up in our table
                                mNumberOfDofsPerElement =
                                        this->number_of_dofs_per_element(
                                        reinterpret_cast< SideSet * > ( aGroup ) );
                                break ;
                            }
                            default:
                            {
                                BELFEM_ERROR( false, "Invalid group type for Group %u : %u",
                                             ( unsigned int ) aGroup->id(), ( unsigned int ) aGroup->type() );
                            }
                        }
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid Dof Mode" );
                    }
                }

                if( aGroup->type() == GroupType::CUT )
                {
                    mNumberOfIntegrationPoints = 0 ;
                }
                else
                {
                    mNumberOfIntegrationPoints = aGroup->integration_weights().length();
                }


                // RHS dofs must be set if we have an L2 projection
                if( mNumberOfRhsEdgeDofsPerElement == 0 )
                {
                    mNumberOfRhsEdgeDofsPerElement = mNumberOfEdgeDofsPerElement ;
                }

                this->allocate_work_matrices( aGroup );
            }
        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_dofs_per_element( Block * aBlock ) const
        {
            // get block index from map
            index_t tBlockIndex = mBlockIndices( aBlock->id() );

            return aBlock->number_of_nodes_per_element()
                   * mBlockDofs( tBlockIndex )->Node.length()
                   + aBlock->number_of_edges_per_element()
                   * this->number_of_dofs_per_edge( aBlock->id() )
                   + aBlock->number_of_faces_per_element()
                   * this->number_of_dofs_per_face( aBlock->id() ) ;

        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_nodes_per_element( SideSet * aSideSet ) const
        {
            // check if this is an empty sideset
            if ( aSideSet->number_of_elements() == 0 )
            {
                return 0 ;
            }
            else
            {
                switch ( this->sideset_dof_link_mode() )
                {
                    case( SideSetDofLinkMode::FacetOnly ) :
                    {
                        return aSideSet->number_of_nodes_per_element() ;
                    }
                    case( SideSetDofLinkMode::FacetAndMaster ) :
                    {
                        return aSideSet->parent()->enforce_linear_interpolation() ?
                            mesh::number_of_corner_nodes( aSideSet->master_type() ) :
                            mesh::number_of_nodes( aSideSet->master_type() );
                    }
                    case( SideSetDofLinkMode::MasterAndSlave ) :
                    {
                        if( this->has_edge_dofs() )
                        {
                            return aSideSet->parent()->enforce_linear_interpolation() ?
                                   mesh::number_of_corner_nodes(
                                           aSideSet->slave_type() == ElementType::EMPTY ? aSideSet->master_type() : aSideSet->slave_type() ) :
                                   mesh::number_of_nodes( aSideSet->slave_type() == ElementType::EMPTY ? aSideSet->master_type() : aSideSet->slave_type() );
                        }
                        else
                        {
                            return aSideSet->parent()->enforce_linear_interpolation() ?
                                   mesh::number_of_corner_nodes( aSideSet->master_type() ) :
                                   mesh::number_of_nodes( aSideSet->master_type() );
                        }
                    }
                    case( SideSetDofLinkMode::FacetAndSlave ) :
                    {
                        return aSideSet->parent()->enforce_linear_interpolation() ?
                               mesh::number_of_corner_nodes( aSideSet->slave_type() ) :
                               mesh::number_of_nodes( aSideSet->slave_type() );
                    }
                    default :
                    {
                        BELFEM_ERROR( false , "Unsupported mode" );
                        return 0 ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_edges_per_element( SideSet * aSideSet ) const
        {
            // check if this is an empty sideset
            if ( aSideSet->number_of_elements() == 0 )
            {
                return 0 ;
            }
            else
            {
                switch ( this->sideset_dof_link_mode() )
                {
                    case( SideSetDofLinkMode::FacetOnly ) :
                    case( SideSetDofLinkMode::FacetAndSlave ) :
                    {
                        return aSideSet->number_of_edges_per_element() ;
                    }
                    case( SideSetDofLinkMode::FacetAndMaster ) :
                    case( SideSetDofLinkMode::MasterAndSlave ) :
                    {
                        return mesh::number_of_edges( aSideSet->master_type() );
                    }
                    default :
                    {
                        BELFEM_ERROR( false , "Unsupported mode" );
                        return 0 ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_faces_per_element( SideSet * aSideSet ) const
        {
            // check if this is an empty sideset
            if ( aSideSet->number_of_elements() == 0 )
            {
                return 0 ;
            }
            else
            {
                switch ( this->sideset_dof_link_mode() )
                {
                    case( SideSetDofLinkMode::FacetOnly ) :
                    case( SideSetDofLinkMode::FacetAndSlave ) :
                    {
                        return aSideSet->number_of_faces_per_element() ;
                    }
                    case( SideSetDofLinkMode::FacetAndMaster ) :
                    case( SideSetDofLinkMode::MasterAndSlave ) :
                    {
                        return mesh::number_of_faces( aSideSet->master_type() );
                    }
                    default :
                    {
                        BELFEM_ERROR( false , "Unsupported mode" );
                        return 0 ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_dofs_per_element( SideSet * aSideSet ) const
        {
            if( aSideSet->number_of_elements() > 0 )
            {
                uint tIndex = mSideSetIndices( aSideSet->id() ) ;

                if( ! this->has_edge_dofs() )
                {
                    return  this->number_of_nodes_per_element( aSideSet )
                            * mSideSetDofs( tIndex )->Node.length()
                            + mSideSetDofs( tIndex )->Lambda.length() ;
                }
                else if( aSideSet->type() == GroupType::CUT )
                {
                    return 3 * mSideSetDofs( tIndex )->Node.length() ;
                }
                else if ( aSideSet->domain_type() == DomainType::ThinShell )
                {
                    uint tNumberOfDofsPerEdge = mSideSetDofs( tIndex )->Edge.length() * mEdgeDofMultiplicity ;
                    uint tNumberOfDofsPerFace = mSideSetDofs( tIndex )->Face.length() * mFaceDofMultiplicity ;
                    uint tNumberOfVolumeDofs  = ( mesh::number_of_nodes( aSideSet->slave_type() ) +
                                                   mesh::number_of_nodes( aSideSet->master_type() ) ) *
                                                   mSideSetDofs( tIndex )->Node.length() ;
                    uint tNumberOfDofsPerFacet =  mSideSetDofs( tIndex )->Lambda.length() ;

                    if ( tNumberOfDofsPerEdge == 0 && tNumberOfDofsPerFace == 0 )
                    {
                        return tNumberOfVolumeDofs + tNumberOfDofsPerFacet ;
                    }
                    else
                    {
                        return    tNumberOfVolumeDofs + tNumberOfDofsPerFacet
                                  + aSideSet->number_of_thin_shell_layers()  *
                                    ( mesh::number_of_edges( aSideSet->element_type() )
                                    * tNumberOfDofsPerEdge
                                    + mesh::number_of_faces( aSideSet->element_type() )
                                    * tNumberOfDofsPerFace ) ;
                    }
                }
                else
                {
                    // update number of dofs
                    uint tNumberOfDofsPerEdge = mSideSetDofs( tIndex )->Edge.length() * mEdgeDofMultiplicity ;
                    uint tNumberOfDofsPerFace = mSideSetDofs( tIndex )->Face.length() * mFaceDofMultiplicity ;

                    if ( tNumberOfDofsPerEdge == 0 && tNumberOfDofsPerFace == 0 )
                    {
                        return   this->number_of_nodes_per_element( aSideSet )
                                 * mSideSetDofs( tIndex )->Node.length()
                                 + mSideSetDofs( tIndex )->Lambda.length() ;
                    }
                    else
                    {
                        return    mesh::number_of_nodes( aSideSet->slave_type() )
                                  * mSideSetDofs( tIndex )->Node.length()
                                  + mesh::number_of_edges( aSideSet->master_type() )
                                    * tNumberOfDofsPerEdge
                                  + mesh::number_of_faces( aSideSet->master_type() )
                                    * tNumberOfDofsPerFace
                                  + mSideSetDofs( tIndex )->Lambda.length() ;
                    }

                }
            }
            else
            {
                return 0 ;
            }
        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_lambda_dofs( const id_t aSideSetID ) const
        {

            BELFEM_ASSERT( mRank == 0,
                          "this function should only be called by the master proc" );

            return mSideSetDofs(  mSideSetIndices( aSideSetID ) )->Lambda.length() ;
        }

//------------------------------------------------------------------------------

        uint
        IWG::number_of_dofs_per_node() const
        {
            return mNumberOfDofsPerNode;
        }

//------------------------------------------------------------------------------

        uint
        IWG::num_rhs_cols() const
        {
            return mNumberOfRhsCols;
        }

//------------------------------------------------------------------------------

        void
        IWG::set_field( DofManagerBase * aField )
        {
            mField  = aField;
            mMesh   = aField->mesh();
        }

//------------------------------------------------------------------------------

        void
        IWG::collect_node_data(
                Element        * aElement,
                Cell< string > & aFieldLabels,
                Matrix< real > & aData )
        {
            BELFEM_ASSERT(
                    aData.n_rows() == mNumberOfNodesPerElement,
                    "Number of rows for matrix does not fit ( is %u, but need %u )",
                    ( unsigned int ) aData.n_rows(),
                    ( unsigned int ) mNumberOfNodesPerElement );

            BELFEM_ASSERT(
                    aData.n_cols() >= aFieldLabels.size(),
                    "Number of columns for matrix does not fit ( is %u, but need at least %u )",
                    ( unsigned int ) aData.n_cols(),
                    ( unsigned int ) aFieldLabels.size() );

            uint tN = aFieldLabels.size();
            for( uint j=0; j<tN; ++j )
            {
                // get ref to field on mesh
                Vector< real > & tField = mMesh->field_data( aFieldLabels( j ) );

                BELFEM_ASSERT(
                        mMesh->field( aFieldLabels( j ) )->entity_type() == EntityType::NODE,
                        "Field '%s' is not a node field", aFieldLabels( j ).c_str() );

                // loop over all nodes
                for( uint i=0; i< mNumberOfNodesPerElement; ++i )
                {
                    aData( i, j ) = tField( aElement->element()->node( i )->index() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::collect_node_data(
                Element        * aElement,
                const string   & aFieldLabel,
                Vector< real > & aData )
        {
            BELFEM_ASSERT(
                    aData.length() == mNumberOfNodesPerElement,
                    "Number of rows for matrix does not fit ( is %u, but need %u )",
                    ( unsigned int ) aData.length(),
                    ( unsigned int ) mNumberOfNodesPerElement );


            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( aFieldLabel );

            BELFEM_ASSERT(
                    mMesh->field( aFieldLabel )->entity_type() == EntityType::NODE,
                    "Field '%s' is not a node field", aFieldLabel.c_str() );

            // loop over all nodes
            for( uint i=0; i< mNumberOfNodesPerElement; ++i )
            {
                aData( i  ) = tField( aElement->element()->node( i )->index() );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::collect_node_data(
                Element        * aElement,
                const string   & aFieldLabel,
                Vector< real > & aData,
                          uint & aOffset )
        {
            BELFEM_ASSERT(
                    aData.length() >= mNumberOfNodesPerElement,
                    "Number of rows for matrix does not fit ( is %u, but need at least %u )",
                    ( unsigned int ) aData.length(),
                    ( unsigned int ) mNumberOfNodesPerElement );


            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( aFieldLabel );

            BELFEM_ASSERT(
                    mMesh->field( aFieldLabel )->entity_type() == EntityType::NODE,
                    "Field '%s' is not a node field", aFieldLabel.c_str() );

            // loop over all nodes
            for( uint i=0; i< mNumberOfNodesPerElement; ++i )
            {
                aData( aOffset++  ) = tField( aElement->element()->node( i )->index() );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::collect_edge_data(
                Element        * aElement,
                const string   & aEdgeFieldLabel,
                Vector< real > & aData )
        {
            BELFEM_ASSERT( mNumberOfRhsDofsPerEdge == 1,
                          "this function can be used for linear interpolation only ( one dof per edge )" );

            BELFEM_ASSERT( mNumberOfRhsDofsPerFace == 0,
                          "this function can be used for quadratic interpolation only ( zero dofs per face )" );

            BELFEM_ASSERT(
                    aData.length() >= mNumberOfEdgesPerElement,
                    "Length of vector does not fit does not fit ( is %u, but need at least %u )",
                    ( unsigned int ) aData.length(),
                    ( unsigned int ) mNumberOfEdgesPerElement );

            BELFEM_ASSERT(
                    mMesh->field( aEdgeFieldLabel )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeFieldLabel.c_str() );

            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( aEdgeFieldLabel );

            // loop over all edges
            for( uint e=0; e< mNumberOfEdgesPerElement; ++e )
            {
                aData( e ) = tField( aElement->element()->edge( e )->index() );
            }
         }

//------------------------------------------------------------------------------

        void
        IWG::collect_edge_data(
                 Element        * aElement,
                 const string   & aEdgeFieldLabel,
                 const string   & aFaceFieldLabel,
                 Vector< real > & aData )
        {
            BELFEM_ASSERT( mNumberOfRhsDofsPerEdge == 2,
                          "this function can be used for quadratic interpolation only ( two dofs per edge )" );

            BELFEM_ASSERT( mNumberOfRhsDofsPerFace == 2,
                          "this function can be used for quadratic interpolation only ( two dofs per face )" );

            // check length of memory container
            BELFEM_ASSERT(
                    aData.length() >= 2 * ( mNumberOfEdgesPerElement + mNumberOfFacesPerElement ),
                    "Length of vector does not fit does not fit ( is %u, but need at least %u )",
                    ( unsigned int ) aData.length(),
                    ( unsigned int ) 2 * ( mNumberOfEdgesPerElement + mNumberOfFacesPerElement ) );

            // get ref to edge field on mesh
            Vector< real > & tEdgeField = mMesh->field_data( aEdgeFieldLabel );

            // get ref to face field on mesh
            Vector< real > & tFaceField = mMesh->field_data( aFaceFieldLabel );

            // initialize counter
            uint tCount = 0 ;

            for( uint e=0; e< mNumberOfEdgesPerElement; ++e )
            {

                // get index of edge
                index_t tIndex = aElement->element()->edge( e )->index() ;

                // check direction of edge
                if( aElement->edge_direction( e ) )
                {
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex );
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex + 1 );
                }
                else
                {
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex + 1 );
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex  );
                }
            }

            // face dofs orientation are handeled differently in 3D, we can just populate here
            for( uint f=0; f<mNumberOfFacesPerElement; ++f )
            {
                // get index of face
                index_t tIndex = aElement->element()->face( f )->index() ;

                // write data into container
                aData( tCount++ ) = tFaceField( tIndex + tIndex );
                aData( tCount++ ) = tFaceField( tIndex + tIndex + 1 );
            }

        }

//------------------------------------------------------------------------------

        void
        IWG::collect_lambda_data(
             Element        * aElement,
             const string   & aFieldLabel,
             real & aData )
        {
            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( aFieldLabel );

            BELFEM_ASSERT(
                    mMesh->field( aFieldLabel )->entity_type() == EntityType::FACET,
                    "Field '%s' is not a lambda field", aFieldLabel.c_str() );

            aData = tField( aElement->element()->index() );
        }

//------------------------------------------------------------------------------

        void
        IWG::allocate_work_matrices( Group * aGroup )
        {
            // get number of dimensions per element
            mNumberOfSpatialDimensions = mesh::dimension( aGroup->element_type() );

            if ( aGroup->type() == GroupType::SIDESET )
            {
                // get number of nodes per element
                mNumberOfNodesPerElement = this->number_of_nodes_per_element(
                        reinterpret_cast< SideSet * >( aGroup ) );

                // get number of edges per element
                mNumberOfEdgesPerElement = mesh::number_of_edges( aGroup->element_type() );

                aGroup->node_coords().set_size(
                        mNumberOfNodesPerElement,
                        mNumberOfSpatialDimensions + 1 );
            }
            else if ( aGroup->type() == GroupType::BLOCK )
            {
                // get number of nodes per element
                mNumberOfNodesPerElement = mesh::number_of_nodes( aGroup->element_type() );

                // get number of edges per element
                mNumberOfEdgesPerElement = mesh::number_of_edges( aGroup->element_type() );

                // node coordinate vector
                aGroup->node_coords().set_size(
                        mNumberOfNodesPerElement,
                        mNumberOfSpatialDimensions );
            }

            // Geometry Jacobian
            aGroup->work_J().set_size( mNumberOfSpatialDimensions,
                                        mNumberOfSpatialDimensions,
                                        0.0 );

            // matrices for second derivative
            aGroup->work_K().set_size(
                    mNumberOfDerivativeDimensions,
                    mNumberOfSpatialDimensions,
                       0.0 );

            aGroup->work_L().set_size( mNumberOfDerivativeDimensions,
                    mNumberOfDerivativeDimensions,
                    0.0 );

            // for elasticity
            aGroup->work_C().set_size(
                    mNumberOfDerivativeDimensions,
                    mNumberOfDerivativeDimensions,
                    0.0 );

            // gradient matrix
            aGroup->work_dNdX().set_size(
                    mNumberOfSpatialDimensions,
                    mNumberOfNodesPerElement );

            // field data
            aGroup->work_phi().set_size(
                    mNumberOfNodesPerElement,
                    0.0 );

            if ( mNumberOfDofsPerNode > 0 )
            {
                aGroup->work_Phi().set_size(
                        mNumberOfNodesPerElement,
                        mNumberOfDofsPerNode,
                        0.0 );

                // N-Matrix for element assembly
                aGroup->work_N().set_size( mNumberOfSpatialDimensions,
                        mNumberOfDofsPerNode * mNumberOfNodesPerElement,
                        0.0 );

                // B-Matrix for element assembly
                aGroup->work_B().set_size( mNumberOfDerivativeDimensions,
                                           mNumberOfDofsPerNode * mNumberOfNodesPerElement,
                                           0.0 );
            }

            if( aGroup->type() == GroupType::SIDESET )
            {
                // get facets of this sideset
                Cell< mesh::Facet * > & tFacets = mMesh->sideset( aGroup->id() )->facets() ;

                // make sure that this is not an empty sideset
                if( tFacets.size() > 0 )
                {
                    // get element type of master
                    ElementType tType = tFacets( 0 )->master()->type() ;

                    // get the number of nodes for this element
                    uint tNumNodes = mesh::number_of_nodes( tType );

                    // allocate memory
                    aGroup->work_Xi().set_size( tNumNodes,
                                               mNumberOfSpatialDimensions+1 );

                    // allocate memory
                    aGroup->work_X().set_size( mNumberOfNodesPerElement,
                                                mNumberOfSpatialDimensions+1 );
                }
            }
            else
            {
                aGroup->work_X().set_size( mNumberOfNodesPerElement,
                                           mNumberOfSpatialDimensions  );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::set_omega( const real & aOmega )
        {
            if( mField->parent()->master() == mField->rank() )
            {
                mOmega = aOmega;

                Vector< real > tOmega(
                        mField->parent()->comm_table().length(),
                        aOmega );

                send( mField->parent()->comm_table(), tOmega );
            }
            else
            {
                receive( mField->parent()->master(), mOmega );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::set_psi( const real & aPsi )
        {
            if( mField->parent()->master() == mField->rank() )
            {
                mPsi = aPsi;

                Vector< real > tPsi(
                        mField->parent()->comm_table().length(),
                        aPsi );

                send( mField->parent()->comm_table(), tPsi );
            }
            else
            {
                receive( mField->parent()->master(), mPsi );
            }
        }


//------------------------------------------------------------------------------

        const real &
        IWG::timestep() const
        {
            return mDeltaTime;
        }

//------------------------------------------------------------------------------

        void
        IWG::add_fields( const Cell< string > & aFieldLabels )
        {
            // get number of fields
            uint tNumFields = aFieldLabels.size() ;

            // loop over all fields
            for( uint k=0; k<tNumFields; ++k )
            {
                // get string
                const string & tLabel = aFieldLabels( k );

                // set the alpha flag if this is an alpha-heatflux-boundary condition
                if( tLabel == "alpha" )
                {
                    mHasAlpha = true ;
                }

                // check if field exists on field list
                bool tFlag = false ;
                for( string & tField : mOtherFields )
                {
                    if( tField == tLabel )
                    {
                        tFlag = true ;
                        break ;
                    }
                }

                // add field to others and all
                if( ! tFlag )
                {
                    mOtherFields.push( tLabel );
                    mAllFields.push( tLabel );
                }
            }
        }

//------------------------------------------------------------------------------

        // convection term with alpha boundary condition
        void
        IWG::compute_alpha_boundary_condition(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false,
                         "this boundary condition is not implemented for this IWG");
        }

//------------------------------------------------------------------------------

        void
        IWG::initialize()
        {
            // collect mAllFields
            this->concatenate_field_lists();

            // check if dofmap has been built
            if( mDofMap.size() == 0 )
            {
                this->unique_and_rearrange( mDofFields, true );
            }
            this->create_doftype_map();
            this->count_dofs_per_block();
            this->count_dofs_per_sideset();
            mIsInitialized = true ;
        }

//------------------------------------------------------------------------------

        void
        IWG::create_doftype_map()
        {
            // reset the map
            mDofTypeMap.clear() ;

            // initialize counter
            uint tCount = 0 ;
            uint tEdgeCount = 0 ;
            uint tFaceCount = 0 ;

            mEdgeFieldIndices.set_size( BELFEM_MAX_DOFTYPES, gNoIndex );
            mFaceFieldIndices.set_size( BELFEM_MAX_DOFTYPES, gNoIndex );

            // remember #shell here
            for( string tDofLabel : mDofFields )
            {
                mDofTypeMap[ tDofLabel ] = tCount ;
                switch( entity_type( tDofLabel ) )
                {
                    case( EntityType::NODE ) :
                    {
                        // node values don't have a multiplicity
                        ++tCount ;
                        break ;
                    }
                    case( EntityType::EDGE ) :
                    {
                        for( uint i=0; i<mEdgeDofMultiplicity ; ++i )
                        {
                            mEdgeFieldIndices( tCount++ ) = tEdgeCount++;
                        }
                        break ;
                    }
                    case( EntityType::FACE ) :
                    {
                        for( uint i=0; i<mFaceDofMultiplicity ; ++i )
                        {
                            mFaceFieldIndices( tCount++ ) = tFaceCount++;
                        }
                        break ;
                    }
                    case( EntityType::CELL ) :
                    {
                        for( uint i=0; i<mCellDofMultiplicity ; ++i )
                        {
                            // mCellFieldIndices( tCount++ ) = ...
                            ++tCount ;
                        }
                        break ;
                    }
                    case( EntityType::FACET ) :
                    {
                        for( uint i=0; i<mLambdaDofMultiplicity ; ++i )
                        {
                            // mLambdaFieldIndices( tCount++ ) = ...
                            ++tCount ;
                        }
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid entity type.");
                    }
                }

            }

            BELFEM_ERROR( tEdgeCount <= mEdgeDofMultiplicity, "can't have more than one edge field" );
            BELFEM_ERROR( tFaceCount <= mFaceDofMultiplicity, "can't have more than one face field" );

            mDefaultDofTypes.set_size( tCount );

            // this is going to be a consecutive array
            for( uint k=0; k<tCount; ++k )
            {
                mDefaultDofTypes( k ) = k ;
            }

            mDofEntityTypes.set_size( tCount , 0 );

            // reset counter
            tCount = 0 ;

            for( string tDofLabel : mDofFields )
            {
                EntityType tType = entity_type( tDofLabel ) ;
                uint       tUType = ( uint ) tType ;


                switch( tType )
                {
                    case( EntityType::NODE ) :
                    {
                        // nodes don't have a multiplicity
                        mDofEntityTypes( tCount++ ) = tUType ;
                        break ;
                    }
                    case( EntityType::EDGE ) :
                    {
                        for( uint k=0; k<mEdgeDofMultiplicity; ++k )
                        {
                            mDofEntityTypes( tCount++ ) = tUType ;
                        }
                        break ;
                    }
                    case( EntityType::FACE ) :
                    {
                        for( uint k=0; k<mFaceDofMultiplicity; ++k )
                        {
                            mDofEntityTypes( tCount++ ) = tUType ;
                        }
                        break ;
                    }
                    case( EntityType::CELL ) :
                    {
                        for( uint k=0; k<mCellDofMultiplicity; ++k )
                        {
                            mDofEntityTypes( tCount++ ) = tUType ;
                        }
                        break ;
                    }
                    case( EntityType::FACET ) :
                    {
                        for( uint k=0; k<mLambdaDofMultiplicity; ++k )
                        {
                            mDofEntityTypes( tCount++ ) = tUType ;
                        }
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid entity type.");
                    }
                }
            }

            // we don't support more than 32 dofs per entity
            BELFEM_ERROR( tCount< BELFEM_MAX_DOFTYPES,
                         "This program is not supposed to have more than %u doftypes per problem. \nRedefine BELFEM_MAX_DOFTYPES in cl_IWG.hpp and compile again",
                         ( unsigned int ) BELFEM_MAX_DOFTYPES );
        }

//------------------------------------------------------------------------------

        void
        IWG::concatenate_field_lists()
        {
            mAllFields.clear() ;
            if ( mTensorFields.size() > 0 )
            {
                for ( string tLabel : mTensorFields )
                {
                    mAllFields.push( tLabel );
                }
            }
            else
            {
                for ( string tLabel : mDofFields )
                {
                    mAllFields.push( tLabel );
                }
            }
            for ( string tLabel : mFluxFields )
            {
                mAllFields.push( tLabel );
            }
            for ( string tLabel : mOtherFields )
            {
                mAllFields.push( tLabel );
            }
        }

//---------------------------------------------------------------------------------

        void
        IWG::create_block_dof_tables( const uint aNumBlocks )
        {
            // reset
            this->delete_block_dof_tables() ;

            mBlockDofs.set_size( aNumBlocks, nullptr );
            for( uint b=0; b<aNumBlocks; ++b )
            {
                mBlockDofs( b ) = new DofTable();
            }
        }

//---------------------------------------------------------------------------------

        void
        IWG::create_sideset_dof_tables(const uint aNumSideSets, const uint aNumGhostSidesets )
        {
            // reset
            this->delete_sideset_dof_tables() ;

            mSideSetDofs.set_size( aNumSideSets + aNumGhostSidesets, nullptr );
            uint tCount = 0 ;
            for( uint s=0; s<aNumSideSets; ++s )
            {
                mSideSetDofs( tCount++ ) = new DofTable();
            }

            // create doppelgangers for ghosts
            if( mGhostSideSetMaster > 0 )
            {
                for ( uint g = 0; g < aNumGhostSidesets; ++g )
                {
                    mSideSetDofs( tCount++ ) = mSideSetDofs( mSideSetIndices( mGhostSideSetMaster ) );
                }
            }
        }

//------------------------------------------------------------------------------


        void
        IWG::assign_dofs_per_block( const Vector< id_t > & aBlockIDs )
        {
            // assume that all dofs sit on all blocks
            mDofsPerBlock.set_size( aBlockIDs.length(), mDefaultDofTypes );
        }

//------------------------------------------------------------------------------


        void
        IWG::assign_dofs_per_sideset( const Vector< id_t > & aSideSetIDs )
        {
            // assume that all dofs sit on all sidesets
            mDofsPerSideSet.set_size( aSideSetIDs.length(), mDefaultDofTypes );
        }


//------------------------------------------------------------------------------

        void
        IWG::count_dofs_per_block()
        {
            // #reindex
            // counter for all types
            Vector< uint > tCount( 4 );

            // get number of blocks
            uint tNumBlocks = mDofsPerBlock.size() ;

            // loop over all blocks
            for( uint b=0; b<tNumBlocks; ++b )
            {
                // reset counter
                tCount.fill( 0 );

                Vector< index_t > & tDofsPerBlock = mDofsPerBlock( b );


                // loop over all dofs
                for( uint tDofNum : tDofsPerBlock )
                {
                    // get the Entity Type
                    EntityType tEntityType = entity_type( mDofLabels( tDofNum ) );

                    // check what it is
                    switch ( tEntityType )
                    {
                        case( EntityType::NODE ) :
                        {
                            ++tCount( 0 );
                            break ;
                        }
                        case( EntityType::EDGE ) :
                        {
                            ++tCount( 1 );
                            break ;
                        }
                        case( EntityType::FACE ) :
                        {
                            ++tCount( 2 );
                            break ;
                        }
                        case( EntityType::CELL ) :
                        {
                            ++tCount( 3 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Unknown Entity Type");
                        }
                    }
                }

                // grab vectors
                Vector< index_t > & tDofsPerNode = mBlockDofs( b )->Node ;
                Vector< index_t > & tDofsPerEdge = mBlockDofs( b )->Edge ;
                Vector< index_t > & tDofsPerFace = mBlockDofs( b )->Face ;
                Vector< index_t > & tDofsPerCell = mBlockDofs( b )->Cell ;

                // allocate memory
                tDofsPerNode.set_size( tCount( 0 ) );
                tDofsPerEdge.set_size( tCount( 1 ) );
                tDofsPerFace.set_size( tCount( 2 ) );
                tDofsPerCell.set_size( tCount( 3 ) );

                // reset counter
                tCount.fill( 0 );

                // loop over all dofs
                for( uint tDofNum : tDofsPerBlock )
                {
                    // get the label of the dof
                    const string & tDofLabel = mDofLabels( tDofNum );

                    // get the Entity Type
                    EntityType tEntityType =  entity_type( tDofLabel );

                    // get the type number of the dof
                    uint tDofType = mDofMap( tDofLabel );

                    // check what it is
                    switch ( tEntityType )
                    {
                        case( EntityType::NODE ) :
                        {
                            tDofsPerNode( tCount( 0 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::EDGE ) :
                        {
                            tDofsPerEdge( tCount( 1 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::FACE ) :
                        {
                            tDofsPerFace( tCount( 2 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::CELL ) :
                        {
                            tDofsPerCell( tCount( 3 )++ ) = tDofType ;
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Unknown Entity Type");
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::count_dofs_per_sideset()
        {
            // #reindex

            // counter for all types
            Vector< uint > tCount( 5 );

            // get number of blocks
            uint tNumSideSets = mDofsPerSideSet.size() ;

            // loop over all blocks
            for( uint s=0; s<tNumSideSets; ++s )
            {
                // reset counter
                tCount.fill( 0 );

                Vector< index_t > & tDofsPerSideSet = mDofsPerSideSet( s );

                // loop over all dofs
                for( uint tDofNum : tDofsPerSideSet )
                {
                    // get the Entity Type
                    EntityType tEntityType = ( EntityType ) mDofEntityTypes( tDofNum ) ;

                    // check what it is
                    switch ( tEntityType )
                    {
                        case( EntityType::NODE ) :
                        {
                            ++tCount( 0 );
                            break ;
                        }
                        case( EntityType::EDGE ) :
                        {
                            ++tCount( 1 );
                            break ;
                        }
                        case( EntityType::FACE ) :
                        {
                            ++tCount( 2 );
                            break ;
                        }
                        case( EntityType::CELL ) :
                        {
                            ++tCount( 3 );
                            break ;
                        }
                        case( EntityType::FACET ) : // lambda dof
                        {
                            ++tCount( 4 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Unknown Entity Type");
                        }
                    }
                }

                // grab vectors
                Vector< index_t > & tDofsPerNode = mSideSetDofs( s )->Node ;
                Vector< index_t > & tDofsPerEdge = mSideSetDofs( s )->Edge ;
                Vector< index_t > & tDofsPerFace = mSideSetDofs( s )->Face ;
                Vector< index_t > & tDofsPerCell = mSideSetDofs( s )->Cell ;
                Vector< index_t > & tLambdaDofs  = mSideSetDofs( s )->Lambda ;

                // allocate memory
                tDofsPerNode.set_size( tCount( 0 ) );
                tDofsPerEdge.set_size( tCount( 1 ) );
                tDofsPerFace.set_size( tCount( 2 ) );
                tDofsPerCell.set_size( tCount( 3 ) );
                tLambdaDofs.set_size( tCount( 4 ) );

                // reset counter
                tCount.fill( 0 );

                // loop over all dofs
                for( uint tDofNum : tDofsPerSideSet )
                {
                    // get the label of the dof
                    const string & tDofLabel = mDofLabels( tDofNum );

                    // get the Entity Type
                    EntityType tEntityType =  entity_type( tDofLabel );

                    // get the type number of the dof
                    uint tDofType = mDofMap( tDofLabel );

                    // check what it is
                    switch ( tEntityType )
                    {
                        case( EntityType::NODE ) :
                        {
                            tDofsPerNode( tCount( 0 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::EDGE ) :
                        {
                            tDofsPerEdge( tCount( 1 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::FACE ) :
                        {
                            tDofsPerFace( tCount( 2 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::CELL ) :
                        {
                            tDofsPerCell( tCount( 3 )++ ) = tDofType ;
                            break ;
                        }
                        case( EntityType::FACET ) : // lambda dof
                        {
                            tLambdaDofs( tCount( 4 )++ ) = tDofType ;
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Unknown Entity Type");
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_boundary_flux_matrix(
                Element        * aElement,
                const uint       aDirection,
                Matrix< real > & aJacobian )
        {
            BELFEM_ERROR( false,
                         "compute_boundary_flux_matrix() not implemented for this IWG");
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG::N( const uint & aIntegrationPoint )
        {
            Matrix< real > & aN = mGroup->work_N();
            const Vector< real > & tN = mGroup->n( aIntegrationPoint );

            aN.fill( 0.0 );

            uint j=0;

            for( uint k=0; k<mNumberOfNodesPerElement; ++k )
            {
                for( uint i=0; i<mNumberOfDofsPerNode; ++i )
                {
                    aN( i, j++ ) = tN( k );
                }
            }

            return aN ;
        }
//------------------------------------------------------------------------------

        void
        IWG::set_blocks(
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
        IWG::set_sidesets(
                const Vector< id_t > & aSideSetIDs,
                const Cell< DomainType > & aSideSetTypes )
        {
            mSideSetIDs = aSideSetIDs;
            mSideSetTypes.clear();

            index_t tCount = 0 ;
            for ( id_t tID : aSideSetIDs )
            {
                mSideSetTypes[ tID ] = aSideSetTypes( tCount++ );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::set_ghost_sidesets( const id_t                aMasterSideSet,
                                 const Vector< id_t >    & aGhostSideSet )
        {
            mGhostSideSetMaster = aMasterSideSet,
            mGhostSideSetIDs    = aGhostSideSet ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG::node_coords_from_master( Element * aElement )
        {
            BELFEM_ASSERT( mGroup->type() == GroupType::SIDESET,
                          "IWG must be linked to sideset when requesting the master node coordinates of an element" );

            // grab container
            Matrix< real > & aX = mGroup->work_Xi();

            // get master element
            mesh::Element * tMaster = aElement->element()->element( 0 );


            for( uint i=0; i<=mNumberOfSpatialDimensions; ++i )
            {
                for( uint k=0; k<tMaster->number_of_nodes(); ++k )
                {
                    aX( k, i ) = tMaster->node( k )->x( i );
                }
            }

            return aX ;
        }

//------------------------------------------------------------------------------

        void
        IWG::hide_fields_from_exodus( Mesh * aMesh )
        {
            const uint tNumFields = mHiddenFields.size();

            for( uint f=0; f<tNumFields; ++f )
            {
                aMesh->field( mHiddenFields( f ) )
                    ->set_write_to_file_flag( false );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::shift_fields()
        {
            BELFEM_ERROR( false, "shift_fields() not implemented for this IWG");
        }

//------------------------------------------------------------------------------

        void
        IWG::add_boundary_condition( BoundaryCondition * aBoundaryCondition )
        {
            mBoundaryConditions.push( aBoundaryCondition );
        }

//------------------------------------------------------------------------------

        void
        IWG::delete_boundary_conditions()
        {
            for ( BoundaryCondition * tBC : mBoundaryConditions )
            {
                delete tBC ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::delete_block_dof_tables()
        {
            for( DofTable * tTable : mBlockDofs )
            {
                delete tTable ;
            }
            mBlockDofs.clear();
        }

//------------------------------------------------------------------------------

        void
        IWG::delete_sideset_dof_tables()
        {
            // note, this array may also contain doppelgangers for ghosts
            // this is why we do the delete loop this way
            if( mSideSetDofs.size() > 0 )
            {
                for ( uint k = 0; k < mSideSetIDs.length(); ++k )
                {
                    delete mSideSetDofs( k );
                }
            }
            mSideSetDofs.clear();
        }

//------------------------------------------------------------------------------

        void
        IWG::compute_boundary_conditions( const real aTimeStamp )
        {
            for ( BoundaryCondition * tBC : mBoundaryConditions )
            {
                tBC->compute( aTimeStamp );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::unique_and_rearrange( Cell< string > & aDofs, const bool aMakeMap )
        {
            // catch special case
            if ( aDofs.size() == 0 ) return ;

            // make dofs unique
            unique( aDofs );

            // temporary cells
            Cell< string > tNodeDofs ;
            Cell< string > tFaceDofs ;
            Cell< string > tEdgeDofs ;
            Cell< string > tCellDofs ;
            Cell< string > tLambdaDofs ;

            for( const string & tDof : aDofs )
            {
                // get type of entity
                switch( entity_type( tDof ) )
                {
                    case( EntityType::NODE ) :
                    {
                        tNodeDofs.push( tDof );
                        break ;
                    }
                    case( EntityType::EDGE ) :
                    {
                        tEdgeDofs.push( tDof );
                        break ;
                    }
                    case( EntityType::FACE ) :
                    {
                        tFaceDofs.push( tDof );
                        break ;
                    }
                    case( EntityType::CELL ) :
                    {
                        tCellDofs.push( tDof );
                        break ;
                    }
                    case( EntityType::FACET ) :
                    {
                        tLambdaDofs.push( tDof );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid entity type for dof %s", tDof.c_str() );
                    }
                }
            }

            uint tCount = aDofs.size() ;

            // now we rearrange the dofs, beginning with the nodes
            aDofs.set_size( tCount, "" );
            tCount = 0 ;

            for( const string & tDof : tNodeDofs )
            {
                aDofs( tCount++ ) = tDof ;
            }
            for( const string & tDof : tLambdaDofs )
            {
                aDofs( tCount++ ) = tDof ;
            }
            for( const string & tDof : tEdgeDofs )
            {
                aDofs( tCount++ ) = tDof ;
            }
            for( const string & tDof : tFaceDofs )
            {
                aDofs( tCount++ ) = tDof ;
            }
            for( const string & tDof : tCellDofs )
            {
                aDofs( tCount++ ) = tDof ;
            }

            BELFEM_ASSERT( tCount = aDofs.size(), "Error when rearranging dof cell");


            // special function to create labels for printout
            if( aMakeMap )
            {
                // reset counter
                tCount = 0 ;

                // reset map
                mDofMap.clear();

                // reset label cell
                mDofLabels.clear();

                for( const string & tDof : tNodeDofs )
                {
                    mDofLabels.push( tDof );
                    mDofMap[ tDof ] = tCount++ ;
                }

                for( const string & tDof : tLambdaDofs )
                {
                    mDofMap[ tDof ] = tCount ;
                    tCount += this->lambda_multiplicity();
                    BELFEM_ASSERT( this->lambda_multiplicity() > 0, "can't add lambda dof if multilplicity is zero" );

                    for( uint k=0; k<this->lambda_multiplicity(); ++k )
                    {
                        mDofLabels.push( tDof );
                    }
                }
                for( const string & tDof : tEdgeDofs )
                {
                    mDofMap[ tDof ] = tCount ;
                    tCount += this->edge_multiplicity() ;
                    BELFEM_ASSERT( this->edge_multiplicity() > 0, "can't add edge dof if multilplicity is zero" );

                    for( uint k=0; k<this->edge_multiplicity(); ++k )
                    {
                        mDofLabels.push( tDof );
                    }
                }
                for( const string & tDof : tFaceDofs )
                {
                    mDofMap[ tDof ] = tCount ;
                    tCount += this->face_multiplicity();
                    BELFEM_ASSERT( this->face_multiplicity() > 0, "can't add face dof if multilplicity is zero" );

                    for( uint k=0; k<this->face_multiplicity(); ++k )
                    {
                        mDofLabels.push( tDof );
                    }
                }
                for( const string & tDof : tCellDofs )
                {
                    mDofMap[ tDof ] = tCount ;
                    tCount += this->cell_multiplicity();
                    for( uint k=0; k<this->cell_multiplicity(); ++k )
                    {
                        mDofLabels.push( tDof );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG::print_dofs( Element * aElement )
        {
            std::cout << "Element : " << aElement->id() << std::endl ;

            unsigned int tCount = 0 ;
            uint tNumDofs = aElement->number_of_dofs() ;

            for( uint d=0; d<tNumDofs; ++d )
            {
                Dof * tDof = aElement->dof( d );

                // get label
                const string & tLabel = mDofLabels( tDof->type_id() ) ;

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

                fprintf( stdout, "    %2u  %s%s id %8u     index %8u     %s %8u  field %8u\n",
                         tCount++,
                         tLabel.c_str(),
                         tSpace.c_str(),
                         ( unsigned int ) tDof->id(),
                         ( unsigned int ) tDof->index(),
                         tType.c_str(),
                         ( unsigned int ) tDof->mesh_vertex()->id(),
                         tDof->dof_index_on_field()
                );
            }
        }


//---------------------------------------------------------------------------------

        DomainType
        IWG::block_type( const id_t aID ) const
        {
            if( mBlockTypes.key_exists( aID ) )
            {
                return mBlockTypes( aID );
            }
            else
            {
                return DomainType::UNDEFINED ;
            }
        }

//------------------------------------------------------------------------------

        DomainType
        IWG::sideset_type( const id_t aID ) const
        {
            if( mSideSetTypes.key_exists( aID ) )
            {
                return mSideSetTypes( aID );
            }
            else
            {
                return DomainType::UNDEFINED ;
            }
        }

//------------------------------------------------------------------------------
    }
}

