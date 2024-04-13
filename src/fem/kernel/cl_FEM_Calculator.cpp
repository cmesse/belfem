//
// Created by christian on 6/9/23.
//
#include <iostream>

#include "assert.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_Calculator.hpp"
#include "meshtools.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_entity_type.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
        namespace calculator
        {
//------------------------------------------------------------------------------

            VectorData::VectorData( const string & aLabel,
                                    uint aSize,
                                    const EntityType aType )  :
                                    mLabel( aLabel ),
                                    mType( aType )
            {
                mVectorData.set_size( aSize, BELFEM_QUIET_NAN );
            }

//------------------------------------------------------------------------------

            MatrixData::MatrixData( const string & aLabel,
                                    const uint aNumRows,
                                    const uint aNumCols ) :
                                    mLabel( aLabel )
            {
                mMatrixData.set_size( aNumRows, aNumCols, BELFEM_QUIET_NAN );
            }

        }
//------------------------------------------------------------------------------

        Calculator::Calculator( Group * aGroup, const ModelDimensionality aDimensionality ) :
                mGroup( aGroup ),
                mMesh( aGroup->parent()->mesh() ),
                mDimensionality( aDimensionality )
        {

        }

//------------------------------------------------------------------------------

        Calculator::~Calculator()
        {
            // delete calculator matrices
            for( calculator::MatrixData * tMatrix : mMatrices )
            {
                delete tMatrix ;
            }

            // delete calculator vectors
            for( calculator::VectorData * tVector : mVectors )
            {
                delete tVector ;
            }


            if( mDomainIntegration != nullptr )
            {
                delete mDomainIntegration ;
            }
            if ( mLinearIntegration != nullptr )
            {
                delete mLinearIntegration ;
            }
        }

//------------------------------------------------------------------------------


        void
        Calculator::initialize_integration(
                const ElementType       aElementType,
                const InterpolationType aInterpolationType )
        {
            if( mDomainIntegration != nullptr )
            {
                delete mDomainIntegration ;
            }
            mDomainIntegration = new IntegrationData( aElementType,
                                                aInterpolationType );

            if ( mLinearIntegration != nullptr )
            {
                delete mLinearIntegration ;
            }

            mLinearIntegration = new IntegrationData( mesh::linear_element_type( aElementType ),
                                                      aInterpolationType );

        }

//------------------------------------------------------------------------------

        void
        Calculator::set_integration_order( const uint aOrder )
        {
            mNumberOfNodes = mesh::number_of_nodes( mGroup->element_type() );
            mNumberOfCornerNodes = mesh::number_of_corner_nodes( mGroup->element_type() );

            if( mDomainIntegration  != nullptr )
            {
                mDomainIntegration->populate( aOrder, mGroup->parent()->integration_scheme());
            }
            if( mLinearIntegration != nullptr )
            {
                mLinearIntegration->populate( aOrder, mGroup->parent()->integration_scheme());
            }
            this->allocate_memory() ;

            mNumberOfIntegrationPoints = mDomainIntegration->weights().length() ;
        }

//------------------------------------------------------------------------------

        void
        Calculator::allocate_memory()
        {

            if( mGroup->parent()->iwg() == nullptr )
            {
                return  ;
            }

            if(
                    mGroup->type() == GroupType::BLOCK ||
                    mGroup->type() == GroupType::SIDESET ||
                    mGroup->type() == GroupType::SHELL )
            {
                // get pointer to equation object
                IWG * tEquation = mGroup->parent()->iwg();

                // reset the list
                mVectors.clear() ;

                // reset the map
                mVectorMap.clear() ;

                // flags for special cases
                bool tHaveEdgeH = false ;
                bool tHaveFaceH = false ;
                bool tHaveCellH = false ;
                bool tHaveEdgeA = false ;
                bool tHaveFaceA = false ;
                bool tHaveCellA = false ;
                for ( const string & tLabel: tEquation->all_fields())
                {
                    // the size of the vector
                    uint tSize = 0;

                    EntityType tType = entity_type( tLabel );

                    if( tLabel == "edge_a" )
                    {
                        tHaveEdgeA = true ;
                    }
                    if( tLabel == "face_a" )
                    {
                        tHaveFaceA = true ;
                    }
                    if( tLabel == "cell_a" )
                    {
                        tHaveCellA = true ;
                    }

                    if( tLabel == "edge_h" )
                    {
                        tHaveEdgeH = true ;
                    }
                    if( tLabel == "face_h" )
                    {
                        tHaveFaceH = true ;
                    }
                    if( tLabel == "cell_h" )
                    {
                        tHaveCellH = true ;
                    }


                    // determine size
                    switch ( tType )
                    {
                        case ( EntityType::EDGE ) :
                        {
                            tSize = tEquation->edge_multiplicity() * mesh::number_of_edges( mGroup->element_type());
                        }
                        case ( EntityType::FACE ) :
                        {
                            tSize = tEquation->face_multiplicity() * mesh::number_of_faces( mGroup->element_type());
                        }
                        case ( EntityType::CELL ) :
                        {
                            tSize = tEquation->cell_multiplicity();
                        }
                        case ( EntityType::FACET ) :
                        {
                            tSize = tEquation->lambda_multiplicity();
                        }
                        default:
                        {
                            tSize = mesh::number_of_nodes( mGroup->element_type());
                        }
                    }

                    // allocate size
                    this->create_vector( tLabel, tSize, tType );
                }

                // catch special cases for higher order Nedelecs
                // catch special case for higher order Nedelecs
                if( tHaveEdgeA && tHaveFaceA )
                {
                    uint tSize =
                            tEquation->edge_multiplicity() * mesh::number_of_edges( mGroup->element_type())
                            + tEquation->face_multiplicity() * mesh::number_of_faces( mGroup->element_type());

                    if( tHaveCellA )
                    {
                        tSize += tEquation->cell_multiplicity();
                    }

                    this->create_vector( "nedelec_a", tSize, EntityType::UNDEFINED );
                    this->create_vector( "nedelec_a0", tSize, EntityType::UNDEFINED );
                }

                if( tHaveEdgeH && tHaveFaceH )
                {
                    uint tSize =
                             tEquation->edge_multiplicity() * mesh::number_of_edges( mGroup->element_type())
                           + tEquation->face_multiplicity() * mesh::number_of_faces( mGroup->element_type());

                    if( tHaveCellH )
                    {
                        tSize += tEquation->cell_multiplicity();
                    }

                    this->create_vector( "nedelec_h", tSize, EntityType::UNDEFINED );
                    this->create_vector( "nedelec_h0", tSize, EntityType::UNDEFINED );
                }
            }
        }

//------------------------------------------------------------------------------


        IntegrationData *
        Calculator::slave_integration_2d( const Element * aElement )
        {
            return mGroup->slave_integration( aElement->facet()->slave_index() );
        }

        IntegrationData *
        Calculator::slave_integration_tet( const Element * aElement )
        {
            return mGroup->slave_integration(
                    aElement->facet()->slave_index() * 3
                    + aElement->facet()->orientation_on_slave() );
        }

        IntegrationData *
        Calculator::slave_integration_hex( const Element * aElement )
        {
            return mGroup->slave_integration(
                    aElement->facet()->slave_index() * 4
                    + aElement->facet()->orientation_on_slave() );
        }

//------------------------------------------------------------------------------

        void
        Calculator::allocate()
        {

            if(  mGroup->number_of_elements() == 0 || this->integration() == nullptr )
            {
                return;
            }

            // make sure that allocation funciton is only called once
            BELFEM_ERROR( ! mIsAllocated, "calculator is already allocated" );
            mIsAllocated = true ;

            // number of nodes per element
            mNumberOfNodes = mesh::number_of_nodes( mGroup->element_type() );

            // corner nodes per element
            mNumberOfCornerNodes = mesh::number_of_corner_nodes( mGroup->element_type() );

            // number of edges per element
            //uint tNumEdges = mesh::number_of_edges( mGroup->element_type() );

            // number of faces per element
            //uint tNumFaces = mesh::number_of_faces( mGroup->element_type() );

            // uint tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->element_type() );

            // get the number of dimensions
            uint tNumDimensions = mMesh->number_of_dimensions() ;

            // get the number of dofs from first element
            uint tNumDofs = mGroup->elements()(0)->number_of_dofs() ;

            // matrices for X-Coordinates
            mX.set_size( mNumberOfNodes, tNumDimensions, BELFEM_QUIET_NAN );
            mXc.set_size( mNumberOfCornerNodes, tNumDimensions, BELFEM_QUIET_NAN );

            // matrix for Jacobian and its inverse
            mJ = this->create_matrix( "J", tNumDimensions, tNumDimensions );
            mInvJ = this->create_matrix( "InvJ", tNumDimensions, tNumDimensions );

            // special node and matrix interpolators
            switch( mGroup->parent()->iwg()->type() )
            {
                case( IwgType::Gradient2D ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 2, tNumDofs );
                    mN->matrix().fill( 0.0 );

                    // link functions
                    mFunN = & Calculator::N2D ;
                    mFunB = & Calculator::Bscalar ;

                    mB = this->create_matrix( "B", tNumDimensions, mNumberOfNodes );
                    mB->matrix().fill( 0.0 );

                    break;
                }
                case( IwgType::Gradient3D ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 3, tNumDofs );
                    mN->matrix().fill( 0.0 );

                    // link functions
                    mFunN = & Calculator::N3D ;
                    mFunB = & Calculator::Bscalar ;

                    mB = this->create_matrix( "B", tNumDimensions, mNumberOfNodes );
                    mB->matrix().fill( 0.0 );

                    break;
                }
                case( IwgType::PlaneStress ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 2, tNumDofs );
                    mN->matrix().fill( 0.0 );

                    // matrices for gradient operator
                    mdN = this->create_matrix( "dN", 2, 2*mNumberOfNodes );
                    mB = this->create_matrix( "B", 3, tNumDofs );
                    mB->matrix().fill( 0.0 );

                    // link functions
                    mFunN = & Calculator::N2D ;
                    mFunB = & Calculator::Bplanestress ;
                    break ;
                }
                case( IwgType::LinearElasticity ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 3, tNumDofs );
                    mN->matrix().fill( 0.0 );

                    // matrices for gradient operator
                    mdN = this->create_matrix( "dN", 3, 3*mNumberOfNodes );
                    mB = this->create_matrix( "B", 6, tNumDofs );
                    mB->matrix().fill( 0.0 );

                    // link functions
                    mFunN = & Calculator::N3D ;
                    mFunB = & Calculator::Bvoigt ;
                    break ;
                }
                default:
                {
                    mFunN = & Calculator::Nscalar ;
                    mFunB = & Calculator::Bscalar ;
                    mB = this->create_matrix( "B", tNumDimensions, mNumberOfNodes );

                    mB->matrix().fill( 0.0 );

                    break ;
                }
            }

            // stiffness matrix
            mK.set_size( tNumDofs, tNumDofs, BELFEM_QUIET_NAN );

            // load vector
            mf.set_size( tNumDofs, BELFEM_QUIET_NAN );

            // dof vectors
            mq0.set_size( tNumDofs, BELFEM_QUIET_NAN );
            mq1.set_size( tNumDofs, BELFEM_QUIET_NAN );
            mq.set_size( tNumDofs, BELFEM_QUIET_NAN );

            // link function to invert J
            switch( tNumDimensions )
            {
                case( 2 ) :
                {
                    mFunInvertJ = & inv2 ;
                    break ;
                }
                case( 3 ) :
                {
                    mFunInvertJ = & inv3 ;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid number of dimensions") ;
                }
            }

            // done if this is a block
            if ( mGroup->type() == GroupType::BLOCK )
            {
                switch( mGroup->parent()->iwg()->model_dimensionality() )
                {
                    case( ModelDimensionality::TwoD ) :
                    {
                        mFundV = & Calculator::dV_2D_3D ;
                        break ;
                    }
                    case( ModelDimensionality::AxSymmX ) :
                    {
                        mFundV = & Calculator::dV_axsymmx ;
                        break ;
                    }
                    case( ModelDimensionality::AxSymmY ) :
                    {
                        mFundV = & Calculator::dV_axsymmy ;
                        break ;
                    }
                    case( ModelDimensionality::ThreeD ) :
                    {
                        mFundV = & Calculator::dV_2D_3D ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid Model Dimensionality");
                    }
                }
                return ;
            }

            switch( mGroup->parent()->iwg()->model_dimensionality() )
            {
                case( ModelDimensionality::TwoD ) :
                {
                    mFundS = & Calculator::dS_line ;
                    break ;
                }
                case( ModelDimensionality::AxSymmX ) :
                {
                    mFundS = & Calculator::dS_axsymmx ;
                    break ;
                }
                case( ModelDimensionality::AxSymmY ) :
                {
                    mFundS = & Calculator::dS_axsymmy ;
                    break ;
                }
                case( ModelDimensionality::ThreeD ) :
                {
                    BELFEM_ERROR( false, "No dS function assigned");
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Model Dimensionality");
                }
            }

            // check if we are allocating master and slave elements
            if( mGroup->master_type() != ElementType::EMPTY )
            {
                uint tNumNodes = mesh::number_of_nodes( mGroup->master_type() );

                mXm.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

                // allocate normal vector
                mNormal.set_size( tNumDimensions, BELFEM_QUIET_NAN );

                mJm = this->create_matrix( "Jm", tNumDimensions, tNumDimensions );

                // get the interpolation order of the master block
                InterpolationOrder tOrder = mesh::interpolation_order( mGroup->master_type() ) ;

                switch( mesh::geometry_type( mGroup->master_type() ) )
                {
                    case( GeometryType::TRI ) :
                    {
                        if( tOrder == InterpolationOrder::LINEAR )
                        {
                            mFunNormal = & Calculator::normal_tri_straight ;
                        }
                        else
                        {
                            mFunNormal = & Calculator::normal_tri_curved ;
                        }
                        break ;
                    }
                    case( GeometryType::QUAD ) :
                    {
                        if( tOrder == InterpolationOrder::LINEAR )
                        {
                            mFunNormal = & Calculator::normal_quad_straight ;
                        }
                        else
                        {
                            mFunNormal = & Calculator::normal_quad_curved ;
                        }
                        break ;
                    }
                    case( GeometryType::TET ) :
                    {
                        if( tOrder == InterpolationOrder::LINEAR )
                        {
                            mFunNormal = & Calculator::normal_tet_straight ;
                        }
                        else
                        {
                            mFunNormal = & Calculator::normal_tet_curved ;
                        }
                        break ;
                    }
                    case( GeometryType::HEX ) :
                    {
                        mFunNormal = & Calculator::normal_hex ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "No normal function assigned");
                    }
                }
            }

            if( mGroup->slave_type() != ElementType::EMPTY )
            {
                uint tNumNodes = mesh::number_of_nodes( mGroup->slave_type() );

                mXs.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

                //tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->slave_type() );

                // mIndexXs = this->create_matrix( "Xs", tNumNodes, tNumDimensions );
                // mIndexBs = this->create_matrix( "Bs", tNumNodes, tNumDimensions );

                mJs = this->create_matrix( "Js", tNumDimensions, tNumDimensions );

                switch( mesh::geometry_type( mGroup->slave_type() ) )
                {
                    case( GeometryType::TRI ) :
                    case( GeometryType::QUAD ) :
                    {
                        mFunSlaveIntegration = & Calculator::slave_integration_2d  ;
                        break ;
                    }
                    case( GeometryType::TET ) :
                    {
                        mFunSlaveIntegration = & Calculator::slave_integration_tet  ;
                        break ;
                    }
                    case( GeometryType::HEX ) :
                    {
                        mFunSlaveIntegration = & Calculator::slave_integration_hex  ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "unsupported type for slave side of facets");
                    }
                }
            }
            else
            {
                mFunSlaveIntegration = nullptr ;
            }

        }

//--------------------------------------- J inv function---------------------------------------

        void
        Calculator::link( Element * aElement )
        {
            // reset vectors
            for( calculator::VectorData * tVector : mVectors )
            {
                tVector->set_index( BELFEM_UINT_MAX );
            }

            // reset matrices
            for( calculator::MatrixData * tMatrix : mMatrices )
            {
                tMatrix->set_index( BELFEM_UINT_MAX );
            }

            // reset determinant
            mDetJ      = BELFEM_QUIET_NAN ;
            mDetJIndex = BELFEM_UINT_MAX ;
            mSurfaceIncrement = BELFEM_QUIET_NAN ;

            mNormalIndex = BELFEM_UINT_MAX ;

            aElement->get_node_coors( mX );


            mIsCurved = aElement->element()->is_curved() ;

            // for curved elements, we can copy the corner nodes
            if( ! mIsCurved )
            {
                for( uint j=0; j<mX.n_cols(); ++j )
                {
                    for( uint i=0; i<mNumberOfCornerNodes; ++i )
                    {
                        mXc( i, j ) = mX( i, j );
                    }
                }
            }


            mJ->set_index( BELFEM_UINT_MAX );
            if( aElement->master() != nullptr )
            {
                aElement->master()->get_node_coors( mXm );
                mMasterIndex = aElement->facet()->master_index() ;
                mMasterIntegration = mGroup->master_integration( mMasterIndex );
                mJm->set_index( BELFEM_UINT_MAX );
            }

            if( aElement->slave() != nullptr )
            {
                aElement->slave()->get_node_coors( mXs );
                mSlaveIntegration = ( this->*mFunSlaveIntegration )( aElement );
                mJs->set_index( BELFEM_UINT_MAX );
            }

            mElement = aElement ;
        }

//------------------------------------------------------------------------------

        calculator::VectorData *
        Calculator::create_vector( const string & aLabel, const uint aSize, const EntityType aType )
        {
            BELFEM_ASSERT( ! mVectorMap.key_exists( aLabel ),
                           "Vector %s already exists in calculator object.",
                           aLabel.c_str() );

            // create a new vector
            calculator::VectorData * aVector = new calculator::VectorData( aLabel, aSize, aType );

            // add vector to user vectors
            mVectors.push( aVector );

            // add vector to map
            mVectorMap[ aLabel ] = aVector ;

            return aVector ;
        }

//------------------------------------------------------------------------------

        calculator::MatrixData *
        Calculator::create_matrix( const string & aLabel,
                             const uint aNumRows,
                             const uint aNumCols )
        {
            BELFEM_ASSERT( ! mMatrixMap.key_exists( aLabel ),
                           "Matrix %s already exists in calculator object.",
                           aLabel.c_str() );

            // create a new vector
            calculator::MatrixData * aMatrix = new calculator::MatrixData( aLabel,
                                                               aNumRows,
                                                               aNumCols );

            // add vector to user vectors
            mMatrices.push( aMatrix );

            // add vector to map
            mMatrixMap[ aLabel ] = aMatrix ;

            return aMatrix ;
        }


//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::node_data_0( const string & aNodeField )
        {
            string tNodeField = aNodeField + "0";

            BELFEM_ASSERT(
                    mMesh->field( tNodeField )->entity_type() == EntityType::NODE,
                    "Field '%s' is not a node field", tNodeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( tNodeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // grab data object from mesh
            const Vector< real > & tNodeData = mMesh->field_data( tNodeField );

            // loop over all nodes
            for( uint i=0; i< mElement->element()->number_of_nodes(); ++i )
            {
                // copy node data from mesh
                aData( i  ) = tNodeData( mElement->element()->node( i )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::node_data_1( const string & aNodeField )
        {
            BELFEM_ASSERT(
                    mMesh->field( aNodeField )->entity_type() == EntityType::NODE,
                    "Field '%s' is not a node field", aNodeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( aNodeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // grab data object from mesh
            const Vector< real > & tNodeData = mMesh->field_data( aNodeField );

            // loop over all nodes
            for( uint i=0; i< mElement->element()->number_of_nodes(); ++i )
            {
                // copy node data from mesh
                aData( i  ) = tNodeData( mElement->element()->node( i )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::node_data_theta( const string & aNodeField )
        {
            string tNodeField = aNodeField + "0";

            BELFEM_ASSERT(
                    mMesh->field( tNodeField )->entity_type() == EntityType::NODE,
                    "Field '%s' is not a node field", tNodeField.c_str() );

            BELFEM_ASSERT(
                    mMesh->field( aNodeField )->entity_type() == EntityType::NODE,
                    "Field '%s' is not a node field", aNodeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( aNodeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // grab data object from mesh
            const Vector< real > & tNodeData0 = mMesh->field_data( tNodeField );
            const Vector< real > & tNodeData1 = mMesh->field_data( aNodeField );

            // loop over all nodes
            for( uint i=0; i< mElement->element()->number_of_nodes(); ++i )
            {
                // get the node index
                index_t tIndex = mElement->element()->node( i )->index() ;

                // interpolate data at desired timestep
                aData( i  ) =
                          tNodeData0( tIndex ) * mOneMinusTheta
                        + tNodeData1( tIndex ) * mTheta ;
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_linear_0( const string & aEdgeField )
        {
            string tEdgeField = aEdgeField + "0";

            BELFEM_ASSERT(
                    mMesh->field( tEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", tEdgeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( tEdgeField );

            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( tEdgeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // loop over all edges
            for( uint e=0; e< mElement->element()->number_of_edges(); ++e )
            {
                aData( e ) = tField( mElement->element()->edge( e )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_linear_1( const string & aEdgeField )
        {
            BELFEM_ASSERT(
                    mMesh->field( aEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( aEdgeField );

            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( aEdgeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // loop over all edges
            for( uint e=0; e< mElement->element()->number_of_edges(); ++e )
            {
                aData( e ) = tField( mElement->element()->edge( e )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_linear_theta( const string & aEdgeField )
        {
            string tEdgeField = aEdgeField + "0";

            BELFEM_ASSERT(
                    mMesh->field( tEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", tEdgeField.c_str() );

            BELFEM_ASSERT(
                    mMesh->field( aEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( aEdgeField );

            // get ref to field on mesh
            Vector< real > & tField0 = mMesh->field_data( tEdgeField );
            Vector< real > & tField1 = mMesh->field_data( aEdgeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // loop over all edges
            for( uint e=0; e< mElement->element()->number_of_edges(); ++e )
            {
                index_t tIndex = mElement->element()->edge( e )->index() ;
                aData( e ) =  tField0( tIndex ) * mOneMinusTheta
                                   + tField1( tIndex ) * mTheta ;
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_quadratic_2d_0(
                const string & aEdgeField,
                const string & aFaceField,
                const string & aVectorLabel )
        {
            string tEdgeField = aEdgeField + "0";
            string tFaceField = aFaceField + "0";

            BELFEM_ASSERT(
                    mMesh->field( tEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", tEdgeField.c_str());


            BELFEM_ASSERT(
                    mMesh->field( tFaceField )->entity_type() == EntityType::FACE,
                    "Field '%s' is not an edge field", tFaceField.c_str());


            Vector< real > & aData = mMesh->field_data( aVectorLabel );

            Vector< real > & tEdgeData = mMesh->field_data( tEdgeField );
            Vector< real > & tFaceData = mMesh->field_data( tFaceField );

            uint tCount = 0;

            for ( uint e = 0; e < mElement->element()->number_of_edges(); ++e )
            {

                // get index of edge
                index_t tIndex = mElement->element()->edge( e )->index();

                // check direction of edge
                if ( mElement->edge_direction( e ))
                {
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex );
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex + 1 );
                }
                else
                {
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex + 1 );
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex );
                }
            }

            // write data into container
            index_t tIndex = 2 * mElement->element()->index();
            aData( tCount++ ) = tFaceData( tIndex );
            aData( tCount++ ) = tFaceData( tIndex + 1 );

            return aData ;
        }

//-------------------------------------   -----------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_quadratic_2d_1(
                const string & aEdgeField,
                const string & aFaceField,
                const string & aVectorLabel )
        {

            BELFEM_ASSERT(
                    mMesh->field( aEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeField.c_str());


            BELFEM_ASSERT(
                    mMesh->field( aFaceField )->entity_type() == EntityType::FACE,
                    "Field '%s' is not an edge field", aFaceField.c_str());


            Vector< real > & aData = mMesh->field_data( aVectorLabel );

            Vector< real > & tEdgeData = mMesh->field_data( aEdgeField );
            Vector< real > & tFaceData = mMesh->field_data( aFaceField );

            uint tCount = 0;

            for ( uint e = 0; e < mElement->element()->number_of_edges(); ++e )
            {

                // get index of edge
                index_t tIndex = mElement->element()->edge( e )->index();

                // check direction of edge
                if ( mElement->edge_direction( e ))
                {
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex );
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex + 1 );
                }
                else
                {
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex + 1 );
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex );
                }
            }

            // write data into container
            index_t tIndex = 2 * mElement->element()->index();
            aData( tCount++ ) = tFaceData( tIndex );
            aData( tCount++ ) = tFaceData( tIndex + 1 );

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_quadratic_2d_theta(
                const string & aEdgeField,
                const string & aFaceField,
                const string & aVectorLabel )
        {
            string tEdgeField = aEdgeField + "0";
            string tFaceField = aFaceField + "0";

            BELFEM_ASSERT(
                    mMesh->field( tEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", tEdgeField.c_str());


            BELFEM_ASSERT(
                    mMesh->field( tFaceField )->entity_type() == EntityType::FACE,
                    "Field '%s' is not an edge field", tFaceField.c_str());


            BELFEM_ASSERT(
                    mMesh->field( aEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeField.c_str());


            BELFEM_ASSERT(
                    mMesh->field( aFaceField )->entity_type() == EntityType::FACE,
                    "Field '%s' is not an edge field", aFaceField.c_str());

            Vector< real > & aData = mMesh->field_data( aVectorLabel );

            Vector< real > & tEdgeData0 = mMesh->field_data( tEdgeField );
            Vector< real > & tFaceData0 = mMesh->field_data( tFaceField );

            Vector< real > & tEdgeData1 = mMesh->field_data( aEdgeField );
            Vector< real > & tFaceData1 = mMesh->field_data( aFaceField );

            uint tCount = 0;

            for ( uint e = 0; e < mElement->element()->number_of_edges(); ++e )
            {

                // get index of edge
                index_t tIndex = mElement->element()->edge( e )->index();

                // check direction of edge
                if ( mElement->edge_direction( e ))
                {
                    aData( tCount++ ) =   tEdgeData0( tIndex + tIndex ) * mOneMinusTheta
                                              + tEdgeData1( tIndex + tIndex )  * mTheta ;

                    aData( tCount++ ) = tEdgeData0( tIndex + tIndex + 1 ) * mOneMinusTheta
                                              + tEdgeData1( tIndex + tIndex + 1 ) * mTheta ;
                }
                else
                {
                    aData( tCount++ ) = tEdgeData0( tIndex + tIndex + 1 ) * mOneMinusTheta
                                        + tEdgeData1( tIndex + tIndex + 1 ) * mTheta ;

                    aData( tCount++ ) =   tEdgeData0( tIndex + tIndex ) * mOneMinusTheta
                                          + tEdgeData1( tIndex + tIndex )  * mTheta ;
                }
            }

            // write data into container
            index_t tIndex = 2 * mElement->element()->index();
            aData( tCount++ ) =   tFaceData0( tIndex ) * mOneMinusTheta
                                       + tFaceData1( tIndex ) * mTheta ;

            aData( tCount++ ) = tFaceData0( tIndex + 1 )
                                     + tFaceData1( tIndex + 1 ) * mTheta ;

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tri_straight( const uint aIndex  )
        {
            if( mNormalIndex == BELFEM_UINT_MAX )
            {
                // remember the index
                mNormalIndex = aIndex;

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) = mXm( 1, 1 ) - mXm( 0, 1 );
                        mNormal( 1 ) = mXm( 0, 0 ) - mXm( 1, 0 );
                        break;
                    }
                    case ( 1 ) :
                    {

                        mNormal( 0 ) = mXm( 2, 1 ) - mXm( 1, 1 );
                        mNormal( 1 ) = mXm( 1, 0 ) - mXm( 2, 0 );
                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = mXm( 0, 1 ) - mXm( 2, 1 );
                        mNormal( 1 ) = mXm( 2, 0 ) - mXm( 0, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this value will contain the length of the side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // finally, we must adapt this value, since along the edge
                // we integrate from -1 to 1 rather than 0 to 1
                mSurfaceIncrement *= 0.5;
            }

            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tri_curved( const uint aIndex  )
        {
            if( ! mIsCurved )
            {
                return this->normal_tri_straight( aIndex );
            }
            else if( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix
                const Matrix< real > & tJ = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) =
                                tJ( 1, 1 )
                                - tJ( 0, 1 );

                        mNormal( 1 ) =
                                tJ( 0, 0 )
                                - tJ( 1, 0 );

                        break;
                    }
                    case ( 1 ) :
                    {
                        mNormal( 0 ) = -tJ( 1, 1 );
                        mNormal( 1 ) = tJ( 1, 0 );

                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = tJ( 0, 1 );
                        mNormal( 1 ) = -tJ( 0, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // if this edge was straight, this would be the length of this side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // finally, we must adapt this value, since along the edge
                // we integrate from -1 to 1 rather than 0 to 1
                mSurfaceIncrement *= 0.5;
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_quad_straight( const uint aIndex )
        {
            if( mNormalIndex == BELFEM_UINT_MAX )
            {
                // remember the index
                mNormalIndex = aIndex;

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) = mXm( 1, 1 ) - mXm( 0, 1 );
                        mNormal( 1 ) = mXm( 0, 0 ) - mXm( 1, 0 );
                        break;
                    }
                    case ( 1 ) :
                    {
                        mNormal( 0 ) = mXm( 2, 1 ) - mXm( 1, 1 );
                        mNormal( 1 ) = mXm( 1, 0 ) - mXm( 2, 0 );
                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = mXm( 3, 1 ) - mXm( 2, 1 );
                        mNormal( 1 ) = mXm( 2, 0 ) - mXm( 3, 0 );
                        break;
                    }
                    case ( 3 ) :
                    {
                        mNormal( 0 ) = mXm( 0, 1 ) - mXm( 3, 1 );
                        mNormal( 1 ) = mXm( 3, 0 ) - mXm( 0, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this value will contain the length of the side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // finally, we must adapt this value, since along the edge
                // we integrate from -1 to 1 rather than 0 to 1
                mSurfaceIncrement *= 0.5;
            }

            // now we can return the vector
            return mNormal;
        }
//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_quad_curved( const uint aIndex )
        {
            if( ! mIsCurved )
            {
                return this->normal_quad_straight( aIndex );
            }
            else if( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix
                const Matrix< real > & tJ = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) =  tJ( 0, 1 );
                        mNormal( 1 ) = -tJ( 0, 0 );
                        break;
                    }
                    case ( 1 ) :
                    {
                        mNormal( 0 ) =  tJ( 1, 1 );
                        mNormal( 1 ) = -tJ( 1, 0 );
                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = -tJ( 0, 1 );
                        mNormal( 1 ) =  tJ( 0, 0 );
                        break;
                    }
                    case ( 3 ) :
                    {
                        mNormal( 0 ) = -tJ( 1, 1 );
                        mNormal( 1 ) =  tJ( 1, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this value will contain the length of the side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // no multiplication of mSurfaceIncrement with 0.5 here!
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tet_straight( const uint aIndex )
        {
            if( mNormalIndex == BELFEM_UINT_MAX )
            {
                mNormal( 0 ) =
                          ( mX( 0, 1 ) - mX( 1, 1 ) )
                        * ( mX( 0, 2 ) - mX( 2, 2 ) )
                        - ( mX( 0, 1 ) - mX( 2, 1 ) )
                        * ( mX( 0, 2 ) - mX( 1, 2 ) );


                mNormal( 1 ) =
                          ( mX( 0, 0 ) - mX( 2, 0 ) )
                        * ( mX( 0, 2 ) - mX( 1, 2 ) )
                        - ( mX( 0, 0 ) - mX( 1, 0 ) )
                        * ( mX( 0, 2 ) - mX( 2, 2 ) );


                mNormal( 2 ) = ( mX( 0, 0 ) - mX( 1, 0 ) )
                          * ( mX( 0, 1 ) - mX( 2, 1 ) )
                          - ( mX( 0, 0 ) - mX( 2, 0 ) )
                          * ( mX( 0, 1 ) - mX( 1, 1 ) ) ;

                // This is a bit confusing: the sum of the weights for a triangle is 0.5,
                // so this should be multiplied by 2. On the other hand, however,
                // however, this is twice the surface of the triangle since we span a parallelogram.
                // Eventually, we multiply by 0.5 * 2 = 1
                mSurfaceIncrement = norm( mNormal );
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tet_curved( const uint aIndex )
        {
            if( ! mIsCurved )
            {
                return this->normal_tet_straight( aIndex );
            }
            else if ( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix
                const Matrix< real > & tJ = this->Jm( aIndex );

                // normal
                mNormal( 0 ) =   tJ( 0, 1 ) * tJ( 1, 2 )
                            - tJ( 0, 2 ) * tJ( 1, 1 );

                mNormal( 1 ) =   tJ( 0, 2 ) * tJ( 1, 0 )
                            - tJ( 0, 0 ) * tJ( 1, 2 );

                mNormal( 2 ) =   tJ( 0, 0 ) * tJ( 1, 1 )
                            - tJ( 0, 1 ) * tJ( 1, 0 );


                mSurfaceIncrement = norm( mNormal );
            }


            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_hex( const uint aIndex )
        {
            if ( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix for the master element
                const Matrix< real > & tJ = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case( 0 ) :
                    {
                        mNormal( 0 ) = tJ( 0, 1 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 1 );
                        mNormal( 1 ) = tJ( 0, 2 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 2 );
                        mNormal( 2 ) = tJ( 0, 0 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 0 );
                        break ;
                    }
                    case( 1 ) :
                    {
                        mNormal( 0 ) = tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 );
                        mNormal( 1 ) = tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 );
                        mNormal( 2 ) = tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 );
                        break ;
                    }
                    case( 2 ) :
                    {
                        mNormal( 0 ) = tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 );
                        mNormal( 1 ) = tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 );
                        mNormal( 2 ) = tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 );
                        break ;
                    }
                    case( 3 ) :
                    {
                        mNormal( 0 ) = tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 );
                        mNormal( 1 ) = tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 );
                        mNormal( 2 ) = tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 );
                        break ;
                    }
                    case( 4 ) :
                    {
                        mNormal( 0 ) = tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 );
                        mNormal( 1 ) = tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 );
                        mNormal( 2 ) = tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 );
                        break ;
                    }
                    case( 5 ) :
                    {
                        mNormal( 0 ) = tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 );
                        mNormal( 1 ) = tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 );
                        mNormal( 2 ) = tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // note that the sum of the weights along the surface is 4, so we need to
                // multiply by 0.25
                mSurfaceIncrement = norm( mNormal ) * 0.25 ;
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::q0()
        {
            for( uint k=0; k<mElement->number_of_dofs(); ++k )
            {
                Dof * tDof = mElement->dof( k );

                mq0( k ) = mMesh->field( mMesh->field( tDof->field_index() )->label() + "0" )->data()( tDof->dof_index_on_field() );
            }
            return mq0 ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::q1()
        {
            for( uint k=0; k<mElement->number_of_dofs(); ++k )
            {
                Dof * tDof = mElement->dof( k );

                mq1( k ) = mMesh->field( tDof->field_index() )->data()( tDof->dof_index_on_field() );
            }
            return mq1 ;
        }

//------------------------------------------------------------------------------
    }
}
