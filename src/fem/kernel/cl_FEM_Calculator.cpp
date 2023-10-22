//
// Created by christian on 6/9/23.
//
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

        Calculator::Calculator( Group * aGroup ) :
                mGroup( aGroup ),
                mMesh( aGroup->parent()->mesh() )
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
        }

//------------------------------------------------------------------------------


        void
        Calculator::initialize_integration(
                const ElementType       aElementType,
                const InterpolationType aInterpolationType )
        {
            if( mIntegration != nullptr )
            {
                delete mIntegration ;
            }
            mIntegration = new IntegrationData( aElementType,
                                                aInterpolationType );
        }

//------------------------------------------------------------------------------

        void
        Calculator::set_integration_order( const uint aOrder )
        {
            mNumberOfNodes = mesh::number_of_nodes( mGroup->element_type() );

            mIntegration->populate( aOrder, mGroup->parent()->integration_scheme() );

            this->allocate_memory() ;
        }

//------------------------------------------------------------------------------

        void
        Calculator::initialize_master_integration( const ElementType aElementType,
                                       const InterpolationType aInterpolationType )
        {

        }

//------------------------------------------------------------------------------

        void
        initialize_slave_integration( const ElementType aElementType,
                                      const InterpolationType aInterpolationType )
        {

        }

//------------------------------------------------------------------------------

        void
        Calculator::allocate_memory()
        {

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

                for ( const string & tLabel: tEquation->all_fields())
                {
                    // the size of the vector
                    uint tSize = 0;

                    EntityType tType = entity_type( tLabel );

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
            }
        }

//------------------------------------------------------------------------------


        const IntegrationData *
        Calculator::slave_integration_2d( const Element * aElement )
        {
            return mGroup->slave_integration( aElement->facet()->slave_index() );
        }

        const IntegrationData *
        Calculator::slave_integration_tet( const Element * aElement )
        {
            return mGroup->slave_integration(
                    aElement->facet()->slave_index() * 3
                    + aElement->facet()->orientation_on_slave() );
        }

        const IntegrationData *
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

            if(  mGroup->number_of_elements() == 0 )
            {
                return;
            }

            // copy integration weights from group
            const Vector< real > & tW = mGroup->integration_weights();
            mIntegrationWeights.set_size( tW.length() );
            for ( uint k=0; k< tW.length(); ++k )
            {
                mIntegrationWeights( k ) = tW( k );
            }

            // make sure that allocation funciton is only called once
            BELFEM_ERROR( ! mIsAllocated, "calculator is already allocated" );
            mIsAllocated = true ;

            // number of nodes per element
            mNumberOfNodes = mesh::number_of_nodes( mGroup->element_type() );

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
                /*
                // check of we have vector fields

                // grab first element on block
                Element * tElement = mGroup->elements()( 0 ) ;

                // grab equation objectmB->matrix() ;
                IWG * tIWG = mGroup->parent()->iwg() ;

                for( uint k = 0 ; k<tElement->number_of_dofs() ; ++k  )
                {
                    std::cout << k << " " << tIWG->dof_label( tElement->dof( k )->type_id() ) << std::endl ;
                }*/

                return ;
            }

            // check if we are allocating master and slave elements
            if( mGroup->master_type() != ElementType::EMPTY )
            {
                uint tNumNodes = mesh::number_of_nodes( mGroup->master_type() );

                mXm.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

                // allocate normal vector
                mNormal.set_size( tNumDimensions, BELFEM_QUIET_NAN );

            }

            if( mGroup->slave_type() != ElementType::EMPTY )
            {
                uint tNumNodes = mesh::number_of_nodes( mGroup->slave_type() );

                mXs.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

                //tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->slave_type() );

                // mIndexXs = this->create_matrix( "Xs", tNumNodes, tNumDimensions );
                //mIndexJs = this->create_matrix( "Js", tNumDimensions, tNumDimensions );
                // mIndexBs = this->create_matrix( "Bs", tNumNodes, tNumDimensions );

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

            if( aElement->master() != nullptr )
            {
                aElement->master()->get_node_coors( mXm );
                mMasterIndex = aElement->facet()->master_index() ;
                mMasterIntegration = mGroup->master_integration( mMasterIndex );
            }

            if( aElement->slave() != nullptr )
            {
                aElement->slave()->get_node_coors( mXs );
                mSlaveIntegration = ( this->*mFunSlaveIntegration )( aElement );
            }
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
            if( mNormalIndex != aIndex )
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
            if( mNormalIndex != aIndex )
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
            if ( mNormalIndex != aIndex )
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
    }
}
