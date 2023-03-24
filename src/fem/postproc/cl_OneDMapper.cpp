//
// Created by Christian Messe on 30.12.20.
//

#include "cl_OneDMapper.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "fn_trans.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    OneDMapper::OneDMapper(
            const Vector< real > & aTargetCoords,
            const uint aOrder,
            const proc_t aMasterRank ) :
            mNumTargetNodes( aTargetCoords.length() ),
            mOrder( aOrder ),
            mMasterRank( aMasterRank ),
            mMyRank( comm_rank())
    {
        if ( mMyRank == mMasterRank )
        {
            // create the mesh
            this->create_target_nodes( aTargetCoords );
            this->create_elements();
            this->compute_adjency();

            // create the system matrix
            this->create_shape_function();
            this->create_mass_matrix();

            // creat the solver
            mSolver = new Solver( SolverType::UMFPACK, mMasterRank );

            // allocate the RHS
            mRHS.set_size( mNodes.size() );
        }
    }

//------------------------------------------------------------------------------

    OneDMapper::~OneDMapper()
    {
        if ( mMyRank == mMasterRank )
        {
            delete mSolver ;
            delete mMassMatrix;
            delete mShape ;

            for( mesh::Element * tElement : mElements )
            {
                delete tElement ;
            }

            for( mesh::Node * tNode : mNodes )
            {
                delete tNode ;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::project(
            const Vector< real > & aSourceNodes,
            const Vector< real > & aSourceValues,
                  Vector< real > & aTargetValues )
    {
        if( mMyRank == mMasterRank )
        {
            // compute the RHS
            this->compute_rhs( aSourceNodes, aSourceValues, mRHS );

            // allocate size of target vector
            aTargetValues.set_size( mNodes.size() );

            // solve the system
            mSolver->solve( *mMassMatrix, aTargetValues, mRHS );
        }
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::derive( const Vector< real > & aValues,
                              Vector< real > & aDerivatives )
    {
        if( mMyRank == mMasterRank )
        {
            // compute the RHS
            this->compute_derivative_rhs( aValues, mRHS );

            // allocate size of target vector
            aDerivatives.set_size( mNodes.size() );

            // solve the system
            mSolver->solve( *mMassMatrix, aDerivatives, mRHS );
        }
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::create_target_nodes( const Vector< real > & aTargetCoords )
    {
        // allocate node container
        mNodes.set_size( mNumTargetNodes, nullptr );

        for ( index_t k = 0; k < mNumTargetNodes; ++k )
        {
            mNodes( k ) = new mesh::Node( k + 1, aTargetCoords( k ));
            mNodes( k )->set_index( k );
        }

    }

//------------------------------------------------------------------------------

    void
    OneDMapper::create_elements()
    {

        index_t tNumElements = 0;

        switch ( mOrder )
        {
            case ( 1 ) :
            {
                tNumElements = mNumTargetNodes - 1;
                mElementType = ElementType::LINE2;
                break;
            }
            case ( 2 ) :
            {
                BELFEM_ERROR( mNumTargetNodes % 2 == 1,
                             "Need an odd number for second order interpolation" );
                tNumElements = ( mNumTargetNodes - 1 ) / 2;
                mElementType = ElementType::LINE3;
                break;
            }
            case ( 3 ) :
            {
                BELFEM_ERROR( mNumTargetNodes % 3 == 1,
                             "Wrong number of target nodes" );
                tNumElements = ( mNumTargetNodes - 1 ) / 3;
                mElementType = ElementType::LINE4;
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid order : %u", ( unsigned int ) mOrder );
                break;
            }
        }

        // create an element factory
        mesh::ElementFactory tFactory ;

        // allocate memory
        mElements.set_size( tNumElements, nullptr );

        // connect nodes to elements
        index_t off = 0;

        // create the elements
        for ( index_t e = 0; e < tNumElements; ++e )
        {
            mesh::Element * tElement = tFactory.create_element( mElementType, e + 1 );

            tElement->set_index( e );

            // link element to nodes
            tElement->insert_node( mNodes( off ), 0 );
            tElement->insert_node( mNodes( off + mOrder ), 1 );

            for ( index_t i = 1; i < mOrder; ++i )
            {
                tElement->insert_node( mNodes( off + i ), i + 1 );
            }

            // increment offset
            off += mOrder;

            // add element to container
            mElements( e ) = tElement;
        }
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::compute_adjency()
    {
        uint tNumNodesPerElement = mOrder + 1 ;

        // reset node containers
        for( mesh::Node * tNode : mNodes )
        {
            tNode->reset_element_container() ;
        }

        // count elements per node
        for( mesh::Element * tElement : mElements )
        {
            for( uint k=0; k<tNumNodesPerElement; ++k )
            {
                tElement->node( k )->increment_element_counter() ;
            }
        }

        // allocate node containers
        for( mesh::Node * tNode : mNodes )
        {
            tNode->allocate_element_container();
        }

        // add elements to node
        for( mesh::Element * tElement : mElements )
        {
            for( uint k=0; k<tNumNodesPerElement; ++k )
            {
                tElement->node( k )->add_element( tElement );
            }
        }

        Cell< mesh::Node * > tNodes ;

        // create node adjency
        for( mesh::Node * tNode : mNodes )
        {
            // allocate cell
            tNodes.set_size( tNode->number_of_elements() * tNumNodesPerElement, nullptr );
            uint tCount = 0 ;

            // collect nodes
            for( uint e=0; e<tNode->number_of_elements(); ++e )
            {
                mesh::Element * tElement = tNode->element( e );
                for( uint k=0; k<tNumNodesPerElement; ++k )
                {
                    tNodes( tCount++ ) = tElement->node( k );
                }
            }

            // make nodes unique
            unique( tNodes );

            // add nodes to nodes
            tCount = tNodes.size() ;

            tNode->init_vertex_container( tCount );

            for( uint k=0; k<tCount; ++k )
            {
                tNode->insert_vertex( tNodes( k ) );
            }
        }

    }

//------------------------------------------------------------------------------

    void
    OneDMapper::create_shape_function()
    {
        // the factory
        fem::InterpolationFunctionFactory tFactory ;

        // create the shape pointer
        mShape = tFactory.create_lagrange_function( mElementType );

        // integration points for this element
        intpoints( IntegrationScheme::GAUSS,
                   GeometryType::LINE,
                   2 * mOrder + 1 ,
                   mWeights,
                   mPoints );

        uint tNumPoints = mWeights.length() ;

        // initialize memory
        mN.set_size( tNumPoints, {} );
        mdNdXi.set_size( tNumPoints, {} );

        // evaluate the shape function
        for( uint k=0; k<tNumPoints; ++k )
        {
            // get ref to entry
            Matrix< real > & tN = mN( k );

            // get ref to entry
            Matrix< real > & tdNdxi = mdNdXi( k );

            // allocate vector
            tN.set_size( 1, mOrder + 1 );
            tdNdxi.set_size( 1, mOrder + 1 );

            // compute data
            mShape->N( mPoints.col( k ), tN );
            mShape->dNdXi( mPoints.col( k ), tdNdxi );
        }
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::create_mass_matrix()
    {
        mMassMatrix = new SpMatrix(
                reinterpret_cast< Cell< graph::Vertex * > & > ( mNodes ) );

        // element matrix for element of length 2
        Matrix< real > tNN( mOrder+1, mOrder+1, 0.0 );

        // individual element matrix
        uint tNumPoints = mWeights.length() ;
        for( uint k=0; k<tNumPoints; ++k )
        {
            tNN += mWeights( k ) * trans( mN( k ) ) * mN( k ) ;
        }

        // scale matrix to represent an element of length 1
        tNN *= 0.5 ;

        // assemble full matrix
        mMassMatrix->fill( 0.0 );

        // convert matrix to ref
        SpMatrix & tM = *mMassMatrix ;

        // get number of nodes per element
        uint tNumNodes = mOrder + 1 ;

        // loop over all elements
        for( mesh::Element * tElement : mElements )
        {
            // compute the element length
            real tLength = tElement->node( 1 )->x() - tElement->node( 0 )->x() ;

            // assemble matrix
            for( uint j=0; j<tNumNodes; ++j )
            {
                for( uint i=0; i<tNumNodes; ++i )
                {
                    tM( tElement->node( i )->index(),
                                 tElement->node( j )->index() ) +=
                                         tNN( i, j ) * tLength ;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    index_t
    OneDMapper::find_interval( const Vector< real > & aSourceNodes, const real aX )
    {
        index_t i = 0 ;
        index_t k = aSourceNodes.length() - 1 ;
        index_t j = ( i + k +1 ) / 2 ;
        index_t n = k ;

        for( index_t e=0; e<n; ++e )
        {
            if( aSourceNodes( i ) <= aX  && aSourceNodes( j ) >= aX )
            {
                k = j ;
            }
            else
            {
                i = j ;
            }
            if( j-i == 1 )
            {
                BELFEM_ASSERT( aSourceNodes( i ) <= aX  && aSourceNodes( j ) >= aX,
                    "Something went wrong : X( %u ) = %12.3g, X( %u ) = %12.3g, X = %12.3g",
                              ( unsigned int ) i,
                              ( double ) aSourceNodes( i ),
                              ( unsigned  int ) j,
                              ( double ) aSourceNodes( j ),
                              ( double ) aX );
                return i ;
            }
            else
            {
                j = ( i + k +1 ) / 2 ;
            }
        }

        return BELFEM_UINT_MAX ;
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::compute_rhs(
            const Vector< real > & aSourceNodes,
            const Vector< real > & aSourceValues,
                  Vector< real > & aRHS )
    {
        aRHS.fill( 0.0 );

        Vector< real > tRHS( mOrder + 1 );

        uint tNumPoints = mWeights.length() ;

        // loop over all elements
        for( mesh::Element * tElement : mElements )
        {
            // fill element contribution
            tRHS.fill( 0.0 );

            // loop over all integration points
            for( uint k=0; k<tNumPoints; ++k )
            {
                // get shape function
                Matrix< real > & tN = mN( k );

                // compute x-value for this integration point
                real tX = ( tElement->node( 0 )->x() * ( 1. - mPoints( 0, k ) )
                          + tElement->node( 1 )->x() * ( 1. + mPoints( 0, k ) ) ) * 0.5 ;

                // find interval in source
                index_t j = this->find_interval( aSourceNodes, tX );

                real tXi = ( tX - aSourceNodes( j ) ) / ( aSourceNodes( j+1 ) - aSourceNodes( j ) );

                // compute target value
                real tY = aSourceValues( j ) + tXi * ( aSourceValues( j + 1 ) - aSourceValues( j ) );

                // compute contribution to RHS
                for( uint i=0; i<=mOrder; ++i )
                {
                    tRHS( i ) += mWeights( k ) * tN( 0, i ) * tY ;
                }
            }

            // scale RHS
            tRHS *= 0.5 * ( tElement->node( 1 )->x() - tElement->node( 0 )->x() );

            // add element rhs to system rhs
            for( uint k=0; k<= mOrder; ++k )
            {
                aRHS( tElement->node( k )->index() ) += tRHS( k );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    OneDMapper::compute_derivative_rhs( const Vector< real > & aValues,
                                              Vector< real > & aRHS )
    {
        aRHS.fill( 0.0 );

        Vector< real > tRHS( mOrder + 1 );

        // element values
        Vector< real > tValues( mOrder + 1 );

        uint tNumPoints = mWeights.length() ;

        // compute projection matrix
        Matrix< real > tP( mOrder + 1, mOrder + 1, 0.0 ) ;

        for( uint k=0; k<tNumPoints; ++k )
        {
            // get shape functions
            Matrix< real > & tN     = mN( k );
            Matrix< real > & tdNdXi = mdNdXi( k );

            tP += mWeights( k ) * trans( tN ) * tdNdXi ;
        }

        // note: P is not scaled, sind the expression dx/dxi = l/2 cancels itself

        // loop over all elements
        for( mesh::Element * tElement : mElements )
        {
            // collect element values
            for( uint k=0; k<=mOrder; ++k )
            {
                tValues( k ) = aValues( tElement->node( k )->index() );
            }

            // compute RHS
            tRHS = tP * tValues ;

            // add rhs to right hand side
            for( uint k=0; k<=mOrder; ++k )
            {
                aRHS( tElement->node( k )->index() ) += tRHS( k );
            }
        }
    }

//------------------------------------------------------------------------------
}