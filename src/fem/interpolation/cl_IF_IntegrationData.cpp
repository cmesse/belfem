//
// Created by Christian Messe on 18.01.22.
//

#include "cl_IF_IntegrationData.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "meshtools.hpp"
#include "fn_IF_initialize_integration_points_on_facet.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IntegrationData::IntegrationData( const ElementType aElementType,
                                          const InterpolationType aType,
                                          InterpolationFunction * aShapeFunction ) :
            mElementType( aElementType )
        {
            // create a temporary factory
            InterpolationFunctionFactory tFactory ;


            if( aShapeFunction == nullptr )
            {
                // create the shape function
                mShapeFunction = tFactory.create_function( aElementType, aType );

                // set owning flag
                mOwnShapeFunction = true ;
            }
            else
            {
                mShapeFunction = aShapeFunction ;
                mOwnShapeFunction = false ;
            }
        }

//------------------------------------------------------------------------------

        IntegrationData::~IntegrationData()
        {
            // delete the shape function if we own it
            if( mOwnShapeFunction )
            {
                delete mShapeFunction ;
            }
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::populate( const uint aIntegrationOrder, const IntegrationScheme aScheme )
        {
            // set integration order
            uint tIntegrationOrder = aIntegrationOrder == 0 ?
                                auto_integration_order( mElementType ) : aIntegrationOrder ;


            // compute the integration points
            intpoints( aScheme,
                       mesh::geometry_type( mElementType ),
                       tIntegrationOrder,
                       mWeights,
                       mPoints );

            this->evaluate_function();
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::populate_for_master(
                const uint aMasterIndex,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {
            initialize_integration_points_on_facet( mElementType,
                                                    aMasterIndex,
                                                    mWeights,
                                                    mPoints,
                                                    aIntegrationOrder,
                                                    aScheme );

            this->evaluate_function();
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::populate_for_slave(
                const uint aSlaveIndex,
                const uint aOrientation,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {
            switch ( mesh::geometry_type( mElementType ) )
            {
                case ( GeometryType::TRI ) :
                case ( GeometryType::QUAD ) :
                {
                    this->populate_for_slave_tri( aSlaveIndex,
                                                  aIntegrationOrder,
                                                  aScheme );
                    break ;
                }
                case ( GeometryType::TET ) :
                {
                    this->populate_for_slave_tet( aSlaveIndex,
                                                  aOrientation,
                                                  aIntegrationOrder,
                                                  aScheme );
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    this->populate_for_slave_hex( aSlaveIndex,
                                                   aOrientation,
                                                   aIntegrationOrder,
                                                   aScheme );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Geometry Type for slave points" );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::populate_for_slave_tri(
                const uint aFaceIndex,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {

            Vector< real > tWeights ;
            Matrix< real > tPoints ;

            initialize_integration_points_on_facet( mElementType,
                                                    aFaceIndex,
                                                    tWeights,
                                                    tPoints,
                                                    aIntegrationOrder,
                                                    aScheme );

            // for 2D elements, all we need to do is flip the directions of
            // on the edges
            int tNumPoints = tWeights.length() ;
            int tNumDim    = tPoints.n_rows() ;

            mWeights.set_size( tNumPoints );
            mPoints.set_size( tNumDim, tNumPoints );

            uint tCount = 0 ;
            for( int k=tNumPoints-1; k>=0; k-- )
            {
                mPoints.set_col( tCount, tPoints.col( k ) );
                mWeights( tCount++ ) = tWeights( k );
            }

            this->evaluate_function();
        }

//------------------------------------------------------------------------------

        /**
         * special function for popularization if this is a sideset
         */
        void
        IntegrationData::populate_for_slave_tet(
                const uint aSlaveIndex,
                const uint aOrientation,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {

            facetintpoints::intpoints_tet(
                    aSlaveIndex,
                    aOrientation,
                    mWeights,
                    mPoints,
                    aIntegrationOrder,
                    aScheme );

            this->evaluate_function();
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::populate_for_slave_hex(
                const uint aSlaveIndex,
                const uint aOrientation,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {

            facetintpoints::intpoints_hex(
                    aSlaveIndex,
                    aOrientation,
                    mWeights,
                    mPoints,
                    aIntegrationOrder,
                    aScheme );

            this->evaluate_function();
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::evaluate_function()
        {
            mNumberOfIntegrationPoints = mWeights.length() ;
            uint tNumberOfNodes  = mesh::number_of_nodes( mElementType );

            // allocate cells
            mN.set_size( mNumberOfIntegrationPoints,
                         Matrix< real >( 1, tNumberOfNodes ) );


            // evaluate shape function
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mShapeFunction->N( mPoints.col( k ), mN( k ) );
            }


            // create phi vector
            mPhi.set_size( mNumberOfIntegrationPoints,
                           Vector< real >( tNumberOfNodes ) );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                Vector< real > & tPhi = mPhi( k );
                const Matrix< real > & tN = mN( k );
                for( uint i=0; i< tNumberOfNodes; ++i )
                {
                    tPhi( i ) = tN( 0, i );
                }
            }


            uint tNumDim = mesh::dimension( mShapeFunction->element_type() ) ;

            mNvector.set_size(
                    mNumberOfIntegrationPoints,
                    Matrix< real >( tNumDim, tNumDim*tNumberOfNodes,
                                    0.0 ) );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Vector< real > & tPhi = mPhi( k );
                Matrix< real > & tN = mNvector( k );
                uint tCount = 0 ;
                for( uint i=0; i< tNumberOfNodes; ++i )
                {
                    for( uint j=0; j<tNumDim; ++j )
                    {
                        tN( j, tCount++ ) = tPhi( i );
                    }
                }
            }

            // allocate cells for first derivative
            mdNdXi.set_size( mNumberOfIntegrationPoints,
                             Matrix< real >( tNumDim, tNumberOfNodes ) );

            // evaluate shape function
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mShapeFunction->dNdXi( mPoints.col( k ), mdNdXi( k ) );
            }


            // also write derivative into vector
            // ( special function for line elements only)
            if( tNumDim == 1 )
            {
                mdPhidxi.set_size( mNumberOfIntegrationPoints, Vector< real >( tNumberOfNodes ) );
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    Vector< real > & tPhi_xi = mdPhidxi( k );
                    const Matrix< real > & tdNdXi  = mdNdXi( k );

                    for( uint i=0; i<tNumberOfNodes; ++i )
                    {
                        tPhi_xi( i ) = tdNdXi( 0, i );
                    }
                }
            }

            // allocate cells for second derivative
            md2NdXi2.set_size( mNumberOfIntegrationPoints,
                               Matrix< real >( tNumDim, tNumberOfNodes ) );

            // evaluate shape function
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mShapeFunction->d2NdXi2( mPoints.col( k ), md2NdXi2( k ) );
            }
        }

//------------------------------------------------------------------------------
    }
}