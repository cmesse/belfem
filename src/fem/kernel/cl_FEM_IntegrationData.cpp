//
// Created by Christian Messe on 18.01.22.
//

#include "cl_FEM_IntegrationData.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "meshtools.hpp"
#include "fn_FEM_initialize_integration_points_on_facet.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IntegrationData::IntegrationData( const ElementType aElementType, InterpolationFunction * aShapeFunction ) :
            mElementType( aElementType )
        {
            // create a temporary factory
            InterpolationFunctionFactory tFactory ;


            if( aShapeFunction == nullptr )
            {
                // create the shape function
                mShapeFunction = tFactory.create_lagrange_function( aElementType );

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
            // compute the integration points
            intpoints( aScheme,
                       mesh::geometry_type( mElementType ),
                       aIntegrationOrder == 0 ?
                       auto_integration_order( mElementType ) : aIntegrationOrder,
                       mWeights,
                       mPoints );

            this->evaluate_function();
        }

//------------------------------------------------------------------------------

        void
        IntegrationData::populate_for_master(
                const uint aSideSetIndex,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {
            initialize_integration_points_on_facet( mElementType,
                                                    aSideSetIndex,
                                                    mWeights,
                                                    mPoints,
                                                    aIntegrationOrder,
                                                    aScheme );

            this->evaluate_function();

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
                const uint aFaceIndex,
                const uint aOrientation,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {
            // compute integration order
            uint tIntegrationOrder = aIntegrationOrder == 0 ?
                    auto_integration_order( mElementType ) : aIntegrationOrder ;

            facetintpoints::intpoints_tet(
                    aFaceIndex,
                    aOrientation,
                    mWeights,
                    mPoints,
                    tIntegrationOrder,
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