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

            if( aShapeFunction == nullptr )
            {
                // create a temporary factory
                InterpolationFunctionFactory tFactory ;

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

        IntegrationData::IntegrationData( const ElementType aElementType,
                 InterpolationFunction * aShapeFunction,
                 const bool aClaimOwnership ) :
            mElementType( aElementType )
        {
            mShapeFunction    = aShapeFunction ;
            mOwnShapeFunction = aClaimOwnership ;
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
            // Thin-shell elements share their facet topology with the
            // corresponding volume element, so dispatch is by geometry type.
            // Exception: QUAD*TS thin shells need orientation-aware slave
            // integration (2 orientations per LINE facet), while volume
            // TRI/QUAD use the legacy column-reversal path (orientation
            // implicit, count = 1).
            switch ( mesh::geometry_type( mElementType ) )
            {
                case ( GeometryType::TRI ) :
                {
                    this->populate_for_slave_tri( aSlaveIndex,
                                                  aIntegrationOrder,
                                                  aScheme );
                    break ;
                }
                case ( GeometryType::QUAD ) :
                {
                    if ( mElementType == ElementType::QUAD4TS
                      || mElementType == ElementType::QUAD9TS )
                    {
                        this->populate_for_slave_quad( aSlaveIndex,
                                                       aOrientation,
                                                       aIntegrationOrder,
                                                       aScheme );
                    }
                    else
                    {
                        this->populate_for_slave_tri( aSlaveIndex,
                                                      aIntegrationOrder,
                                                      aScheme );
                    }
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
                case( GeometryType::PENTA ) :
                {
                    this->populate_for_slave_penta( aSlaveIndex,
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

        void
        IntegrationData::populate_for_slave_quad(
                const uint aSlaveIndex,
                const uint aOrientation,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {
            facetintpoints::intpoints_quad(
                    aSlaveIndex,
                    aOrientation,
                    mWeights,
                    mPoints,
                    aIntegrationOrder,
                    aScheme );

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
        IntegrationData::populate_for_slave_penta(
                const uint aSlaveIndex,
                const uint aOrientation,
                const uint aIntegrationOrder,
                const IntegrationScheme aScheme )
        {

            facetintpoints::intpoints_penta(
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
            uint tNumBases = mShapeFunction->number_of_bases() ;

            // allocate cells
            mN.set_size( mNumberOfIntegrationPoints,
                         Matrix< real >( 1, tNumBases ) );


            // evaluate shape function
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mShapeFunction->N( mPoints.col( k ), mN( k ) );
            }

            // create phi vector
            mPhi.set_size( mNumberOfIntegrationPoints,
                           Vector< real >( tNumBases ) );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {

                Vector< real > & tPhi = mPhi( k );
                Matrix< real > & tN = mN( k );
                for( uint i=0; i< tNumBases; ++i )
                {
                    if( std::abs( tN( 0, i ) ) < 1e-16 )
                    {
                       tN( 0, i ) = 0.0 ;
                    }

                    tPhi( i ) = tN( 0, i );
                }
            }


            uint tNumDim = mesh::dimension( mShapeFunction->element_type() ) ;

            mNvector.set_size(
                    mNumberOfIntegrationPoints,
                    Matrix< real >( tNumDim, tNumDim*tNumBases,
                                    0.0 ) );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Vector< real > & tPhi = mPhi( k );
                Matrix< real > & tN = mNvector( k );
                uint tCount = 0 ;
                for( uint i=0; i< tNumBases; ++i )
                {
                    for( uint j=0; j<tNumDim; ++j )
                    {
                        tN( j, tCount++ ) = tPhi( i );
                    }
                }
            }

            // allocate cells for first derivative
            mdNdXi.set_size( mNumberOfIntegrationPoints,
                             Matrix< real >( tNumDim, tNumBases ) );

            // evaluate shape function
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mShapeFunction->dNdXi( mPoints.col( k ), mdNdXi( k ) );
            }


            // also write derivative into vector
            // ( special function for line elements only)
            if( tNumDim == 1 )
            {
                mdPhidxi.set_size( mNumberOfIntegrationPoints, Vector< real >( tNumBases ) );
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    Vector< real > & tPhi_xi = mdPhidxi( k );
                    Matrix< real > & tdNdXi  = mdNdXi( k );

                    for( uint i=0; i<tNumBases; ++i )
                    {
                        if( std::abs( tdNdXi( 0, i ) ) < BELFEM_EPSILON )
                        {
                            tdNdXi( 0, i ) = 0.0 ;
                        }
                        tPhi_xi( i ) = tdNdXi( 0, i );
                    }
                }
            }

            // allocate cells for second derivative
            md2NdXi2.set_size( mNumberOfIntegrationPoints,
                               Matrix< real >( tNumDim, tNumBases ) );

            // evaluate shape function
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mShapeFunction->d2NdXi2( mPoints.col( k ), md2NdXi2( k ) );
            }
        }

//------------------------------------------------------------------------------
    }
}