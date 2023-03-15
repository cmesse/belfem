//
// Created by christian on 8/23/22.
//
#include "cl_Pipette.hpp"
#include "fn_det.hpp"
#include "geometrytools.hpp"

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        Pipette::Pipette()
        {

        }

//-----------------------------------------------------------------------------

        Pipette::~Pipette()
        {
            this->reset_containers();
        }

//-----------------------------------------------------------------------------

        void
        Pipette::reset_containers()
        {


            if( mIntegrationData != nullptr )
            {
                delete mIntegrationData ;
                mIntegrationData = nullptr ;
            }

            if( mIntegrationDataLinear != nullptr )
            {
                delete mIntegrationDataLinear ;
                mIntegrationDataLinear = nullptr ;
            }

            if( mHaveW )
            {
                free( mW );
                mW = nullptr ;
            }
            if( mNumNodes != 0 )
            {
                free( mX );
                mX = nullptr ;
            }
            if( mNumNodes != 0 )
            {
                free( mY );
                mY = nullptr ;
            }
            if( mNumNodes != 0 )
            {
                free( mZ );
                mZ = nullptr ;
            }

            mNumDim = 0 ;
            mNumNodes = 0 ;
            mNumIntPoints = 0 ;
            mNumIntPointsLinear = 0 ;
        }

//-----------------------------------------------------------------------------

        void
        Pipette::set_element_type( const ElementType aType )
        {
            this->reset_containers() ;

            mType = aType ;
            mNumDim   = mesh::dimension( aType );
            mNumNodes = mesh::number_of_nodes( aType );

            mX = ( real * ) malloc( mNumNodes * sizeof( real ) );
            mY = ( real * ) malloc( mNumNodes * sizeof( real ) );
            mZ = ( real * ) malloc( mNumNodes * sizeof( real ) );

            mIntegrationData = new belfem::fem::IntegrationData( aType ) ;
            mIntegrationData->populate(  mesh::interpolation_order_numeric( aType ) );
            mNumIntPoints = mIntegrationData->weights().length() ;


            mNodeCoordsFacet.set_size( mNumNodes, 3 );
            mNodeCoords.set_size( mNumNodes, mNumDim );

            // set the linear elemenf function
            switch( mesh::geometry_type( aType ) )
            {
                case( GeometryType::TRI ) :
                {
                    mW = ( real * ) malloc( 3 * sizeof( real ) );
                    mHaveW = true ;
                    mVolumeFunctionLinear = & Pipette::measure_tri3 ;
                    mSurfaceFunction      = & Pipette::measure_higher_order_tri ;
                    break ;
                }
                case( GeometryType::QUAD ) :
                {
                    mVolumeFunctionLinear = & Pipette::measure_quad4 ;
                    mSurfaceFunction = & Pipette::measure_higher_order_quad ;
                    break ;
                }
                case( GeometryType::TET ) :
                {
                    mW = ( real * ) malloc( 9 * sizeof( real ) );
                    mHaveW = true ;
                    mVolumeFunctionLinear = & Pipette::measure_tet4 ;
                    break ;
                }
                case( GeometryType::PENTA ) :
                {
                    mNodeCoordsLinear.set_size( 6, 3 );
                    mIntegrationDataLinear = new belfem::fem::IntegrationData( ElementType::PENTA6 ) ;
                    mIntegrationDataLinear->populate(  2 );
                    mNumIntPointsLinear = mIntegrationDataLinear->weights().length() ;
                    mVolumeFunctionLinear = & Pipette::measure_penta6 ;
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    mNodeCoordsLinear.set_size( 8, 3 );
                    mIntegrationDataLinear = new belfem::fem::IntegrationData( ElementType::HEX8 ) ;
                    mIntegrationDataLinear->populate(  2 );
                    mNumIntPointsLinear = mIntegrationDataLinear->weights().length() ;
                    mVolumeFunctionLinear = & Pipette::measure_hex8 ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unsupported Geometry Type");
                }
            }

            switch( aType )
            {
                case( ElementType::TRI3 ) :
                {
                    mVolumeFunction  = & Pipette::measure_tri3 ;
                    mSurfaceFunction = & Pipette::measure_linear_tri ;
                    break ;
                }
                case( ElementType::QUAD4 ) :
                {
                    mVolumeFunction    = & Pipette::measure_quad4 ;
                    break ;
                }
                case( ElementType::TET4 ) :
                {
                    mVolumeFunction = & Pipette::measure_tet4 ;
                    mSurfaceFunction = nullptr ;
                    break ;
                }
                case( ElementType::PENTA6 ) :
                {
                    mVolumeFunction = & Pipette::measure_penta6 ;
                    mSurfaceFunction = nullptr ;
                    break ;
                }
                case( ElementType::HEX8 ) :
                {
                    mVolumeFunction = & Pipette::measure_hex8 ;
                    mSurfaceFunction = nullptr ;
                    break ;
                }
                default :
                {
                    mVolumeFunction = & Pipette::measure_higher_order ;
                    mSurfaceFunction = nullptr ;
                    break ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Pipette::collect_node_coords( const Element * aElement )
        {
            for( uint j=0; j<mNumDim; ++j )
            {
                for( uint i=0; i<mNumNodes; ++i )
                {
                    mNodeCoords( i, j )
                        = aElement->node( i )->x( j );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Pipette::collect_node_coords( const Facet * aFacet )
        {
            for( uint j=0; j<3; ++j )
            {
                for( uint i=0; i<mNumNodes; ++i )
                {
                    mNodeCoordsFacet( i, j )
                            = aFacet->node( i )->x( j );
                }
            }
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_tri3( const Element * aElement )
        {
            // collect the node coordinates
            mX[ 0 ] = aElement->node( 0 )->x() ;
            mX[ 1 ] = aElement->node( 1 )->x() ;
            mX[ 2 ] = aElement->node( 2 )->x() ;

            mY[ 0 ] = aElement->node( 0 )->y() ;
            mY[ 1 ] = aElement->node( 1 )->y() ;
            mY[ 2 ] = aElement->node( 2 )->y() ;

            // compute the value
            return 0.5*((mX[1]-mX[0])*(mY[2]-mY[1])-(mY[1]-mY[0])*(mX[2]-mX[1]));
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_quad4( const Element * aElement )
        {
            // collect the node coordinates
            mX[ 0 ] = aElement->node( 0 )->x() ;
            mX[ 1 ] = aElement->node( 1 )->x() ;
            mX[ 2 ] = aElement->node( 2 )->x() ;
            mX[ 3 ] = aElement->node( 3 )->x() ;

            mY[ 0 ] = aElement->node( 0 )->y() ;
            mY[ 1 ] = aElement->node( 1 )->y() ;
            mY[ 2 ] = aElement->node( 2 )->y() ;
            mY[ 3 ] = aElement->node( 3 )->y() ;

            // compute the value
            return  0.5*((mY[1]-mY[3])*(mX[0]-mX[2])
                                  +(mY[0]-mY[2])*(mX[3]-mX[1]));
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_tet4( const Element * aElement )
        {
            // collect the node coordinates
            mX[ 0 ] = aElement->node( 0 )->x() ;
            mX[ 1 ] = aElement->node( 1 )->x() ;
            mX[ 2 ] = aElement->node( 2 )->x() ;
            mX[ 3 ] = aElement->node( 3 )->x() ;

            mY[ 0 ] = aElement->node( 0 )->y() ;
            mY[ 1 ] = aElement->node( 1 )->y() ;
            mY[ 2 ] = aElement->node( 2 )->y() ;
            mY[ 3 ] = aElement->node( 3 )->y() ;

            mZ[ 0 ] = aElement->node( 0 )->z() ;
            mZ[ 1 ] = aElement->node( 1 )->z() ;
            mZ[ 2 ] = aElement->node( 2 )->z() ;
            mZ[ 3 ] = aElement->node( 3 )->z() ;

            // Compute Jacobian
            mW[ 0 ] = mX[ 0 ] - mX[ 3 ];
            mW[ 1 ] = mY[ 0 ] - mY[ 3 ];
            mW[ 2 ] = mZ[ 0 ] - mZ[ 3 ];
            mW[ 3 ] = mX[ 2 ] - mX[ 3 ];
            mW[ 4 ] = mY[ 2 ] - mY[ 3 ];
            mW[ 5 ] = mZ[ 2 ] - mZ[ 3 ];
            mW[ 6 ] = mX[ 1 ] - mX[ 3 ];
            mW[ 7 ] = mY[ 1 ] - mY[ 3 ];
            mW[ 8 ] = mZ[ 1 ] - mZ[ 3 ];

            // compute the determinant
            return (
                      mW[ 0 ]*(mW[ 4 ]*mW[ 8 ] - mW[ 5 ]*mW[ 7 ])
                    + mW[ 1 ]*(mW[ 5 ]*mW[ 6 ] - mW[ 3 ]*mW[ 8 ])
                    + mW[ 2 ]*(mW[ 3 ]*mW[ 7 ] - mW[ 4 ]*mW[ 6 ]) )/6. ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_penta6( const Element * aElement )
        {
            // copy the node coordinates of the element
            mNodeCoordsLinear( 0, 0 ) = aElement->node( 0 )->x() ;
            mNodeCoordsLinear( 1, 0 ) = aElement->node( 1 )->x() ;
            mNodeCoordsLinear( 2, 0 ) = aElement->node( 2 )->x() ;
            mNodeCoordsLinear( 3, 0 ) = aElement->node( 3 )->x() ;
            mNodeCoordsLinear( 4, 0 ) = aElement->node( 4 )->x() ;
            mNodeCoordsLinear( 5, 0 ) = aElement->node( 5 )->x() ;

            mNodeCoordsLinear( 0, 1 ) = aElement->node( 0 )->y() ;
            mNodeCoordsLinear( 1, 1 ) = aElement->node( 1 )->y() ;
            mNodeCoordsLinear( 2, 1 ) = aElement->node( 2 )->y() ;
            mNodeCoordsLinear( 3, 1 ) = aElement->node( 3 )->y() ;
            mNodeCoordsLinear( 4, 1 ) = aElement->node( 4 )->y() ;
            mNodeCoordsLinear( 5, 1 ) = aElement->node( 5 )->y() ;

            mNodeCoordsLinear( 0, 2 ) = aElement->node( 0 )->z() ;
            mNodeCoordsLinear( 1, 2 ) = aElement->node( 1 )->z() ;
            mNodeCoordsLinear( 2, 2 ) = aElement->node( 2 )->z() ;
            mNodeCoordsLinear( 3, 2 ) = aElement->node( 3 )->z() ;
            mNodeCoordsLinear( 4, 2 ) = aElement->node( 4 )->z() ;
            mNodeCoordsLinear( 5, 2 ) = aElement->node( 5 )->z() ;
            // get integration weights
            const Vector< real > & tW = mIntegrationDataLinear->weights() ;

            // the domain
            real aOmega = 0.0 ;

            for( uint k=0; k<mNumIntPointsLinear; ++k )
            {
                aOmega += tW( k ) *
                          std::abs( det( mIntegrationDataLinear->dNdXi( k ) * mNodeCoordsLinear ) );

            }

            return aOmega ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_hex8( const Element * aElement )
        {
            // copy the node coordinates of the element
            mNodeCoordsLinear( 0, 0 ) = aElement->node( 0 )->x() ;
            mNodeCoordsLinear( 1, 0 ) = aElement->node( 1 )->x() ;
            mNodeCoordsLinear( 2, 0 ) = aElement->node( 2 )->x() ;
            mNodeCoordsLinear( 3, 0 ) = aElement->node( 3 )->x() ;
            mNodeCoordsLinear( 4, 0 ) = aElement->node( 4 )->x() ;
            mNodeCoordsLinear( 5, 0 ) = aElement->node( 5 )->x() ;
            mNodeCoordsLinear( 6, 0 ) = aElement->node( 6 )->x() ;
            mNodeCoordsLinear( 7, 0 ) = aElement->node( 7 )->x() ;

            mNodeCoordsLinear( 0, 1 ) = aElement->node( 0 )->y() ;
            mNodeCoordsLinear( 1, 1 ) = aElement->node( 1 )->y() ;
            mNodeCoordsLinear( 2, 1 ) = aElement->node( 2 )->y() ;
            mNodeCoordsLinear( 3, 1 ) = aElement->node( 3 )->y() ;
            mNodeCoordsLinear( 4, 1 ) = aElement->node( 4 )->y() ;
            mNodeCoordsLinear( 5, 1 ) = aElement->node( 5 )->y() ;
            mNodeCoordsLinear( 6, 1 ) = aElement->node( 6 )->y() ;
            mNodeCoordsLinear( 7, 1 ) = aElement->node( 7 )->y() ;

            mNodeCoordsLinear( 0, 2 ) = aElement->node( 0 )->z() ;
            mNodeCoordsLinear( 1, 2 ) = aElement->node( 1 )->z() ;
            mNodeCoordsLinear( 2, 2 ) = aElement->node( 2 )->z() ;
            mNodeCoordsLinear( 3, 2 ) = aElement->node( 3 )->z() ;
            mNodeCoordsLinear( 4, 2 ) = aElement->node( 4 )->z() ;
            mNodeCoordsLinear( 5, 2 ) = aElement->node( 5 )->z() ;
            mNodeCoordsLinear( 6, 2 ) = aElement->node( 6 )->z() ;
            mNodeCoordsLinear( 7, 2 ) = aElement->node( 7 )->z() ;

            // get integration weights
            const Vector< real > & tW = mIntegrationDataLinear->weights() ;

            // the domain
            real aOmega = 0.0 ;

            for( uint k=0; k<mNumIntPointsLinear; ++k )
            {
                aOmega += tW( k ) *
                          ( det( mIntegrationDataLinear->dNdXi( k ) * mNodeCoordsLinear ) );

            }

            return aOmega ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_higher_order( const Element * aElement )
        {
            // use simple function if element has straight edges
            if( ! aElement->is_curved() )
            {
                return ( this->*mVolumeFunctionLinear )( aElement );
            }

            // copy the node coordinates of the element
            this->collect_node_coords( aElement );

            // get integration weights
            const Vector< real > & tW = mIntegrationData->weights() ;

            // the domain ( surface in 2D, volume in 3D)
            real aOmega = 0.0 ;

            for( uint k=0; k<mNumIntPoints; ++k )
            {
                aOmega += tW( k ) *
                        ( det( mIntegrationData->dNdXi( k ) * mNodeCoords ) );

            }

            return aOmega ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_higher_order_quad( const Facet * aFacet )
        {
            this->collect_node_coords( aFacet );

            const Vector< real > & tW = mIntegrationData->weights() ;

            real aGamma = 0.0 ;
            for( uint k=0; k<mNumIntPoints; ++k )
            {
                aGamma += tW( k ) * mesh::compute_surface_increment(
                        mIntegrationData->dNdXi( k ),
                        mNodeCoordsFacet );

            }

            return aGamma ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_linear_tri( const Facet * aFacet )
        {
            // collect the node coordinates
            mX[ 0 ] = aFacet->node( 0 )->x() ;
            mX[ 1 ] = aFacet->node( 1 )->x() ;
            mX[ 2 ] = aFacet->node( 2 )->x() ;

            mY[ 0 ] = aFacet->node( 0 )->y() ;
            mY[ 1 ] = aFacet->node( 1 )->y() ;
            mY[ 2 ] = aFacet->node( 2 )->y() ;

            mZ[ 0 ] = aFacet->node( 0 )->z() ;
            mZ[ 1 ] = aFacet->node( 1 )->z() ;
            mZ[ 2 ] = aFacet->node( 2 )->z() ;

            // A = 0.5*norm(cross(p1-p0,p2-p0))

            mW[ 0 ] = (mY[0]-mY[1])*(mZ[0]-mZ[2]) -(mY[0]-mY[2])*(mZ[0]-mZ[1]);
            mW[ 1 ] = (mX[0]-mX[2])*(mZ[0]-mZ[1]) -(mX[0]-mX[1])*(mZ[0]-mZ[2]);
            mW[ 2 ] = (mX[0]-mX[1])*(mY[0]-mY[2]) -(mX[0]-mX[2])*(mY[0]-mY[1]);

            return 0.5 * std::sqrt( mW[0]*mW[0]+mW[1]*mW[1]+mW[2]*mW[2]);
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_higher_order_tri( const Facet * aFacet )
        {
            if( aFacet->is_curved() )
            {
                return this->measure_higher_order_quad( aFacet );
            }
            else
            {
                return this->measure_linear_tri( aFacet );
            }
        }

//------------------------------------------------------------------------------

    }
}