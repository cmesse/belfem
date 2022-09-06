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
                    mW = ( real * ) malloc( 30 * sizeof( real ) );
                    mHaveW = true ;
                    mVolumeFunctionLinear = & Pipette::measure_penta6 ;
                    break ;
                }
                case( GeometryType::HEX ) :
                {
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
                    // overwrite, using simpler function
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
            return 0.5*std::abs((mX[1]-mX[0])*(mY[2]-mY[1])-(mY[1]-mY[0])*(mX[2]-mX[1]));
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
            return  0.5*std::abs((mY[1]-mY[3])*(mX[0]-mX[2])
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
            mW[ 3 ] = mX[ 1 ] - mX[ 3 ];
            mW[ 4 ] = mY[ 1 ] - mY[ 3 ];
            mW[ 5 ] = mZ[ 1 ] - mZ[ 3 ];
            mW[ 6 ] = mX[ 2 ] - mX[ 3 ];
            mW[ 7 ] = mY[ 2 ] - mY[ 3 ];
            mW[ 8 ] = mZ[ 2 ] - mZ[ 3 ];

            // compute the determinant
            return std::abs(
                      mW[ 0 ]*(mW[ 4 ]*mW[ 8 ] - mW[ 5 ]*mW[ 7 ])
                    + mW[ 1 ]*(mW[ 5 ]*mW[ 6 ] - mW[ 3 ]*mW[ 8 ])
                    + mW[ 2 ]*(mW[ 3 ]*mW[ 7 ] - mW[ 4 ]*mW[ 6 ]) )/6. ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_penta6( const Element * aElement )
        {
            // collect the node coordinates
            mX[ 0 ] = aElement->node( 0 )->x() ;
            mX[ 1 ] = aElement->node( 1 )->x() ;
            mX[ 2 ] = aElement->node( 2 )->x() ;
            mX[ 3 ] = aElement->node( 3 )->x() ;
            mX[ 4 ] = aElement->node( 4 )->x() ;
            mX[ 5 ] = aElement->node( 5 )->x() ;

            mY[ 0 ] = aElement->node( 0 )->y() ;
            mY[ 1 ] = aElement->node( 1 )->y() ;
            mY[ 2 ] = aElement->node( 2 )->y() ;
            mY[ 3 ] = aElement->node( 3 )->y() ;
            mY[ 4 ] = aElement->node( 4 )->y() ;
            mY[ 5 ] = aElement->node( 5 )->y() ;

            mZ[ 0 ] = aElement->node( 0 )->z() ;
            mZ[ 1 ] = aElement->node( 1 )->z() ;
            mZ[ 2 ] = aElement->node( 2 )->z() ;
            mZ[ 3 ] = aElement->node( 3 )->z() ;
            mZ[ 4 ] = aElement->node( 4 )->z() ;
            mZ[ 5 ] = aElement->node( 5 )->z() ;

            // Work Vector
            mW[  0 ] = mX[ 1 ] * mY[ 0 ];
            mW[  1 ] = mX[ 2 ] * mY[ 0 ];
            mW[  2 ] = mX[ 3 ] * mY[ 0 ];
            mW[  3 ] = mX[ 4 ] * mY[ 0 ];
            mW[  4 ] = mX[ 5 ] * mY[ 0 ];
            mW[  5 ] = mX[ 0 ] * mY[ 1 ];
            mW[  6 ] = mX[ 2 ] * mY[ 1 ];
            mW[  7 ] = mX[ 3 ] * mY[ 1 ];
            mW[  8 ] = mX[ 4 ] * mY[ 1 ];
            mW[  9 ] = mX[ 5 ] * mY[ 1 ];
            mW[ 10 ] = mX[ 0 ] * mY[ 2 ];
            mW[ 11 ] = mX[ 1 ] * mY[ 2 ];
            mW[ 12 ] = mX[ 3 ] * mY[ 2 ];
            mW[ 13 ] = mX[ 4 ] * mY[ 2 ];
            //mW[ 14 ] = mX[ 5 ] * mY[ 2 ];
            mW[ 15 ] = mX[ 0 ] * mY[ 3 ];
            mW[ 16 ] = mX[ 1 ] * mY[ 3 ];
            mW[ 17 ] = mX[ 2 ] * mY[ 3 ];
            mW[ 18 ] = mX[ 4 ] * mY[ 3 ];
            mW[ 19 ] = mX[ 5 ] * mY[ 3 ];
            mW[ 20 ] = mX[ 0 ] * mY[ 4 ];
            mW[ 21 ] = mX[ 1 ] * mY[ 4 ];
            mW[ 22 ] = mX[ 2 ] * mY[ 4 ];
            mW[ 23 ] = mX[ 3 ] * mY[ 4 ];
            mW[ 24 ] = mX[ 5 ] * mY[ 4 ];
            mW[ 25 ] = mX[ 0 ] * mY[ 5 ];
            mW[ 26 ] = mX[ 1 ] * mY[ 5 ];
            //mW[ 27 ] = mX[ 2 ] * mY[ 5 ];
            mW[ 28 ] = mX[ 3 ] * mY[ 5 ];
            mW[ 29 ] = mX[ 4 ] * mY[ 5 ];

            return std::abs(
                     (mZ[3]-mZ[5])*(mW[5]-mW[0]+2.*(mW[21]-mW[8]))
                    +(mZ[2]-mZ[0])*(mW[23]-mW[18]-2.*(mW[21]-mW[8]))
                    +(mZ[4]-mZ[0])*(mW[17]-mW[12]+2.*(mW[11]-mW[6]))
                    +(mZ[3]-mZ[1])*(mW[4]-mW[25]+2.*(mW[29]-mW[24]))
                    +(mZ[5]-mZ[1])*(mW[20]-mW[3]+2.*(mW[23]-mW[18]))
                    +(mZ[2]-mZ[4])*(mW[20]-mW[3]+2.*(mW[0]-mW[5]))
                    +2.*mZ[1]*(mW[10]-mW[13]+mW[22]-mW[1])
                    +2.*mZ[4]*(mW[26]-mW[16]+mW[7]-mW[28]+mW[19]-mW[9])
                    +(mZ[0]+mZ[3])*(mW[9]-mW[26]+mW[13]-mW[22])
                    +(mZ[2]+mZ[3])*(mW[20]-mW[3])
                    +(mZ[0]+mZ[5])*(mW[16]-mW[7])
                    +(mZ[3]+mZ[4])*(mW[1]-mW[10])
                    +(mZ[1]+mZ[0])*(mW[28]-mW[19])
                    +(mZ[1]-mZ[2]+mZ[4]-mZ[5])*(mW[2]-mW[15])
                    +(mZ[4]-mZ[2]-mZ[1]-mZ[5])*(mW[20]-mW[3]) )/12.0 ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_hex8( const Element * aElement )
        {
            mX[ 0 ] = aElement->node( 0 )->x() ;
            mX[ 1 ] = aElement->node( 1 )->x() ;
            mX[ 2 ] = aElement->node( 2 )->x() ;
            mX[ 3 ] = aElement->node( 3 )->x() ;
            mX[ 4 ] = aElement->node( 4 )->x() ;
            mX[ 5 ] = aElement->node( 5 )->x() ;
            mX[ 6 ] = aElement->node( 6 )->x() ;
            mX[ 7 ] = aElement->node( 7 )->x() ;

            mY[ 0 ] = aElement->node( 0 )->y() ;
            mY[ 1 ] = aElement->node( 1 )->y() ;
            mY[ 2 ] = aElement->node( 2 )->y() ;
            mY[ 3 ] = aElement->node( 3 )->y() ;
            mY[ 4 ] = aElement->node( 4 )->y() ;
            mY[ 5 ] = aElement->node( 5 )->y() ;
            mY[ 6 ] = aElement->node( 6 )->y() ;
            mY[ 7 ] = aElement->node( 7 )->y() ;

            mZ[ 0 ] = aElement->node( 0 )->z() ;
            mZ[ 1 ] = aElement->node( 1 )->z() ;
            mZ[ 2 ] = aElement->node( 2 )->z() ;
            mZ[ 3 ] = aElement->node( 3 )->z() ;
            mZ[ 4 ] = aElement->node( 4 )->z() ;
            mZ[ 5 ] = aElement->node( 5 )->z() ;
            mZ[ 6 ] = aElement->node( 6 )->z() ;
            mZ[ 7 ] = aElement->node( 7 )->z() ;

            // compute the value
            return 0.5 * std::abs(
            (mZ[1]+mZ[2]-mZ[4]-mZ[7])*(mX[0]*mY[3]-mX[3]*mY[0]+mX[6]*mY[5]-mX[5]*mY[6])
            +(mZ[4]+mZ[0]-mZ[2]-mZ[6])*(mX[1]*mY[5]-mX[5]*mY[1]+mX[7]*mY[3]-mX[3]*mY[7])
            +(mZ[5]+mZ[6]-mZ[3]-mZ[0])*(mX[1]*mY[2]-mX[2]*mY[1]+mX[7]*mY[4]-mX[4]*mY[7])
            +(mZ[6]+mZ[7]-mZ[1]-mZ[0])*(mX[2]*mY[3]-mX[3]*mY[2]+mX[4]*mY[5]-mX[5]*mY[4])
            +(mZ[4]+mZ[5]-mZ[2]-mZ[3])*(mX[0]*mY[1]-mX[1]*mY[0])
            +(mZ[3]+mZ[7]-mZ[1]-mZ[5])*(mX[0]*mY[4]-mX[4]*mY[0])
            +(mZ[1]+mZ[5]-mZ[3]-mZ[7])*(mX[2]*mY[6]-mX[6]*mY[2])
            +(mZ[4]+mZ[5]-mZ[2]-mZ[3])*(mX[6]*mY[7]-mX[7]*mY[6])
            +(mZ[1]-mZ[3])*(mX[0]*mY[2]-mX[2]*mY[0])
            +(mZ[4]-mZ[1])*(mX[0]*mY[5]-mX[5]*mY[0])
            +(mZ[3]-mZ[4])*(mX[0]*mY[7]-mX[7]*mY[0])
            +(mZ[0]-mZ[5])*(mX[1]*mY[4]-mX[4]*mY[1])
            +(mZ[5]-mZ[2])*(mX[1]*mY[6]-mX[6]*mY[1])
            +(mZ[1]-mZ[6])*(mX[2]*mY[5]-mX[5]*mY[2])
            +(mZ[6]-mZ[3])*(mX[2]*mY[7]-mX[7]*mY[2])
            +(mZ[7]-mZ[0])*(mX[3]*mY[4]-mX[4]*mY[3])
            +(mZ[2]-mZ[7])*(mX[3]*mY[6]-mX[6]*mY[3])
            +(mZ[7]-mZ[5])*(mX[4]*mY[6]-mX[6]*mY[4])
            +(mZ[4]-mZ[6])*(mX[5]*mY[7]-mX[7]*mY[5])
            +(mZ[0]-mZ[2])*(mX[3]*mY[1]-mY[3]*mX[1]));
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
                        std::abs( det( mIntegrationData->dNdXi( k ) * mNodeCoords ) );

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