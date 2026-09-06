/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */
#include <cmath>

#include "cl_Pipette.hpp"

#include <units.hpp>

#include "fn_det.hpp"
#include "fn_norm.hpp"

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

            if( mW != nullptr )
            {
                free( mW );
                mW = nullptr ;
            }
            if( mX != nullptr )
            {
                free( mX );
                mX = nullptr ;
            }
            if( mY != nullptr )
            {
                free( mY );
                mY = nullptr ;
            }
            if( mZ != nullptr )
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
            mSurfaceFunction = nullptr ;

            this->reset_containers() ;

            mNumDim         = dimension( aType );
            mNumNodes       = number_of_nodes( aType ) ;
            mNumCornerNodes = number_of_corner_nodes( aType ) ;

            GeometryType tGeoType = geometry_type( aType ) ;

            // tri and tet elements use full polynomials
            // while for quad, hex, penta and pyra, we have incomplete polynomials as well
            // to account for these, we must use a higher order integration
            uint tOrder = tGeoType == GeometryType::TRI or tGeoType == GeometryType::TET ?
                          interpolation_order_numeric( aType ) :
                          2 * interpolation_order_numeric( aType ) ;

            // allocate general interpolation data
            mIntegrationData = new belfem::fem::IntegrationData( aType ) ;
            mIntegrationData->populate(  tOrder );
            mNumIntPoints = mIntegrationData->weights().length() ;

            // for the linear data, we have to do the same
            tOrder = tGeoType == GeometryType::TRI or tGeoType == GeometryType::TET ? 1 : 2 ;
            mIntegrationDataLinear = new belfem::fem::IntegrationData( aType ) ;
            mIntegrationDataLinear->populate(  tOrder );
            mNumIntPointsLinear = mIntegrationDataLinear->weights().length() ;

            mNodeCoordsLinear.set_size( mNumCornerNodes, mNumDim ) ;
            mNodeCoords.set_size( mNumNodes, mNumDim ) ;

            // now we link the integration functions
            switch( tGeoType )
            {
                case GeometryType::TRI :
                {
                    mX = ( real * ) malloc( 3 * sizeof( real ) );
                    mY = ( real * ) malloc( 3 * sizeof( real ) );
                    mVolumeFunctionLinear = & Pipette::measure_tri3 ;
                    if ( aType == ElementType::TRI3 )
                    {
                        mVolumeFunction = & Pipette::measure_tri3 ;
                    }
                    else
                    {
                        mVolumeFunction = & Pipette::measure_higher_order ;
                    }
                    break ;
                }
                case( GeometryType::QUAD ) :
                {
                    mX = ( real * ) malloc( 4 * sizeof( real ) );
                    mY = ( real * ) malloc( 4 * sizeof( real ) );
                    mVolumeFunctionLinear = & Pipette::measure_quad4 ;
                    if ( aType == ElementType::QUAD4 )
                    {
                        mVolumeFunction = & Pipette::measure_quad4 ;
                    }
                    else if (aType == ElementType::QUAD4TS)
                    {
                        // | signed area | : thin-shell layer elements are
                        // wound clockwise by construction, see measure_quad4ts
                        mVolumeFunctionLinear = & Pipette::measure_quad4ts ;
                        mVolumeFunction = & Pipette::measure_quad4ts ;
                    }
                    else
                    {
                        mVolumeFunction = & Pipette::measure_higher_order ;
                    }
                    break;
                }
                case GeometryType::TET :
                {
                    mX = ( real * ) malloc( 4 * sizeof( real ) );
                    mY = ( real * ) malloc( 4 * sizeof( real ) );
                    mZ = ( real * ) malloc( 4 * sizeof( real ) );

                    mW = ( real * ) malloc( 10 * sizeof( real ) );
                    mW[ 9 ] = 1.0/6.0 ;

                    mVolumeFunctionLinear = & Pipette::measure_tet4 ;

                    if ( aType == ElementType::TET4 )
                    {
                        mVolumeFunction = & Pipette::measure_tet4 ;
                    }
                    else
                    {
                        mVolumeFunction = & Pipette::measure_higher_order ;
                    }
                    break ;
                }
                default:
                {
                    uint n = number_of_corner_nodes( aType );

                    mX = ( real * ) malloc( n * sizeof( real ) );
                    mY = ( real * ) malloc( n * sizeof( real ) );
                    mZ = ( real * ) malloc( n * sizeof( real ) );

                    mVolumeFunction = & Pipette::measure_linear ;
                    if ( interpolation_order_numeric( aType ) == 1 )
                    {
                        mVolumeFunctionLinear = & Pipette::measure_linear ;
                    }
                    else
                    {
                        mVolumeFunction = & Pipette::measure_higher_order ;
                    }
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
                    + mW[ 2 ]*(mW[ 3 ]*mW[ 7 ] - mW[ 4 ]*mW[ 6 ]) ) * mW[ 9 ];
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_linear( const Element * aElement )
        {
            // copy the node coordinates of the element
            for ( uint i=0; i<mNumDim; ++i )
            {
                for ( uint k=0; k<mNumCornerNodes; ++k )
                {
                    mNodeCoordsLinear( k, i ) = aElement->node( k )->x( i ) ;
                }
            }

            // get integration weights
            const Vector< real > & tW = mIntegrationDataLinear->weights() ;

            // the domain
            real aOmega = 0.0 ;

            for( uint k=0; k<mNumIntPointsLinear; ++k )
            {
                aOmega += tW( k ) * det( mIntegrationDataLinear->dNdXi( k ) * mNodeCoordsLinear );

            }

            return aOmega ;
        }

//------------------------------------------------------------------------------

        real
        Pipette::measure_quad4ts( const Element * aElement )
        {
            // 2D thin-shell layer elements are wound clockwise by
            // construction ( see Calculator::dV_quad4ts for why that winding
            // is forced ), so the signed integral above comes out negative
            // for a perfectly healthy element. Report the magnitude.
            //
            // Taking the absolute value of the integral equals the integral
            // of the absolute value only while the sign of det J is uniform
            // over the element. That holds for a factory-generated QUAD4TS,
            // which is a non-inverting extrusion of its facet; it would NOT
            // hold for a folded bilinear quad, so this function must not be
            // reused for one.
            return std::abs( this->measure_linear( aElement ) );
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
            for ( uint i=0; i<mNumDim; ++i )
            {
                for ( uint k=0; k<mNumNodes; ++k )
                {
                    mNodeCoords( k, i ) = aElement->node( k )->x( i );
                }
            }

            // get integration weights
            const Vector< real > & tW = mIntegrationData->weights() ;

            // the domain ( surface in 2D, volume in 3D)
            real aOmega = 0.0 ;

            for( uint k=0; k<mNumIntPoints; ++k )
            {
                aOmega += tW( k ) *
                        det( mIntegrationData->dNdXi( k ) * mNodeCoords ) ;

            }

            return aOmega ;
        }

//------------------------------------------------------------------------------


        void
        Pipette::set_facet_type( const ElementType aType )
        {
            mVolumeFunction = nullptr ;

            this->reset_containers() ;

            uint d = dimension( aType ) + 1 ;
            uint n = number_of_nodes( aType );

            mN.set_size( d );

            mJ.set_size( d, n );
            mNodeCoords.set_size( n, d );

            // allocate general interpolation data
            if ( aType == ElementType::TRI3 )
            {
                mSurfaceFunction = & Pipette::measure_surface_tri3 ;
            }
            else if (aType == ElementType::LINE2)
            {
                mSurfaceFunction = & Pipette::measure_surface_line2 ;
            }
            else
            {
                mSurfaceFunction = & Pipette::measure_surface_higher_order ;

                mIntegrationData = new fem::IntegrationData( aType ) ;

                uint tOrder = interpolation_order_numeric( aType ) ;

                if ( geometry_type( aType ) == GeometryType::QUAD )
                {
                    tOrder *= 2 ;
                }
                mIntegrationData->populate( tOrder  );
                mNumIntPoints = mIntegrationData->weights().length() ;
            }

        }

        real
        Pipette::measure_surface_tri3( const Facet * aFacet )
        {

            // grab nodes
            const Node * tNode0 = aFacet->node( 0 );
            const Node * tNode1 = aFacet->node( 1 );
            const Node * tNode2 = aFacet->node( 2 );

            // compute jacobian
            mJ( 0, 0 ) = tNode0->x() - tNode2->x();
            mJ( 1, 0 ) = tNode1->x() - tNode2->x();

            mJ( 0, 1 ) = tNode0->y() - tNode2->y();
            mJ( 1, 1 ) = tNode1->y() - tNode2->y();

            mJ( 0, 2 ) = tNode0->z() - tNode2->z();
            mJ( 1, 2 ) = tNode1->z() - tNode2->z();

            // compute the normal
            mN( 0 ) = mJ( 0, 1 ) * mJ( 1, 2 ) - mJ( 0, 2 ) * mJ( 1, 1 );
            mN( 1 ) = mJ( 0, 2 ) * mJ( 1, 0 ) - mJ( 0, 0 ) * mJ( 1, 2 );
            mN( 2 ) = mJ( 0, 0 ) * mJ( 1, 1 ) - mJ( 0, 1 ) * mJ( 1, 0 );

            real tNorm = norm( mN );

            mN /= tNorm ;

            return 0.5 * tNorm ;
        }

        real
        Pipette::measure_surface_line2( const Facet * aFacet )
        {

            // grab nodes
            const Node * tNode0 = aFacet->node( 0 );
            const Node * tNode1 = aFacet->node( 1 );

            // compute jacobian
            mJ( 0, 0 ) = tNode0->x() - tNode1->x();

            mJ( 0, 1 ) = tNode0->y() - tNode1->y();

            // compute the normal
            mN( 0 ) = mJ( 0, 1 ) ;
            mN( 1 ) = - mJ( 0, 0 ) ;

            real aLength = norm( mN );

            mN /= aLength ;

            return aLength ;
        }

        real
        Pipette::measure_surface_higher_order( const Facet * aFacet )
        {
            // copy the node coordinates of the element
            // copy the node coordinates of the element
            for ( uint i=0; i<mNumDim; ++i )
            {
                for ( uint k=0; k<mNumNodes; ++k )
                {
                    mNodeCoords( k, i ) = aFacet->node( k )->x( i );
                }
            }

            // get integration weights
            const Vector< real > & tW = mIntegrationData->weights() ;

            real aOmega = 0.0 ;

            for( uint k=0; k<mNumIntPoints; ++k )
            {
                // compute the jacobian (note: it is not quadratic here!)
                mJ = mIntegrationData->dNdXi( k ) * mNodeCoords ;

                // compute the normal
                mN( 0 ) = mJ( 0, 1 ) * mJ( 1, 2 ) - mJ( 0, 2 ) * mJ( 1, 1 );
                mN( 1 ) = mJ( 0, 2 ) * mJ( 1, 0 ) - mJ( 0, 0 ) * mJ( 1, 2 );
                mN( 2 ) = mJ( 0, 0 ) * mJ( 1, 1 ) - mJ( 0, 1 ) * mJ( 1, 0 );

                aOmega += tW( k ) * norm( mN );
            }

            return aOmega ;
        }

//------------------------------------------------------------------------------

    }
}