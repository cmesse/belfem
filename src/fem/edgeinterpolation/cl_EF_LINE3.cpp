//
// Created by christian on 8/11/22.
//

#include "cl_EF_LINE3.hpp"
#include "fn_dot.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "intpoints.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_LINE3::EF_LINE3()
        {
            mXi.set_size( 2 );

            mX.set_size(  3, 1 );
            mY.set_size( 3, 1 );
            mJ.set_size( 2, 2 );

            InterpolationFunctionFactory tFactory ;
            mVolumeLagrange = tFactory.create_lagrange_function( ElementType::TRI6 );
        }

//------------------------------------------------------------------------------

        EF_LINE3::~EF_LINE3()
        {
            delete mVolumeLagrange ;
        }

//------------------------------------------------------------------------------

        void
        EF_LINE3::precompute( const Matrix< real > & aXi )
        {
            // number of integration points
            uint tNumIntPoints = aXi.n_cols() ;

            mXi( 0 ) = aXi( 0, 0 );
            mXi( 1 ) = aXi( 0, tNumIntPoints - 1 );

            // evaluate the node function
            mNxi.set_size( tNumIntPoints+2, { 0,0,0 } );
            for( uint k=0; k<tNumIntPoints; ++k )
            {
                real xi = aXi( 0, k );
                Vector< real > & tdNdxi = mNxi( k );
                tdNdxi( 0 ) =  xi - 0.5;
                tdNdxi( 1 ) =  xi + 0.5;
                tdNdxi( 2 ) = -2.0 * xi;
            }
            for( uint k=0; k<2; ++k )
            {
                real xi = mXi( k );
                Vector< real > & tdNdxi = mNxi( k+tNumIntPoints );
                tdNdxi( 0 ) =  xi - 0.5;
                tdNdxi( 1 ) =  xi + 0.5;
                tdNdxi( 2 ) = -2.0 * xi;
            }

            // we compute the coefficients for the edge function here
            // since we set them in the include file

            // evaluate the edge function
            mE.set_size( tNumIntPoints+2, Matrix< real >(1,2) );


            // the vandermonde matrix
            Matrix< real > tV( 2, 2);
            tV( 0, 0 )  =  1.0 ;
            tV( 1, 0 ) = -mXi( 1 );
            tV( 0, 1 ) = -1.0 ;
            tV( 1, 1 ) =  mXi( 0 );
            tV /= mXi( 0 ) - mXi( 1 );

            Vector< real > tF( 2 );
            tF( 0 ) = 1.0 ;
            tF( 1 ) = 0.0 ;
            mC0 = tV * tF ;
            tF( 0 ) = 0.0 ;
            tF( 1 ) = 1.0 ;
            mC1 = tV * tF ;

            for( uint k=0; k<tNumIntPoints; ++k )
            {
                real xi = aXi( 0, k );
                Matrix< real > & tE = mE( k );
                tE( 0, 0 ) = mC0( 0 ) * xi + mC0( 1 );
                tE( 0, 1 ) = mC1( 0 ) * xi + mC1( 1 );
            }

            for( uint k=0; k<2; ++k )
            {
                real xi = mX( k );
                Matrix< real > & tE = mE( k=tNumIntPoints );
                tE( 0, 0 ) = mC0( 0 ) * xi + mC0( 1 );
                tE( 0, 1 ) = mC1( 0 ) * xi + mC1( 1 );
            }

            // now we precompute the data for the volume element
            mVolumeNxi.set_size( 6, Matrix< real >( 2, 6 ) );

            // local values for xi and ena
            Matrix< real > tXi( 2, 6 );

            // coordinates edge dof 0 side 0
            tXi( 0, 0 ) = 0.5 - 0.5*mXi( 0 );
            tXi( 1, 0 ) = 0.5 + 0.5*mXi( 0 );

            // coordinates edge dof 1 side 0
            tXi( 0, 1 ) = 0.5 - 0.5*mXi( 1 );
            tXi( 1, 1 ) = 0.5 + 0.5*mXi( 1 );

            // coordinates edge dof 0 side 1
            tXi( 0, 2 ) = 0.0 ;
            tXi( 1, 2 ) = 0.5 - 0.5*mXi( 0 );

            // coordinates edge dof 1 side 1
            tXi( 0, 3 ) = 0.0 ;
            tXi( 1, 3 ) = 0.5 - 0.5*mXi( 1 );

            // coordinates edge dof 0 side 2
            tXi( 0, 4 ) = 0.5 + 0.5*mXi( 0 );
            tXi( 1, 4 ) = 0.0 ;

            // coordinates edge dof 1 side 2
            tXi( 0, 5 ) = 0.5 + 0.5*mXi( 1 );
            tXi( 1, 5 ) = 0.0 ;

            for( uint k=0; k<6; ++k )
            {
                mVolumeLagrange->dNdXi( tXi.col( k ),  mVolumeNxi( k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        EF_LINE3::link(
                Element * aElement,
                const bool aOperatorCurlH,
                const bool aOperatorGradPhi,
                const bool aOperatorCurlA )
        {
            // grab facet element
            mesh::Element * tElement = aElement->facet()->element() ;

            BELFEM_ASSERT( tElement->type() == ElementType::LINE3 ,
                "Can only link to LINE3 element" );

            // collect node coords
            mX( 0 ) = tElement->node( 0 )->x() ;
            mX( 1 ) = tElement->node( 1 )->x() ;
            mX( 2 ) = tElement->node( 2 )->x() ;
            mY( 0 ) = tElement->node( 0 )->y() ;
            mY( 1 ) = tElement->node( 1 )->y() ;
            mY( 2 ) = tElement->node( 2 )->y() ;

        }

//------------------------------------------------------------------------------

        real
        EF_LINE3::compute_det_J( const uint aIndex )
        {
            real tdX =  dot( mNxi( aIndex ), mX ) ;
            real tdY =  dot( mNxi( aIndex ), mY ) ;
            return std::sqrt( tdX*tdX + tdY * tdY );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_LINE3::C( const uint aIndex )
        {
            BELFEM_ERROR( false, "function not implemented");
            return mC ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_LINE3::B( const uint aIndex )
        {
            BELFEM_ERROR( false, "function not implemented");
            return mB ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_LINE3::CA( const uint aIndex )
        {
            BELFEM_ERROR( false, "function not implemented");
            return mC ;
        }

//------------------------------------------------------------------------------
    }
}