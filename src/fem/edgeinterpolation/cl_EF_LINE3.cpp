//
// Created by christian on 8/11/22.
//

#include "cl_EF_LINE3.hpp"
#include "fn_dot.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_LINE3::EF_LINE3()
        {
            mX.set_size(  3, 1 );
            mY.set_size( 3, 1 );
            mJ.set_size( 2, 2 );


        }
//------------------------------------------------------------------------------

        EF_LINE3::~EF_LINE3()
        {

        }

//------------------------------------------------------------------------------

        void
        EF_LINE3::precompute( const Matrix< real > & aXi )
        {
            // number of integration points
            uint tNumIntPoints = aXi.n_cols() ;

            // evaluate the node function
            mNxi.set_size( tNumIntPoints, { 0,0,0 } );
            for( uint k=0; k<tNumIntPoints; ++k )
            {
                real xi = aXi( 0, k );
                Vector< real > & tdNdxi = mNxi( k );
                tdNdxi( 0 ) =  xi - 0.5;
                tdNdxi( 1 ) =  xi + 0.5;
                tdNdxi( 2 ) = -2.0 * xi;
            }

            // evaluate the edge function
            mE.set_size( tNumIntPoints, Matrix< real >(1,2) );

            for( uint k=0; k<tNumIntPoints; ++k )
            {
                real xi = aXi( 0, k );
                Matrix< real > & tE = mE( k );
                tE( 0, 0 ) =  1.5 * xi + 0.5 ;
                tE( 0, 1 ) = -1.5 * xi + 0.5 ;
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