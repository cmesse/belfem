//
// Created by christian on 12/3/21.
//

#include "cl_EF_TRI3.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_dot.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_TRI3::EF_TRI3()
        {
            mJ.set_size( 2, 2 );
            mInvJ.set_size( 2, 2 );
            mE.set_size( 2, 3 );
            mC.set_size( 1, 3 );
            mB.set_size( 2, 3 );
            mCurlA.set_size( 2, 3, 0.0 );

        }

//------------------------------------------------------------------------------

        void
        EF_TRI3::link( Element * aElement,
                       const bool aOperatorCurlH,
                       const bool aOperatorGradPhi,
                       const bool aOperatorCurlA )
        {
            // grab nodes
            mesh::Node * tNode0 = aElement->element()->node( 0 );
            mesh::Node * tNode1 = aElement->element()->node( 1 );
            mesh::Node * tNode2 = aElement->element()->node( 2 );

            // compute jacobian
            mJ( 0, 0 ) = tNode0->x() - tNode2->x();
            mJ( 1, 0 ) = tNode1->x() - tNode2->x();
            mJ( 0, 1 ) = tNode0->y() - tNode2->y();
            mJ( 1, 1 ) = tNode1->y() - tNode2->y();

            // compute the determinant
            mDetJ = det( mJ );

            // compute the absolute value
            mAbsDetJ = std::abs( mDetJ );

            // inverse components
            mInvJ = inv( mJ );

            if( aOperatorCurlH )
            {
                // grab signs
                aElement->edge_directions( mS );

                // compute nabla values
                mNablaXi[ 0 ]   = mInvJ( 0, 0 ) ;
                mNablaXi[ 1 ]   = mInvJ( 1, 0 ) ;
                mNablaEta[ 0 ]  = mInvJ( 0, 1 ) ;
                mNablaEta[ 1 ]  = mInvJ( 1, 1 ) ;
                mNablaZeta[ 0 ] = - mNablaXi[ 0 ] - mNablaEta[ 0 ] ;
                mNablaZeta[ 1 ] = - mNablaXi[ 1 ] - mNablaEta[ 1 ] ;

                // compute the curl matrix (it is constant for this element)
                mC( 0, 0 ) = mS[ 0 ];
                mC( 0, 1 ) = mS[ 1 ];
                mC( 0, 2 ) = mS[ 2 ];
                mC *= 2./mDetJ;
            }

            if( aOperatorGradPhi )
            {
                mB( 0, 0 ) =  mInvJ( 0, 0 );
                mB( 1, 0 ) =  mInvJ( 1, 0 );
                mB( 0, 1 ) =  mInvJ( 0, 1 );
                mB( 1, 1 ) =  mInvJ( 1, 1 );
                mB( 0, 2 ) = -mInvJ( 0, 0 )-mInvJ( 0, 1 );
                mB( 1, 2 ) = -mInvJ( 1, 0 )-mInvJ( 1, 1 );
            }

            if( aOperatorCurlA )
            {
                mCurlA( 0, 0 ) =  mInvJ( 1, 0 );
                mCurlA( 1, 0 ) = -mInvJ( 0, 0 );
                mCurlA( 0, 1 ) =  mInvJ( 1, 1 );
                mCurlA( 1, 1 ) = -mInvJ( 0, 1 );
                mCurlA( 0, 2 ) = -mInvJ( 1, 0 )-mInvJ( 1, 1 );
                mCurlA( 1, 2 ) =  mInvJ( 0, 0 )+mInvJ( 0, 1 ) ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI3::precompute( const Matrix< real > & aXi )
        {
            // get number of integration points
            uint tN = aXi.n_cols();

            // compute first factor
            mG.set_size( 3, tN );
            for( uint k=0; k<tN; ++k )
            {
                mG( 0, k ) = aXi( 0, k );
                mG( 1, k ) = aXi( 1, k );
                mG( 2, k ) = 1.0 - aXi( 0, k ) - aXi( 1, k );
            }

            // compute second factor
            mH.set_size( 3, tN );
            for( uint k=0; k<tN; ++k )
            {
                mH( 0, k ) = aXi( 1, k );
                mH( 1, k ) = 1.0 - aXi( 0, k ) - aXi( 1, k );
                mH( 2, k ) =  aXi( 0, k );
            }
        }

//------------------------------------------------------------------------------

        const Matrix <real> &
        EF_TRI3::E( const uint aIndex )
        {
            // edge 0
            mE( 0, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaEta[ 0 ] - mH( 0, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaEta[ 1 ] - mH( 0, aIndex ) * mNablaXi[ 1 ] );

            // edge 1
            mE( 0, 1 ) = mS[ 1 ] * (  mG( 1, aIndex ) * mNablaZeta[ 0 ] - mH( 1, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 1 ) = mS[ 1 ] * (  mG( 1, aIndex )* mNablaZeta[ 1 ] - mH( 1, aIndex ) * mNablaEta[ 1 ] );

            // edge 2
            mE( 0, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 0 ] - mH( 2, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 1 ] - mH( 2, aIndex ) * mNablaZeta[ 1 ] );

            return mE ;
        }

//------------------------------------------------------------------------------
    }
}