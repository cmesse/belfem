//
// Created by christian on 12/3/21.
//

#include "cl_EF_TET4.hpp"
#include "fn_inv.hpp"
#include "fn_det.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_TET4::EF_TET4()
        {
            mJ.set_size( 3, 3 );
            mInvJ.set_size( 3, 3 );
            mE.set_size( 3, 6 );
            mC.set_size( 3, 6 );
            mB.set_size( 3, 4 );
            mCurlA.set_size( 3, 12, 0.0 );
        }

//------------------------------------------------------------------------------

        void
        EF_TET4::link( Element * aElement,
                       const bool aOperatorCurlH,
                       const bool aOperatorGradPhi,
                       const bool aOperatorCurlA )
        {
            // grab nodes
            mesh::Node * tNode0 = aElement->element()->node( 0 );
            mesh::Node * tNode1 = aElement->element()->node( 1 );
            mesh::Node * tNode2 = aElement->element()->node( 2 );
            mesh::Node * tNode3 = aElement->element()->node( 3 );

            // components from the Jacobian
            real & a = mJ( 0, 0 );
            real & d = mJ( 1, 0 );
            real & g = mJ( 2, 0 );
            real & b = mJ( 0, 1 );
            real & e = mJ( 1, 1 );
            real & h = mJ( 2, 1 );
            real & c = mJ( 0, 2 );
            real & f = mJ( 1, 2 );
            real & i = mJ( 2, 2 );

            // components from the inversed Jacobian
            real & l = mInvJ( 0, 0 );
            real & p = mInvJ( 1, 0 );
            real & u = mInvJ( 2, 0 );
            real & m = mInvJ( 0, 1 );
            real & q = mInvJ( 1, 1 );
            real & v = mInvJ( 2, 1 );
            real & n = mInvJ( 0, 2 );
            real & r = mInvJ( 1, 2 );
            real & w = mInvJ( 2, 2 );

            // compute the values for the jacobian
            a = tNode0->x() - tNode3->x();
            d = tNode1->x() - tNode3->x();
            g = tNode2->x() - tNode3->x();
            b = tNode0->y() - tNode3->y();
            e = tNode1->y() - tNode3->y();
            h = tNode2->y() - tNode3->y();
            c = tNode0->z() - tNode3->z();
            f = tNode1->z() - tNode3->z();
            i = tNode2->z() - tNode3->z();

            // compute the determinant
            mDetJ = a*(e*i - h*f) + b*(f*g - d*i) + c*(d*h - e*g);

            // compute the absolute value
            mAbsDetJ = std::abs( mDetJ );

            // inverse components of J
            l = e*i-f*h;
            p = f*g-d*i;
            u = d*h-e*g;
            m = c*h-b*i;
            q = a*i-c*g;
            v = b*g-a*h;
            n = b*f-c*e;
            r = c*d-a*f;
            w = a*e-b*d;

            mInvJ /= mDetJ ;

            if( aOperatorGradPhi )
            {
                mB( 0, 0 ) = l;
                mB( 1, 0 ) = p;
                mB( 2, 0 ) = u;
                mB( 0, 1 ) = m;
                mB( 1, 1 ) = q;
                mB( 2, 1 ) = v;
                mB( 0, 2 ) = n;
                mB( 1, 2 ) = r;
                mB( 2, 2 ) = w;
                mB( 0, 3 ) = -l-m-n;
                mB( 1, 3 ) = -p-q-r;
                mB( 2, 3 ) = -u-v-w;
            }

            if( aOperatorCurlH )
            {
                // grab signs
                aElement->edge_directions( mS );

                // extract nabla variables
                mNablaXi[ 0 ] = l ;
                mNablaXi[ 1 ] = p ;
                mNablaXi[ 2 ] = u ;

                mNablaEta[ 0 ] = m ;
                mNablaEta[ 1 ] = q ;
                mNablaEta[ 2 ] = v ;

                mNablaZeta[ 0 ] = n ;
                mNablaZeta[ 1 ] = r ;
                mNablaZeta[ 2 ] = w ;

                mNablaTau[ 0 ] = -(l+m+n);
                mNablaTau[ 1 ] = -(p+q+r);
                mNablaTau[ 2 ] = -(u+v+w);

                // compute curl
                mC( 0, 0 ) = mS[ 0 ] * ( p*v - q*u );
                mC( 1, 0 ) = mS[ 0 ] * ( m*u - l*v );
                mC( 2, 0 ) = mS[ 0 ] * ( l*q - m*p );

                mC( 0, 1 ) = mS[ 1 ] * ( q*w - r*v );
                mC( 1, 1 ) = mS[ 1 ] * ( n*v - m*w );
                mC( 2, 1 ) = mS[ 1 ] * ( m*r - n*q );

                mC( 0, 2 ) = mS[ 2 ] * ( r*u - p*w );
                mC( 1, 2 ) = mS[ 2 ] * ( l*w - n*u );
                mC( 2, 2 ) = mS[ 2 ] * ( n*p - l*r );

                mC( 0, 3 ) = mS[ 3 ] * ( u*(q+r) - p*(v+w) );
                mC( 1, 3 ) = mS[ 3 ] * ( l*(v+w) - u*(m+n) );
                mC( 2, 3 ) = mS[ 3 ] * ( p*(n+m) - l*(q+r) );

                mC( 0, 4 ) = mS[ 4 ] * ( v*(p+r) - q*(u+w) );
                mC( 1, 4 ) = mS[ 4 ] * ( m*(u+w) - v*(l+n) );
                mC( 2, 4 ) = mS[ 4 ] * ( q*(l+n) - m*(p+r) );

                mC( 0, 5 ) = mS[ 5 ] * ( w*(p+q) - r*(u+v) );
                mC( 1, 5 ) = mS[ 5 ] * ( n*(u+v) - w*(l+m) );
                mC( 2, 5 ) = mS[ 5 ] * ( r*(l+m) - n*(p+q) );

                mC *= 2.0 ;

            }

            if( aOperatorCurlA )
            {
                mCurlA( 1,  0 ) =  u;
                mCurlA( 2,  0 ) = -p;
                mCurlA( 0,  1 ) = -u;
                mCurlA( 2,  1 ) =  l;
                mCurlA( 0,  2 ) =  p;
                mCurlA( 1,  2 ) = -l;
                mCurlA( 1,  3 ) =  v;
                mCurlA( 2,  3 ) = -q;
                mCurlA( 0,  4 ) = -v;
                mCurlA( 2,  4 ) =  m;
                mCurlA( 0,  5 ) =  q;
                mCurlA( 1,  5 ) = -m;
                mCurlA( 1,  6 ) =  w;
                mCurlA( 2,  6 ) = -r;
                mCurlA( 0,  7 ) = -w;
                mCurlA( 2,  7 ) =  n;
                mCurlA( 0,  8 ) =  r;
                mCurlA( 1,  8 ) = -n;
                mCurlA( 1,  9 ) = -u-v-w;
                mCurlA( 2,  9 ) =  p+q+r;
                mCurlA( 0, 10 ) =  u+v+w;
                mCurlA( 2, 10 ) = -l-m-n;
                mCurlA( 0, 11 ) = -p-q-r;
                mCurlA( 1, 11 ) =  l+m+n ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET4::precompute( const Matrix< real > & aXi )
        {
            // get number of integration points
            uint tN = aXi.n_cols();

            // compute first factor
            mG.set_size( 6, tN );
            for( uint k=0; k<tN; ++k )
            {
                mG( 0, k ) = aXi( 0, k );
                mG( 1, k ) = aXi( 1, k );
                mG( 2, k ) = aXi( 2, k );
                mG( 3, k ) = aXi( 0, k );
                mG( 4, k ) = aXi( 1, k );
                mG( 5, k ) = aXi( 2, k );
            }

            // compute second factor
            mH.set_size( 6, tN );
            for( uint k=0; k<tN; ++k )
            {
                mH( 0, k ) = aXi( 1, k );
                mH( 1, k ) = aXi( 2, k );
                mH( 2, k ) = aXi( 0, k );
                mH( 3, k ) = 1. - aXi( 0, k ) - aXi( 1, k ) - aXi( 2, k );
                mH( 4, k ) = 1. - aXi( 0, k ) - aXi( 1, k ) - aXi( 2, k );
                mH( 5, k ) = 1. - aXi( 0, k ) - aXi( 1, k ) - aXi( 2, k );
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TET4::E( const uint aIndex )
        {
            // edge 0 : xi -> zeta
            mE( 0, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaEta[ 0 ] - mH( 0, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaEta[ 1 ] - mH( 0, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaEta[ 2 ] - mH( 0, aIndex ) * mNablaXi[ 2 ] );

            // edge 1 : zeta -> eta
            mE( 0, 1 ) = mS[ 1 ] * ( mG( 1, aIndex ) * mNablaZeta[ 0 ] - mH( 1, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 1 ) = mS[ 1 ] * ( mG( 1, aIndex ) * mNablaZeta[ 1 ] - mH( 1, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 1 ) = mS[ 1 ] * ( mG( 1, aIndex ) * mNablaZeta[ 2 ] - mH( 1, aIndex ) * mNablaEta[ 2 ] );

            // edge 2 : eta -> xi
            mE( 0, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 0 ] - mH( 2, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 1 ] - mH( 2, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 2 ] - mH( 2, aIndex ) * mNablaZeta[ 2 ] );

            // edge 3 : xi -> tau
            mE( 0, 3 ) = mS[ 3 ] * ( mG( 3, aIndex ) * mNablaTau[ 0 ] - mH( 3, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 3 ) = mS[ 3 ] * ( mG( 3, aIndex ) * mNablaTau[ 1 ] - mH( 3, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 3 ) = mS[ 3 ] * ( mG( 3, aIndex ) * mNablaTau[ 2 ] - mH( 3, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 : zeta -> tau
            mE( 0, 4 ) = mS[ 4 ] * ( mG( 4, aIndex ) * mNablaTau[ 0 ] - mH( 4, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 4 ) = mS[ 4 ] * ( mG( 4, aIndex ) * mNablaTau[ 1 ] - mH( 4, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 4 ) = mS[ 4 ] * ( mG( 4, aIndex ) * mNablaTau[ 2 ] - mH( 4, aIndex ) * mNablaEta[ 2 ] );

            // edge 5 : eta -> tau
            mE( 0, 5 ) = mS[ 5 ] * ( mG( 5, aIndex ) * mNablaTau[ 0 ] - mH( 5, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 5 ) = mS[ 5 ] * ( mG( 5, aIndex ) * mNablaTau[ 1 ] - mH( 5, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 5 ) = mS[ 5 ] * ( mG( 5, aIndex ) * mNablaTau[ 2 ] - mH( 5, aIndex ) * mNablaZeta[ 2 ] );
            return mE ;

        }

//------------------------------------------------------------------------------
    }
}