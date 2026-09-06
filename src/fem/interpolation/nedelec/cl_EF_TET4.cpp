//
// Created by christian on 12/3/21.
//

#include "nedelec/cl_EF_TET4.hpp"
#include "fn_inv.hpp"
#include "cl_FEM_Group.hpp"

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
            mGrad.set_size( 9, 6 );

            mNumDofs = 6 ;
            mSumW = 1.0/6.0 ;
        }

//------------------------------------------------------------------------------

        void
        EF_TET4::link( Element * aElement )
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
            d = tNode2->x() - tNode3->x();
            g = tNode1->x() - tNode3->x();

            b = tNode0->y() - tNode3->y();
            e = tNode2->y() - tNode3->y();
            h = tNode1->y() - tNode3->y();

            c = tNode0->z() - tNode3->z();
            f = tNode2->z() - tNode3->z();
            i = tNode1->z() - tNode3->z();

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
            mC( 0, 0 ) = mS[ 0 ] * ( p*w-r*u );
            mC( 1, 0 ) = mS[ 0 ] * ( n*u-l*w );
            mC( 2, 0 ) = mS[ 0 ] * ( l*r-n*p );

            mC( 0, 1 ) = mS[ 1 ] * ( r*v-q*w );
            mC( 1, 1 ) = mS[ 1 ] * ( m*w-n*v );
            mC( 2, 1 ) = mS[ 1 ] * ( n*q-m*r );

            mC( 0, 2 ) = mS[ 2 ] * ( q*u-p*v );
            mC( 1, 2 ) = mS[ 2 ] * ( l*v-m*u );
            mC( 2, 2 ) = mS[ 2 ] * ( m*p-l*q );

            mC( 0, 3 ) = mS[ 3 ] * ( u*(q+r)-p*(v+w) );
            mC( 1, 3 ) = mS[ 3 ] * ( l*(v+w)-u*(m+n) );
            mC( 2, 3 ) = mS[ 3 ] * ( p*(n+m)-l*(q+r) );

            mC( 0, 4 ) = mS[ 4 ] * ( w*(p+q)-r*(u+v) );
            mC( 1, 4 ) = mS[ 4 ] * ( n*(u+v)-w*(l+m) );
            mC( 2, 4 ) = mS[ 4 ] * ( r*(l+m)-n*(p+q) );

            mC( 0, 5 ) = mS[ 5 ] * ( v*(p+r)-q*(u+w) );
            mC( 1, 5 ) = mS[ 5 ] * ( m*(u+w)-v*(l+n) );
            mC( 2, 5 ) = mS[ 5 ] * ( q*(l+n)-m*(p+r) );

            mC *= 2.0 ;

            // compute the gradient operator ( constant for this element,
            // layout contract in EdgeFunction::mGrad ); pair table identical
            // to E(): e0 (xi,zeta), e1 (zeta,eta), e2 (eta,xi),
            //         e3 (xi,tau),  e4 (zeta,tau), e5 (eta,tau)
            // mGrad is sized in the constructor, never resized here
            const real * tA[ 6 ] = { mNablaXi, mNablaZeta, mNablaEta,
                                     mNablaXi, mNablaZeta, mNablaEta };
            const real * tB[ 6 ] = { mNablaZeta, mNablaEta, mNablaXi,
                                     mNablaTau,  mNablaTau, mNablaTau };

            // tEdge/tRow/tCol instead of e/i/j: the Jacobian reference
            // aliases e and i are live in this scope and must not be
            // shadowed ( a later cleanup dropping the uint would smash mJ )
            for ( uint tEdge = 0; tEdge < 6; ++tEdge )
            {
                for ( uint tCol = 0; tCol < 3; ++tCol )        // field component
                {
                    for ( uint tRow = 0; tRow < 3; ++tRow )    // derivative direction
                    {
                        mGrad( tRow + 3 * tCol, tEdge ) = mS[ tEdge ] *
                              ( tA[ tEdge ][ tRow ] * tB[ tEdge ][ tCol ]
                              - tB[ tEdge ][ tRow ] * tA[ tEdge ][ tCol ] );
                    }
                }
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
                mG( 1, k ) = aXi( 2, k );
                mG( 2, k ) = aXi( 1, k );
                mG( 3, k ) = aXi( 0, k );
                mG( 4, k ) = aXi( 2, k );
                mG( 5, k ) = aXi( 1, k );
            }

            // compute second factor
            mH.set_size( 6, tN );
            for( uint k=0; k<tN; ++k )
            {
                mH( 0, k ) = aXi( 2, k );
                mH( 1, k ) = aXi( 1, k );
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
            // xi* ∇*zeta  - zeta * ∇*xi
            mE( 0, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaZeta[ 0 ] - mH( 0, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaZeta[ 1 ] - mH( 0, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaZeta[ 2 ] - mH( 0, aIndex ) * mNablaXi[ 2 ] );

            // edge 1 : zeta -> eta
            // zeta* ∇*eta - eta * ∇*zeta
            mE( 0, 1 ) = mS[ 1 ] * ( mG( 1, aIndex ) * mNablaEta[ 0 ] - mH( 1, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 1 ) = mS[ 1 ] * ( mG( 1, aIndex ) * mNablaEta[ 1 ] - mH( 1, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 1 ) = mS[ 1 ] * ( mG( 1, aIndex ) * mNablaEta[ 2 ] - mH( 1, aIndex ) * mNablaZeta[ 2 ] );

            // edge 2 : eta -> xi
            // 
            // eta* ∇*xi  - xi ∇*eta
            mE( 0, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 0 ] - mH( 2, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 1 ] - mH( 2, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 2 ) = mS[ 2 ] * ( mG( 2, aIndex ) * mNablaXi[ 2 ] - mH( 2, aIndex ) * mNablaEta[ 2 ] );

            // edge 3 : xi -> tau
            // xi * ∇*tau - tau * ∇*xi
            mE( 0, 3 ) = mS[ 3 ] * ( mG( 3, aIndex ) * mNablaTau[ 0 ] - mH( 3, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 3 ) = mS[ 3 ] * ( mG( 3, aIndex ) * mNablaTau[ 1 ] - mH( 3, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 3 ) = mS[ 3 ] * ( mG( 3, aIndex ) * mNablaTau[ 2 ] - mH( 3, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 : zeta -> tau
            // zeta * ∇*tau - tau *  ∇*zeta
            mE( 0, 4 ) = mS[ 4 ] * ( mG( 4, aIndex ) * mNablaTau[ 0 ] - mH( 4, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 4 ) = mS[ 4 ] * ( mG( 4, aIndex ) * mNablaTau[ 1 ] - mH( 4, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 4 ) = mS[ 4 ] * ( mG( 4, aIndex ) * mNablaTau[ 2 ] - mH( 4, aIndex ) * mNablaZeta[ 2 ] );

            // edge 5 : eta -> tau
            // eta *  ∇*tau - tau * ∇*eta
            mE( 0, 5 ) = mS[ 5 ] * ( mG( 5, aIndex ) * mNablaTau[ 0 ] - mH( 5, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 5 ) = mS[ 5 ] * ( mG( 5, aIndex ) * mNablaTau[ 1 ] - mH( 5, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 5 ) = mS[ 5 ] * ( mG( 5, aIndex ) * mNablaTau[ 2 ] - mH( 5, aIndex ) * mNablaEta[ 2 ] );
            return mE ;

        }

//------------------------------------------------------------------------------
    }
}