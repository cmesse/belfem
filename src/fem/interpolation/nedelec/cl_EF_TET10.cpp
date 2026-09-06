//
// Created by christian on 12/1/21.
//

#include "nedelec/cl_EF_TET10.hpp"
#include "assert.hpp"
#include "fn_dot.hpp"
#include "cl_Element.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_TET10::EF_TET10()
        {
            mJ.set_size( 3, 3 );
            mInvJ.set_size( 3, 3 );
            mX.set_size( 10 );
            mY.set_size( 10 );
            mZ.set_size( 10 );


            mE.set_size( 3, 20, 0.0 );
            mC.set_size( 3, 20, 0.0 );
            mGrad.set_size( 9, 20, 0.0 );

            mF.set_size( 3, 12, 0.0 );

            mNodeCoords.set_size( 10, 3 );

            mExi.set_size( 3, 20, 0.0 );
            mEeta.set_size( 3, 20, 0.0 );
            mEzeta.set_size( 3, 20, 0.0 );

            mFxi.set_size( 3, 12, 0.0 );
            mFeta.set_size( 3, 12, 0.0 );
            mFzeta.set_size( 3, 12, 0.0 );

            mNumDofs = 20 ;
            mSumW = 1.0/6.0 ;

            // N_xi2, N_eta2, N_zeta2, N_etazeta, N_xizeta, N_xieta
            mD = {{ 4., 0., 0., 4., 0., 0., 0., -8., 0., 0. },
                    { 0., 0., 4., 4., 0., 0., 0., 0., 0., -8. },
                    { 0., 4., 0., 4., 0., 0., 0., 0., -8., 0. },
                    { 0., 0., 0., 4., 0., 4., 0., 0., -4., -4.},
                    { 0., 0., 0., 4., 4., 0., 0., -4., -4., 0.},
                    { 0., 0., 0., 4., 0., 0., 4., -4., 0., -4.}};
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::link( Element * aElement )
        {
            // make sure that this is the correct element type
            BELFEM_ASSERT( aElement->element()->type() == ElementType::TET10,
                          "Element %lu is not of type TET10",
                          ( long unsigned int ) aElement->element()->id() );

            // grab mesh element
            mesh::Element * tElement = aElement->element() ;

            // get x-coordinates
            mX( 0 ) = tElement->node( 0 )->x();
            mX( 1 ) = tElement->node( 1 )->x();
            mX( 2 ) = tElement->node( 2 )->x();
            mX( 3 ) = tElement->node( 3 )->x();
            mX( 4 ) = tElement->node( 4 )->x();
            mX( 5 ) = tElement->node( 5 )->x();
            mX( 6 ) = tElement->node( 6 )->x();
            mX( 7 ) = tElement->node( 7 )->x();
            mX( 8 ) = tElement->node( 8 )->x();
            mX( 9 ) = tElement->node( 9 )->x();

            // get y-coordinates
            mY( 0 ) = tElement->node( 0 )->y();
            mY( 1 ) = tElement->node( 1 )->y();
            mY( 2 ) = tElement->node( 2 )->y();
            mY( 3 ) = tElement->node( 3 )->y();
            mY( 4 ) = tElement->node( 4 )->y();
            mY( 5 ) = tElement->node( 5 )->y();
            mY( 6 ) = tElement->node( 6 )->y();
            mY( 7 ) = tElement->node( 7 )->y();
            mY( 8 ) = tElement->node( 8 )->y();
            mY( 9 ) = tElement->node( 9 )->y();

            // get z-coordinates
            mZ( 0 ) = tElement->node( 0 )->z();
            mZ( 1 ) = tElement->node( 1 )->z();
            mZ( 2 ) = tElement->node( 2 )->z();
            mZ( 3 ) = tElement->node( 3 )->z();
            mZ( 4 ) = tElement->node( 4 )->z();
            mZ( 5 ) = tElement->node( 5 )->z();
            mZ( 6 ) = tElement->node( 6 )->z();
            mZ( 7 ) = tElement->node( 7 )->z();
            mZ( 8 ) = tElement->node( 8 )->z();
            mZ( 9 ) = tElement->node( 9 )->z();

            // reset parameters
            mLastJ = BELFEM_UINT_MAX ;
            mLastNabla = BELFEM_UINT_MAX ;

            aElement->edge_directions( mS );

            // get edge directions
            if(  tElement->is_curved() )
            {
                // link functions
                mFunInterpolation = & EF_TET10::E_curved ;
                mFunDerivatives   = & EF_TET10::E_xi_curved ;
                mFunGrad          = & EF_TET10::G_curved ;

                // contracted second derivatives of the geometry map,
                // constant because the map is quadratic; mD IS the d2N
                // table ( ctor, IF-Voigt rows xx, yy, zz, yz, xz, xy )
                for ( uint v = 0; v < 6; ++v )
                {
                    mCurv[ v ][ 0 ] = 0.0 ;
                    mCurv[ v ][ 1 ] = 0.0 ;
                    mCurv[ v ][ 2 ] = 0.0 ;
                    for ( uint n = 0; n < 10; ++n )
                    {
                        mCurv[ v ][ 0 ] += mD( v, n ) * mX( n );
                        mCurv[ v ][ 1 ] += mD( v, n ) * mY( n );
                        mCurv[ v ][ 2 ] += mD( v, n ) * mZ( n );
                    }
                }
            }
            else // we do this if this element has straight edges
            {
                // get the directions of the edges
                aElement->edge_directions( mS );

                this->compute_nabla( 0 );

                mFunInterpolation  = & EF_TET10::E_straight ;
                mFunDerivatives    = & EF_TET10::E_xi_straight ;
                mFunGrad           = & EF_TET10::G_straight ;

                // the hess channel must never be read on this branch
                // ( S = 0 on the affine map ) — NaN rather than
                // uninitialized so a missed skip screams
                for ( uint v = 0; v < 6; ++v )
                {
                    mCurv[ v ][ 0 ] = BELFEM_QUIET_NAN ;
                    mCurv[ v ][ 1 ] = BELFEM_QUIET_NAN ;
                    mCurv[ v ][ 2 ] = BELFEM_QUIET_NAN ;
                }
            }

            // check face orientation
            for( uint f=0; f<4; ++f )
            {
                if( tElement->face( f )->master()->id() == tElement->id() )
                {
                    mT[ f ] = 0 ;
                }
                else
                {
                    // 1, 2 or 3
                    mT[ f ] = tElement->face( f )->orientation_on_slave() ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::precompute( const Matrix< real > & aXi )
        {
            uint tN = aXi.n_cols() ;

            mNxi.set_size( 10, tN );
            mNeta.set_size( 10, tN );
            mNzeta.set_size( 10, tN );

            mG.set_size( 12, tN );
            mH.set_size( 12, tN );

            mU.set_size( 12, tN );
            mV.set_size( 12, tN );
            mW.set_size( 12, tN );

            mGxi.set_size( 12, tN, 0.0 );
            mGeta.set_size( 12, tN, 0.0  );
            mGzeta.set_size( 12, tN, 0.0  );

            mHxi.set_size( 12, tN, 0.0  );
            mHeta.set_size( 12, tN, 0.0  );
            mHzeta.set_size( 12, tN, 0.0  );

            mUxi.set_size( 12, tN );
            mUeta.set_size( 12, tN );
            mUzeta.set_size( 12, tN  );

            mVxi.set_size( 12, tN );
            mVeta.set_size( 12, tN  );
            mVzeta.set_size( 12, tN  );

            mWxi.set_size( 12, tN );
            mWeta.set_size( 12, tN  );
            mWzeta.set_size( 12, tN );

            for( uint k=0; k<tN; ++k )
            {
                real   xi = aXi( 0, k ) ;
                real  eta = aXi( 1, k ) ;
                real zeta = aXi( 2, k ) ;

                real  tau = 1.-xi-eta-zeta ;

                real   xi2 = xi + xi ;
                real  eta2 = eta + eta ;
                real zeta2 = zeta + zeta ;
                real  tau2 = tau + tau ;

                real   xi4 = xi2 + xi2 ;
                real  eta4 = eta2 + eta2 ;
                real zeta4 = zeta2 + zeta2 ;
                real  tau4 = tau2 + tau2 ;

                real   xi8 = xi4 + xi4 ;
                real  eta8 = eta4 + eta4 ;
                real zeta8 = zeta4 + zeta4 ;
                real tau8 = tau4 + tau4 ;

                real   xi16 = xi8 + xi8 ;
                real  eta16 = eta8 + eta8 ;
                real zeta16 = zeta8 + zeta8 ;
                real tau16 = tau8 + tau8 ;

                real xi32 = xi16 + xi16 ;
                real eta32 = eta16 + eta16 ;
                real zeta32 = zeta16 + zeta16 ;

                mNxi( 0, k ) = xi4-1.0;
                mNxi( 1, k ) = 0.0;
                mNxi( 2, k ) = 0.0;
                mNxi( 3, k ) = (xi4+eta4+zeta4)-3.0;
                mNxi( 4, k ) = zeta4 ;
                mNxi( 5, k ) = 0.0;
                mNxi( 6, k ) = eta4 ;
                mNxi( 7, k ) = 4. - xi4 - xi4 - eta4 - zeta4 ;
                mNxi( 8, k ) = -zeta4 ;
                mNxi( 9, k ) = -eta4 ;

                mNeta( 0, k ) = 0.0;
                mNeta( 1, k ) = 0.0 ;
                mNeta( 2, k ) = eta4-1.0;
                mNeta( 3, k ) = (xi4+eta4+zeta4)-3.0;
                mNeta( 4, k ) = 0.0 ;
                mNeta( 5, k ) = zeta4 ;
                mNeta( 6, k ) = xi4 ;
                mNeta( 7, k ) = -xi4 ;
                mNeta( 8, k ) = -zeta4 ;
                mNeta( 9, k ) = 4. - xi4 - eta4 - eta4 - zeta4 ;

                mNzeta( 0, k ) = 0.0;
                mNzeta( 1, k ) = zeta4-1.0;
                mNzeta( 2, k ) = 0.0 ;
                mNzeta( 3, k ) = (xi4+eta4+zeta4)-3.0;
                mNzeta( 4, k ) = xi4 ;
                mNzeta( 5, k ) = eta4 ;
                mNzeta( 6, k ) = 0.0 ;
                mNzeta( 7, k ) = -xi4 ;
                mNzeta( 8, k ) = 4. - xi4 - eta4 - zeta4 - zeta4 ;
                mNzeta( 9, k ) = -eta4 ;

                // tables corrected 2026-08-14: eta and zeta exchanged relative
                // to the original transcription, which followed the naive node
                // map (EXODUS trap: node 1 carries zeta, node 2 carries eta).
                // Source of truth: tmp/tet10/tet10_function.m / tet10_derivatives.m
                mG(  0, k ) = xi4*(xi2-1.);
                mG(  1, k ) = xi2*(zeta4-1.);
                mG(  2, k ) = zeta4*(zeta2-1.);
                mG(  3, k ) = zeta2*(eta4-1.);
                mG(  4, k ) = eta4*(eta2-1.);
                mG(  5, k ) = eta2*(xi4-1.);
                mG(  6, k ) = xi4*(xi2-1.);
                mG(  7, k ) = xi2*(tau4-1.);
                mG(  8, k ) = zeta4*(zeta2-1.);
                mG(  9, k ) = zeta2*(tau4-1.);
                mG( 10, k ) = eta4*(eta2-1.);
                mG( 11, k ) = eta2*(tau4-1.);

                mH(  0, k ) = zeta2*(xi4-1.);
                mH(  1, k ) = zeta4*(zeta2-1.);
                mH(  2, k ) = eta2*(zeta4-1.);
                mH(  3, k ) = eta4*(eta2-1.);
                mH(  4, k ) = xi2*(eta4-1.);
                mH(  5, k ) = xi4*(xi2-1.);
                mH(  6, k ) = -(1.-xi4)*tau2;
                mH(  7, k ) = tau4*(tau2-1.);
                mH(  8, k ) = -(1.-zeta4)*tau2;
                mH(  9, k ) = tau4*(tau2-1.);
                mH( 10, k ) = -(1.-eta4)*tau2;
                mH( 11, k ) = tau4*(tau2-1.);

                mU(  0, k ) = zeta16*tau;
                mU(  1, k ) = -zeta8*tau;
                mU(  2, k ) = -zeta8*tau;
                mU(  3, k ) = eta16*tau;
                mU(  4, k ) = -eta8*tau;
                mU(  5, k ) = -eta8*tau;
                mU(  6, k ) = eta*tau16 ;
                mU(  7, k ) = -eta*tau8;
                mU(  8, k ) = -eta*tau8;
                mU(  9, k ) = eta*zeta16;
                mU( 10, k ) = -eta*zeta8;
                mU( 11, k ) = -eta*zeta8;

                mV(  0, k ) = -xi*tau8;
                mV(  1, k ) = xi*tau16;
                mV(  2, k ) = -xi*tau8;
                mV(  3, k ) = -zeta*tau8;
                mV(  4, k ) = zeta*tau16;
                mV(  5, k ) = -zeta*tau8;
                mV(  6, k ) = -xi*eta8;
                mV(  7, k ) = xi*eta16;
                mV(  8, k ) = -xi*eta8;
                mV(  9, k ) = -zeta8*xi;
                mV( 10, k ) = zeta16*xi;
                mV( 11, k ) = -zeta8*xi;

                mW(  0, k ) = -zeta8*xi;
                mW(  1, k ) = -zeta8*xi;
                mW(  2, k ) = zeta16*xi;
                mW(  3, k ) = -eta*zeta8;
                mW(  4, k ) = -eta*zeta8;
                mW(  5, k ) = eta*zeta16;
                mW(  6, k ) = -xi8*tau;
                mW(  7, k ) = -xi8*tau;
                mW(  8, k ) = xi16*tau;
                mW(  9, k ) = -xi*eta8;
                mW( 10, k ) = -xi*eta8;
                mW( 11, k ) = xi*eta16;

                mGxi(  0, k ) = xi16-4.;
                mGxi(  1, k ) = zeta8-2.;
                mGxi(  2, k ) = 0.;
                mGxi(  3, k ) = 0.;
                mGxi(  4, k ) = 0.;
                mGxi(  5, k ) = eta8;
                mGxi(  6, k ) = xi16-4.;
                mGxi(  7, k ) = 6.-xi16-zeta8-eta8;
                mGxi(  8, k ) = 0.;
                mGxi(  9, k ) = -zeta8;
                mGxi( 10, k ) = 0.;
                mGxi( 11, k ) = -eta8;

                mGeta(  0, k ) = 0.;
                mGeta(  1, k ) = 0.;
                mGeta(  2, k ) = 0.;
                mGeta(  3, k ) = zeta8;
                mGeta(  4, k ) = eta16-4.;
                mGeta(  5, k ) = xi8-2.;
                mGeta(  6, k ) = 0.;
                mGeta(  7, k ) = -xi8;
                mGeta(  8, k ) = 0.;
                mGeta(  9, k ) = -zeta8;
                mGeta( 10, k ) = eta16-4.;
                mGeta( 11, k ) = 6.-xi8-zeta8-eta16;

                mGzeta(  0, k ) = 0.;
                mGzeta(  1, k ) = xi8;
                mGzeta(  2, k ) = zeta16-4.;
                mGzeta(  3, k ) = eta8-2.;
                mGzeta(  4, k ) = 0.;
                mGzeta(  5, k ) = 0.;
                mGzeta(  6, k ) = 0.;
                mGzeta(  7, k ) = -xi8;
                mGzeta(  8, k ) = zeta16-4.;
                mGzeta(  9, k ) = 6.-xi8-zeta16-eta8;
                mGzeta( 10, k ) = 0.;
                mGzeta( 11, k ) = -eta8;

                mHxi(  0, k ) = zeta8;
                mHxi(  1, k ) = 0.;
                mHxi(  2, k ) = 0.;
                mHxi(  3, k ) = 0.;
                mHxi(  4, k ) = eta8-2.;
                mHxi(  5, k ) = xi16-4.;
                mHxi(  6, k ) = 10-xi16-zeta8-eta8;
                mHxi(  7, k ) = eta16+xi16+zeta16-12.;
                mHxi(  8, k ) = 2.-zeta8;
                mHxi(  9, k ) = eta16+xi16+zeta16-12.;
                mHxi( 10, k ) = 2.-eta8;
                mHxi( 11, k ) = eta16+xi16+zeta16-12.;

                mHeta(  0, k ) = 0.;
                mHeta(  1, k ) = 0.;
                mHeta(  2, k ) = zeta8-2.;
                mHeta(  3, k ) = eta16-4.;
                mHeta(  4, k ) = xi8;
                mHeta(  5, k ) = 0.;
                mHeta(  6, k ) = 2.-xi8;
                mHeta(  7, k ) = eta16+xi16+zeta16-12.;
                mHeta(  8, k ) = 2.-zeta8;
                mHeta(  9, k ) = eta16+xi16+zeta16-12.;
                mHeta( 10, k ) = 10-xi8-zeta8-eta16;
                mHeta( 11, k ) = eta16+xi16+zeta16-12.;

                mHzeta(  0, k ) = xi8-2.;
                mHzeta(  1, k ) = zeta16-4.;
                mHzeta(  2, k ) = eta8;
                mHzeta(  3, k ) = 0.;
                mHzeta(  4, k ) = 0.;
                mHzeta(  5, k ) = 0.;
                mHzeta(  6, k ) = 2.-xi8;
                mHzeta(  7, k ) = eta16+xi16+zeta16-12.;
                mHzeta(  8, k ) = 10-xi8-zeta16-eta8;
                mHzeta(  9, k ) = eta16+xi16+zeta16-12.;
                mHzeta( 10, k ) = 2.-eta8;
                mHzeta( 11, k ) = eta16+xi16+zeta16-12.;

                mUxi(  0, k ) = -zeta16;
                mUxi(  1, k ) = zeta8;
                mUxi(  2, k ) = zeta8;
                mUxi(  3, k ) = -eta16;
                mUxi(  4, k ) = eta8;
                mUxi(  5, k ) = eta8;
                mUxi(  6, k ) = -eta16;
                mUxi(  7, k ) = eta8;
                mUxi(  8, k ) = eta8;
                mUxi(  9, k ) = 0.;
                mUxi( 10, k ) = 0.;
                mUxi( 11, k ) = 0.;

                mUeta(  0, k ) = -zeta16;
                mUeta(  1, k ) = zeta8;
                mUeta(  2, k ) = zeta8;
                mUeta(  3, k ) = 16.-xi16-zeta16-eta32;
                mUeta(  4, k ) = eta16+xi8+zeta8-8.;
                mUeta(  5, k ) = eta16+xi8+zeta8-8.;
                mUeta(  6, k ) = 16.-xi16-zeta16-eta32;
                mUeta(  7, k ) = eta16+xi8+zeta8-8.;
                mUeta(  8, k ) = eta16+xi8+zeta8-8.;
                mUeta(  9, k ) = zeta16;
                mUeta( 10, k ) = -zeta8;
                mUeta( 11, k ) = -zeta8;

                mUzeta(  0, k ) = 16.-xi16-zeta32-eta16;
                mUzeta(  1, k ) = eta8+xi8+zeta16-8.;
                mUzeta(  2, k ) = eta8+xi8+zeta16-8.;
                mUzeta(  3, k ) = -eta16;
                mUzeta(  4, k ) = eta8;
                mUzeta(  5, k ) = eta8;
                mUzeta(  6, k ) = -eta16;
                mUzeta(  7, k ) = eta8;
                mUzeta(  8, k ) = eta8;
                mUzeta(  9, k ) = eta16;
                mUzeta( 10, k ) = -eta8;
                mUzeta( 11, k ) = -eta8;

                mVxi(  0, k ) = eta8+xi16+zeta8-8.;
                mVxi(  1, k ) = 16.-xi32-zeta16-eta16;
                mVxi(  2, k ) = eta8+xi16+zeta8-8.;
                mVxi(  3, k ) = zeta8;
                mVxi(  4, k ) = -zeta16;
                mVxi(  5, k ) = zeta8;
                mVxi(  6, k ) = -eta8;
                mVxi(  7, k ) = eta16;
                mVxi(  8, k ) = -eta8;
                mVxi(  9, k ) = -zeta8;
                mVxi( 10, k ) = zeta16;
                mVxi( 11, k ) = -zeta8;

                mVeta(  0, k ) = xi8;
                mVeta(  1, k ) = -xi16;
                mVeta(  2, k ) = xi8;
                mVeta(  3, k ) = zeta8;
                mVeta(  4, k ) = -zeta16;
                mVeta(  5, k ) = zeta8;
                mVeta(  6, k ) = -xi8;
                mVeta(  7, k ) = xi16;
                mVeta(  8, k ) = -xi8;
                mVeta(  9, k ) = 0.;
                mVeta( 10, k ) = 0.;
                mVeta( 11, k ) = 0.;

                mVzeta(  0, k ) = xi8;
                mVzeta(  1, k ) = -xi16;
                mVzeta(  2, k ) = xi8;
                mVzeta(  3, k ) = eta8+xi8+zeta16-8.;
                mVzeta(  4, k ) = 16.-xi16-zeta32-eta16;
                mVzeta(  5, k ) = eta8+xi8+zeta16-8.;
                mVzeta(  6, k ) = 0.;
                mVzeta(  7, k ) = 0.;
                mVzeta(  8, k ) = 0.;
                mVzeta(  9, k ) = -xi8;
                mVzeta( 10, k ) = xi16;
                mVzeta( 11, k ) = -xi8;

                mWxi(  0, k ) = -zeta8;
                mWxi(  1, k ) = -zeta8;
                mWxi(  2, k ) = zeta16;
                mWxi(  3, k ) = 0.;
                mWxi(  4, k ) = 0.;
                mWxi(  5, k ) = 0.;
                mWxi(  6, k ) = eta8+xi16+zeta8-8.;
                mWxi(  7, k ) = eta8+xi16+zeta8-8.;
                mWxi(  8, k ) = 16.-xi32-zeta16-eta16;
                mWxi(  9, k ) = -eta8;
                mWxi( 10, k ) = -eta8;
                mWxi( 11, k ) = eta16;

                mWeta(  0, k ) = 0.;
                mWeta(  1, k ) = 0.;
                mWeta(  2, k ) = 0.;
                mWeta(  3, k ) = -zeta8;
                mWeta(  4, k ) = -zeta8;
                mWeta(  5, k ) = zeta16;
                mWeta(  6, k ) = xi8;
                mWeta(  7, k ) = xi8;
                mWeta(  8, k ) = -xi16;
                mWeta(  9, k ) = -xi8;
                mWeta( 10, k ) = -xi8;
                mWeta( 11, k ) = xi16;

                mWzeta(  0, k ) = -xi8;
                mWzeta(  1, k ) = -xi8;
                mWzeta(  2, k ) = xi16;
                mWzeta(  3, k ) = -eta8;
                mWzeta(  4, k ) = -eta8;
                mWzeta(  5, k ) = eta16;
                mWzeta(  6, k ) = xi8;
                mWzeta(  7, k ) = xi8;
                mWzeta(  8, k ) = -xi16;
                mWzeta(  9, k ) = 0.;
                mWzeta( 10, k ) = 0.;
                mWzeta( 11, k ) = 0.;

            }
        }

 //------------------------------------------------------------------------------

        void
        EF_TET10::E_xi_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );

            this->compute_edge_derivatives( aIndex );
            this->compute_face_derivatives( aIndex );

            this->combine_functions( mExi, mFxi );
            this->combine_functions( mEeta, mFeta );
            this->combine_functions( mEzeta, mFzeta );
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::E_xi_straight( const uint aIndex )
        {
            this->compute_edge_derivatives( aIndex );
            this->compute_face_derivatives( aIndex );

            this->combine_functions( mExi, mFxi );
            this->combine_functions( mEeta, mFeta );
            this->combine_functions( mEzeta, mFzeta );
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_nabla( const uint aIndex )
        {
            if( aIndex != mLastNabla )
            {
                // J = [ a d g; b e h; c f i ] ;
                mM[ 27 ] = dot( mNxi.col( aIndex ),   mX );  // a
                mM[ 28 ] = dot( mNeta.col( aIndex ),  mX );  // b
                mM[ 29 ] = dot( mNzeta.col( aIndex ), mX );  // c
                mM[ 30 ] = dot( mNxi.col( aIndex ),   mY );  // d
                mM[ 31 ] = dot( mNeta.col( aIndex ),  mY );  // e
                mM[ 32 ] = dot( mNzeta.col( aIndex ), mY );  // f
                mM[ 33 ] = dot( mNxi.col( aIndex ),   mZ );  // d
                mM[ 34 ] = dot( mNeta.col( aIndex ),  mZ );  // e
                mM[ 35 ] = dot( mNzeta.col( aIndex ), mZ );  // f


                // inv(J)*det(J) = [ l p u ; m q v ; n r w ];
                mM[ 36 ] = mM[ 31 ] * mM[ 35 ]  - mM[ 32 ] * mM[ 34 ]; // l = e*i - f*h ;
                mM[ 37 ] = mM[ 29 ] * mM[ 34 ]  - mM[ 28 ] * mM[ 35 ]; // m = c*h - b*i ;
                mM[ 38 ] = mM[ 28 ] * mM[ 32 ]  - mM[ 29 ] * mM[ 31 ]; // n = b*f - c*e ;
                mM[ 39 ] = mM[ 32 ] * mM[ 33 ]  - mM[ 30 ] * mM[ 35 ]; // p = f*g - d*i ;
                mM[ 40 ] = mM[ 27 ] * mM[ 35 ]  - mM[ 29 ] * mM[ 33 ]; // q = a*i - c*g ;
                mM[ 41 ] = mM[ 29 ] * mM[ 30 ]  - mM[ 27 ] * mM[ 32 ]; // r = c*d - a*f ;
                mM[ 42 ] = mM[ 30 ] * mM[ 34 ]  - mM[ 31 ] * mM[ 33 ]; // u = d*h - e*g ;
                mM[ 43 ] = mM[ 28 ] * mM[ 33 ]  - mM[ 27 ] * mM[ 34 ]; // v = b*g - a*h ;
                mM[ 44 ] = mM[ 27 ] * mM[ 31 ]  - mM[ 28 ] * mM[ 30 ]; // w = a*e - b*d ;

                // det( J ) = a*l + b*p + c*u ;
                mDetJ = mM[ 27 ] * mM[ 36 ] + mM[ 28 ] * mM[ 39 ] + mM[ 29 ] * mM[ 42 ] ;

                mAbsDetJ = std::abs( mDetJ ) ;

                real tInvDetJ = 1./mDetJ ;

                // Nablas
                mNablaXi[ 0 ]   = mM[ 36 ] * tInvDetJ ;
                mNablaXi[ 1 ]   = mM[ 37 ] * tInvDetJ ;
                mNablaXi[ 2 ]   = mM[ 38 ] * tInvDetJ ;

                mNablaEta[ 0 ]  = mM[ 39 ] * tInvDetJ ;
                mNablaEta[ 1 ]  = mM[ 40 ] * tInvDetJ ;
                mNablaEta[ 2 ]  = mM[ 41 ] * tInvDetJ ;

                mNablaZeta[ 0 ] = mM[ 42 ] * tInvDetJ ;
                mNablaZeta[ 1 ] = mM[ 43 ] * tInvDetJ ;
                mNablaZeta[ 2 ] = mM[ 44 ] * tInvDetJ ;

                mNablaTau[ 0 ] = -( mNablaXi[ 0 ] + mNablaEta[ 0 ] + mNablaZeta[ 0 ] );
                mNablaTau[ 1 ] = -( mNablaXi[ 1 ] + mNablaEta[ 1 ] + mNablaZeta[ 1 ] );
                mNablaTau[ 2 ] = -( mNablaXi[ 2 ] + mNablaEta[ 2 ] + mNablaZeta[ 2 ] );

                mLastNabla = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_edge_functions( const uint aIndex )
        {
            // edge 0 from xi to zeta
            mE( 0, 0 ) = mS[ 0 ] * ( mG( 0, aIndex )  * mNablaZeta[ 0 ] - mH( 0, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 0 ) = mS[ 0 ] * ( mG( 0, aIndex )  * mNablaZeta[ 1 ] - mH( 0, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 0 ) = mS[ 0 ] * ( mG( 0, aIndex )  * mNablaZeta[ 2 ] - mH( 0, aIndex ) * mNablaXi[ 2 ] );
            mE( 0, 1 ) = mS[ 0 ] * ( mG( 1, aIndex )  * mNablaZeta[ 0 ] - mH( 1, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 1 ) = mS[ 0 ] * ( mG( 1, aIndex )  * mNablaZeta[ 1 ] - mH( 1, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 1 ) = mS[ 0 ] * ( mG( 1, aIndex )  * mNablaZeta[ 2 ] - mH( 1, aIndex ) * mNablaXi[ 2 ] );


            // edge 1 from zeta to eta
            mE( 0, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaEta[ 0 ] -  mH( 2, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaEta[ 1 ] -  mH( 2, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaEta[ 2 ] -  mH( 2, aIndex ) * mNablaZeta[ 2 ] );
            mE( 0, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaEta[ 0 ] -  mH( 3, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaEta[ 1 ] -  mH( 3, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaEta[ 2 ] -  mH( 3, aIndex ) * mNablaZeta[ 2 ] );

            // edge 2 from eta to xi
            mE( 0, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 0 ] - mH( 4, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 1 ] - mH( 4, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 2 ] - mH( 4, aIndex ) * mNablaEta[ 2 ] );
            mE( 0, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 0 ] - mH( 5, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 1 ] - mH( 5, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 2 ] - mH( 5, aIndex ) * mNablaEta[ 2 ] );

            // edge 3 from xi to tau
            mE( 0, 6 ) = mS[ 3 ] * ( mG( 6, aIndex ) * mNablaTau[ 0 ] - mH( 6, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 6 ) = mS[ 3 ] * ( mG( 6, aIndex ) * mNablaTau[ 1 ] - mH( 6, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 6 ) = mS[ 3 ] * ( mG( 6, aIndex ) * mNablaTau[ 2 ] - mH( 6, aIndex ) * mNablaXi[ 2 ] );
            mE( 0, 7 ) = mS[ 3 ] * ( mG( 7, aIndex ) * mNablaTau[ 0 ] - mH( 7, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 7 ) = mS[ 3 ] * ( mG( 7, aIndex ) * mNablaTau[ 1 ] - mH( 7, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 7 ) = mS[ 3 ] * ( mG( 7, aIndex ) * mNablaTau[ 2 ] - mH( 7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 from zeta to tau
            mE( 0, 8 ) = mS[ 4 ] * (  mG( 8, aIndex ) * mNablaTau[ 0 ] -  mH( 8, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 8 ) = mS[ 4 ] * (  mG( 8, aIndex ) * mNablaTau[ 1 ] -  mH( 8, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 8 ) = mS[ 4 ] * (  mG( 8, aIndex ) * mNablaTau[ 2 ] -  mH( 8, aIndex ) * mNablaZeta[ 2 ] );
            mE( 0, 9 ) = mS[ 4 ] * (  mG( 9, aIndex ) * mNablaTau[ 0 ] -  mH( 9, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 9 ) = mS[ 4 ] * (  mG( 9, aIndex ) * mNablaTau[ 1 ] -  mH( 9, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 9 ) = mS[ 4 ] * (  mG( 9, aIndex ) * mNablaTau[ 2 ] -  mH( 9, aIndex ) * mNablaZeta[ 2 ] );

            // edge 5 from  eta to tau
            mE( 0, 10 ) = mS[ 5 ] * ( mG( 10, aIndex ) * mNablaTau[ 0 ] - mH( 10, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 10 ) = mS[ 5 ] * ( mG( 10, aIndex ) * mNablaTau[ 1 ] - mH( 10, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 10 ) = mS[ 5 ] * ( mG( 10, aIndex ) * mNablaTau[ 2 ] - mH( 10, aIndex ) * mNablaEta[ 2 ] );
            mE( 0, 11 ) = mS[ 5 ] * ( mG( 11, aIndex ) * mNablaTau[ 0 ] - mH( 11, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 11 ) = mS[ 5 ] * ( mG( 11, aIndex ) * mNablaTau[ 1 ] - mH( 11, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 11 ) = mS[ 5 ] * ( mG( 11, aIndex ) * mNablaTau[ 2 ] - mH( 11, aIndex ) * mNablaEta[ 2 ] );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_face_functions( const uint aIndex )
        {
            // - - - - - - - - - - - - - - - - - - -
            // FACE 0  : xi->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mF( 0, 0 ) = mU(  0, aIndex ) * mNablaXi[ 0 ] + mV(  0, aIndex ) * mNablaZeta[ 0 ] + mW( 0, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 0 ) = mU(  0, aIndex ) * mNablaXi[ 1 ] + mV(  0, aIndex ) * mNablaZeta[ 1 ] + mW( 0, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 0 ) = mU(  0, aIndex ) * mNablaXi[ 2 ] + mV(  0, aIndex ) * mNablaZeta[ 2 ] + mW( 0, aIndex ) * mNablaTau[ 2 ];

            // zeta
            mF( 0, 1 ) = mU(  1, aIndex ) * mNablaXi[ 0 ] + mV(  1, aIndex ) * mNablaZeta[ 0 ] + mW( 1, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 1 ) = mU(  1, aIndex ) * mNablaXi[ 1 ] + mV(  1, aIndex ) * mNablaZeta[ 1 ] + mW( 1, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 1 ) = mU(  1, aIndex ) * mNablaXi[ 2 ] + mV(  1, aIndex ) * mNablaZeta[ 2 ] + mW( 1, aIndex ) * mNablaTau[ 2 ];

            // tau = -xi-zeta
            mF( 0, 2 ) = -mF( 0, 0 ) - mF( 0, 1 );
            mF( 1, 2 ) = -mF( 1, 0 ) - mF( 1, 1 );
            mF( 2, 2 ) = -mF( 2, 0 ) - mF( 2, 1 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 zeta->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mF( 0, 3 ) = mU( 3, aIndex ) * mNablaZeta[ 0 ] + mV( 3, aIndex ) * mNablaEta[ 0 ] + mW( 3, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 3 ) = mU( 3, aIndex ) * mNablaZeta[ 1 ] + mV( 3, aIndex ) * mNablaEta[ 1 ] + mW( 3, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 3 ) = mU( 3, aIndex ) * mNablaZeta[ 2 ] + mV( 3, aIndex ) * mNablaEta[ 2 ] + mW( 3, aIndex ) * mNablaTau[ 2 ];

            // zeta
            mF( 0, 4 ) = mU( 4, aIndex ) * mNablaZeta[ 0 ] + mV( 4, aIndex ) * mNablaEta[ 0 ] + mW( 4, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 4 ) = mU( 4, aIndex ) * mNablaZeta[ 1 ] + mV( 4, aIndex ) * mNablaEta[ 1 ] + mW( 4, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 4 ) = mU( 4, aIndex ) * mNablaZeta[ 2 ] + mV( 4, aIndex ) * mNablaEta[ 2 ] + mW( 4, aIndex ) * mNablaTau[ 2 ];

            // tau = -eta-zeta
            mF( 0, 5 ) = -mF( 0, 3 ) - mF( 0, 4 );
            mF( 1, 5 ) = -mF( 1, 3 ) - mF( 1, 4 );
            mF( 2, 5 ) = -mF( 2, 3 ) - mF( 2, 4 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 xi->tau->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mF( 0, 6 ) = mU( 6, aIndex ) * mNablaXi[ 0 ] + mV( 6, aIndex ) * mNablaTau[ 0 ] + mW( 6, aIndex ) * mNablaEta[ 0 ];
            mF( 1, 6 ) = mU( 6, aIndex ) * mNablaXi[ 1 ] + mV( 6, aIndex ) * mNablaTau[ 1 ] + mW( 6, aIndex ) * mNablaEta[ 1 ];
            mF( 2, 6 ) = mU( 6, aIndex ) * mNablaXi[ 2 ] + mV( 6, aIndex ) * mNablaTau[ 2 ] + mW( 6, aIndex ) * mNablaEta[ 2 ];

            // tau
            mF( 0, 7 ) = mU( 7, aIndex ) * mNablaXi[ 0 ] + mV( 7, aIndex ) * mNablaTau[ 0 ] + mW( 7, aIndex ) * mNablaEta[ 0 ];
            mF( 1, 7 ) = mU( 7, aIndex ) * mNablaXi[ 1 ] + mV( 7, aIndex ) * mNablaTau[ 1 ] + mW( 7, aIndex ) * mNablaEta[ 1 ];
            mF( 2, 7 ) = mU( 7, aIndex ) * mNablaXi[ 2 ] + mV( 7, aIndex ) * mNablaTau[ 2 ] + mW( 7, aIndex ) * mNablaEta[ 2 ];

            // eta = -xi-tau
            mF( 0, 8 ) = -mF( 0, 6 ) - mF( 0, 7 );
            mF( 1, 8 ) = -mF( 1, 6 ) - mF( 1, 7 );
            mF( 2, 8 ) = -mF( 2, 6 ) - mF( 2, 7 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->eta->zeta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mF( 0, 9 ) = mU( 9, aIndex ) * mNablaXi[ 0 ] + mV( 9, aIndex ) * mNablaEta[ 0 ] + mW( 9, aIndex ) * mNablaZeta[ 0 ];
            mF( 1, 9 ) = mU( 9, aIndex ) * mNablaXi[ 1 ] + mV( 9, aIndex ) * mNablaEta[ 1 ] + mW( 9, aIndex ) * mNablaZeta[ 1 ];
            mF( 2, 9 ) = mU( 9, aIndex ) * mNablaXi[ 2 ] + mV( 9, aIndex ) * mNablaEta[ 2 ] + mW( 9, aIndex ) * mNablaZeta[ 2 ];

            // eta
            mF( 0, 10 ) = mU( 10, aIndex ) * mNablaXi[ 0 ] + mV( 10, aIndex ) * mNablaEta[ 0 ] + mW( 10, aIndex ) * mNablaZeta[ 0 ];
            mF( 1, 10 ) = mU( 10, aIndex ) * mNablaXi[ 1 ] + mV( 10, aIndex ) * mNablaEta[ 1 ] + mW( 10, aIndex ) * mNablaZeta[ 1 ];
            mF( 2, 10 ) = mU( 10, aIndex ) * mNablaXi[ 2 ] + mV( 10, aIndex ) * mNablaEta[ 2 ] + mW( 10, aIndex ) * mNablaZeta[ 2 ];

            // zeta = -xi-eta
            mF( 0, 11 ) = -mF( 0, 9 ) - mF( 0, 10 );
            mF( 1, 11 ) = -mF( 1, 9 ) - mF( 1, 10 );
            mF( 2, 11 ) = -mF( 2, 9 ) - mF( 2, 10 );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_edge_derivatives( const uint aIndex )
        {
           	// edge 0 : xi->zeta
            mExi( 0, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaZeta[ 0 ] - mHxi(  0, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaZeta[ 1 ] - mHxi(  0, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaZeta[ 2 ] - mHxi(  0, aIndex ) * mNablaXi[ 2 ] );
            mExi( 0, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaZeta[ 0 ] - mHxi(  1, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaZeta[ 1 ] - mHxi(  1, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaZeta[ 2 ] - mHxi(  1, aIndex ) * mNablaXi[ 2 ] );

            // edge 1  : zeta->eta
            mExi( 0, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaEta[ 0 ] -  mHxi(  2, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaEta[ 1 ] -  mHxi(  2, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaEta[ 2 ] -  mHxi(  2, aIndex ) * mNablaZeta[ 2 ] );
            mExi( 0, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaEta[ 0 ] -  mHxi(  3, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaEta[ 1 ] -  mHxi(  3, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaEta[ 2 ] -  mHxi(  3, aIndex ) * mNablaZeta[ 2 ] );

            // edge 2 : eta->xi
            mExi( 0, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 0 ] - mHxi(  4, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 1 ] - mHxi(  4, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 2 ] - mHxi(  4, aIndex ) * mNablaEta[ 2 ] );
            mExi( 0, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 0 ] - mHxi(  5, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 1 ] - mHxi(  5, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 2 ] - mHxi(  5, aIndex ) * mNablaEta[ 2 ] );

            // edge 3 : xi->tau
            mExi( 0, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 0 ] - mHxi(  6, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 1 ] - mHxi(  6, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 2 ] - mHxi(  6, aIndex ) * mNablaXi[ 2 ] );
            mExi( 0, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 0 ] - mHxi(  7, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 1 ] - mHxi(  7, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 2 ] - mHxi(  7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 : zeta->tau
            mExi( 0, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 0 ] - mHxi(  8, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 1 ] - mHxi(  8, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 2 ] - mHxi(  8, aIndex ) * mNablaZeta[ 2 ] );
            mExi( 0, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 0 ] - mHxi(  9, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 1 ] - mHxi(  9, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 2 ] - mHxi(  9, aIndex ) * mNablaZeta[ 2 ] );

            // edge 5 : eta->tau
            mExi( 0, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 0 ] - mHxi( 10, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 1 ] - mHxi( 10, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 2 ] - mHxi( 10, aIndex ) * mNablaEta[ 2 ] );
            mExi( 0, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 0 ] - mHxi( 11, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 1 ] - mHxi( 11, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 2 ] - mHxi( 11, aIndex ) * mNablaEta[ 2 ] );

            // edge 0 : xi->zeta
            mEeta( 0, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaZeta[ 0 ] - mHeta(  0, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaZeta[ 1 ] - mHeta(  0, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaZeta[ 2 ] - mHeta(  0, aIndex ) * mNablaXi[ 2 ] );
            mEeta( 0, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaZeta[ 0 ] - mHeta(  1, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaZeta[ 1 ] - mHeta(  1, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaZeta[ 2 ] - mHeta(  1, aIndex ) * mNablaXi[ 2 ] );

            // edge 1  : zeta->eta
            mEeta( 0, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaEta[ 0 ] -  mHeta(  2, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaEta[ 1 ] -  mHeta(  2, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaEta[ 2 ] -  mHeta(  2, aIndex ) * mNablaZeta[ 2 ] );
            mEeta( 0, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaEta[ 0 ] -  mHeta(  3, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaEta[ 1 ] -  mHeta(  3, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaEta[ 2 ] -  mHeta(  3, aIndex ) * mNablaZeta[ 2 ] );

            // edge 2 : eta->xi
            mEeta( 0, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 0 ] - mHeta(  4, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 1 ] - mHeta(  4, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 2 ] - mHeta(  4, aIndex ) * mNablaEta[ 2 ] );
            mEeta( 0, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 0 ] - mHeta(  5, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 1 ] - mHeta(  5, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 2 ] - mHeta(  5, aIndex ) * mNablaEta[ 2 ] );

            // edge 3 : xi->tau
            mEeta( 0, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 0 ] - mHeta(  6, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 1 ] - mHeta(  6, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 2 ] - mHeta(  6, aIndex ) * mNablaXi[ 2 ] );
            mEeta( 0, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 0 ] - mHeta(  7, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 1 ] - mHeta(  7, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 2 ] - mHeta(  7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 : zeta->tau
            mEeta( 0, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 0 ] - mHeta(  8, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 1 ] - mHeta(  8, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 2 ] - mHeta(  8, aIndex ) * mNablaZeta[ 2 ] );
            mEeta( 0, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 0 ] - mHeta(  9, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 1 ] - mHeta(  9, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 2 ] - mHeta(  9, aIndex ) * mNablaZeta[ 2 ] );

            // edge 5 : eta->tau
            mEeta( 0, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 0 ] - mHeta( 10, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 1 ] - mHeta( 10, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 2 ] - mHeta( 10, aIndex ) * mNablaEta[ 2 ] );
            mEeta( 0, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 0 ] - mHeta( 11, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 1 ] - mHeta( 11, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 2 ] - mHeta( 11, aIndex ) * mNablaEta[ 2 ] );

             // edge 0 : xi->zeta
            mEzeta( 0, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaZeta[ 0 ] - mHzeta(  0, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaZeta[ 1 ] - mHzeta(  0, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaZeta[ 2 ] - mHzeta(  0, aIndex ) * mNablaXi[ 2 ] );
            mEzeta( 0, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaZeta[ 0 ] - mHzeta(  1, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaZeta[ 1 ] - mHzeta(  1, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaZeta[ 2 ] - mHzeta(  1, aIndex ) * mNablaXi[ 2 ] );

            // edge 1  : zeta->eta
            mEzeta( 0, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaEta[ 0 ] -  mHzeta(  2, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaEta[ 1 ] -  mHzeta(  2, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaEta[ 2 ] -  mHzeta(  2, aIndex ) * mNablaZeta[ 2 ] );
            mEzeta( 0, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaEta[ 0 ] -  mHzeta(  3, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaEta[ 1 ] -  mHzeta(  3, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaEta[ 2 ] -  mHzeta(  3, aIndex ) * mNablaZeta[ 2 ] );

            // edge 2 : eta->xi
            mEzeta( 0, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 0 ] - mHzeta(  4, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 1 ] - mHzeta(  4, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 2 ] - mHzeta(  4, aIndex ) * mNablaEta[ 2 ] );
            mEzeta( 0, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 0 ] - mHzeta(  5, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 1 ] - mHzeta(  5, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 2 ] - mHzeta(  5, aIndex ) * mNablaEta[ 2 ] );

            // edge 3 : xi->tau
            mEzeta( 0, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 0 ] - mHzeta(  6, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 1 ] - mHzeta(  6, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 2 ] - mHzeta(  6, aIndex ) * mNablaXi[ 2 ] );
            mEzeta( 0, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 0 ] - mHzeta(  7, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 1 ] - mHzeta(  7, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 2 ] - mHzeta(  7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 : zeta->tau
            mEzeta( 0, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 0 ] - mHzeta(  8, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 1 ] - mHzeta(  8, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 2 ] - mHzeta(  8, aIndex ) * mNablaZeta[ 2 ] );
            mEzeta( 0, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 0 ] - mHzeta(  9, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 1 ] - mHzeta(  9, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 2 ] - mHzeta(  9, aIndex ) * mNablaZeta[ 2 ] );

            // edge 5 : eta->tau
            mEzeta( 0, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 0 ] - mHzeta( 10, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 1 ] - mHzeta( 10, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 2 ] - mHzeta( 10, aIndex ) * mNablaEta[ 2 ] );
            mEzeta( 0, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 0 ] - mHzeta( 11, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 1 ] - mHzeta( 11, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 2 ] - mHzeta( 11, aIndex ) * mNablaEta[ 2 ] );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_face_derivatives( const uint aIndex )
        {
            // - - - - - - - - - - - - - - - - - - -
            // FACE 0  : xi->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

			// xi
            mFxi( 0, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 0 ] + mVxi(  0, aIndex ) * mNablaZeta[ 0 ] + mWxi( 0, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 1 ] + mVxi(  0, aIndex ) * mNablaZeta[ 1 ] + mWxi( 0, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 2 ] + mVxi(  0, aIndex ) * mNablaZeta[ 2 ] + mWxi( 0, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFxi( 0, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 0 ] + mVxi(  1, aIndex ) * mNablaZeta[ 0 ] + mWxi( 1, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 1 ] + mVxi(  1, aIndex ) * mNablaZeta[ 1 ] + mWxi( 1, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 2 ] + mVxi(  1, aIndex ) * mNablaZeta[ 2 ] + mWxi( 1, aIndex ) * mNablaTau[ 2 ] ;

			// tau = -xi-zeta
            mFxi( 0, 2 ) = -mFxi( 0, 0 ) - mFxi( 0, 1 );
            mFxi( 1, 2 ) = -mFxi( 1, 0 ) - mFxi( 1, 1 );
            mFxi( 2, 2 ) = -mFxi( 2, 0 ) - mFxi( 2, 1 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 zeta->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFxi( 0, 3 ) = mUxi( 3, aIndex ) * mNablaZeta[ 0 ] + mVxi( 3, aIndex ) * mNablaEta[ 0 ] + mWxi( 3, aIndex ) * mNablaTau[ 0 ]  ;
            mFxi( 1, 3 ) = mUxi( 3, aIndex ) * mNablaZeta[ 1 ] + mVxi( 3, aIndex ) * mNablaEta[ 1 ] + mWxi( 3, aIndex ) * mNablaTau[ 1 ]  ;
            mFxi( 2, 3 ) = mUxi( 3, aIndex ) * mNablaZeta[ 2 ] + mVxi( 3, aIndex ) * mNablaEta[ 2 ] + mWxi( 3, aIndex ) * mNablaTau[ 2 ]  ;

            // zeta
            mFxi( 0, 4 ) = mUxi( 4, aIndex ) * mNablaZeta[ 0 ] + mVxi( 4, aIndex ) * mNablaEta[ 0 ] + mWxi( 4, aIndex ) * mNablaTau[ 0 ]  ;
            mFxi( 1, 4 ) = mUxi( 4, aIndex ) * mNablaZeta[ 1 ] + mVxi( 4, aIndex ) * mNablaEta[ 1 ] + mWxi( 4, aIndex ) * mNablaTau[ 1 ]  ;
            mFxi( 2, 4 ) = mUxi( 4, aIndex ) * mNablaZeta[ 2 ] + mVxi( 4, aIndex ) * mNablaEta[ 2 ] + mWxi( 4, aIndex ) * mNablaTau[ 2 ]  ;

            // tau = -eta-zeta
            mFxi( 0, 5 ) = -mFxi( 0, 3 ) - mFxi( 0, 4 );
            mFxi( 1, 5 ) = -mFxi( 1, 3 ) - mFxi( 1, 4 );
            mFxi( 2, 5 ) = -mFxi( 2, 3 ) - mFxi( 2, 4 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 xi->tau->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFxi( 0, 6 ) = mUxi( 6, aIndex ) * mNablaXi[ 0 ] + mVxi( 6, aIndex ) * mNablaTau[ 0 ] + mWxi( 6, aIndex ) * mNablaEta[ 0 ] ;
            mFxi( 1, 6 ) = mUxi( 6, aIndex ) * mNablaXi[ 1 ] + mVxi( 6, aIndex ) * mNablaTau[ 1 ] + mWxi( 6, aIndex ) * mNablaEta[ 1 ] ;
            mFxi( 2, 6 ) = mUxi( 6, aIndex ) * mNablaXi[ 2 ] + mVxi( 6, aIndex ) * mNablaTau[ 2 ] + mWxi( 6, aIndex ) * mNablaEta[ 2 ] ;

            // tau
            mFxi( 0, 7 ) = mUxi( 7, aIndex ) * mNablaXi[ 0 ] + mVxi( 7, aIndex ) * mNablaTau[ 0 ] + mWxi( 7, aIndex ) * mNablaEta[ 0 ] ;
            mFxi( 1, 7 ) = mUxi( 7, aIndex ) * mNablaXi[ 1 ] + mVxi( 7, aIndex ) * mNablaTau[ 1 ] + mWxi( 7, aIndex ) * mNablaEta[ 1 ] ;
            mFxi( 2, 7 ) = mUxi( 7, aIndex ) * mNablaXi[ 2 ] + mVxi( 7, aIndex ) * mNablaTau[ 2 ] + mWxi( 7, aIndex ) * mNablaEta[ 2 ] ;

            // eta = -xi-tau
            mFxi( 0, 8 ) = -mFxi( 0, 6 ) - mFxi( 0, 7 );
            mFxi( 1, 8 ) = -mFxi( 1, 6 ) - mFxi( 1, 7 );
            mFxi( 2, 8 ) = -mFxi( 2, 6 ) - mFxi( 2, 7 );

			// - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->eta->zeta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFxi( 0, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 0 ] + mVxi( 9, aIndex ) * mNablaEta[ 0 ] + mWxi( 9, aIndex ) * mNablaZeta[ 0 ] ;
            mFxi( 1, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 1 ] + mVxi( 9, aIndex ) * mNablaEta[ 1 ] + mWxi( 9, aIndex ) * mNablaZeta[ 1 ] ;
            mFxi( 2, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 2 ] + mVxi( 9, aIndex ) * mNablaEta[ 2 ] + mWxi( 9, aIndex ) * mNablaZeta[ 2 ] ;

            // eta
            mFxi( 0, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 0 ] + mVxi( 10, aIndex ) * mNablaEta[ 0 ] + mWxi( 10, aIndex ) * mNablaZeta[ 0 ] ;
            mFxi( 1, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 1 ] + mVxi( 10, aIndex ) * mNablaEta[ 1 ] + mWxi( 10, aIndex ) * mNablaZeta[ 1 ] ;
            mFxi( 2, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 2 ] + mVxi( 10, aIndex ) * mNablaEta[ 2 ] + mWxi( 10, aIndex ) * mNablaZeta[ 2 ] ;

            // zeta= = -xi-eta
            mFxi( 0, 11 ) = -mFxi( 0, 9 ) - mFxi( 0, 10 );
            mFxi( 1, 11 ) = -mFxi( 1, 9 ) - mFxi( 1, 10 );
            mFxi( 2, 11 ) = -mFxi( 2, 9 ) - mFxi( 2, 10 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 0  : xi->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

			// xi
            mFeta( 0, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 0 ] + mVeta(  0, aIndex ) * mNablaZeta[ 0 ] + mWeta( 0, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 1 ] + mVeta(  0, aIndex ) * mNablaZeta[ 1 ] + mWeta( 0, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 2 ] + mVeta(  0, aIndex ) * mNablaZeta[ 2 ] + mWeta( 0, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFeta( 0, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 0 ] + mVeta(  1, aIndex ) * mNablaZeta[ 0 ] + mWeta( 1, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 1 ] + mVeta(  1, aIndex ) * mNablaZeta[ 1 ] + mWeta( 1, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 2 ] + mVeta(  1, aIndex ) * mNablaZeta[ 2 ] + mWeta( 1, aIndex ) * mNablaTau[ 2 ] ;

			// tau = -xi-zeta
            mFeta( 0, 2 ) = -mFeta( 0, 0 ) - mFeta( 0, 1 );
            mFeta( 1, 2 ) = -mFeta( 1, 0 ) - mFeta( 1, 1 );
            mFeta( 2, 2 ) = -mFeta( 2, 0 ) - mFeta( 2, 1 );

             // - - - - - - - - - - - - - - - - - - -
            // FACE 1 zeta->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFeta( 0, 3 ) = mUeta( 3, aIndex ) * mNablaZeta[ 0 ] + mVeta( 3, aIndex ) * mNablaEta[ 0 ] + mWeta( 3, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 3 ) = mUeta( 3, aIndex ) * mNablaZeta[ 1 ] + mVeta( 3, aIndex ) * mNablaEta[ 1 ] + mWeta( 3, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 3 ) = mUeta( 3, aIndex ) * mNablaZeta[ 2 ] + mVeta( 3, aIndex ) * mNablaEta[ 2 ] + mWeta( 3, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFeta( 0, 4 ) = mUeta( 4, aIndex ) * mNablaZeta[ 0 ] + mVeta( 4, aIndex ) * mNablaEta[ 0 ] + mWeta( 4, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 4 ) = mUeta( 4, aIndex ) * mNablaZeta[ 1 ] + mVeta( 4, aIndex ) * mNablaEta[ 1 ] + mWeta( 4, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 4 ) = mUeta( 4, aIndex ) * mNablaZeta[ 2 ] + mVeta( 4, aIndex ) * mNablaEta[ 2 ] + mWeta( 4, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -eta-zeta
            mFeta( 0, 5 ) = -mFeta( 0, 3 ) - mFeta( 0, 4 );
            mFeta( 1, 5 ) = -mFeta( 1, 3 ) - mFeta( 1, 4 );
            mFeta( 2, 5 ) = -mFeta( 2, 3 ) - mFeta( 2, 4 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 xi->tau->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFeta( 0, 6 ) = mUeta( 6, aIndex ) * mNablaXi[ 0 ] + mVeta( 6, aIndex ) * mNablaTau[ 0 ] + mWeta( 6, aIndex ) * mNablaEta[ 0 ] ;
            mFeta( 1, 6 ) = mUeta( 6, aIndex ) * mNablaXi[ 1 ] + mVeta( 6, aIndex ) * mNablaTau[ 1 ] + mWeta( 6, aIndex ) * mNablaEta[ 1 ] ;
            mFeta( 2, 6 ) = mUeta( 6, aIndex ) * mNablaXi[ 2 ] + mVeta( 6, aIndex ) * mNablaTau[ 2 ] + mWeta( 6, aIndex ) * mNablaEta[ 2 ] ;

            // tau
            mFeta( 0, 7 ) = mUeta( 7, aIndex ) * mNablaXi[ 0 ] + mVeta( 7, aIndex ) * mNablaTau[ 0 ] + mWeta( 7, aIndex ) * mNablaEta[ 0 ] ;
            mFeta( 1, 7 ) = mUeta( 7, aIndex ) * mNablaXi[ 1 ] + mVeta( 7, aIndex ) * mNablaTau[ 1 ] + mWeta( 7, aIndex ) * mNablaEta[ 1 ] ;
            mFeta( 2, 7 ) = mUeta( 7, aIndex ) * mNablaXi[ 2 ] + mVeta( 7, aIndex ) * mNablaTau[ 2 ] + mWeta( 7, aIndex ) * mNablaEta[ 2 ] ;

            // eta = -xi-tau
            mFeta( 0, 8 ) = -mFeta( 0, 6 ) - mFeta( 0, 7 );
            mFeta( 1, 8 ) = -mFeta( 1, 6 ) - mFeta( 1, 7 );
            mFeta( 2, 8 ) = -mFeta( 2, 6 ) - mFeta( 2, 7 );

			// - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->eta->zeta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFeta( 0, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 0 ] + mVeta( 9, aIndex ) * mNablaEta[ 0 ] + mWeta( 9, aIndex ) * mNablaZeta[ 0 ] ;
            mFeta( 1, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 1 ] + mVeta( 9, aIndex ) * mNablaEta[ 1 ] + mWeta( 9, aIndex ) * mNablaZeta[ 1 ] ;
            mFeta( 2, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 2 ] + mVeta( 9, aIndex ) * mNablaEta[ 2 ] + mWeta( 9, aIndex ) * mNablaZeta[ 2 ] ;

            // eta
            mFeta( 0, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 0 ] + mVeta( 10, aIndex ) * mNablaEta[ 0 ] + mWeta( 10, aIndex ) * mNablaZeta[ 0 ] ;
            mFeta( 1, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 1 ] + mVeta( 10, aIndex ) * mNablaEta[ 1 ] + mWeta( 10, aIndex ) * mNablaZeta[ 1 ] ;
            mFeta( 2, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 2 ] + mVeta( 10, aIndex ) * mNablaEta[ 2 ] + mWeta( 10, aIndex ) * mNablaZeta[ 2 ] ;

            // zeta= = -xi-eta
            mFeta( 0, 11 ) = -mFeta( 0, 9 ) - mFeta( 0, 10 );
            mFeta( 1, 11 ) = -mFeta( 1, 9 ) - mFeta( 1, 10 );
            mFeta( 2, 11 ) = -mFeta( 2, 9 ) - mFeta( 2, 10 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 0  : xi->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

			// xi
            mFzeta( 0, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 0 ] + mVzeta(  0, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 0, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 1 ] + mVzeta(  0, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 0, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 2 ] + mVzeta(  0, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 0, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFzeta( 0, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 0 ] + mVzeta(  1, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 1, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 1 ] + mVzeta(  1, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 1, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 2 ] + mVzeta(  1, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 1, aIndex ) * mNablaTau[ 2 ] ;

			// tau = -xi-zeta
            mFzeta( 0, 2 ) = -mFzeta( 0, 0 ) - mFzeta( 0, 1 );
            mFzeta( 1, 2 ) = -mFzeta( 1, 0 ) - mFzeta( 1, 1 );
            mFzeta( 2, 2 ) = -mFzeta( 2, 0 ) - mFzeta( 2, 1 );

             // - - - - - - - - - - - - - - - - - - -
            // FACE 1 zeta->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFzeta( 0, 3 ) = mUzeta( 3, aIndex ) * mNablaZeta[ 0 ] + mVzeta( 3, aIndex ) * mNablaEta[ 0 ] + mWzeta( 3, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 3 ) = mUzeta( 3, aIndex ) * mNablaZeta[ 1 ] + mVzeta( 3, aIndex ) * mNablaEta[ 1 ] + mWzeta( 3, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 3 ) = mUzeta( 3, aIndex ) * mNablaZeta[ 2 ] + mVzeta( 3, aIndex ) * mNablaEta[ 2 ] + mWzeta( 3, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFzeta( 0, 4 ) = mUzeta( 4, aIndex ) * mNablaZeta[ 0 ] + mVzeta( 4, aIndex ) * mNablaEta[ 0 ] + mWzeta( 4, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 4 ) = mUzeta( 4, aIndex ) * mNablaZeta[ 1 ] + mVzeta( 4, aIndex ) * mNablaEta[ 1 ] + mWzeta( 4, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 4 ) = mUzeta( 4, aIndex ) * mNablaZeta[ 2 ] + mVzeta( 4, aIndex ) * mNablaEta[ 2 ] + mWzeta( 4, aIndex ) * mNablaTau[ 2 ];

            // tau = -eta-zeta
            mFzeta( 0, 5 ) = -mFzeta( 0, 3 ) - mFzeta( 0, 4 );
            mFzeta( 1, 5 ) = -mFzeta( 1, 3 ) - mFzeta( 1, 4 );
            mFzeta( 2, 5 ) = -mFzeta( 2, 3 ) - mFzeta( 2, 4 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 xi->tau->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFzeta( 0, 6 ) = mUzeta( 6, aIndex ) * mNablaXi[ 0 ] + mVzeta( 6, aIndex ) * mNablaTau[ 0 ] + mWzeta( 6, aIndex ) * mNablaEta[ 0 ] ;
            mFzeta( 1, 6 ) = mUzeta( 6, aIndex ) * mNablaXi[ 1 ] + mVzeta( 6, aIndex ) * mNablaTau[ 1 ] + mWzeta( 6, aIndex ) * mNablaEta[ 1 ] ;
            mFzeta( 2, 6 ) = mUzeta( 6, aIndex ) * mNablaXi[ 2 ] + mVzeta( 6, aIndex ) * mNablaTau[ 2 ] + mWzeta( 6, aIndex ) * mNablaEta[ 2 ] ;

            // tau
            mFzeta( 0, 7 ) = mUzeta( 7, aIndex ) * mNablaXi[ 0 ] + mVzeta( 7, aIndex ) * mNablaTau[ 0 ] + mWzeta( 7, aIndex ) * mNablaEta[ 0 ] ;
            mFzeta( 1, 7 ) = mUzeta( 7, aIndex ) * mNablaXi[ 1 ] + mVzeta( 7, aIndex ) * mNablaTau[ 1 ] + mWzeta( 7, aIndex ) * mNablaEta[ 1 ] ;
            mFzeta( 2, 7 ) = mUzeta( 7, aIndex ) * mNablaXi[ 2 ] + mVzeta( 7, aIndex ) * mNablaTau[ 2 ] + mWzeta( 7, aIndex ) * mNablaEta[ 2 ] ;

            // eta = -xi-tau
            mFzeta( 0, 8 ) = -mFzeta( 0, 6 ) - mFzeta( 0, 7 );
            mFzeta( 1, 8 ) = -mFzeta( 1, 6 ) - mFzeta( 1, 7 );
            mFzeta( 2, 8 ) = -mFzeta( 2, 6 ) - mFzeta( 2, 7 );

			// - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->eta->zeta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFzeta( 0, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 0 ] + mVzeta( 9, aIndex ) * mNablaEta[ 0 ] + mWzeta( 9, aIndex ) * mNablaZeta[ 0 ] ;
            mFzeta( 1, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 1 ] + mVzeta( 9, aIndex ) * mNablaEta[ 1 ] + mWzeta( 9, aIndex ) * mNablaZeta[ 1 ] ;
            mFzeta( 2, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 2 ] + mVzeta( 9, aIndex ) * mNablaEta[ 2 ] + mWzeta( 9, aIndex ) * mNablaZeta[ 2 ] ;

            // eta
            mFzeta( 0, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 0 ] + mVzeta( 10, aIndex ) * mNablaEta[ 0 ] + mWzeta( 10, aIndex ) * mNablaZeta[ 0 ] ;
            mFzeta( 1, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 1 ] + mVzeta( 10, aIndex ) * mNablaEta[ 1 ] + mWzeta( 10, aIndex ) * mNablaZeta[ 1 ] ;
            mFzeta( 2, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 2 ] + mVzeta( 10, aIndex ) * mNablaEta[ 2 ] + mWzeta( 10, aIndex ) * mNablaZeta[ 2 ] ;

            // zeta= = -xi-eta
            mFzeta( 0, 11 ) = -mFzeta( 0, 9 ) - mFzeta( 0, 10 );
            mFzeta( 1, 11 ) = -mFzeta( 1, 9 ) - mFzeta( 1, 10 );
            mFzeta( 2, 11 ) = -mFzeta( 2, 9 ) - mFzeta( 2, 10 );
        }


//------------------------------------------------------------------------------

        void
        EF_TET10::E_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );
            this->compute_edge_functions( aIndex );
            this->compute_face_functions( aIndex );
            this->combine_functions( mE, mF );
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::E_straight( const uint aIndex )
        {
            this->compute_edge_functions( aIndex );
            this->compute_face_functions( aIndex );
            this->combine_functions( mE, mF );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TET10::C( const uint aIndex )
        {
            // do the math
            ( this->*mFunDerivatives )( aIndex );

            // merge and consider orientation of faces
            this->combine_functions( mExi, mFxi );
            this->combine_functions( mEeta, mFeta );
            this->combine_functions( mEzeta, mFzeta );

            for ( uint k = 0; k < 20; ++k )
            {

                // Ez,y - Ey,z
                mC( 0, k ) =
                          mExi( 2, k ) * mNablaXi[ 1 ]       // Ez,xi*xi,y
                        + mEeta( 2, k ) * mNablaEta[ 1 ]     // Ez,eta*eta,y
                        + mEzeta( 2, k ) * mNablaZeta[ 1 ]   // Ez,zeta*zeta,y
                        - mExi( 1, k ) * mNablaXi[ 2 ]       // Ey,xi * xi,z
                        - mEeta( 1, k ) * mNablaEta[ 2 ]     // Ey,eta * eta,z
                        - mEzeta( 1, k ) * mNablaZeta[ 2 ] ; // Ey,eta * eta,z

                // Ex,z - Ez,x
                mC( 1, k ) =
                          mExi( 0, k ) * mNablaXi[ 2 ]   // Ex,xi*xi,z
                        + mEeta( 0, k ) * mNablaEta[ 2 ]    // Ex,eta*eta,z
                        + mEzeta( 0, k ) * mNablaZeta[ 2 ]    // Ex,zeta*zeta,z
                        - mExi( 2, k ) * mNablaXi[ 0 ]    // Ez,xi * xi,x
                        - mEeta( 2, k ) * mNablaEta[ 0 ]    // Ez,eta * eta,x
                        - mEzeta( 2, k ) * mNablaZeta[ 0 ] ; // Ez,eta * eta,x

                // Ey,x - Ex,y
                mC( 2, k ) =
                          mExi( 1, k ) * mNablaXi[ 0 ]    // Ey,xi*xi,x
                        + mEeta( 1, k ) * mNablaEta[ 0 ]    // Ey,eta*eta,x
                        + mEzeta( 1, k ) * mNablaZeta[ 0 ]    // Ey,zeta*zeta,x
                        - mExi( 0, k ) * mNablaXi[ 1 ]    // Ex,xi * xi,y
                        - mEeta( 0, k ) * mNablaEta[ 1 ] // Ex,eta * eta,y
                        - mEzeta( 0, k ) * mNablaZeta[ 1 ] ; // Ex,eta * eta,y
            }
            return mC ;
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::combine_functions( Matrix< real > & aE, const Matrix< real > & aF )
        {
            // write data into edge function
            uint tCount = 12 ;
            uint tOff = 0 ;

            // loop over all faces
            for( uint f=0; f<4; ++f )
            {
                // check orientation
                switch( mT[ f ] )
                {
                    case( 0 ) : // the element owns this face
                    {
                        aE( 0, tCount ) = aF( 0, tOff );
                        aE( 1, tCount ) = aF( 1, tOff );
                        aE( 2, tCount ) = aF( 2, tOff );
                        ++tCount ;

                        aE( 0, tCount ) = aF( 0, tOff+1 );
                        aE( 1, tCount ) = aF( 1, tOff+1 );
                        aE( 2, tCount ) = aF( 2, tOff+1 );
                        ++tCount ;

                        break ;
                    }
                    case( 1 ) : // node 0 of face is identical M0->S0, M1->S2
                    {
                        aE( 0, tCount ) = -aF( 0, tOff );
                        aE( 1, tCount ) = -aF( 1, tOff );
                        aE( 2, tCount ) = -aF( 2, tOff );
                        ++tCount ;

                        aE( 0, tCount ) = -aF( 0, tOff+2 );
                        aE( 1, tCount ) = -aF( 1, tOff+2 );
                        aE( 2, tCount ) = -aF( 2, tOff+2 );
                        ++tCount ;

                        break ;
                    }
                    case( 2 ) : // node 1 of face is identical, M0->S1, M1->S0
                    {
                        aE( 0, tCount ) = -aF( 0, tOff+1 );
                        aE( 1, tCount ) = -aF( 1, tOff+1 );
                        aE( 2, tCount ) = -aF( 2, tOff+1 );
                        ++tCount ;

                        aE( 0, tCount ) = -aF( 0, tOff );
                        aE( 1, tCount ) = -aF( 1, tOff );
                        aE( 2, tCount ) = -aF( 2, tOff );
                        ++tCount ;

                        break ;
                    }
                    case( 3 ) : // node 2 of face is identical, M0->S2, M1->S1
                    {
                        aE( 0, tCount ) = -aF( 0, tOff+2 );
                        aE( 1, tCount ) = -aF( 1, tOff+2 );
                        aE( 2, tCount ) = -aF( 2, tOff+2 );
                        ++tCount ;

                        aE( 0, tCount ) = -aF( 0, tOff+1 );
                        aE( 1, tCount ) = -aF( 1, tOff+1 );
                        aE( 2, tCount ) = -aF( 2, tOff+1 );
                        ++tCount ;

                        break ;
                    }
                    default : // this should never happen!
                    {
                        BELFEM_ERROR( false, "Invalid face orientation");
                    }
                }

                tOff += 3 ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::G_term1( const uint aIndex )
        {
            // reference partials with the nablas held fixed, in EXACTLY
            // C()'s sequence: the fork, then the three combines ( the
            // combines are idempotent — E_xi_* already combined — but the
            // sequence is bound to C()'s to stay immune to refactors; the
            // face partial columns are ZERO before combine ). The tables
            // are PRE-SIGNED and orientation-combined for this class.
            ( this->*mFunDerivatives )( aIndex );
            this->combine_functions( mExi,   mFxi );
            this->combine_functions( mEeta,  mFeta );
            this->combine_functions( mEzeta, mFzeta );

            for ( uint k = 0; k < 20; ++k )
            {
                for ( uint j = 0; j < 3; ++j )        // field component
                {
                    for ( uint i = 0; i < 3; ++i )    // derivative direction
                    {
                        mGrad( i + 3 * j, k ) =
                              mExi(   j, k ) * mNablaXi[ i ]
                            + mEeta(  j, k ) * mNablaEta[ i ]
                            + mEzeta( j, k ) * mNablaZeta[ i ];
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::G_straight( const uint aIndex )
        {
            // frozen nablas from link(); the affine map has S = 0, so the
            // hess channel does not exist here and MUST not be touched:
            // mCurv is NaN on this branch by design ( and mJ is never
            // written by this class at all )
            this->G_term1( aIndex );
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::G_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );
            this->G_term1( aIndex );

            // E at this point for the Q-folded hess channel ( refreshes
            // mE — the documented curved-path contract ); compute_nabla
            // above is cached, E_curved will not redo it
            this->E_curved( aIndex );

            // the three inverse-map Hessians H^(d) = -A * Cm_d * A^T,
            // dof-independent at this point; A columns are the nabla
            // vectors — NEVER the mJ / mInvJ members ( never written by
            // this class )
            const real * tA[ 3 ] = { mNablaXi, mNablaEta, mNablaZeta };
            real tH[ 3 ][ 9 ];
            for ( uint d = 0; d < 3; ++d )
            {
                // Cm from the IF-Voigt rows ( xx, yy, zz, yz, xz, xy )
                // contracted with nabla_d — the diagonal pairs are LIVE
                // for the quadratic map
                real tCm[ 3 ][ 3 ];
                real tS[ 6 ];
                for ( uint v = 0; v < 6; ++v )
                {
                    tS[ v ] = mCurv[ v ][ 0 ] * tA[ d ][ 0 ]
                            + mCurv[ v ][ 1 ] * tA[ d ][ 1 ]
                            + mCurv[ v ][ 2 ] * tA[ d ][ 2 ];
                }
                tCm[ 0 ][ 0 ] = tS[ 0 ];
                tCm[ 1 ][ 1 ] = tS[ 1 ];
                tCm[ 2 ][ 2 ] = tS[ 2 ];
                tCm[ 1 ][ 2 ] = tCm[ 2 ][ 1 ] = tS[ 3 ];
                tCm[ 0 ][ 2 ] = tCm[ 2 ][ 0 ] = tS[ 4 ];
                tCm[ 0 ][ 1 ] = tCm[ 1 ][ 0 ] = tS[ 5 ];

                for ( uint j = 0; j < 3; ++j )        // field component
                {
                    for ( uint i = 0; i < 3; ++i )    // derivative direction
                    {
                        real tSum = 0.0 ;
                        for ( uint p = 0; p < 3; ++p )
                        {
                            for ( uint q = 0; q < 3; ++q )
                            {
                                tSum += tA[ p ][ i ] * tCm[ p ][ q ]
                                      * tA[ q ][ j ];
                            }
                        }
                        tH[ d ][ i + 3 * j ] = -tSum ;
                    }
                }
            }

            // Q-folded channel: Q^d = sum_i mM[ 27 + 3*i + d ] * e_i —
            // compute_nabla's own Jacobian dump, fresh for aIndex
            for ( uint k = 0; k < 20; ++k )
            {
                real tQ[ 3 ];
                for ( uint d = 0; d < 3; ++d )
                {
                    tQ[ d ] = mM[ 27 + d ]     * mE( 0, k )
                            + mM[ 27 + 3 + d ] * mE( 1, k )
                            + mM[ 27 + 6 + d ] * mE( 2, k );
                }
                for ( uint r = 0; r < 9; ++r )
                {
                    mGrad( r, k ) += tQ[ 0 ] * tH[ 0 ][ r ]
                                   + tQ[ 1 ] * tH[ 1 ][ r ]
                                   + tQ[ 2 ] * tH[ 2 ][ r ];
                }
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TET10::G( const uint aIndex )
        {
            ( this->*mFunGrad )( aIndex );
            return mGrad ;
        }

    }
}
