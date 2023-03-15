//
// Created by christian on 12/1/21.
//

#include "cl_EF_TET10.hpp"
#include "assert.hpp"
#include "fn_dot.hpp"
#include "cl_Element.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_det.hpp"

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

            mF.set_size( 3, 12, 0.0 );

            mB.set_size( 3, 10 );
            mCurlA.set_size( 3, 30, 0.0 );

            mNodeCoords.set_size( 10, 3 );

            mExi.set_size( 3, 20, 0.0 );
            mEeta.set_size( 3, 20, 0.0 );
            mEzeta.set_size( 3, 20, 0.0 );

            mFxi.set_size( 3, 12, 0.0 );
            mFeta.set_size( 3, 12, 0.0 );
            mFzeta.set_size( 3, 12, 0.0 );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::link( Element * aElement,
                        const bool aOperatorCurlH,
                        const bool aOperatorGradPhi,
                        const bool aOperatorCurlA )
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
            mLastNablaDeriv = BELFEM_UINT_MAX ;

            // get edge directions
            if(  tElement->is_curved() )
            {
                if( aOperatorCurlH )
                {
                    aElement->edge_directions( mS );

                    // derivatives of jacobian matrix to xi
                    mM[  0 ] = 4. * ( mX( 0 ) + mX( 3 ) - mX( 7 ) - mX( 7 ) ); // a_xi
                    mM[  1 ] = 4. * ( mY( 0 ) + mY( 3 ) - mY( 7 ) - mY( 7 ) ); // b_xi
                    mM[  2 ] = 4. * ( mZ( 0 ) + mZ( 3 ) - mZ( 7 ) - mZ( 7 ) ); // c_xi
                    mM[  3 ] = 4. * ( mX( 3 ) + mX( 4 ) - mX( 7 ) - mX( 8 ) ); // d_xi
                    mM[  4 ] = 4. * ( mY( 3 ) + mY( 4 ) - mY( 7 ) - mY( 8 ) ); // e_xi
                    mM[  5 ] = 4. * ( mZ( 3 ) + mZ( 4 ) - mZ( 7 ) - mZ( 8 ) ); // f_xi
                    mM[  6 ] = 4. * ( mX( 0 ) + mX( 3 ) - mX( 7 ) - mX( 7 ) ); // g_xi
                    mM[  7 ] = 4. * ( mY( 0 ) + mY( 3 ) - mY( 7 ) - mY( 7 ) ); // h_xi
                    mM[  8 ] = 4. * ( mZ( 0 ) + mZ( 3 ) - mZ( 7 ) - mZ( 7 ) ); // i_xi

                    // derivatives of jacobian matrix to eta
                    mM[  9 ] = 4. * ( mX( 3 ) + mX( 4 ) - mX( 7 ) - mX( 8 ) ); // a_eta
                    mM[ 10 ] = 4. * ( mY( 3 ) + mY( 4 ) - mY( 7 ) - mY( 8 ) ); // b_eta
                    mM[ 11 ] = 4. * ( mZ( 3 ) + mZ( 4 ) - mZ( 7 ) - mZ( 8 ) ); // c_eta
                    mM[ 12 ] = 4. * ( mX( 1 ) + mX( 3 ) - mX( 8 ) - mX( 8 ) ); // d_eta
                    mM[ 13 ] = 4. * ( mY( 1 ) + mY( 3 ) - mY( 8 ) - mY( 8 ) ); // e_eta
                    mM[ 14 ] = 4. * ( mZ( 1 ) + mZ( 3 ) - mZ( 8 ) - mZ( 8 ) ); // f_eta
                    mM[ 15 ] = 4. * ( mX( 3 ) + mX( 5 ) - mX( 8 ) - mX( 9 ) ); // g_eta
                    mM[ 16 ] = 4. * ( mY( 3 ) + mY( 5 ) - mY( 8 ) - mY( 9 ) ); // h_eta
                    mM[ 17 ] = 4. * ( mZ( 3 ) + mZ( 5 ) - mZ( 8 ) - mZ( 9 ) ); // i_eta

                    // derivatives of jacobian matrix to zeta
                    mM[ 18 ] = 4. * ( mX( 0 ) + mX( 3 ) - mX( 7 ) - mX( 7 ) ); // a_zeta
                    mM[ 19 ] = 4. * ( mY( 0 ) + mY( 3 ) - mY( 7 ) - mY( 7 ) ); // b_zeta
                    mM[ 20 ] = 4. * ( mZ( 0 ) + mZ( 3 ) - mZ( 7 ) - mZ( 7 ) ); // c_zeta
                    mM[ 21 ] = 4. * ( mX( 3 ) + mX( 5 ) - mX( 8 ) - mX( 9 ) ); // d_zeta
                    mM[ 22 ] = 4. * ( mY( 3 ) + mY( 5 ) - mY( 8 ) - mY( 9 ) ); // e_zeta
                    mM[ 23 ] = 4. * ( mZ( 3 ) + mZ( 5 ) - mZ( 8 ) - mZ( 9 ) ); // f_zeta
                    mM[ 24 ] = 4. * ( mX( 2 ) + mX( 3 ) - mX( 9 ) - mX( 9 ) ); // g_zeta
                    mM[ 25 ] = 4. * ( mY( 2 ) + mY( 3 ) - mY( 9 ) - mY( 9 ) ); // h_zeta
                    mM[ 26 ] = 4. * ( mZ( 2 ) + mZ( 3 ) - mZ( 9 ) - mZ( 9 ) ); // i_zeta

                    // link functions
                    mFunInterpolation = & EF_TET10::E_curved ;
                    mFunDerivatives   = & EF_TET10::E_xi_curved ;
                }
                else
                {
                    mFunInterpolation = nullptr ;
                    mFunDerivatives   = nullptr ;
                }

                if( aOperatorGradPhi || aOperatorCurlA )
                {
                    // get the node coordinates
                    mNodeCoords.set_col( 0, mX );
                    mNodeCoords.set_col( 1, mY );
                    mNodeCoords.set_col( 2, mZ );

                    // link curved grad function
                    mFunGrad = & EF_TET10::B_curved ;
                }
                else
                {
                    mFunGrad = nullptr ;
                    mJ.fill( BELFEM_QUIET_NAN );
                    mInvJ.fill( BELFEM_QUIET_NAN );
                    mB.fill( BELFEM_QUIET_NAN );

                }
            }
            else // we do if this element has straight edges
            {
                if( aOperatorCurlH )
                {
                    // get the directions of the edges
                    aElement->edge_directions( mS );

                    this->compute_nabla( 0 );

                    mFunInterpolation  = & EF_TET10::E_straight ;
                    mFunDerivatives    = & EF_TET10::E_xi_straight ;

                }
                else
                {
                    mFunInterpolation = nullptr ;
                    mFunDerivatives   = nullptr ;
                }

                if( aOperatorGradPhi || aOperatorCurlA )
                {
                    // get the node coordinates
                    mNodeCoords.set_col( 0, mX );
                    mNodeCoords.set_col( 1, mY );
                    mNodeCoords.set_col( 2, mZ );

                    // compute the jacobian
                    this->compute_jacobian( 0 );

                    // link straight grad function
                    mFunGrad = & EF_TET10::B_straight ;
                }
                else
                {
                    mFunGrad = nullptr ;
                    mJ.fill( BELFEM_QUIET_NAN );
                    mInvJ.fill( BELFEM_QUIET_NAN );
                    mB.fill( BELFEM_QUIET_NAN );
                }
            }

            if( aOperatorCurlH )
            {
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

                real  xi16 = xi8 + xi8 ;
                real  eta16 = eta8 + eta8 ;
                real zeta16 = zeta8 + zeta8 ;

                mNxi( 0 , k ) = xi4-1.;
                mNxi( 1 , k ) = 0.;
                mNxi( 2 , k ) = 0.;
                mNxi( 3 , k ) = xi4+eta4+zeta4-3.;
                mNxi( 4 , k ) = eta4;
                mNxi( 5 , k ) = 0.;
                mNxi( 6 , k ) = zeta4;
                mNxi( 7 , k ) = 4.-xi4-xi4-eta4-zeta4;
                mNxi( 8 , k ) = -eta4;
                mNxi( 9 , k ) = -zeta4;

                mNeta( 0 , k ) = 0.0;
                mNeta( 1 , k ) = eta4-1.;
                mNeta( 2 , k ) = 0.;
                mNeta( 3 , k ) = xi4+eta4+zeta4-3.;
                mNeta( 4 , k ) = xi4;
                mNeta( 5 , k ) = zeta4;
                mNeta( 6 , k ) = 0.;
                mNeta( 7 , k ) = -xi4;
                mNeta( 8 , k ) = 4.-xi4-eta4-eta4-zeta4;
                mNeta( 9 , k ) = -zeta4;

                mNzeta( 0 , k ) = 0.;
                mNzeta( 1 , k ) = 0.;
                mNzeta( 2 , k ) = zeta4-1.;
                mNzeta( 3 , k ) = xi4+eta4+zeta4-3.;
                mNzeta( 4 , k ) = 0.;
                mNzeta( 5 , k ) = eta4;
                mNzeta( 6 , k ) = xi4;
                mNzeta( 7 , k ) = -xi4;
                mNzeta( 8 , k ) = -eta4;
                mNzeta( 9 , k ) = 4.-xi4-eta4-zeta4-zeta4;

                mG(  0, k )	= xi4*(xi2 - 1.);
                mG(  1, k )	= xi2*(eta4 - 1.);
                mG(  2, k )	= eta4*(eta2 - 1.);
                mG(  3, k )	= eta2*(zeta4 - 1.);
                mG(  4, k )	= zeta4*(zeta2 - 1.);
                mG(  5, k )	= zeta2*(xi4 - 1.);
                mG(  6, k )	= xi4*(xi2 - 1.);
                mG(  7, k )	= xi2*(tau4 - 1.);
                mG(  8, k )	= eta4*(eta2 - 1.);
                mG(  9, k )	= eta2*(tau4 - 1.);
                mG( 10, k )	= zeta4*(zeta2 - 1.);
                mG( 11, k )	= zeta2*(tau4 - 1.);

                mH(  0, k )	= eta2*(1.- xi4 );
                mH(  1, k )	= eta4*(1.- eta2 );
                mH(  2, k )	= zeta2*(1.- eta4 );
                mH(  3, k )	= zeta4*(1.- zeta2 );
                mH(  4, k )	= xi2*(1.- zeta4 );
                mH(  5, k )	= xi4*(1.- xi2 );
                mH(  6, k )	= tau2*(1.- xi4 );
                mH(  7, k )	= tau4*(1.- tau2 );
                mH(  8, k )	= tau2*(1.- eta4 );
                mH(  9, k )	= tau4*(1.- tau2 );
                mH( 10, k )	= tau2*(1.- zeta4 );
                mH( 11, k )	= tau4*(1.- tau2 );

                mGxi( 0, k ) = xi16-4.0;
                mGxi( 1, k ) = eta8-2.0;
                mGxi( 5, k ) = zeta8;
                mGxi( 6, k ) = xi16-4.0;
                mGxi( 7, k ) = -eta8-xi16-zeta8+6.0;
                mGxi( 9, k ) = -eta8;
                mGxi( 11, k ) = -zeta8;

                mGeta( 1, k ) = xi8;
                mGeta( 2, k ) = eta16-4.0;
                mGeta( 3, k ) = zeta8-2.0;
                mGeta( 7, k ) = -xi8;
                mGeta( 8, k ) = eta16-4.0;
                mGeta( 9, k ) = -eta16-xi8-zeta8+6.0;
                mGeta( 11, k ) = -zeta8;

                mGzeta( 3, k ) = eta8;
                mGzeta( 4, k ) = zeta16-4.0;
                mGzeta( 5, k ) = xi8-2.0;
                mGzeta( 7, k ) = -xi8;
                mGzeta( 9, k ) = -eta8;
                mGzeta( 10, k ) = zeta16-4.0;
                mGzeta( 11, k ) = -eta8-xi8-zeta16+6.0;

                mHxi( 0, k ) = -eta8;
                mHxi( 4, k ) = -zeta8+2.0;
                mHxi( 5, k ) = -xi16+4.0;
                mHxi( 6, k ) = eta8+xi16+zeta8-10.;
                mHxi( 7, k ) = -eta16-xi16-zeta16+12.;
                mHxi( 8, k ) = eta8-2.0;
                mHxi( 9, k ) = -eta16-xi16-zeta16+12.;
                mHxi( 10, k ) = zeta8-2.0;
                mHxi( 11, k ) = 12.-eta16-xi16-zeta16;

                mHeta( 0, k ) = 2.0 -xi8 ;
                mHeta( 1, k ) = 4.0 - eta16 ;
                mHeta( 2, k ) = -zeta8 ;
                mHeta( 6, k ) = xi8 - 2.0 ;
                mHeta( 7, k ) = 12.-eta16-xi16-zeta16;
                mHeta( 8, k ) = xi8 + eta16 + zeta8 - 10. ;
                mHeta( 9, k ) =  12.-eta16-xi16-zeta16;
                mHeta( 10, k ) = zeta8 - 2. ;
                mHeta( 11, k ) =  12.-eta16-xi16-zeta16;


                mHzeta( 2, k ) = -eta8+2.0;
                mHzeta( 3, k ) = -zeta16+4.0;
                mHzeta( 4, k ) = -xi8;
                mHzeta( 6, k ) = xi8-2.0;
                mHzeta( 7, k ) = -eta16-xi16-zeta16+12.;
                mHzeta( 8, k ) = eta8-2.0;
                mHzeta( 9, k ) = -eta16-xi16-zeta16+12.;
                mHzeta( 10, k ) = eta8+xi8+zeta16-10.;
                mHzeta( 11, k ) = -eta16-xi16-zeta16+12.;

                // face coeffs
                real tA = xi*eta8 ;
                real tB = eta*zeta8 ;
                real tC = zeta*tau8 ;
                real tD = eta*tau8 ;
                real tE = tau*xi8 ;
                real tF = zeta*xi8 ;

                mU(  0, k ) = tD+tD;
                mU(  1, k ) = -tD;
                //mU(  2, k ) = -tD;
                mU(  3, k ) = tC+tC;
                mU(  4, k ) = -tC;
                //mU(  5, k ) = -tC;
                mU(  6, k ) = tE+tE;
                mU(  7, k ) = -tE;
               // mU(  8, k ) = -tE;
                mU(  9, k ) = tB+tB;
                mU( 10, k ) = -tB;
                //mU( 11, k ) = -tB;

                mV(  0, k ) = -tE;
                mV(  1, k ) = tE+tE;
                //mV(  2, k ) = -tE;
                mV(  3, k ) = -tD;
                mV(  4, k ) = tD+tD ;
                //mV(  5, k ) = -tD;
                mV(  6, k ) = -tC;
                mV(  7, k ) = tC+tC;
                //mV(  8, k ) = -tC;
                mV(  9, k ) = -tA;
                mV( 10, k ) = tA+tA;
                //mV( 11, k ) = -tA;

                mW(  0, k ) = -tA;
                mW(  1, k ) = -tA;
                //mW(  2, k ) = tA+tA;
                mW(  3, k ) = -tB;
                mW(  4, k ) = -tB;
                //mW(  5, k ) = tB+tB;
                mW(  6, k ) = -tF;
                mW(  7, k ) = -tF;
                //mW(  8, k ) = tF+tF;
                mW(  9, k ) = -tF;
                mW( 10, k ) = -tF;
                //mW( 11, k ) = tF+tF;

                // derivatives to xi
                tA = eta8;
                tB = 0.0;
                tC = -zeta8;
                tD = -eta8;
                tE = -eta8-xi16-zeta8+8.0;
                tF = zeta8;
                mUxi(  0, k ) = tD+tD;
                mUxi(  1, k ) = -tD;
                //mUxi(  2, k ) = -tD;
                mUxi(  3, k ) = tC+tC;
                mUxi(  4, k ) = -tC;
                //mUxi(  5, k ) = -tC;
                mUxi(  6, k ) = tE+tE;
                mUxi(  7, k ) = -tE;
                //mUxi(  8, k ) = -tE;
                mUxi(  9, k ) = tB+tB;
                mUxi( 10, k ) = -tB;
                //mUxi( 11, k ) = -tB;

                mVxi(  0, k ) = -tE;
                mVxi(  1, k ) = tE+tE;
                //mVxi(  2, k ) = -tE;
                mVxi(  3, k ) = -tD;
                mVxi(  4, k ) = tD+tD;
                //mVxi(  5, k ) = -tD;
                mVxi(  6, k ) = -tC;
                mVxi(  7, k ) = tC+tC;
                //mVxi(  8, k ) = -tC;
                mVxi(  9, k ) = -tA;
                mVxi( 10, k ) = tA+tA;
                //mVxi( 11, k ) = -tA;

                mWxi(  0, k ) = -tA;
                mWxi(  1, k ) = -tA;
                //mWxi(  2, k ) = tA+tA;
                mWxi(  3, k ) = -tB;
                mWxi(  4, k ) = -tB;
                //mWxi(  5, k ) = tB+tB;
                mWxi(  6, k ) = -tF;
                mWxi(  7, k ) = -tF;
                //mWxi(  8, k ) = tF+tF;
                mWxi(  9, k ) = -tF;
                mWxi( 10, k ) = -tF;
                //mWxi( 11, k ) = tF+tF;

                // derivatives to eta
                tA = xi8;
                tB = zeta8;
                tC = -zeta8;
                tD = -eta16-xi8-zeta8+8.0;
                tE = -xi8;
                tF = 0.0 ;
                mUeta(  0, k ) = tD+tD;
                mUeta(  1, k ) = -tD;
                //mUeta(  2, k ) = -tD;
                mUeta(  3, k ) = tC+tC;
                mUeta(  4, k ) = -tC;
                //mUeta(  5, k ) = -tC;
                mUeta(  6, k ) = tE+tE;
                mUeta(  7, k ) = -tE;
                //mUeta(  8, k ) = -tE;
                mUeta(  9, k ) = tB+tB;
                mUeta( 10, k ) = -tB;
                //mUeta( 11, k ) = -tB;

                mVeta(  0, k ) = -tE;
                mVeta(  1, k ) = tE+tE;
                //mVeta(  2, k ) = -tE;
                mVeta(  3, k ) = -tD;
                mVeta(  4, k ) = tD+tD;
                //mVeta(  5, k ) = -tD;
                mVeta(  6, k ) = -tC;
                mVeta(  7, k ) = tC+tC;
                //mVeta(  8, k ) = -tC;
                mVeta(  9, k ) = -tA;
                mVeta( 10, k ) = tA+tA;
                //mVeta( 11, k ) = -tA;

                mWeta(  0, k ) = -tA;
                mWeta(  1, k ) = -tA;
                //mWeta(  2, k ) = tA+tA;
                mWeta(  3, k ) = -tB;
                mWeta(  4, k ) = -tB;
                //mWeta(  5, k ) = tB+tB;
                mWeta(  6, k ) = -tF;
                mWeta(  7, k ) = -tF;
                //mWeta(  8, k ) = tF+tF;
                mWeta(  9, k ) = -tF;
                mWeta( 10, k ) = -tF;
                //mWeta( 11, k ) = tF+tF;

                // derivatives to zeta
                tA = 0.0;
                tB = eta8;
                tC = -eta8-xi8-zeta16+8.0;
                tD = -eta8;
                tE = -xi8;
                tF = xi8;
                mUzeta(  0, k ) = tD+tD;
                mUzeta(  1, k ) = -tD;
                //mUzeta(  2, k ) = -tD;
                mUzeta(  3, k ) = tC+tC;
                mUzeta(  4, k ) = -tC;
                //mUzeta(  5, k ) = -tC;
                mUzeta(  6, k ) = tE+tE;
                mUzeta(  7, k ) = -tE;
                //mUzeta(  8, k ) = -tE;
                mUzeta(  9, k ) = tB+tB;
                mUzeta( 10, k ) = -tB;
                //mUzeta( 11, k ) = -tB;

                mVzeta(  0, k ) = -tE;
                mVzeta(  1, k ) = tE+tE;
                //mVzeta(  2, k ) = -tE;
                mVzeta(  3, k ) = -tD;
                mVzeta(  4, k ) = tD+tD;
                //mVzeta(  5, k ) = -tD;
                mVzeta(  6, k ) = -tC;
                mVzeta(  7, k ) = tC+tC;
                //mVzeta(  8, k ) = -tC;
                mVzeta(  9, k ) = -tA;
                mVzeta( 10, k ) = tA+tA;
                //mVzeta( 11, k ) = -tA;

                mWzeta(  0, k ) = -tA;
                mWzeta(  1, k ) = -tA;
                //mWzeta(  2, k ) = tA+tA;
                mWzeta(  3, k ) = -tB;
                mWzeta(  4, k ) = -tB;
                //mWzeta(  5, k ) = tB+tB;
                mWzeta(  6, k ) = -tF;
                mWzeta(  7, k ) = -tF;
                //mWzeta(  8, k ) = tF+tF;
                mWzeta(  9, k ) = -tF;
                mWzeta( 10, k ) = -tF;
                //mWzeta( 11, k ) = tF+tF;

            }
        }

 //------------------------------------------------------------------------------

        void
        EF_TET10::E_xi_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );
            this->compute_nabla_derivatives( aIndex );

            this->compute_edge_derivatives_curved( aIndex );
            this->compute_face_derivatives_curved( aIndex );

            this->combine_functions( mExi, mFxi );
            this->combine_functions( mEeta, mFeta );
            this->combine_functions( mEzeta, mFzeta );
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::E_xi_straight( const uint aIndex )
        {
            this->compute_edge_derivatives_straight( aIndex );
            this->compute_face_derivatives_straight( aIndex );

            this->combine_functions( mExi, mFxi );
            this->combine_functions( mEeta, mFeta );
            this->combine_functions( mEzeta, mFzeta );
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_jacobian( const uint aIndex )
        {
            // check if Jacobian is updated
            if( aIndex != mLastJ )
            {
                mJ = mGroup->dNdXi( aIndex ) * mNodeCoords ;

                mDetJ =   mJ( 0, 2 )*( mJ( 1, 0 )*mJ( 2, 1 ) - mJ( 1, 1 )*mJ( 2, 0 ) )
                        - mJ( 0, 1 )*( mJ( 1, 0 )*mJ( 2, 2 ) - mJ( 1, 2 )*mJ( 2, 0 ) )
                        + mJ( 0, 0 )*( mJ( 1, 1 )*mJ( 2, 2 ) - mJ( 1, 2 )*mJ( 2, 1 ) );

                mAbsDetJ = std::abs( mDetJ ) ;

                mInvJ( 0,0 ) = mJ( 1, 1 ) * mJ( 2, 2 ) - mJ( 1, 2 ) * mJ( 2, 1 );
                mInvJ( 1,0 ) = mJ( 1, 2 ) * mJ( 2, 0 ) - mJ( 1, 0 ) * mJ( 2, 2 );
                mInvJ( 2,0 ) = mJ( 1, 0 ) * mJ( 2, 1 ) - mJ( 1, 1 ) * mJ( 2, 0 );

                mInvJ( 0,1 ) = mJ( 0, 2 ) * mJ( 2, 1 ) - mJ( 0, 1 ) * mJ( 2, 2 );
                mInvJ( 1,1 ) = mJ( 0, 0 ) * mJ( 2, 2 ) - mJ( 0, 2 ) * mJ( 2, 0 );
                mInvJ( 2,1 ) = mJ( 0, 1 ) * mJ( 2, 0 ) - mJ( 0, 0 ) * mJ( 2, 1 );

                mInvJ( 0,2 ) = mJ( 0, 1 ) * mJ( 1, 2 ) - mJ( 0, 2 ) * mJ( 1, 1 );
                mInvJ( 1,2 ) = mJ( 0, 2 ) * mJ( 1, 0 ) - mJ( 0, 0 ) * mJ( 1, 2 );
                mInvJ( 2,2 ) = mJ( 0, 0 ) * mJ( 1, 1 ) - mJ( 0, 1 ) * mJ( 1, 0 );

                mInvJ /= mDetJ ;

                mLastJ = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_nabla( const uint aIndex )
        {
            if( aIndex != mLastNabla )
            {
                // J = [ a b c ; d e f ; g h i ] ;
                mM[ 27 ] = dot( mNxi.col( aIndex ), mX );   // a
                mM[ 28 ] = dot( mNxi.col( aIndex ), mY );   // b
                mM[ 29 ] = dot( mNxi.col( aIndex ), mZ );   // c
                mM[ 30 ] = dot( mNeta.col( aIndex ), mX );  // d
                mM[ 31 ] = dot( mNeta.col( aIndex ), mY );  // e
                mM[ 32 ] = dot( mNeta.col( aIndex ), mZ );  // f
                mM[ 33 ] = dot( mNzeta.col( aIndex ), mX ); // g
                mM[ 34 ] = dot( mNzeta.col( aIndex ), mY ); // h
                mM[ 35 ] = dot( mNzeta.col( aIndex ), mZ ); // i


                // t = a * ( e * i - h * f ) + b * ( f * g - d * i ) + c * ( d * h - e * g );
                mDetJ =    mM[ 29 ]*( mM[ 30 ]*mM[ 34 ] - mM[ 31 ]*mM[ 33 ] )
                         - mM[ 28 ]*( mM[ 30 ]*mM[ 35 ] - mM[ 32 ]*mM[ 33 ] )
                         + mM[ 27 ]*( mM[ 31 ]*mM[ 35 ] - mM[ 32 ]*mM[ 34 ] ) ;

                mAbsDetJ = std::abs( mDetJ ) ;

                mM[ 36 ] = 1.0 / ( mDetJ * mDetJ );

                // temporary variable, will overwrite later
                mM[ 37 ] = 1.0 / mDetJ ;

                // Nablas
                mNablaXi[ 0 ]   = ( mM[ 31 ]*mM[ 35 ] - mM[ 32 ]*mM[ 34 ] ) * mM[ 37 ] ;
                mNablaXi[ 1 ]   = ( mM[ 32 ]*mM[ 33 ] - mM[ 30 ]*mM[ 35 ] ) * mM[ 37 ] ;
                mNablaXi[ 2 ]   = ( mM[ 30 ]*mM[ 34 ] - mM[ 31 ]*mM[ 33 ] ) * mM[ 37 ] ;

                mNablaEta[ 0 ]  = ( mM[ 29 ]*mM[ 34 ] - mM[ 28 ]*mM[ 35 ] ) * mM[ 37 ] ;
                mNablaEta[ 1 ]  = ( mM[ 27 ]*mM[ 35 ] - mM[ 29 ]*mM[ 33 ] ) * mM[ 37 ] ;
                mNablaEta[ 2 ]  = ( mM[ 28 ]*mM[ 33 ] - mM[ 27 ]*mM[ 34 ] ) * mM[ 37 ] ;

                mNablaZeta[ 0 ] = ( mM[ 28 ]*mM[ 32 ] - mM[ 29 ]*mM[ 31 ] ) * mM[ 37 ] ;
                mNablaZeta[ 1 ] = ( mM[ 29 ]*mM[ 30 ] - mM[ 27 ]*mM[ 32 ] ) * mM[ 37 ] ;
                mNablaZeta[ 2 ] = ( mM[ 27 ]*mM[ 31 ] - mM[ 28 ]*mM[ 30 ] ) * mM[ 37 ] ;

                mNablaTau[ 0 ] = -( mNablaXi[ 0 ] + mNablaEta[ 0 ] + mNablaZeta[ 0 ] );
                mNablaTau[ 1 ] = -( mNablaXi[ 1 ] + mNablaEta[ 1 ] + mNablaZeta[ 1 ] );
                mNablaTau[ 2 ] = -( mNablaXi[ 2 ] + mNablaEta[ 2 ] + mNablaZeta[ 2 ] );

                mLastNabla = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_nabla_derivatives( const uint aIndex )
        {
            if( aIndex != mLastNablaDeriv )
            {

                /* derivative of det(J) to xi
                t_xi = a * ( e_xi * i + e * i_xi - h_xi * f )
                        + b * ( f_xi * g + f * g_xi - i_xi * d )
                        + c * ( d_xi * h + d * h_xi - g_xi * e )
                        + g * ( b_xi * f - e_xi * c - e * c_xi )
                        + h * ( c_xi * d - f_xi * a - f * a_xi )
                        + i * ( a_xi * e - d_xi * b - d * b_xi ); */
                mM[ 37 ] = mM[ 29 ] * ( mM[  3 ] * mM[ 34 ] - mM[  6 ] * mM[ 31 ] + mM[  7 ] * mM[ 30 ] )
                         - mM[ 34 ] * ( mM[  0 ] * mM[ 32 ] - mM[  2 ] * mM[ 30 ] + mM[  5 ] * mM[ 27 ] )
                         - mM[ 33 ] * ( mM[  2 ] * mM[ 31 ] - mM[  1 ] * mM[ 32 ] + mM[  4 ] * mM[ 29 ] )
                         - mM[ 35 ] * ( mM[  1 ] * mM[ 30 ] - mM[  0 ] * mM[ 31 ] + mM[  3 ] * mM[ 28 ] )
                         + mM[ 28 ] * ( mM[  5 ] * mM[ 33 ] + mM[  6 ] * mM[ 32 ] - mM[  8 ] * mM[ 30 ] )
                         + mM[ 27 ] * ( mM[  4 ] * mM[ 35 ] - mM[  7 ] * mM[ 32 ] + mM[  8 ] * mM[ 31 ] );



                /* derivative of det(J) to eta
                t_eta = a * ( e_eta * i + e * i_eta - h_eta * f )
                        + b * ( f_eta * g + f * g_eta - i_eta * d )
                        + c * ( d_eta * h + d * h_eta - g_eta * e )
                        + g * ( b_eta * f - e_eta * c - e * c_eta )
                        + h * ( c_eta * d - f_eta * a - f * a_eta )
                        + i * ( a_eta * e - d_eta * b - d * b_eta ); */
                mM[ 38 ] = mM[ 29 ] * ( mM[ 12 ] * mM[ 34 ] - mM[ 15 ] * mM[ 31 ] + mM[ 16 ] * mM[ 30 ] )
                         - mM[ 34 ] * ( mM[  9 ] * mM[ 32 ] - mM[ 11 ] * mM[ 30 ] + mM[ 14 ] * mM[ 27 ] )
                         - mM[ 33 ] * ( mM[ 11 ] * mM[ 31 ] - mM[ 10 ] * mM[ 32 ] + mM[ 13 ] * mM[ 29 ] )
                         - mM[ 35 ] * ( mM[ 10 ] * mM[ 30 ] - mM[  9 ] * mM[ 31 ] + mM[ 12 ] * mM[ 28 ] )
                         + mM[ 28 ] * ( mM[ 14 ] * mM[ 33 ] + mM[ 15 ] * mM[ 32 ] - mM[ 17 ] * mM[ 30 ] )
                         + mM[ 27 ] * ( mM[ 13 ] * mM[ 35 ] - mM[ 16 ] * mM[ 32 ] + mM[ 17 ] * mM[ 31 ] ) ;

                /* derivative of det(J) to zeta
                t_zeta = a * ( e_zeta * i + e * i_zeta - h_zeta * f )
                        + b * ( f_zeta * g + f * g_zeta - i_zeta * d )
                        + c * ( d_zeta * h + d * h_zeta - g_zeta * e )
                        + g * ( b_zeta * f - e_zeta * c - e * c_zeta )
                        + h * ( c_zeta * d - f_zeta * a - f * a_zeta )
                        + i * ( a_zeta * e - d_zeta * b - d * b_zeta ); */
                mM[ 39 ] = mM[ 29 ] * ( mM[ 21 ] * mM[ 34 ] - mM[ 24 ] * mM[ 31 ] + mM[ 25 ] * mM[ 30 ] )
                         - mM[ 34 ] * ( mM[ 18 ] * mM[ 32 ] - mM[ 20 ] * mM[ 30 ] + mM[ 23 ] * mM[ 27 ] )
                         - mM[ 33 ] * ( mM[ 20 ] * mM[ 31 ] - mM[ 19 ] * mM[ 32 ] + mM[ 22 ] * mM[ 29 ] )
                         - mM[ 35 ] * ( mM[ 19 ] * mM[ 30 ] - mM[ 18 ] * mM[ 31 ] + mM[ 21 ] * mM[ 28 ] )
                         + mM[ 28 ] * ( mM[ 23 ] * mM[ 33 ] + mM[ 24 ] * mM[ 32 ] - mM[ 26 ] * mM[ 30 ] )
                         + mM[ 27 ] * ( mM[ 22 ] * mM[ 35 ] - mM[ 25 ] * mM[ 32 ] + mM[ 26 ] * mM[ 31 ] ) ;

                // componens [ l m n ; p q r ; u v w ]  = inv(J) * det(J)

                // l = e * i - f * h;
                mM[ 40 ] = mM[ 31 ] * mM[ 35 ] - mM[ 32 ] * mM[ 34 ];

                // m = c * h - b * i;
                mM[ 41 ] = mM[ 29 ] * mM[ 34 ] - mM[ 28 ] * mM[ 35 ];

                // n = b * f - c * e;
                mM[ 42 ] = mM[ 28 ] * mM[ 32 ] - mM[ 29 ] * mM[ 31 ];

                // p = f * g - d * i;
                mM[ 43 ] = mM[ 32 ] * mM[ 33 ] - mM[ 30 ] * mM[ 35 ];

                // q = a * i - c * g;
                mM[ 44 ] = mM[ 27 ] * mM[ 35 ] - mM[ 29 ] * mM[ 33 ];

                // r = c * d - a * f;
                mM[ 45 ] = mM[ 29 ] * mM[ 30 ] - mM[ 27 ] * mM[ 32 ];

                // u = d * h - e * g;
                mM[ 46 ] = mM[ 30 ] * mM[ 34 ] - mM[ 31 ] * mM[ 33 ];

                // v = b * g - a * h;
                mM[ 47 ] = mM[ 28 ] * mM[ 33 ] - mM[ 27 ] * mM[ 34 ];

                // w = a * e - b * d;
                mM[ 48 ] = mM[ 27 ] * mM[ 31 ] - mM[ 28 ] * mM[ 30 ];

                // derivatives to xi

                // l_xi = e_xi * i + e * i_xi - f_xi * h - f * h_xi;
                mM[ 49 ] = mM[  4 ] * mM[ 35 ] - mM[  5 ] * mM[ 34 ] - mM[  7 ] * mM[ 32 ] + mM[  8 ] * mM[ 31 ];

                // m_xi = c_xi * h + c * h_xi - b_xi * i - b * i_xi;
                mM[ 50 ] = mM[  2 ] * mM[ 34 ] - mM[  1 ] * mM[ 35 ] + mM[  7 ] * mM[ 29 ] - mM[  8 ] * mM[ 28 ];

                // n_xi = b_xi * f + b * f_xi - c_xi * e - c * e_xi;
                mM[ 51 ] = mM[  1 ] * mM[ 32 ] - mM[  2 ] * mM[ 31 ] - mM[  4 ] * mM[ 29 ] + mM[  5 ] * mM[ 28 ];

                // p_xi = f_xi * g + f * g_xi - d_xi * i - d * i_xi;
                mM[ 52 ] = mM[  5 ] * mM[ 33 ] - mM[  3 ] * mM[ 35 ] + mM[  6 ] * mM[ 32 ] - mM[  8 ] * mM[ 30 ];

                // q_xi = a_xi * i + a * i_xi - c_xi * g - c * g_xi;
                mM[ 53 ] = mM[  0 ] * mM[ 35 ] - mM[  2 ] * mM[ 33 ] - mM[  6 ] * mM[ 29 ] + mM[  8 ] * mM[ 27 ];

                // r_xi = c_xi * d + c * d_xi - a_xi * f - a * f_xi;
                mM[ 54 ] = mM[  2 ] * mM[ 30 ] - mM[  0 ] * mM[ 32 ] + mM[  3 ] * mM[ 29 ] - mM[  5 ] * mM[ 27 ];

                // u_xi = d_xi * h + d * h_xi - e_xi * g - e * g_xi;
                mM[ 55 ] = mM[  3 ] * mM[ 34 ] - mM[  4 ] * mM[ 33 ] - mM[  6 ] * mM[ 31 ] + mM[  7 ] * mM[ 30 ];

                // v_xi = b_xi * g + b * g_xi - a_xi * h - a * h_xi;
                mM[ 56 ] = mM[  1 ] * mM[ 33 ] - mM[  0 ] * mM[ 34 ] + mM[  6 ] * mM[ 28 ] - mM[  7 ] * mM[ 27 ];

                // l_eta = e_eta * i + e * i_eta - f_eta * h - f * h_eta;
                mM[ 57 ] = mM[  0 ] * mM[ 31 ] - mM[  1 ] * mM[ 30 ] - mM[  3 ] * mM[ 28 ] + mM[  4 ] * mM[ 27 ];

                // derivatives to eta

                // l_eta = e_eta * i + e * i_eta - f_eta * h - f * h_eta;
                mM[ 58 ] = mM[ 13 ] * mM[ 35 ] - mM[ 14 ] * mM[ 34 ] - mM[ 16 ] * mM[ 32 ] + mM[ 17 ] * mM[ 31 ];

                // m_eta = c_eta * h + c * h_eta - b_eta * i - b * i_eta;
                mM[ 59 ] = mM[ 11 ] * mM[ 34 ] - mM[ 10 ] * mM[ 35 ] + mM[ 16 ] * mM[ 29 ] - mM[ 17 ] * mM[ 28 ];

                // n_eta = b_eta * f + b * f_eta - c_eta * e - c * e_eta;
                mM[ 60 ] = mM[ 10 ] * mM[ 32 ] - mM[ 11 ] * mM[ 31 ] - mM[ 13 ] * mM[ 29 ] + mM[ 14 ] * mM[ 28 ];

                // p_eta = f_eta * g + f * g_eta - d_eta * i - d * i_eta;
                mM[ 61 ] = mM[ 14 ] * mM[ 33 ] - mM[ 12 ] * mM[ 35 ] + mM[ 15 ] * mM[ 32 ] - mM[ 17 ] * mM[ 30 ];

                // q_eta = a_eta * i + a * i_eta - c_eta * g - c * g_eta;
                mM[ 62 ] = mM[  9 ] * mM[ 35 ] - mM[ 11 ] * mM[ 33 ] - mM[ 15 ] * mM[ 29 ] + mM[ 17 ] * mM[ 27 ];

                // r_eta = c_eta * d + c * d_eta - a_eta * f - a * f_eta;
                mM[ 63 ] = mM[ 11 ] * mM[ 30 ] - mM[  9 ] * mM[ 32 ] + mM[ 12 ] * mM[ 29 ] - mM[ 14 ] * mM[ 27 ];

                // u_eta = d_eta * h + d * h_eta - e_eta * g - e * g_eta;
                mM[ 64 ] = mM[ 12 ] * mM[ 34 ] - mM[ 13 ] * mM[ 33 ] - mM[ 15 ] * mM[ 31 ] + mM[ 16 ] * mM[ 30 ];



                // v_eta = b_eta * g + b * g_eta - a_eta * h - a * h_eta;
                mM[ 65 ] = mM[ 10 ] * mM[ 33 ] - mM[  9 ] * mM[ 34 ] + mM[ 15 ] * mM[ 28 ] - mM[ 16 ] * mM[ 27 ];

                // w_eta = a_eta * e + a * e_eta - b_eta * d - b * d_eta;
                mM[ 66 ] = mM[  9 ] * mM[ 31 ] - mM[ 10 ] * mM[ 30 ] - mM[ 12 ] * mM[ 28 ] + mM[ 13 ] * mM[ 27 ];

                // derivatives to zeta

                // l_zeta = e_zeta * i + e * i_zeta - f_zeta * h - f * h_zeta;
                mM[ 67 ] = mM[ 22 ] * mM[ 35 ] - mM[ 23 ] * mM[ 34 ] - mM[ 25 ] * mM[ 32 ] + mM[ 26 ] * mM[ 31 ];

                // m_zeta = c_zeta * h + c * h_zeta - b_zeta * i - b * i_zeta;
                mM[ 68 ] = mM[ 20 ] * mM[ 34 ] - mM[ 19 ] * mM[ 35 ] + mM[ 25 ] * mM[ 29 ] - mM[ 26 ] * mM[ 28 ];

                // n_zeta = b_zeta * f + b * f_zeta - c_zeta * e - c * e_zeta;
                mM[ 69 ] = mM[ 19 ] * mM[ 32 ] - mM[ 20 ] * mM[ 31 ] - mM[ 22 ] * mM[ 29 ] + mM[ 23 ] * mM[ 28 ];

                // p_zeta = f_zeta * g + f * g_zeta - d_zeta * i - d * i_zeta;
                mM[ 70 ] = mM[ 23 ] * mM[ 33 ] - mM[ 21 ] * mM[ 35 ] + mM[ 24 ] * mM[ 32 ] - mM[ 26 ] * mM[ 30 ];

                // q_zeta = a_zeta * i + a * i_zeta - c_zeta * g - c * g_zeta;
                mM[ 71 ] = mM[ 18 ] * mM[ 35 ] - mM[ 20 ] * mM[ 33 ] - mM[ 24 ] * mM[ 29 ] + mM[ 26 ] * mM[ 27 ];

                // r_zeta = c_zeta * d + c * d_zeta - a_zeta * f - a * f_zeta;
                mM[ 72 ] = mM[ 20 ] * mM[ 30 ] - mM[ 18 ] * mM[ 32 ] + mM[ 21 ] * mM[ 29 ] - mM[ 23 ] * mM[ 27 ];

                // u_zeta = d_zeta * h + d * h_zeta - e_zeta * g - e * g_zeta;
                mM[ 73 ] = mM[ 21 ] * mM[ 34 ] - mM[ 22 ] * mM[ 33 ] - mM[ 24 ] * mM[ 31 ] + mM[ 25 ] * mM[ 30 ];

                // v_zeta = b_zeta * g + b * g_zeta - a_zeta * h - a * h_zeta;
                mM[ 74 ] = mM[ 19 ] * mM[ 33 ] - mM[ 18 ] * mM[ 34 ] + mM[ 24 ] * mM[ 28 ] - mM[ 25 ] * mM[ 27 ];

                // w_zeta = a_zeta * e + a * e_zeta - b_zeta * d - b * d_zeta;
                mM[ 75 ] = mM[ 18 ] * mM[ 31 ] - mM[ 19 ] * mM[ 30 ] - mM[ 21 ] * mM[ 28 ] + mM[ 22 ] * mM[ 27 ];

                // t * l_xi - l * t_xi
                mdNablaXidXi[ 0 ] = ( mDetJ * mM[ 49 ] - mM[ 37 ] * mM[ 40 ] ) * mM[ 36 ];

                // t * p_xi - p * t_xi
                mdNablaXidXi[ 1 ] = ( mDetJ * mM[ 52 ] - mM[ 37 ] * mM[ 43 ] ) * mM[ 36 ];

                // t * u_xi - u * t_xi
                mdNablaXidXi[ 2 ] = ( mDetJ * mM[ 55 ] - mM[ 37 ] * mM[ 46 ] ) * mM[ 36 ];

                // t * l_eta - l * t_eta
                mdNablaXidEta[ 0 ] = ( mDetJ * mM[ 58 ] - mM[ 38 ] * mM[ 40 ] ) * mM[ 36 ];

                // t * p_eta - p * t_eta
                mdNablaXidEta[ 1 ] = ( mDetJ * mM[ 61 ] - mM[ 38 ] * mM[ 43 ] ) * mM[ 36 ];

                // t * u_eta - u * t_eta
                mdNablaXidEta[ 2 ] = ( mDetJ * mM[ 64 ] - mM[ 38 ] * mM[ 46 ] ) * mM[ 36 ];

                // t * l_zeta - l * t_zeta
                mdNablaXidZeta[ 0 ] = ( mDetJ * mM[ 67 ] - mM[ 39 ] * mM[ 40 ] ) * mM[ 36 ];

                // t * p_zeta - p * t_zeta
                mdNablaXidZeta[ 1 ] = ( mDetJ * mM[ 70 ] - mM[ 39 ] * mM[ 43 ] ) * mM[ 36 ];

                //  t * u_zeta - u * t_zeta
                mdNablaXidZeta[ 2 ] = ( mDetJ * mM[ 73 ] - mM[ 39 ] * mM[ 46 ] ) * mM[ 36 ];

                // t * m_xi - m * t_xi
                mdNablaEtadXi[ 0 ] = ( mDetJ * mM[ 50 ] - mM[ 37 ] * mM[ 41 ] ) * mM[ 36 ];

                // t * q_xi - q * t_xi
                mdNablaEtadXi[ 1 ] = ( mDetJ * mM[ 53 ] - mM[ 37 ] * mM[ 44 ] ) * mM[ 36 ];

                // t * v_xi - v * t_xi
                mdNablaEtadXi[ 2 ] = ( mDetJ * mM[ 56 ] - mM[ 37 ] * mM[ 47 ] ) * mM[ 36 ];

                // t * m_eta - m * t_eta
                mdNablaEtadEta[ 0 ] = ( mDetJ * mM[ 59 ] - mM[ 38 ] * mM[ 41 ] ) * mM[ 36 ];

                // t * q_eta - q * t_eta
                mdNablaEtadEta[ 1 ] = ( mDetJ * mM[ 62 ] - mM[ 38 ] * mM[ 44 ] ) * mM[ 36 ];

                // t * v_eta - v * t_eta
                mdNablaEtadEta[ 2 ] = ( mDetJ * mM[ 65 ] - mM[ 38 ] * mM[ 47 ] ) * mM[ 36 ];

                // t * m_zeta - m * t_zeta
                mdNablaEtadZeta[ 0 ] = ( mDetJ * mM[ 68 ] - mM[ 39 ] * mM[ 41 ] ) * mM[ 36 ];

                // t * q_zeta - q * t_zeta
                mdNablaEtadZeta[ 1 ] = ( mDetJ * mM[ 71 ] - mM[ 39 ] * mM[ 44 ] ) * mM[ 36 ];

                // t * v_zeta - v * t_zeta
                mdNablaEtadZeta[ 2 ] = ( mDetJ * mM[ 74 ] - mM[ 39 ] * mM[ 47 ] ) * mM[ 36 ];

                // t * n_xi - n * t_xi
                mdNablaZetadXi[ 0 ] = ( mDetJ * mM[ 51 ] - mM[ 37 ] * mM[ 42 ] ) * mM[ 36 ];

                // t * r_xi - r * t_xi
                mdNablaZetadXi[ 1 ] = ( mDetJ * mM[ 54 ] - mM[ 37 ] * mM[ 45 ] ) * mM[ 36 ];

                // t * w_xi - w * t_xi
                mdNablaZetadXi[ 2 ] = ( mDetJ * mM[ 57 ] - mM[ 37 ] * mM[ 48 ] ) * mM[ 36 ];

                // t * n_eta - n * t_eta
                mdNablaZetadEta[ 0 ] = ( mDetJ * mM[ 60 ] - mM[ 38 ] * mM[ 42 ] ) * mM[ 36 ];

                //  t * r_eta - r * t_eta
                mdNablaZetadEta[ 1 ] = ( mDetJ * mM[ 63 ] - mM[ 38 ] * mM[ 45 ] ) * mM[ 36 ];

                // t * w_eta - w * t_eta
                mdNablaZetadEta[ 2 ] = ( mDetJ * mM[ 66 ] - mM[ 38 ] * mM[ 48 ] ) * mM[ 36 ];

                // t * n_zeta - n * t_zeta
                mdNablaZetadZeta[ 0 ] = ( mDetJ * mM[ 69 ] - mM[ 39 ] * mM[ 42 ] ) * mM[ 36 ];

                // t * r_zeta - r * t_zeta
                mdNablaZetadZeta[ 1 ] = ( mDetJ * mM[ 72 ] - mM[ 39 ] * mM[ 45 ] ) * mM[ 36 ];

                // t * w_zeta - w * t_zeta
                mdNablaZetadZeta[ 2 ] = ( mDetJ * mM[ 75 ] - mM[ 39 ] * mM[ 48 ] ) * mM[ 36 ];

                mdNablaTaudXi [ 0 ] = - ( mdNablaXidXi[ 0 ] + mdNablaEtadXi[ 0 ] + mdNablaZetadXi[ 0 ] );
                mdNablaTaudXi [ 1 ] = - ( mdNablaXidXi[ 1 ] + mdNablaEtadXi[ 1 ] + mdNablaZetadXi[ 1 ] );
                mdNablaTaudXi [ 2 ] = - ( mdNablaXidXi[ 2 ] + mdNablaEtadXi[ 2 ] + mdNablaZetadXi[ 2 ] );

                mdNablaTaudEta [ 0 ] = - ( mdNablaXidEta[ 0 ] + mdNablaEtadEta[ 0 ] + mdNablaZetadEta[ 0 ] );
                mdNablaTaudEta [ 1 ] = - ( mdNablaXidEta[ 1 ] + mdNablaEtadEta[ 1 ] + mdNablaZetadEta[ 1 ] );
                mdNablaTaudEta [ 2 ] = - ( mdNablaXidEta[ 2 ] + mdNablaEtadEta[ 2 ] + mdNablaZetadEta[ 2 ] );

                mdNablaTaudZeta [ 0 ] = - ( mdNablaXidZeta[ 0 ] + mdNablaEtadZeta[ 0 ] + mdNablaZetadZeta[ 0 ] );
                mdNablaTaudZeta [ 1 ] = - ( mdNablaXidZeta[ 1 ] + mdNablaEtadZeta[ 1 ] + mdNablaZetadZeta[ 1 ] );
                mdNablaTaudZeta [ 2 ] = - ( mdNablaXidZeta[ 2 ] + mdNablaEtadZeta[ 2 ] + mdNablaZetadZeta[ 2 ] );
                
                mLastNablaDeriv = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_edge_functions( const uint aIndex )
        {
            // edge 0 from xi to eta
            mE( 0, 0 ) = mS[ 0 ] * ( mG( 0, aIndex )  * mNablaEta[ 0 ] + mH( 0, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 0 ) = mS[ 0 ] * ( mG( 0, aIndex )  * mNablaEta[ 1 ] + mH( 0, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 0 ) = mS[ 0 ] * ( mG( 0, aIndex )  * mNablaEta[ 2 ] + mH( 0, aIndex )  * mNablaXi[ 2 ] );
            mE( 0, 1 ) = mS[ 0 ] * ( mG( 1, aIndex )  * mNablaEta[ 0 ] + mH( 1, aIndex )  * mNablaXi[ 0 ] );
            mE( 1, 1 ) = mS[ 0 ] * ( mG( 1, aIndex )  * mNablaEta[ 1 ] + mH( 1, aIndex )  * mNablaXi[ 1 ] );
            mE( 2, 1 ) = mS[ 0 ] * ( mG( 1, aIndex )  * mNablaEta[ 2 ] + mH( 1, aIndex )  * mNablaXi[ 2 ] );


            // edge 1 from eta to zeta
            mE( 0, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaZeta[ 0 ] +  mH( 2, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaZeta[ 1 ] +  mH( 2, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaZeta[ 2 ] +  mH( 2, aIndex ) * mNablaEta[ 2 ] );
            mE( 0, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaZeta[ 0 ] +  mH( 3, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaZeta[ 1 ] +  mH( 3, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaZeta[ 2 ] +  mH( 3, aIndex ) * mNablaEta[ 2 ] );

            // edge 2 from zeta to xi
            mE( 0, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 0 ] + mH( 4, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 1 ] + mH( 4, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 2 ] + mH( 4, aIndex ) * mNablaZeta[ 2 ] );
            mE( 0, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 0 ] + mH( 5, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 1 ] + mH( 5, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 2 ] + mH( 5, aIndex ) * mNablaZeta[ 2 ] );

            // edge 3 from xi to tau
            mE( 0, 6 ) = mS[ 3 ] * ( mG( 6, aIndex ) * mNablaTau[ 0 ] + mH( 6, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 6 ) = mS[ 3 ] * ( mG( 6, aIndex ) * mNablaTau[ 1 ] + mH( 6, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 6 ) = mS[ 3 ] * ( mG( 6, aIndex ) * mNablaTau[ 2 ] + mH( 6, aIndex ) * mNablaXi[ 2 ] );
            mE( 0, 7 ) = mS[ 3 ] * ( mG( 7, aIndex ) * mNablaTau[ 0 ] + mH( 7, aIndex ) * mNablaXi[ 0 ] );
            mE( 1, 7 ) = mS[ 3 ] * ( mG( 7, aIndex ) * mNablaTau[ 1 ] + mH( 7, aIndex ) * mNablaXi[ 1 ] );
            mE( 2, 7 ) = mS[ 3 ] * ( mG( 7, aIndex ) * mNablaTau[ 2 ] + mH( 7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4 from eta to tau
            mE( 0, 8 ) = mS[ 4 ] * (  mG( 8, aIndex ) * mNablaTau[ 0 ] +  mH( 8, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 8 ) = mS[ 4 ] * (  mG( 8, aIndex ) * mNablaTau[ 1 ] +  mH( 8, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 8 ) = mS[ 4 ] * (  mG( 8, aIndex ) * mNablaTau[ 2 ] +  mH( 8, aIndex ) * mNablaEta[ 2 ] );
            mE( 0, 9 ) = mS[ 4 ] * (  mG( 9, aIndex ) * mNablaTau[ 0 ] +  mH( 9, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 9 ) = mS[ 4 ] * (  mG( 9, aIndex ) * mNablaTau[ 1 ] +  mH( 9, aIndex ) * mNablaEta[ 1 ] );
            mE( 2, 9 ) = mS[ 4 ] * (  mG( 9, aIndex ) * mNablaTau[ 2 ] +  mH( 9, aIndex ) * mNablaEta[ 2 ] );

            // edge 5 from zeta to tau
            mE( 0, 10 ) = mS[ 5 ] * ( mG( 10, aIndex ) * mNablaTau[ 0 ] + mH( 10, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 10 ) = mS[ 5 ] * ( mG( 10, aIndex ) * mNablaTau[ 1 ] + mH( 10, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 10 ) = mS[ 5 ] * ( mG( 10, aIndex ) * mNablaTau[ 2 ] + mH( 10, aIndex ) * mNablaZeta[ 2 ] );
            mE( 0, 11 ) = mS[ 5 ] * ( mG( 11, aIndex ) * mNablaTau[ 0 ] + mH( 11, aIndex ) * mNablaZeta[ 0 ] );
            mE( 1, 11 ) = mS[ 5 ] * ( mG( 11, aIndex ) * mNablaTau[ 1 ] + mH( 11, aIndex ) * mNablaZeta[ 1 ] );
            mE( 2, 11 ) = mS[ 5 ] * ( mG( 11, aIndex ) * mNablaTau[ 2 ] + mH( 11, aIndex ) * mNablaZeta[ 2 ] );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_face_functions( const uint aIndex )
        {
            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau # change to eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mF( 0, 0 ) = mU(  0, aIndex ) * mNablaXi[ 0 ] + mV(  0, aIndex ) * mNablaEta[ 0 ] + mW( 0, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 0 ) = mU(  0, aIndex ) * mNablaXi[ 1 ] + mV(  0, aIndex ) * mNablaEta[ 1 ] + mW( 0, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 0 ) = mU(  0, aIndex ) * mNablaXi[ 2 ] + mV(  0, aIndex ) * mNablaEta[ 2 ] + mW( 0, aIndex ) * mNablaTau[ 2 ];

            // eta
            mF( 0, 1 ) = mU(  1, aIndex ) * mNablaXi[ 0 ] + mV(  1, aIndex ) * mNablaEta[ 0 ] + mW( 1, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 1 ) = mU(  1, aIndex ) * mNablaXi[ 1 ] + mV(  1, aIndex ) * mNablaEta[ 1 ] + mW( 1, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 1 ) = mU(  1, aIndex ) * mNablaXi[ 2 ] + mV(  1, aIndex ) * mNablaEta[ 2 ] + mW( 1, aIndex ) * mNablaTau[ 2 ];

            // tau = -xi-eta
            mF( 0, 2 ) = -mF( 0, 0 ) - mF( 0, 1 );
            mF( 1, 2 ) = -mF( 1, 0 ) - mF( 1, 1 );
            mF( 2, 2 ) = -mF( 2, 0 ) - mF( 2, 1 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau # change xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mF( 0, 3 ) = mU( 3, aIndex ) * mNablaEta[ 0 ] + mV( 3, aIndex ) * mNablaZeta[ 0 ] + mW( 3, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 3 ) = mU( 3, aIndex ) * mNablaEta[ 1 ] + mV( 3, aIndex ) * mNablaZeta[ 1 ] + mW( 3, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 3 ) = mU( 3, aIndex ) * mNablaEta[ 2 ] + mV( 3, aIndex ) * mNablaZeta[ 2 ] + mW( 3, aIndex ) * mNablaTau[ 2 ];

            // zeta
            mF( 0, 4 ) = mU( 4, aIndex ) * mNablaEta[ 0 ] + mV( 4, aIndex ) * mNablaZeta[ 0 ] + mW( 4, aIndex ) * mNablaTau[ 0 ]; 
            mF( 1, 4 ) = mU( 4, aIndex ) * mNablaEta[ 1 ] + mV( 4, aIndex ) * mNablaZeta[ 1 ] + mW( 4, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 4 ) = mU( 4, aIndex ) * mNablaEta[ 2 ] + mV( 4, aIndex ) * mNablaZeta[ 2 ] + mW( 4, aIndex ) * mNablaTau[ 2 ];

            // tau = -eta-zeta
            mF( 0, 5 ) = -mF( 0, 3 ) - mF( 0, 4 );
            mF( 1, 5 ) = -mF( 1, 3 ) - mF( 1, 4 );
            mF( 2, 5 ) = -mF( 2, 3 ) - mF( 2, 4 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau # change to eta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mF( 0, 6 ) = mU( 6, aIndex ) * mNablaZeta[ 0 ] + mV( 6, aIndex ) * mNablaXi[ 0 ] + mW( 6, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 6 ) = mU( 6, aIndex ) * mNablaZeta[ 1 ] + mV( 6, aIndex ) * mNablaXi[ 1 ] + mW( 6, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 6 ) = mU( 6, aIndex ) * mNablaZeta[ 2 ] + mV( 6, aIndex ) * mNablaXi[ 2 ] + mW( 6, aIndex ) * mNablaTau[ 2 ];

            // xi
            mF( 0, 7 ) = mU( 7, aIndex ) * mNablaZeta[ 0 ] + mV( 7, aIndex ) * mNablaXi[ 0 ] + mW( 7, aIndex ) * mNablaTau[ 0 ];
            mF( 1, 7 ) = mU( 7, aIndex ) * mNablaZeta[ 1 ] + mV( 7, aIndex ) * mNablaXi[ 1 ] + mW( 7, aIndex ) * mNablaTau[ 1 ];
            mF( 2, 7 ) = mU( 7, aIndex ) * mNablaZeta[ 2 ] + mV( 7, aIndex ) * mNablaXi[ 2 ] + mW( 7, aIndex ) * mNablaTau[ 2 ];

            // tau = -zeta-xi
            mF( 0, 8 ) = -mF( 0, 6 ) - mF( 0, 7 );
            mF( 1, 8 ) = -mF( 1, 6 ) - mF( 1, 7 );
            mF( 2, 8 ) = -mF( 2, 6 ) - mF( 2, 7 );

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta # change to xi->eta->zeta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mF( 0, 9 ) = mU( 9, aIndex ) * mNablaXi[ 0 ] + mV( 9, aIndex ) * mNablaZeta[ 0 ] + mW( 9, aIndex ) * mNablaEta[ 0 ];
            mF( 1, 9 ) = mU( 9, aIndex ) * mNablaXi[ 1 ] + mV( 9, aIndex ) * mNablaZeta[ 1 ] + mW( 9, aIndex ) * mNablaEta[ 1 ];
            mF( 2, 9 ) = mU( 9, aIndex ) * mNablaXi[ 2 ] + mV( 9, aIndex ) * mNablaZeta[ 2 ] + mW( 9, aIndex ) * mNablaEta[ 2 ];

            // zeta
            mF( 0, 10 ) = mU( 10, aIndex ) * mNablaXi[ 0 ] + mV( 10, aIndex ) * mNablaZeta[ 0 ] + mW( 10, aIndex ) * mNablaEta[ 0 ];
            mF( 1, 10 ) = mU( 10, aIndex ) * mNablaXi[ 1 ] + mV( 10, aIndex ) * mNablaZeta[ 1 ] + mW( 10, aIndex ) * mNablaEta[ 1 ];
            mF( 2, 10 ) = mU( 10, aIndex ) * mNablaXi[ 2 ] + mV( 10, aIndex ) * mNablaZeta[ 2 ] + mW( 10, aIndex ) * mNablaEta[ 2 ];

            // eta= = -xi-zeta
            mF( 0, 11 ) = -mF( 0, 9 ) - mF( 0, 10 );
            mF( 1, 11 ) = -mF( 1, 9 ) - mF( 1, 10 );
            mF( 2, 11 ) = -mF( 2, 9 ) - mF( 2, 10 );

        }


//------------------------------------------------------------------------------

        void
        EF_TET10::compute_edge_derivatives_curved( const uint aIndex )
        {
            // edge 0
            mExi( 0, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaEta[ 0 ] + mHxi(  0, aIndex ) * mNablaXi[ 0 ]  + mG(  0, aIndex )  * mdNablaEtadXi[ 0 ] + mH(  0, aIndex ) * mdNablaXidXi[ 0 ] );
            mExi( 1, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaEta[ 1 ] + mHxi(  0, aIndex ) * mNablaXi[ 1 ]  + mG(  0, aIndex )  * mdNablaEtadXi[ 1 ] + mH(  0, aIndex ) * mdNablaXidXi[ 1 ] );
            mExi( 2, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaEta[ 2 ] + mHxi(  0, aIndex )  * mNablaXi[ 2 ] + mG(  0, aIndex )  * mdNablaEtadXi[ 2 ] + mH(  0, aIndex )  * mdNablaXidXi[ 2 ] );
            mExi( 0, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaEta[ 0 ] + mHxi(  1, aIndex )  * mNablaXi[ 0 ] + mG(  1, aIndex )  * mdNablaEtadXi[ 0 ] + mH(  1, aIndex )  * mdNablaXidXi[ 0 ] );
            mExi( 1, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaEta[ 1 ] + mHxi(  1, aIndex )  * mNablaXi[ 1 ] + mG(  1, aIndex )  * mdNablaEtadXi[ 1 ] + mH(  1, aIndex )  * mdNablaXidXi[ 1 ] );
            mExi( 2, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaEta[ 2 ] + mHxi(  1, aIndex )  * mNablaXi[ 2 ] + mG(  1, aIndex )  * mdNablaEtadXi[ 2 ] + mH(  1, aIndex )  * mdNablaXidXi[ 2 ] );


            // edge 1
            mExi( 0, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaZeta[ 0 ] +  mHxi(  2, aIndex ) * mNablaEta[ 0 ] + mG(  2, aIndex ) * mdNablaZetadXi[ 0 ] +  mH(  2, aIndex ) * mdNablaEtadXi[ 0 ] );
            mExi( 1, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaZeta[ 1 ] +  mHxi(  2, aIndex ) * mNablaEta[ 1 ] + mG(  2, aIndex ) * mdNablaZetadXi[ 1 ] +  mH(  2, aIndex ) * mdNablaEtadXi[ 1 ] );
            mExi( 2, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaZeta[ 2 ] +  mHxi(  2, aIndex ) * mNablaEta[ 2 ] + mG(  2, aIndex ) * mdNablaZetadXi[ 2 ] +  mH(  2, aIndex ) * mdNablaEtadXi[ 2 ] );
            mExi( 0, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaZeta[ 0 ] +  mHxi(  3, aIndex ) * mNablaEta[ 0 ] + mG(  3, aIndex ) * mdNablaZetadXi[ 0 ] +  mH(  3, aIndex ) * mdNablaEtadXi[ 0 ] );
            mExi( 1, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaZeta[ 1 ] +  mHxi(  3, aIndex ) * mNablaEta[ 1 ] + mG(  3, aIndex ) * mdNablaZetadXi[ 1 ] +  mH(  3, aIndex ) * mdNablaEtadXi[ 1 ] );
            mExi( 2, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaZeta[ 2 ] +  mHxi(  3, aIndex ) * mNablaEta[ 2 ] + mG(  3, aIndex ) * mdNablaZetadXi[ 2 ] +  mH(  3, aIndex ) * mdNablaEtadXi[ 2 ] );

            // edge 2
            mExi( 0, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 0 ] + mHxi(  4, aIndex ) * mNablaZeta[ 0 ] + mG(  4, aIndex ) * mdNablaXidXi[ 0 ] + mH(  4, aIndex ) * mdNablaZetadXi[ 0 ] );
            mExi( 1, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 1 ] + mHxi(  4, aIndex ) * mNablaZeta[ 1 ] + mG(  4, aIndex ) * mdNablaXidXi[ 1 ] + mH(  4, aIndex ) * mdNablaZetadXi[ 1 ] );
            mExi( 2, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 2 ] + mHxi(  4, aIndex ) * mNablaZeta[ 2 ] + mG(  4, aIndex ) * mdNablaXidXi[ 2 ] + mH(  4, aIndex ) * mdNablaZetadXi[ 2 ] );
            mExi( 0, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 0 ] + mHxi(  5, aIndex ) * mNablaZeta[ 0 ] + mG(  5, aIndex ) * mdNablaXidXi[ 0 ] + mH(  5, aIndex ) * mdNablaZetadXi[ 0 ] );
            mExi( 1, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 1 ] + mHxi(  5, aIndex ) * mNablaZeta[ 1 ] + mG(  5, aIndex ) * mdNablaXidXi[ 1 ] + mH(  5, aIndex ) * mdNablaZetadXi[ 1 ] );
            mExi( 2, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 2 ] + mHxi(  5, aIndex ) * mNablaZeta[ 2 ] + mG(  5, aIndex ) * mdNablaXidXi[ 2 ] + mH(  5, aIndex ) * mdNablaZetadXi[ 2 ] );

            // edge 3
            mExi( 0, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 0 ] + mHxi(  6, aIndex ) * mNablaXi[ 0 ] + mG(  6, aIndex ) * mdNablaTaudXi[ 0 ] + mH(  6, aIndex ) * mdNablaXidXi[ 0 ] );
            mExi( 1, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 1 ] + mHxi(  6, aIndex ) * mNablaXi[ 1 ] + mG(  6, aIndex ) * mdNablaTaudXi[ 1 ] + mH(  6, aIndex ) * mdNablaXidXi[ 1 ] );
            mExi( 2, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 2 ] + mHxi(  6, aIndex ) * mNablaXi[ 2 ] + mG(  6, aIndex ) * mdNablaTaudXi[ 2 ] + mH(  6, aIndex ) * mdNablaXidXi[ 2 ] );
            mExi( 0, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 0 ] + mHxi(  7, aIndex ) * mNablaXi[ 0 ] + mG(  7, aIndex ) * mdNablaTaudXi[ 0 ] + mH(  7, aIndex ) * mdNablaXidXi[ 0 ] );
            mExi( 1, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 1 ] + mHxi(  7, aIndex ) * mNablaXi[ 1 ] + mG(  7, aIndex ) * mdNablaTaudXi[ 1 ] + mH(  7, aIndex ) * mdNablaXidXi[ 1 ] );
            mExi( 2, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 2 ] + mHxi(  7, aIndex ) * mNablaXi[ 2 ] + mG(  7, aIndex ) * mdNablaTaudXi[ 2 ] + mH(  7, aIndex ) * mdNablaXidXi[ 2 ] );

            // edge 4
            mExi( 0, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 0 ] + mHxi(  8, aIndex ) * mNablaEta[ 0 ] + mG(  8, aIndex ) * mdNablaTaudXi[ 0 ] +  mH(  8, aIndex ) * mdNablaEtadXi[ 0 ] );
            mExi( 1, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 1 ] + mHxi(  8, aIndex ) * mNablaEta[ 1 ] + mG(  8, aIndex ) * mdNablaTaudXi[ 1 ] +  mH(  8, aIndex ) * mdNablaEtadXi[ 1 ] );
            mExi( 2, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 2 ] + mHxi(  8, aIndex ) * mNablaEta[ 2 ] + mG(  8, aIndex ) * mdNablaTaudXi[ 2 ] +  mH(  8, aIndex ) * mdNablaEtadXi[ 2 ] );
            mExi( 0, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 0 ] + mHxi(  9, aIndex ) * mNablaEta[ 0 ] + mG(  9, aIndex ) * mdNablaTaudXi[ 0 ] +  mH(  9, aIndex ) * mdNablaEtadXi[ 0 ] );
            mExi( 1, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 1 ] + mHxi(  9, aIndex ) * mNablaEta[ 1 ] + mG(  9, aIndex ) * mdNablaTaudXi[ 1 ] +  mH(  9, aIndex ) * mdNablaEtadXi[ 1 ] );
            mExi( 2, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 2 ] + mHxi(  9, aIndex ) * mNablaEta[ 2 ] + mG(  9, aIndex ) * mdNablaTaudXi[ 2 ] +  mH(  9, aIndex ) * mdNablaEtadXi[ 2 ] );

            // edge 5
            mExi( 0, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 0 ] + mHxi( 10, aIndex ) * mNablaZeta[ 0 ] + mG( 10, aIndex ) * mdNablaTaudXi[ 0 ] + mH( 10, aIndex ) * mdNablaZetadXi[ 0 ] );
            mExi( 1, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 1 ] + mHxi( 10, aIndex ) * mNablaZeta[ 1 ] + mG( 10, aIndex ) * mdNablaTaudXi[ 1 ] + mH( 10, aIndex ) * mdNablaZetadXi[ 1 ] );
            mExi( 2, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 2 ] + mHxi( 10, aIndex ) * mNablaZeta[ 2 ] + mG( 10, aIndex ) * mdNablaTaudXi[ 2 ] + mH( 10, aIndex ) * mdNablaZetadXi[ 2 ] );
            mExi( 0, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 0 ] + mHxi( 11, aIndex ) * mNablaZeta[ 0 ] + mG( 11, aIndex ) * mdNablaTaudXi[ 0 ] + mH( 11, aIndex ) * mdNablaZetadXi[ 0 ] );
            mExi( 1, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 1 ] + mHxi( 11, aIndex ) * mNablaZeta[ 1 ] + mG( 11, aIndex ) * mdNablaTaudXi[ 1 ] + mH( 11, aIndex ) * mdNablaZetadXi[ 1 ] );
            mExi( 2, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 2 ] + mHxi( 11, aIndex ) * mNablaZeta[ 2 ] + mG( 11, aIndex ) * mdNablaTaudXi[ 2 ] + mH( 11, aIndex ) * mdNablaZetadXi[ 2 ] );

            // edge 0
            mEeta( 0, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaEta[ 0 ] + mHeta(  0, aIndex ) * mNablaXi[ 0 ]  + mG(  0, aIndex )  * mdNablaEtadEta[ 0 ] + mH(  0, aIndex ) * mdNablaXidEta[ 0 ] );
            mEeta( 1, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaEta[ 1 ] + mHeta(  0, aIndex ) * mNablaXi[ 1 ]  + mG(  0, aIndex )  * mdNablaEtadEta[ 1 ] + mH(  0, aIndex ) * mdNablaXidEta[ 1 ] );
            mEeta( 2, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaEta[ 2 ] + mHeta(  0, aIndex )  * mNablaXi[ 2 ] + mG(  0, aIndex )  * mdNablaEtadEta[ 2 ] + mH(  0, aIndex )  * mdNablaXidEta[ 2 ] );
            mEeta( 0, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaEta[ 0 ] + mHeta(  1, aIndex )  * mNablaXi[ 0 ] + mG(  1, aIndex )  * mdNablaEtadEta[ 0 ] + mH(  1, aIndex )  * mdNablaXidEta[ 0 ] );
            mEeta( 1, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaEta[ 1 ] + mHeta(  1, aIndex )  * mNablaXi[ 1 ] + mG(  1, aIndex )  * mdNablaEtadEta[ 1 ] + mH(  1, aIndex )  * mdNablaXidEta[ 1 ] );
            mEeta( 2, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaEta[ 2 ] + mHeta(  1, aIndex )  * mNablaXi[ 2 ] + mG(  1, aIndex )  * mdNablaEtadEta[ 2 ] + mH(  1, aIndex )  * mdNablaXidEta[ 2 ] );


            // edge 1
            mEeta( 0, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaZeta[ 0 ] +  mHeta(  2, aIndex ) * mNablaEta[ 0 ] + mG(  2, aIndex ) * mdNablaZetadEta[ 0 ] +  mH(  2, aIndex ) * mdNablaEtadEta[ 0 ] );
            mEeta( 1, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaZeta[ 1 ] +  mHeta(  2, aIndex ) * mNablaEta[ 1 ] + mG(  2, aIndex ) * mdNablaZetadEta[ 1 ] +  mH(  2, aIndex ) * mdNablaEtadEta[ 1 ] );
            mEeta( 2, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaZeta[ 2 ] +  mHeta(  2, aIndex ) * mNablaEta[ 2 ] + mG(  2, aIndex ) * mdNablaZetadEta[ 2 ] +  mH(  2, aIndex ) * mdNablaEtadEta[ 2 ] );
            mEeta( 0, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaZeta[ 0 ] +  mHeta(  3, aIndex ) * mNablaEta[ 0 ] + mG(  3, aIndex ) * mdNablaZetadEta[ 0 ] +  mH(  3, aIndex ) * mdNablaEtadEta[ 0 ] );
            mEeta( 1, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaZeta[ 1 ] +  mHeta(  3, aIndex ) * mNablaEta[ 1 ] + mG(  3, aIndex ) * mdNablaZetadEta[ 1 ] +  mH(  3, aIndex ) * mdNablaEtadEta[ 1 ] );
            mEeta( 2, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaZeta[ 2 ] +  mHeta(  3, aIndex ) * mNablaEta[ 2 ] + mG(  3, aIndex ) * mdNablaZetadEta[ 2 ] +  mH(  3, aIndex ) * mdNablaEtadEta[ 2 ] );

            // edge 2
            mEeta( 0, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 0 ] + mHeta(  4, aIndex ) * mNablaZeta[ 0 ] + mG(  4, aIndex ) * mdNablaXidEta[ 0 ] + mH(  4, aIndex ) * mdNablaZetadEta[ 0 ] );
            mEeta( 1, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 1 ] + mHeta(  4, aIndex ) * mNablaZeta[ 1 ] + mG(  4, aIndex ) * mdNablaXidEta[ 1 ] + mH(  4, aIndex ) * mdNablaZetadEta[ 1 ] );
            mEeta( 2, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 2 ] + mHeta(  4, aIndex ) * mNablaZeta[ 2 ] + mG(  4, aIndex ) * mdNablaXidEta[ 2 ] + mH(  4, aIndex ) * mdNablaZetadEta[ 2 ] );
            mEeta( 0, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 0 ] + mHeta(  5, aIndex ) * mNablaZeta[ 0 ] + mG(  5, aIndex ) * mdNablaXidEta[ 0 ] + mH(  5, aIndex ) * mdNablaZetadEta[ 0 ] );
            mEeta( 1, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 1 ] + mHeta(  5, aIndex ) * mNablaZeta[ 1 ] + mG(  5, aIndex ) * mdNablaXidEta[ 1 ] + mH(  5, aIndex ) * mdNablaZetadEta[ 1 ] );
            mEeta( 2, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 2 ] + mHeta(  5, aIndex ) * mNablaZeta[ 2 ] + mG(  5, aIndex ) * mdNablaXidEta[ 2 ] + mH(  5, aIndex ) * mdNablaZetadEta[ 2 ] );

            // edge 3
            mEeta( 0, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 0 ] + mHeta(  6, aIndex ) * mNablaXi[ 0 ] + mG(  6, aIndex ) * mdNablaTaudEta[ 0 ] + mH(  6, aIndex ) * mdNablaXidEta[ 0 ] );
            mEeta( 1, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 1 ] + mHeta(  6, aIndex ) * mNablaXi[ 1 ] + mG(  6, aIndex ) * mdNablaTaudEta[ 1 ] + mH(  6, aIndex ) * mdNablaXidEta[ 1 ] );
            mEeta( 2, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 2 ] + mHeta(  6, aIndex ) * mNablaXi[ 2 ] + mG(  6, aIndex ) * mdNablaTaudEta[ 2 ] + mH(  6, aIndex ) * mdNablaXidEta[ 2 ] );
            mEeta( 0, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 0 ] + mHeta(  7, aIndex ) * mNablaXi[ 0 ] + mG(  7, aIndex ) * mdNablaTaudEta[ 0 ] + mH(  7, aIndex ) * mdNablaXidEta[ 0 ] );
            mEeta( 1, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 1 ] + mHeta(  7, aIndex ) * mNablaXi[ 1 ] + mG(  7, aIndex ) * mdNablaTaudEta[ 1 ] + mH(  7, aIndex ) * mdNablaXidEta[ 1 ] );
            mEeta( 2, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 2 ] + mHeta(  7, aIndex ) * mNablaXi[ 2 ] + mG(  7, aIndex ) * mdNablaTaudEta[ 2 ] + mH(  7, aIndex ) * mdNablaXidEta[ 2 ] );

            // edge 4
            mEeta( 0, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 0 ] + mHeta(  8, aIndex ) * mNablaEta[ 0 ] + mG(  8, aIndex ) * mdNablaTaudEta[ 0 ] +  mH(  8, aIndex ) * mdNablaEtadEta[ 0 ] );
            mEeta( 1, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 1 ] + mHeta(  8, aIndex ) * mNablaEta[ 1 ] + mG(  8, aIndex ) * mdNablaTaudEta[ 1 ] +  mH(  8, aIndex ) * mdNablaEtadEta[ 1 ] );
            mEeta( 2, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 2 ] + mHeta(  8, aIndex ) * mNablaEta[ 2 ] + mG(  8, aIndex ) * mdNablaTaudEta[ 2 ] +  mH(  8, aIndex ) * mdNablaEtadEta[ 2 ] );
            mEeta( 0, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 0 ] + mHeta(  9, aIndex ) * mNablaEta[ 0 ] + mG(  9, aIndex ) * mdNablaTaudEta[ 0 ] +  mH(  9, aIndex ) * mdNablaEtadEta[ 0 ] );
            mEeta( 1, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 1 ] + mHeta(  9, aIndex ) * mNablaEta[ 1 ] + mG(  9, aIndex ) * mdNablaTaudEta[ 1 ] +  mH(  9, aIndex ) * mdNablaEtadEta[ 1 ] );
            mEeta( 2, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 2 ] + mHeta(  9, aIndex ) * mNablaEta[ 2 ] + mG(  9, aIndex ) * mdNablaTaudEta[ 2 ] +  mH(  9, aIndex ) * mdNablaEtadEta[ 2 ] );

            // edge 5
            mEeta( 0, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 0 ] + mHeta( 10, aIndex ) * mNablaZeta[ 0 ] + mG( 10, aIndex ) * mdNablaTaudEta[ 0 ] + mH( 10, aIndex ) * mdNablaZetadEta[ 0 ] );
            mEeta( 1, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 1 ] + mHeta( 10, aIndex ) * mNablaZeta[ 1 ] + mG( 10, aIndex ) * mdNablaTaudEta[ 1 ] + mH( 10, aIndex ) * mdNablaZetadEta[ 1 ] );
            mEeta( 2, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 2 ] + mHeta( 10, aIndex ) * mNablaZeta[ 2 ] + mG( 10, aIndex ) * mdNablaTaudEta[ 2 ] + mH( 10, aIndex ) * mdNablaZetadEta[ 2 ] );
            mEeta( 0, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 0 ] + mHeta( 11, aIndex ) * mNablaZeta[ 0 ] + mG( 11, aIndex ) * mdNablaTaudEta[ 0 ] + mH( 11, aIndex ) * mdNablaZetadEta[ 0 ] );
            mEeta( 1, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 1 ] + mHeta( 11, aIndex ) * mNablaZeta[ 1 ] + mG( 11, aIndex ) * mdNablaTaudEta[ 1 ] + mH( 11, aIndex ) * mdNablaZetadEta[ 1 ] );
            mEeta( 2, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 2 ] + mHeta( 11, aIndex ) * mNablaZeta[ 2 ] + mG( 11, aIndex ) * mdNablaTaudEta[ 2 ] + mH( 11, aIndex ) * mdNablaZetadEta[ 2 ] );

            // edge 0
            mEzeta( 0, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaEta[ 0 ] + mHzeta(  0, aIndex ) * mNablaXi[ 0 ]  + mG(  0, aIndex )  * mdNablaEtadZeta[ 0 ] + mH(  0, aIndex ) * mdNablaXidZeta[ 0 ] );
            mEzeta( 1, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaEta[ 1 ] + mHzeta(  0, aIndex ) * mNablaXi[ 1 ]  + mG(  0, aIndex )  * mdNablaEtadZeta[ 1 ] + mH(  0, aIndex ) * mdNablaXidZeta[ 1 ] );
            mEzeta( 2, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaEta[ 2 ] + mHzeta(  0, aIndex )  * mNablaXi[ 2 ] + mG(  0, aIndex )  * mdNablaEtadZeta[ 2 ] + mH(  0, aIndex )  * mdNablaXidZeta[ 2 ] );
            mEzeta( 0, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaEta[ 0 ] + mHzeta(  1, aIndex )  * mNablaXi[ 0 ] + mG(  1, aIndex )  * mdNablaEtadZeta[ 0 ] + mH(  1, aIndex )  * mdNablaXidZeta[ 0 ] );
            mEzeta( 1, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaEta[ 1 ] + mHzeta(  1, aIndex )  * mNablaXi[ 1 ] + mG(  1, aIndex )  * mdNablaEtadZeta[ 1 ] + mH(  1, aIndex )  * mdNablaXidZeta[ 1 ] );
            mEzeta( 2, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaEta[ 2 ] + mHzeta(  1, aIndex )  * mNablaXi[ 2 ] + mG(  1, aIndex )  * mdNablaEtadZeta[ 2 ] + mH(  1, aIndex )  * mdNablaXidZeta[ 2 ] );


            // edge 1
            mEzeta( 0, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaZeta[ 0 ] +  mHzeta(  2, aIndex ) * mNablaEta[ 0 ] + mG(  2, aIndex ) * mdNablaZetadZeta[ 0 ] +  mH(  2, aIndex ) * mdNablaEtadZeta[ 0 ] );
            mEzeta( 1, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaZeta[ 1 ] +  mHzeta(  2, aIndex ) * mNablaEta[ 1 ] + mG(  2, aIndex ) * mdNablaZetadZeta[ 1 ] +  mH(  2, aIndex ) * mdNablaEtadZeta[ 1 ] );
            mEzeta( 2, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaZeta[ 2 ] +  mHzeta(  2, aIndex ) * mNablaEta[ 2 ] + mG(  2, aIndex ) * mdNablaZetadZeta[ 2 ] +  mH(  2, aIndex ) * mdNablaEtadZeta[ 2 ] );
            mEzeta( 0, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaZeta[ 0 ] +  mHzeta(  3, aIndex ) * mNablaEta[ 0 ] + mG(  3, aIndex ) * mdNablaZetadZeta[ 0 ] +  mH(  3, aIndex ) * mdNablaEtadZeta[ 0 ] );
            mEzeta( 1, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaZeta[ 1 ] +  mHzeta(  3, aIndex ) * mNablaEta[ 1 ] + mG(  3, aIndex ) * mdNablaZetadZeta[ 1 ] +  mH(  3, aIndex ) * mdNablaEtadZeta[ 1 ] );
            mEzeta( 2, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaZeta[ 2 ] +  mHzeta(  3, aIndex ) * mNablaEta[ 2 ] + mG(  3, aIndex ) * mdNablaZetadZeta[ 2 ] +  mH(  3, aIndex ) * mdNablaEtadZeta[ 2 ] );

            // edge 2
            mEzeta( 0, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 0 ] + mHzeta(  4, aIndex ) * mNablaZeta[ 0 ] + mG(  4, aIndex ) * mdNablaXidZeta[ 0 ] + mH(  4, aIndex ) * mdNablaZetadZeta[ 0 ] );
            mEzeta( 1, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 1 ] + mHzeta(  4, aIndex ) * mNablaZeta[ 1 ] + mG(  4, aIndex ) * mdNablaXidZeta[ 1 ] + mH(  4, aIndex ) * mdNablaZetadZeta[ 1 ] );
            mEzeta( 2, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 2 ] + mHzeta(  4, aIndex ) * mNablaZeta[ 2 ] + mG(  4, aIndex ) * mdNablaXidZeta[ 2 ] + mH(  4, aIndex ) * mdNablaZetadZeta[ 2 ] );
            mEzeta( 0, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 0 ] + mHzeta(  5, aIndex ) * mNablaZeta[ 0 ] + mG(  5, aIndex ) * mdNablaXidZeta[ 0 ] + mH(  5, aIndex ) * mdNablaZetadZeta[ 0 ] );
            mEzeta( 1, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 1 ] + mHzeta(  5, aIndex ) * mNablaZeta[ 1 ] + mG(  5, aIndex ) * mdNablaXidZeta[ 1 ] + mH(  5, aIndex ) * mdNablaZetadZeta[ 1 ] );
            mEzeta( 2, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 2 ] + mHzeta(  5, aIndex ) * mNablaZeta[ 2 ] + mG(  5, aIndex ) * mdNablaXidZeta[ 2 ] + mH(  5, aIndex ) * mdNablaZetadZeta[ 2 ] );

            // edge 3
            mEzeta( 0, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 0 ] + mHzeta(  6, aIndex ) * mNablaXi[ 0 ] + mG(  6, aIndex ) * mdNablaTaudZeta[ 0 ] + mH(  6, aIndex ) * mdNablaXidZeta[ 0 ] );
            mEzeta( 1, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 1 ] + mHzeta(  6, aIndex ) * mNablaXi[ 1 ] + mG(  6, aIndex ) * mdNablaTaudZeta[ 1 ] + mH(  6, aIndex ) * mdNablaXidZeta[ 1 ] );
            mEzeta( 2, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 2 ] + mHzeta(  6, aIndex ) * mNablaXi[ 2 ] + mG(  6, aIndex ) * mdNablaTaudZeta[ 2 ] + mH(  6, aIndex ) * mdNablaXidZeta[ 2 ] );
            mEzeta( 0, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 0 ] + mHzeta(  7, aIndex ) * mNablaXi[ 0 ] + mG(  7, aIndex ) * mdNablaTaudZeta[ 0 ] + mH(  7, aIndex ) * mdNablaXidZeta[ 0 ] );
            mEzeta( 1, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 1 ] + mHzeta(  7, aIndex ) * mNablaXi[ 1 ] + mG(  7, aIndex ) * mdNablaTaudZeta[ 1 ] + mH(  7, aIndex ) * mdNablaXidZeta[ 1 ] );
            mEzeta( 2, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 2 ] + mHzeta(  7, aIndex ) * mNablaXi[ 2 ] + mG(  7, aIndex ) * mdNablaTaudZeta[ 2 ] + mH(  7, aIndex ) * mdNablaXidZeta[ 2 ] );

            // edge 4
            mEzeta( 0, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 0 ] + mHzeta(  8, aIndex ) * mNablaEta[ 0 ] + mG(  8, aIndex ) * mdNablaTaudZeta[ 0 ] +  mH(  8, aIndex ) * mdNablaEtadZeta[ 0 ] );
            mEzeta( 1, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 1 ] + mHzeta(  8, aIndex ) * mNablaEta[ 1 ] + mG(  8, aIndex ) * mdNablaTaudZeta[ 1 ] +  mH(  8, aIndex ) * mdNablaEtadZeta[ 1 ] );
            mEzeta( 2, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 2 ] + mHzeta(  8, aIndex ) * mNablaEta[ 2 ] + mG(  8, aIndex ) * mdNablaTaudZeta[ 2 ] +  mH(  8, aIndex ) * mdNablaEtadZeta[ 2 ] );
            mEzeta( 0, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 0 ] + mHzeta(  9, aIndex ) * mNablaEta[ 0 ] + mG(  9, aIndex ) * mdNablaTaudZeta[ 0 ] +  mH(  9, aIndex ) * mdNablaEtadZeta[ 0 ] );
            mEzeta( 1, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 1 ] + mHzeta(  9, aIndex ) * mNablaEta[ 1 ] + mG(  9, aIndex ) * mdNablaTaudZeta[ 1 ] +  mH(  9, aIndex ) * mdNablaEtadZeta[ 1 ] );
            mEzeta( 2, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 2 ] + mHzeta(  9, aIndex ) * mNablaEta[ 2 ] + mG(  9, aIndex ) * mdNablaTaudZeta[ 2 ] +  mH(  9, aIndex ) * mdNablaEtadZeta[ 2 ] );

            // edge 5
            mEzeta( 0, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 0 ] + mHzeta( 10, aIndex ) * mNablaZeta[ 0 ] + mG( 10, aIndex ) * mdNablaTaudZeta[ 0 ] + mH( 10, aIndex ) * mdNablaZetadZeta[ 0 ] );
            mEzeta( 1, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 1 ] + mHzeta( 10, aIndex ) * mNablaZeta[ 1 ] + mG( 10, aIndex ) * mdNablaTaudZeta[ 1 ] + mH( 10, aIndex ) * mdNablaZetadZeta[ 1 ] );
            mEzeta( 2, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 2 ] + mHzeta( 10, aIndex ) * mNablaZeta[ 2 ] + mG( 10, aIndex ) * mdNablaTaudZeta[ 2 ] + mH( 10, aIndex ) * mdNablaZetadZeta[ 2 ] );
            mEzeta( 0, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 0 ] + mHzeta( 11, aIndex ) * mNablaZeta[ 0 ] + mG( 11, aIndex ) * mdNablaTaudZeta[ 0 ] + mH( 11, aIndex ) * mdNablaZetadZeta[ 0 ] );
            mEzeta( 1, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 1 ] + mHzeta( 11, aIndex ) * mNablaZeta[ 1 ] + mG( 11, aIndex ) * mdNablaTaudZeta[ 1 ] + mH( 11, aIndex ) * mdNablaZetadZeta[ 1 ] );
            mEzeta( 2, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 2 ] + mHzeta( 11, aIndex ) * mNablaZeta[ 2 ] + mG( 11, aIndex ) * mdNablaTaudZeta[ 2 ] + mH( 11, aIndex ) * mdNablaZetadZeta[ 2 ] );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_edge_derivatives_straight( const uint aIndex )
        {
            // edge 0
            mExi( 0, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaEta[ 0 ] + mHxi(  0, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaEta[ 1 ] + mHxi(  0, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 0 ) = mS[ 0 ] * ( mGxi(  0, aIndex )  * mNablaEta[ 2 ] + mHxi(  0, aIndex )  * mNablaXi[ 2 ] );
            mExi( 0, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaEta[ 0 ] + mHxi(  1, aIndex )  * mNablaXi[ 0 ] );
            mExi( 1, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaEta[ 1 ] + mHxi(  1, aIndex )  * mNablaXi[ 1 ] );
            mExi( 2, 1 ) = mS[ 0 ] * ( mGxi(  1, aIndex )  * mNablaEta[ 2 ] + mHxi(  1, aIndex )  * mNablaXi[ 2 ] );


            // edge 1
            mExi( 0, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaZeta[ 0 ] +  mHxi(  2, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaZeta[ 1 ] +  mHxi(  2, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 2 ) = mS[ 1 ] * ( mGxi(  2, aIndex ) * mNablaZeta[ 2 ] +  mHxi(  2, aIndex ) * mNablaEta[ 2 ] );
            mExi( 0, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaZeta[ 0 ] +  mHxi(  3, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaZeta[ 1 ] +  mHxi(  3, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 3 ) = mS[ 1 ] * ( mGxi(  3, aIndex ) * mNablaZeta[ 2 ] +  mHxi(  3, aIndex ) * mNablaEta[ 2 ] );

            // edge 2
            mExi( 0, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 0 ] + mHxi(  4, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 1 ] + mHxi(  4, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 4 ) = mS[ 2 ] * ( mGxi(  4, aIndex ) * mNablaXi[ 2 ] + mHxi(  4, aIndex ) * mNablaZeta[ 2 ] );
            mExi( 0, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 0 ] + mHxi(  5, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 1 ] + mHxi(  5, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 5 ) = mS[ 2 ] * ( mGxi(  5, aIndex ) * mNablaXi[ 2 ] + mHxi(  5, aIndex ) * mNablaZeta[ 2 ] );

            // edge 3
            mExi( 0, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 0 ] + mHxi(  6, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 1 ] + mHxi(  6, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 6 ) = mS[ 3 ] * ( mGxi(  6, aIndex ) * mNablaTau[ 2 ] + mHxi(  6, aIndex ) * mNablaXi[ 2 ] );
            mExi( 0, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 0 ] + mHxi(  7, aIndex ) * mNablaXi[ 0 ] );
            mExi( 1, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 1 ] + mHxi(  7, aIndex ) * mNablaXi[ 1 ] );
            mExi( 2, 7 ) = mS[ 3 ] * ( mGxi(  7, aIndex ) * mNablaTau[ 2 ] + mHxi(  7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4
            mExi( 0, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 0 ] + mHxi(  8, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 1 ] + mHxi(  8, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 8 ) = mS[ 4 ] * (  mGxi(  8, aIndex ) * mNablaTau[ 2 ] + mHxi(  8, aIndex ) * mNablaEta[ 2 ] );
            mExi( 0, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 0 ] + mHxi(  9, aIndex ) * mNablaEta[ 0 ] );
            mExi( 1, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 1 ] + mHxi(  9, aIndex ) * mNablaEta[ 1 ] );
            mExi( 2, 9 ) = mS[ 4 ] * (  mGxi(  9, aIndex ) * mNablaTau[ 2 ] + mHxi(  9, aIndex ) * mNablaEta[ 2 ] );

            // edge 5
            mExi( 0, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 0 ] + mHxi( 10, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 1 ] + mHxi( 10, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 10 ) = mS[ 5 ] * ( mGxi( 10, aIndex ) * mNablaTau[ 2 ] + mHxi( 10, aIndex ) * mNablaZeta[ 2 ] );
            mExi( 0, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 0 ] + mHxi( 11, aIndex ) * mNablaZeta[ 0 ] );
            mExi( 1, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 1 ] + mHxi( 11, aIndex ) * mNablaZeta[ 1 ] );
            mExi( 2, 11 ) = mS[ 5 ] * ( mGxi( 11, aIndex ) * mNablaTau[ 2 ] + mHxi( 11, aIndex ) * mNablaZeta[ 2 ] );

            // edge 0
            mEeta( 0, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaEta[ 0 ] + mHeta(  0, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaEta[ 1 ] + mHeta(  0, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 0 ) = mS[ 0 ] * ( mGeta(  0, aIndex )  * mNablaEta[ 2 ] + mHeta(  0, aIndex )  * mNablaXi[ 2 ] );
            mEeta( 0, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaEta[ 0 ] + mHeta(  1, aIndex )  * mNablaXi[ 0 ] );
            mEeta( 1, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaEta[ 1 ] + mHeta(  1, aIndex )  * mNablaXi[ 1 ] );
            mEeta( 2, 1 ) = mS[ 0 ] * ( mGeta(  1, aIndex )  * mNablaEta[ 2 ] + mHeta(  1, aIndex )  * mNablaXi[ 2 ] );


            // edge 1
            mEeta( 0, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaZeta[ 0 ] +  mHeta(  2, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaZeta[ 1 ] +  mHeta(  2, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 2 ) = mS[ 1 ] * ( mGeta(  2, aIndex ) * mNablaZeta[ 2 ] +  mHeta(  2, aIndex ) * mNablaEta[ 2 ] );
            mEeta( 0, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaZeta[ 0 ] +  mHeta(  3, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaZeta[ 1 ] +  mHeta(  3, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 3 ) = mS[ 1 ] * ( mGeta(  3, aIndex ) * mNablaZeta[ 2 ] +  mHeta(  3, aIndex ) * mNablaEta[ 2 ] );

            // edge 2
            mEeta( 0, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 0 ] + mHeta(  4, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 1 ] + mHeta(  4, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 4 ) = mS[ 2 ] * ( mGeta(  4, aIndex ) * mNablaXi[ 2 ] + mHeta(  4, aIndex ) * mNablaZeta[ 2 ] );
            mEeta( 0, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 0 ] + mHeta(  5, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 1 ] + mHeta(  5, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 5 ) = mS[ 2 ] * ( mGeta(  5, aIndex ) * mNablaXi[ 2 ] + mHeta(  5, aIndex ) * mNablaZeta[ 2 ] );

            // edge 3
            mEeta( 0, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 0 ] + mHeta(  6, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 1 ] + mHeta(  6, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 6 ) = mS[ 3 ] * ( mGeta(  6, aIndex ) * mNablaTau[ 2 ] + mHeta(  6, aIndex ) * mNablaXi[ 2 ] );
            mEeta( 0, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 0 ] + mHeta(  7, aIndex ) * mNablaXi[ 0 ] );
            mEeta( 1, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 1 ] + mHeta(  7, aIndex ) * mNablaXi[ 1 ] );
            mEeta( 2, 7 ) = mS[ 3 ] * ( mGeta(  7, aIndex ) * mNablaTau[ 2 ] + mHeta(  7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4
            mEeta( 0, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 0 ] + mHeta(  8, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 1 ] + mHeta(  8, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 8 ) = mS[ 4 ] * (  mGeta(  8, aIndex ) * mNablaTau[ 2 ] + mHeta(  8, aIndex ) * mNablaEta[ 2 ] );
            mEeta( 0, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 0 ] + mHeta(  9, aIndex ) * mNablaEta[ 0 ] );
            mEeta( 1, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 1 ] + mHeta(  9, aIndex ) * mNablaEta[ 1 ] );
            mEeta( 2, 9 ) = mS[ 4 ] * (  mGeta(  9, aIndex ) * mNablaTau[ 2 ] + mHeta(  9, aIndex ) * mNablaEta[ 2 ] );

            // edge 5
            mEeta( 0, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 0 ] + mHeta( 10, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 1 ] + mHeta( 10, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 10 ) = mS[ 5 ] * ( mGeta( 10, aIndex ) * mNablaTau[ 2 ] + mHeta( 10, aIndex ) * mNablaZeta[ 2 ] );
            mEeta( 0, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 0 ] + mHeta( 11, aIndex ) * mNablaZeta[ 0 ] );
            mEeta( 1, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 1 ] + mHeta( 11, aIndex ) * mNablaZeta[ 1 ] );
            mEeta( 2, 11 ) = mS[ 5 ] * ( mGeta( 11, aIndex ) * mNablaTau[ 2 ] + mHeta( 11, aIndex ) * mNablaZeta[ 2 ] );

            // edge 0
            mEzeta( 0, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaEta[ 0 ] + mHzeta(  0, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaEta[ 1 ] + mHzeta(  0, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 0 ) = mS[ 0 ] * ( mGzeta(  0, aIndex )  * mNablaEta[ 2 ] + mHzeta(  0, aIndex )  * mNablaXi[ 2 ] );
            mEzeta( 0, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaEta[ 0 ] + mHzeta(  1, aIndex )  * mNablaXi[ 0 ] );
            mEzeta( 1, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaEta[ 1 ] + mHzeta(  1, aIndex )  * mNablaXi[ 1 ] );
            mEzeta( 2, 1 ) = mS[ 0 ] * ( mGzeta(  1, aIndex )  * mNablaEta[ 2 ] + mHzeta(  1, aIndex )  * mNablaXi[ 2 ] );


            // edge 1
            mEzeta( 0, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaZeta[ 0 ] +  mHzeta(  2, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaZeta[ 1 ] +  mHzeta(  2, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 2 ) = mS[ 1 ] * ( mGzeta(  2, aIndex ) * mNablaZeta[ 2 ] +  mHzeta(  2, aIndex ) * mNablaEta[ 2 ] );
            mEzeta( 0, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaZeta[ 0 ] +  mHzeta(  3, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaZeta[ 1 ] +  mHzeta(  3, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 3 ) = mS[ 1 ] * ( mGzeta(  3, aIndex ) * mNablaZeta[ 2 ] +  mHzeta(  3, aIndex ) * mNablaEta[ 2 ] );

            // edge 2
            mEzeta( 0, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 0 ] + mHzeta(  4, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 1 ] + mHzeta(  4, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 4 ) = mS[ 2 ] * ( mGzeta(  4, aIndex ) * mNablaXi[ 2 ] + mHzeta(  4, aIndex ) * mNablaZeta[ 2 ] );
            mEzeta( 0, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 0 ] + mHzeta(  5, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 1 ] + mHzeta(  5, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 5 ) = mS[ 2 ] * ( mGzeta(  5, aIndex ) * mNablaXi[ 2 ] + mHzeta(  5, aIndex ) * mNablaZeta[ 2 ] );

            // edge 3
            mEzeta( 0, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 0 ] + mHzeta(  6, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 1 ] + mHzeta(  6, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 6 ) = mS[ 3 ] * ( mGzeta(  6, aIndex ) * mNablaTau[ 2 ] + mHzeta(  6, aIndex ) * mNablaXi[ 2 ] );
            mEzeta( 0, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 0 ] + mHzeta(  7, aIndex ) * mNablaXi[ 0 ] );
            mEzeta( 1, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 1 ] + mHzeta(  7, aIndex ) * mNablaXi[ 1 ] );
            mEzeta( 2, 7 ) = mS[ 3 ] * ( mGzeta(  7, aIndex ) * mNablaTau[ 2 ] + mHzeta(  7, aIndex ) * mNablaXi[ 2 ] );

            // edge 4
            mEzeta( 0, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 0 ] + mHzeta(  8, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 1 ] + mHzeta(  8, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 8 ) = mS[ 4 ] * (  mGzeta(  8, aIndex ) * mNablaTau[ 2 ] + mHzeta(  8, aIndex ) * mNablaEta[ 2 ] );
            mEzeta( 0, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 0 ] + mHzeta(  9, aIndex ) * mNablaEta[ 0 ] );
            mEzeta( 1, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 1 ] + mHzeta(  9, aIndex ) * mNablaEta[ 1 ] );
            mEzeta( 2, 9 ) = mS[ 4 ] * (  mGzeta(  9, aIndex ) * mNablaTau[ 2 ] + mHzeta(  9, aIndex ) * mNablaEta[ 2 ] );

            // edge 5
            mEzeta( 0, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 0 ] + mHzeta( 10, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 1 ] + mHzeta( 10, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 10 ) = mS[ 5 ] * ( mGzeta( 10, aIndex ) * mNablaTau[ 2 ] + mHzeta( 10, aIndex ) * mNablaZeta[ 2 ] );
            mEzeta( 0, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 0 ] + mHzeta( 11, aIndex ) * mNablaZeta[ 0 ] );
            mEzeta( 1, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 1 ] + mHzeta( 11, aIndex ) * mNablaZeta[ 1 ] );
            mEzeta( 2, 11 ) = mS[ 5 ] * ( mGzeta( 11, aIndex ) * mNablaTau[ 2 ] + mHzeta( 11, aIndex ) * mNablaZeta[ 2 ] );

        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_face_derivatives_curved( const uint aIndex )
        {
            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFxi( 0, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 0 ] + mVxi(  0, aIndex ) * mNablaEta[ 0 ] + mWxi( 0, aIndex ) * mNablaTau[ 0 ] + mU(  0, aIndex ) * mdNablaXidXi[ 0 ] + mV(  0, aIndex ) * mdNablaEtadXi[ 0 ] + mW( 0, aIndex ) * mdNablaTaudXi[ 0 ] ;
            mFxi( 1, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 1 ] + mVxi(  0, aIndex ) * mNablaEta[ 1 ] + mWxi( 0, aIndex ) * mNablaTau[ 1 ] + mU(  0, aIndex ) * mdNablaXidXi[ 1 ] + mV(  0, aIndex ) * mdNablaEtadXi[ 1 ] + mW( 0, aIndex ) * mdNablaTaudXi[ 1 ] ;
            mFxi( 2, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 2 ] + mVxi(  0, aIndex ) * mNablaEta[ 2 ] + mWxi( 0, aIndex ) * mNablaTau[ 2 ] + mU(  0, aIndex ) * mdNablaXidXi[ 2 ] + mV(  0, aIndex ) * mdNablaEtadXi[ 2 ] + mW( 0, aIndex ) * mdNablaTaudXi[ 2 ] ;

            // eta
            mFxi( 0, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 0 ] + mVxi(  1, aIndex ) * mNablaEta[ 0 ] + mWxi( 1, aIndex ) * mNablaTau[ 0 ] + mU(  1, aIndex ) * mdNablaXidXi[ 0 ] + mV(  1, aIndex ) * mdNablaEtadXi[ 0 ] + mW( 1, aIndex ) * mdNablaTaudXi[ 0 ] ;
            mFxi( 1, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 1 ] + mVxi(  1, aIndex ) * mNablaEta[ 1 ] + mWxi( 1, aIndex ) * mNablaTau[ 1 ] + mU(  1, aIndex ) * mdNablaXidXi[ 1 ] + mV(  1, aIndex ) * mdNablaEtadXi[ 1 ] + mW( 1, aIndex ) * mdNablaTaudXi[ 1 ] ;
            mFxi( 2, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 2 ] + mVxi(  1, aIndex ) * mNablaEta[ 2 ] + mWxi( 1, aIndex ) * mNablaTau[ 2 ] + mU(  1, aIndex ) * mdNablaXidXi[ 2 ] + mV(  1, aIndex ) * mdNablaEtadXi[ 2 ] + mW( 1, aIndex ) * mdNablaTaudXi[ 2 ] ;

            // tau = -xi-eta
            mFxi( 0, 2 ) = -mFxi( 0, 0 ) - mFxi( 0, 1 ) ;
            mFxi( 1, 2 ) = -mFxi( 1, 0 ) - mFxi( 1, 1 ) ;
            mFxi( 2, 2 ) = -mFxi( 2, 0 ) - mFxi( 2, 1 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFxi( 0, 3 ) = mUxi( 3, aIndex ) * mNablaEta[ 0 ] + mVxi( 3, aIndex ) * mNablaZeta[ 0 ] + mWxi( 3, aIndex ) * mNablaTau[ 0 ] +  mU( 3, aIndex ) * mdNablaEtadXi[ 0 ] + mV( 3, aIndex ) * mdNablaZetadXi[ 0 ] + mW( 3, aIndex ) * mdNablaTaudXi[ 0 ] ;
            mFxi( 1, 3 ) = mUxi( 3, aIndex ) * mNablaEta[ 1 ] + mVxi( 3, aIndex ) * mNablaZeta[ 1 ] + mWxi( 3, aIndex ) * mNablaTau[ 1 ] +  mU( 3, aIndex ) * mdNablaEtadXi[ 1 ] + mV( 3, aIndex ) * mdNablaZetadXi[ 1 ] + mW( 3, aIndex ) * mdNablaTaudXi[ 1 ] ;
            mFxi( 2, 3 ) = mUxi( 3, aIndex ) * mNablaEta[ 2 ] + mVxi( 3, aIndex ) * mNablaZeta[ 2 ] + mWxi( 3, aIndex ) * mNablaTau[ 2 ] +  mU( 3, aIndex ) * mdNablaEtadXi[ 2 ] + mV( 3, aIndex ) * mdNablaZetadXi[ 2 ] + mW( 3, aIndex ) * mdNablaTaudXi[ 2 ] ;

            // zeta
            mFxi( 0, 4 ) = mUxi( 4, aIndex ) * mNablaEta[ 0 ] + mVxi( 4, aIndex ) * mNablaZeta[ 0 ] + mWxi( 4, aIndex ) * mNablaTau[ 0 ] +  mU( 4, aIndex ) * mdNablaEtadXi[ 0 ] + mV( 4, aIndex ) * mdNablaZetadXi[ 0 ] + mW( 4, aIndex ) * mdNablaTaudXi[ 0 ] ;
            mFxi( 1, 4 ) = mUxi( 4, aIndex ) * mNablaEta[ 1 ] + mVxi( 4, aIndex ) * mNablaZeta[ 1 ] + mWxi( 4, aIndex ) * mNablaTau[ 1 ] +  mU( 4, aIndex ) * mdNablaEtadXi[ 1 ] + mV( 4, aIndex ) * mdNablaZetadXi[ 1 ] + mW( 4, aIndex ) * mdNablaTaudXi[ 1 ] ;
            mFxi( 2, 4 ) = mUxi( 4, aIndex ) * mNablaEta[ 2 ] + mVxi( 4, aIndex ) * mNablaZeta[ 2 ] + mWxi( 4, aIndex ) * mNablaTau[ 2 ] +  mU( 4, aIndex ) * mdNablaEtadXi[ 2 ] + mV( 4, aIndex ) * mdNablaZetadXi[ 2 ] + mW( 4, aIndex ) * mdNablaTaudXi[ 2 ] ;

            // tau = -eta-zeta
            mFxi( 0, 5 ) = -mFxi( 0, 3 ) - mFxi( 0, 4 ) ;
            mFxi( 1, 5 ) = -mFxi( 1, 3 ) - mFxi( 1, 4 ) ;
            mFxi( 2, 5 ) = -mFxi( 2, 3 ) - mFxi( 2, 4 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mFxi( 0, 6 ) = mUxi( 6, aIndex ) * mNablaZeta[ 0 ] + mVxi( 6, aIndex ) * mNablaXi[ 0 ] + mWxi( 6, aIndex ) * mNablaTau[ 0 ] + mU( 6, aIndex ) * mdNablaZetadXi[ 0 ] + mV( 6, aIndex ) * mdNablaXidXi[ 0 ] + mW( 6, aIndex ) * mdNablaTaudXi[ 0 ] ;
            mFxi( 1, 6 ) = mUxi( 6, aIndex ) * mNablaZeta[ 1 ] + mVxi( 6, aIndex ) * mNablaXi[ 1 ] + mWxi( 6, aIndex ) * mNablaTau[ 1 ] + mU( 6, aIndex ) * mdNablaZetadXi[ 1 ] + mV( 6, aIndex ) * mdNablaXidXi[ 1 ] + mW( 6, aIndex ) * mdNablaTaudXi[ 1 ] ;
            mFxi( 2, 6 ) = mUxi( 6, aIndex ) * mNablaZeta[ 2 ] + mVxi( 6, aIndex ) * mNablaXi[ 2 ] + mWxi( 6, aIndex ) * mNablaTau[ 2 ] + mU( 6, aIndex ) * mdNablaZetadXi[ 2 ] + mV( 6, aIndex ) * mdNablaXidXi[ 2 ] + mW( 6, aIndex ) * mdNablaTaudXi[ 2 ] ;

            // xi
            mFxi( 0, 7 ) = mUxi( 7, aIndex ) * mNablaZeta[ 0 ] + mVxi( 7, aIndex ) * mNablaXi[ 0 ] + mWxi( 7, aIndex ) * mNablaTau[ 0 ] + mU( 7, aIndex ) * mdNablaZetadXi[ 0 ] + mV( 7, aIndex ) * mdNablaXidXi[ 0 ] + mW( 7, aIndex ) * mdNablaTaudXi[ 0 ] ;
            mFxi( 1, 7 ) = mUxi( 7, aIndex ) * mNablaZeta[ 1 ] + mVxi( 7, aIndex ) * mNablaXi[ 1 ] + mWxi( 7, aIndex ) * mNablaTau[ 1 ] + mU( 7, aIndex ) * mdNablaZetadXi[ 1 ] + mV( 7, aIndex ) * mdNablaXidXi[ 1 ] + mW( 7, aIndex ) * mdNablaTaudXi[ 1 ] ;
            mFxi( 2, 7 ) = mUxi( 7, aIndex ) * mNablaZeta[ 2 ] + mVxi( 7, aIndex ) * mNablaXi[ 2 ] + mWxi( 7, aIndex ) * mNablaTau[ 2 ] + mU( 7, aIndex ) * mdNablaZetadXi[ 2 ] + mV( 7, aIndex ) * mdNablaXidXi[ 2 ] + mW( 7, aIndex ) * mdNablaTaudXi[ 2 ] ;

            // tau = -zeta-xi
            mFxi( 0, 8 ) = -mFxi( 0, 6 ) - mFxi( 0, 7 ) ;
            mFxi( 1, 8 ) = -mFxi( 1, 6 ) - mFxi( 1, 7 ) ;
            mFxi( 2, 8 ) = -mFxi( 2, 6 ) - mFxi( 2, 7 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFxi( 0, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 0 ] + mVxi( 9, aIndex ) * mNablaZeta[ 0 ] + mWxi( 9, aIndex ) * mNablaEta[ 0 ] + mU( 9, aIndex ) * mdNablaXidXi[ 0 ] + mV( 9, aIndex ) * mdNablaZetadXi[ 0 ] + mW( 9, aIndex ) * mdNablaEtadXi[ 0 ] ;
            mFxi( 1, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 1 ] + mVxi( 9, aIndex ) * mNablaZeta[ 1 ] + mWxi( 9, aIndex ) * mNablaEta[ 1 ] + mU( 9, aIndex ) * mdNablaXidXi[ 1 ] + mV( 9, aIndex ) * mdNablaZetadXi[ 1 ] + mW( 9, aIndex ) * mdNablaEtadXi[ 1 ] ;
            mFxi( 2, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 2 ] + mVxi( 9, aIndex ) * mNablaZeta[ 2 ] + mWxi( 9, aIndex ) * mNablaEta[ 2 ] + mU( 9, aIndex ) * mdNablaXidXi[ 2 ] + mV( 9, aIndex ) * mdNablaZetadXi[ 2 ] + mW( 9, aIndex ) * mdNablaEtadXi[ 2 ] ;

            // zeta
            mFxi( 0, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 0 ] + mVxi( 10, aIndex ) * mNablaZeta[ 0 ] + mWxi( 10, aIndex ) * mNablaEta[ 0 ] + mU( 10, aIndex ) * mdNablaXidXi[ 0 ] + mV( 10, aIndex ) * mdNablaZetadXi[ 0 ] + mW( 10, aIndex ) * mdNablaEtadXi[ 0 ] ;
            mFxi( 1, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 1 ] + mVxi( 10, aIndex ) * mNablaZeta[ 1 ] + mWxi( 10, aIndex ) * mNablaEta[ 1 ] + mU( 10, aIndex ) * mdNablaXidXi[ 1 ] + mV( 10, aIndex ) * mdNablaZetadXi[ 1 ] + mW( 10, aIndex ) * mdNablaEtadXi[ 1 ] ;
            mFxi( 2, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 2 ] + mVxi( 10, aIndex ) * mNablaZeta[ 2 ] + mWxi( 10, aIndex ) * mNablaEta[ 2 ] + mU( 10, aIndex ) * mdNablaXidXi[ 2 ] + mV( 10, aIndex ) * mdNablaZetadXi[ 2 ] + mW( 10, aIndex ) * mdNablaEtadXi[ 2 ] ;
            // eta= = -xi-zeta
            mFxi( 0, 11 ) = -mFxi( 0, 9 ) - mFxi( 0, 10 ) ;
            mFxi( 1, 11 ) = -mFxi( 1, 9 ) - mFxi( 1, 10 ) ;
            mFxi( 2, 11 ) = -mFxi( 2, 9 ) - mFxi( 2, 10 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFeta( 0, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 0 ] + mVeta(  0, aIndex ) * mNablaEta[ 0 ] + mWeta( 0, aIndex ) * mNablaTau[ 0 ] + mU(  0, aIndex ) * mdNablaXidEta[ 0 ] + mV(  0, aIndex ) * mdNablaEtadEta[ 0 ] + mW( 0, aIndex ) * mdNablaTaudEta[ 0 ] ;
            mFeta( 1, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 1 ] + mVeta(  0, aIndex ) * mNablaEta[ 1 ] + mWeta( 0, aIndex ) * mNablaTau[ 1 ] + mU(  0, aIndex ) * mdNablaXidEta[ 1 ] + mV(  0, aIndex ) * mdNablaEtadEta[ 1 ] + mW( 0, aIndex ) * mdNablaTaudEta[ 1 ] ;
            mFeta( 2, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 2 ] + mVeta(  0, aIndex ) * mNablaEta[ 2 ] + mWeta( 0, aIndex ) * mNablaTau[ 2 ] + mU(  0, aIndex ) * mdNablaXidEta[ 2 ] + mV(  0, aIndex ) * mdNablaEtadEta[ 2 ] + mW( 0, aIndex ) * mdNablaTaudEta[ 2 ] ;

            // eta
            mFeta( 0, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 0 ] + mVeta(  1, aIndex ) * mNablaEta[ 0 ] + mWeta( 1, aIndex ) * mNablaTau[ 0 ] + mU(  1, aIndex ) * mdNablaXidEta[ 0 ] + mV(  1, aIndex ) * mdNablaEtadEta[ 0 ] + mW( 1, aIndex ) * mdNablaTaudEta[ 0 ] ;
            mFeta( 1, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 1 ] + mVeta(  1, aIndex ) * mNablaEta[ 1 ] + mWeta( 1, aIndex ) * mNablaTau[ 1 ] + mU(  1, aIndex ) * mdNablaXidEta[ 1 ] + mV(  1, aIndex ) * mdNablaEtadEta[ 1 ] + mW( 1, aIndex ) * mdNablaTaudEta[ 1 ] ;
            mFeta( 2, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 2 ] + mVeta(  1, aIndex ) * mNablaEta[ 2 ] + mWeta( 1, aIndex ) * mNablaTau[ 2 ] + mU(  1, aIndex ) * mdNablaXidEta[ 2 ] + mV(  1, aIndex ) * mdNablaEtadEta[ 2 ] + mW( 1, aIndex ) * mdNablaTaudEta[ 2 ] ;

            // tau = -xi-eta
            mFeta( 0, 2 ) = -mFeta( 0, 0 ) - mFeta( 0, 1 ) ;
            mFeta( 1, 2 ) = -mFeta( 1, 0 ) - mFeta( 1, 1 ) ;
            mFeta( 2, 2 ) = -mFeta( 2, 0 ) - mFeta( 2, 1 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFeta( 0, 3 ) = mUeta( 3, aIndex ) * mNablaEta[ 0 ] + mVeta( 3, aIndex ) * mNablaZeta[ 0 ] + mWeta( 3, aIndex ) * mNablaTau[ 0 ] +  mU( 3, aIndex ) * mdNablaEtadEta[ 0 ] + mV( 3, aIndex ) * mdNablaZetadEta[ 0 ] + mW( 3, aIndex ) * mdNablaTaudEta[ 0 ] ;
            mFeta( 1, 3 ) = mUeta( 3, aIndex ) * mNablaEta[ 1 ] + mVeta( 3, aIndex ) * mNablaZeta[ 1 ] + mWeta( 3, aIndex ) * mNablaTau[ 1 ] +  mU( 3, aIndex ) * mdNablaEtadEta[ 1 ] + mV( 3, aIndex ) * mdNablaZetadEta[ 1 ] + mW( 3, aIndex ) * mdNablaTaudEta[ 1 ] ;
            mFeta( 2, 3 ) = mUeta( 3, aIndex ) * mNablaEta[ 2 ] + mVeta( 3, aIndex ) * mNablaZeta[ 2 ] + mWeta( 3, aIndex ) * mNablaTau[ 2 ] +  mU( 3, aIndex ) * mdNablaEtadEta[ 2 ] + mV( 3, aIndex ) * mdNablaZetadEta[ 2 ] + mW( 3, aIndex ) * mdNablaTaudEta[ 2 ] ;

            // zeta
            mFeta( 0, 4 ) = mUeta( 4, aIndex ) * mNablaEta[ 0 ] + mVeta( 4, aIndex ) * mNablaZeta[ 0 ] + mWeta( 4, aIndex ) * mNablaTau[ 0 ] +  mU( 4, aIndex ) * mdNablaEtadEta[ 0 ] + mV( 4, aIndex ) * mdNablaZetadEta[ 0 ] + mW( 4, aIndex ) * mdNablaTaudEta[ 0 ] ;
            mFeta( 1, 4 ) = mUeta( 4, aIndex ) * mNablaEta[ 1 ] + mVeta( 4, aIndex ) * mNablaZeta[ 1 ] + mWeta( 4, aIndex ) * mNablaTau[ 1 ] +  mU( 4, aIndex ) * mdNablaEtadEta[ 1 ] + mV( 4, aIndex ) * mdNablaZetadEta[ 1 ] + mW( 4, aIndex ) * mdNablaTaudEta[ 1 ] ;
            mFeta( 2, 4 ) = mUeta( 4, aIndex ) * mNablaEta[ 2 ] + mVeta( 4, aIndex ) * mNablaZeta[ 2 ] + mWeta( 4, aIndex ) * mNablaTau[ 2 ] +  mU( 4, aIndex ) * mdNablaEtadEta[ 2 ] + mV( 4, aIndex ) * mdNablaZetadEta[ 2 ] + mW( 4, aIndex ) * mdNablaTaudEta[ 2 ] ;

            // tau = -eta-zeta
            mFeta( 0, 5 ) = -mFeta( 0, 3 ) - mFeta( 0, 4 ) ;
            mFeta( 1, 5 ) = -mFeta( 1, 3 ) - mFeta( 1, 4 ) ;
            mFeta( 2, 5 ) = -mFeta( 2, 3 ) - mFeta( 2, 4 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mFeta( 0, 6 ) = mUeta( 6, aIndex ) * mNablaZeta[ 0 ] + mVeta( 6, aIndex ) * mNablaXi[ 0 ] + mWeta( 6, aIndex ) * mNablaTau[ 0 ] + mU( 6, aIndex ) * mdNablaZetadEta[ 0 ] + mV( 6, aIndex ) * mdNablaXidEta[ 0 ] + mW( 6, aIndex ) * mdNablaTaudEta[ 0 ] ;
            mFeta( 1, 6 ) = mUeta( 6, aIndex ) * mNablaZeta[ 1 ] + mVeta( 6, aIndex ) * mNablaXi[ 1 ] + mWeta( 6, aIndex ) * mNablaTau[ 1 ] + mU( 6, aIndex ) * mdNablaZetadEta[ 1 ] + mV( 6, aIndex ) * mdNablaXidEta[ 1 ] + mW( 6, aIndex ) * mdNablaTaudEta[ 1 ] ;
            mFeta( 2, 6 ) = mUeta( 6, aIndex ) * mNablaZeta[ 2 ] + mVeta( 6, aIndex ) * mNablaXi[ 2 ] + mWeta( 6, aIndex ) * mNablaTau[ 2 ] + mU( 6, aIndex ) * mdNablaZetadEta[ 2 ] + mV( 6, aIndex ) * mdNablaXidEta[ 2 ] + mW( 6, aIndex ) * mdNablaTaudEta[ 2 ] ;

            // xi
            mFeta( 0, 7 ) = mUeta( 7, aIndex ) * mNablaZeta[ 0 ] + mVeta( 7, aIndex ) * mNablaXi[ 0 ] + mWeta( 7, aIndex ) * mNablaTau[ 0 ] + mU( 7, aIndex ) * mdNablaZetadEta[ 0 ] + mV( 7, aIndex ) * mdNablaXidEta[ 0 ] + mW( 7, aIndex ) * mdNablaTaudEta[ 0 ] ;
            mFeta( 1, 7 ) = mUeta( 7, aIndex ) * mNablaZeta[ 1 ] + mVeta( 7, aIndex ) * mNablaXi[ 1 ] + mWeta( 7, aIndex ) * mNablaTau[ 1 ] + mU( 7, aIndex ) * mdNablaZetadEta[ 1 ] + mV( 7, aIndex ) * mdNablaXidEta[ 1 ] + mW( 7, aIndex ) * mdNablaTaudEta[ 1 ] ;
            mFeta( 2, 7 ) = mUeta( 7, aIndex ) * mNablaZeta[ 2 ] + mVeta( 7, aIndex ) * mNablaXi[ 2 ] + mWeta( 7, aIndex ) * mNablaTau[ 2 ] + mU( 7, aIndex ) * mdNablaZetadEta[ 2 ] + mV( 7, aIndex ) * mdNablaXidEta[ 2 ] + mW( 7, aIndex ) * mdNablaTaudEta[ 2 ] ;

            // tau = -zeta-xi
            mFeta( 0, 8 ) = -mFeta( 0, 6 ) - mFeta( 0, 7 ) ;
            mFeta( 1, 8 ) = -mFeta( 1, 6 ) - mFeta( 1, 7 ) ;
            mFeta( 2, 8 ) = -mFeta( 2, 6 ) - mFeta( 2, 7 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFeta( 0, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 0 ] + mVeta( 9, aIndex ) * mNablaZeta[ 0 ] + mWeta( 9, aIndex ) * mNablaEta[ 0 ] + mU( 9, aIndex ) * mdNablaXidEta[ 0 ] + mV( 9, aIndex ) * mdNablaZetadEta[ 0 ] + mW( 9, aIndex ) * mdNablaEtadEta[ 0 ] ;
            mFeta( 1, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 1 ] + mVeta( 9, aIndex ) * mNablaZeta[ 1 ] + mWeta( 9, aIndex ) * mNablaEta[ 1 ] + mU( 9, aIndex ) * mdNablaXidEta[ 1 ] + mV( 9, aIndex ) * mdNablaZetadEta[ 1 ] + mW( 9, aIndex ) * mdNablaEtadEta[ 1 ] ;
            mFeta( 2, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 2 ] + mVeta( 9, aIndex ) * mNablaZeta[ 2 ] + mWeta( 9, aIndex ) * mNablaEta[ 2 ] + mU( 9, aIndex ) * mdNablaXidEta[ 2 ] + mV( 9, aIndex ) * mdNablaZetadEta[ 2 ] + mW( 9, aIndex ) * mdNablaEtadEta[ 2 ] ;

            // zeta
            mFeta( 0, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 0 ] + mVeta( 10, aIndex ) * mNablaZeta[ 0 ] + mWeta( 10, aIndex ) * mNablaEta[ 0 ] + mU( 10, aIndex ) * mdNablaXidEta[ 0 ] + mV( 10, aIndex ) * mdNablaZetadEta[ 0 ] + mW( 10, aIndex ) * mdNablaEtadEta[ 0 ] ;
            mFeta( 1, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 1 ] + mVeta( 10, aIndex ) * mNablaZeta[ 1 ] + mWeta( 10, aIndex ) * mNablaEta[ 1 ] + mU( 10, aIndex ) * mdNablaXidEta[ 1 ] + mV( 10, aIndex ) * mdNablaZetadEta[ 1 ] + mW( 10, aIndex ) * mdNablaEtadEta[ 1 ] ;
            mFeta( 2, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 2 ] + mVeta( 10, aIndex ) * mNablaZeta[ 2 ] + mWeta( 10, aIndex ) * mNablaEta[ 2 ] + mU( 10, aIndex ) * mdNablaXidEta[ 2 ] + mV( 10, aIndex ) * mdNablaZetadEta[ 2 ] + mW( 10, aIndex ) * mdNablaEtadEta[ 2 ] ;
            // eta= = -xi-zeta
            mFeta( 0, 11 ) = -mFeta( 0, 9 ) - mFeta( 0, 10 ) ;
            mFeta( 1, 11 ) = -mFeta( 1, 9 ) - mFeta( 1, 10 ) ;
            mFeta( 2, 11 ) = -mFeta( 2, 9 ) - mFeta( 2, 10 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFzeta( 0, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 0 ] + mVzeta(  0, aIndex ) * mNablaEta[ 0 ] + mWzeta( 0, aIndex ) * mNablaTau[ 0 ] + mU(  0, aIndex ) * mdNablaXidZeta[ 0 ] + mV(  0, aIndex ) * mdNablaEtadZeta[ 0 ] + mW( 0, aIndex ) * mdNablaTaudZeta[ 0 ] ;
            mFzeta( 1, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 1 ] + mVzeta(  0, aIndex ) * mNablaEta[ 1 ] + mWzeta( 0, aIndex ) * mNablaTau[ 1 ] + mU(  0, aIndex ) * mdNablaXidZeta[ 1 ] + mV(  0, aIndex ) * mdNablaEtadZeta[ 1 ] + mW( 0, aIndex ) * mdNablaTaudZeta[ 1 ] ;
            mFzeta( 2, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 2 ] + mVzeta(  0, aIndex ) * mNablaEta[ 2 ] + mWzeta( 0, aIndex ) * mNablaTau[ 2 ] + mU(  0, aIndex ) * mdNablaXidZeta[ 2 ] + mV(  0, aIndex ) * mdNablaEtadZeta[ 2 ] + mW( 0, aIndex ) * mdNablaTaudZeta[ 2 ] ;

            // eta
            mFzeta( 0, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 0 ] + mVzeta(  1, aIndex ) * mNablaEta[ 0 ] + mWzeta( 1, aIndex ) * mNablaTau[ 0 ] + mU(  1, aIndex ) * mdNablaXidZeta[ 0 ] + mV(  1, aIndex ) * mdNablaEtadZeta[ 0 ] + mW( 1, aIndex ) * mdNablaTaudZeta[ 0 ] ;
            mFzeta( 1, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 1 ] + mVzeta(  1, aIndex ) * mNablaEta[ 1 ] + mWzeta( 1, aIndex ) * mNablaTau[ 1 ] + mU(  1, aIndex ) * mdNablaXidZeta[ 1 ] + mV(  1, aIndex ) * mdNablaEtadZeta[ 1 ] + mW( 1, aIndex ) * mdNablaTaudZeta[ 1 ] ;
            mFzeta( 2, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 2 ] + mVzeta(  1, aIndex ) * mNablaEta[ 2 ] + mWzeta( 1, aIndex ) * mNablaTau[ 2 ] + mU(  1, aIndex ) * mdNablaXidZeta[ 2 ] + mV(  1, aIndex ) * mdNablaEtadZeta[ 2 ] + mW( 1, aIndex ) * mdNablaTaudZeta[ 2 ] ;

            // tau = -xi-eta
            mFzeta( 0, 2 ) = -mFzeta( 0, 0 ) - mFzeta( 0, 1 ) ;
            mFzeta( 1, 2 ) = -mFzeta( 1, 0 ) - mFzeta( 1, 1 ) ;
            mFzeta( 2, 2 ) = -mFzeta( 2, 0 ) - mFzeta( 2, 1 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFzeta( 0, 3 ) = mUzeta( 3, aIndex ) * mNablaEta[ 0 ] + mVzeta( 3, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 3, aIndex ) * mNablaTau[ 0 ] +  mU( 3, aIndex ) * mdNablaEtadZeta[ 0 ] + mV( 3, aIndex ) * mdNablaZetadZeta[ 0 ] + mW( 3, aIndex ) * mdNablaTaudZeta[ 0 ] ;
            mFzeta( 1, 3 ) = mUzeta( 3, aIndex ) * mNablaEta[ 1 ] + mVzeta( 3, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 3, aIndex ) * mNablaTau[ 1 ] +  mU( 3, aIndex ) * mdNablaEtadZeta[ 1 ] + mV( 3, aIndex ) * mdNablaZetadZeta[ 1 ] + mW( 3, aIndex ) * mdNablaTaudZeta[ 1 ] ;
            mFzeta( 2, 3 ) = mUzeta( 3, aIndex ) * mNablaEta[ 2 ] + mVzeta( 3, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 3, aIndex ) * mNablaTau[ 2 ] +  mU( 3, aIndex ) * mdNablaEtadZeta[ 2 ] + mV( 3, aIndex ) * mdNablaZetadZeta[ 2 ] + mW( 3, aIndex ) * mdNablaTaudZeta[ 2 ] ;

            // zeta
            mFzeta( 0, 4 ) = mUzeta( 4, aIndex ) * mNablaEta[ 0 ] + mVzeta( 4, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 4, aIndex ) * mNablaTau[ 0 ] +  mU( 4, aIndex ) * mdNablaEtadZeta[ 0 ] + mV( 4, aIndex ) * mdNablaZetadZeta[ 0 ] + mW( 4, aIndex ) * mdNablaTaudZeta[ 0 ] ;
            mFzeta( 1, 4 ) = mUzeta( 4, aIndex ) * mNablaEta[ 1 ] + mVzeta( 4, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 4, aIndex ) * mNablaTau[ 1 ] +  mU( 4, aIndex ) * mdNablaEtadZeta[ 1 ] + mV( 4, aIndex ) * mdNablaZetadZeta[ 1 ] + mW( 4, aIndex ) * mdNablaTaudZeta[ 1 ] ;
            mFzeta( 2, 4 ) = mUzeta( 4, aIndex ) * mNablaEta[ 2 ] + mVzeta( 4, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 4, aIndex ) * mNablaTau[ 2 ] +  mU( 4, aIndex ) * mdNablaEtadZeta[ 2 ] + mV( 4, aIndex ) * mdNablaZetadZeta[ 2 ] + mW( 4, aIndex ) * mdNablaTaudZeta[ 2 ] ;

            // tau = -eta-zeta
            mFzeta( 0, 5 ) = -mFzeta( 0, 3 ) - mFzeta( 0, 4 ) ;
            mFzeta( 1, 5 ) = -mFzeta( 1, 3 ) - mFzeta( 1, 4 ) ;
            mFzeta( 2, 5 ) = -mFzeta( 2, 3 ) - mFzeta( 2, 4 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mFzeta( 0, 6 ) = mUzeta( 6, aIndex ) * mNablaZeta[ 0 ] + mVzeta( 6, aIndex ) * mNablaXi[ 0 ] + mWzeta( 6, aIndex ) * mNablaTau[ 0 ] + mU( 6, aIndex ) * mdNablaZetadZeta[ 0 ] + mV( 6, aIndex ) * mdNablaXidZeta[ 0 ] + mW( 6, aIndex ) * mdNablaTaudZeta[ 0 ] ;
            mFzeta( 1, 6 ) = mUzeta( 6, aIndex ) * mNablaZeta[ 1 ] + mVzeta( 6, aIndex ) * mNablaXi[ 1 ] + mWzeta( 6, aIndex ) * mNablaTau[ 1 ] + mU( 6, aIndex ) * mdNablaZetadZeta[ 1 ] + mV( 6, aIndex ) * mdNablaXidZeta[ 1 ] + mW( 6, aIndex ) * mdNablaTaudZeta[ 1 ] ;
            mFzeta( 2, 6 ) = mUzeta( 6, aIndex ) * mNablaZeta[ 2 ] + mVzeta( 6, aIndex ) * mNablaXi[ 2 ] + mWzeta( 6, aIndex ) * mNablaTau[ 2 ] + mU( 6, aIndex ) * mdNablaZetadZeta[ 2 ] + mV( 6, aIndex ) * mdNablaXidZeta[ 2 ] + mW( 6, aIndex ) * mdNablaTaudZeta[ 2 ] ;

            // xi
            mFzeta( 0, 7 ) = mUzeta( 7, aIndex ) * mNablaZeta[ 0 ] + mVzeta( 7, aIndex ) * mNablaXi[ 0 ] + mWzeta( 7, aIndex ) * mNablaTau[ 0 ] + mU( 7, aIndex ) * mdNablaZetadZeta[ 0 ] + mV( 7, aIndex ) * mdNablaXidZeta[ 0 ] + mW( 7, aIndex ) * mdNablaTaudZeta[ 0 ] ;
            mFzeta( 1, 7 ) = mUzeta( 7, aIndex ) * mNablaZeta[ 1 ] + mVzeta( 7, aIndex ) * mNablaXi[ 1 ] + mWzeta( 7, aIndex ) * mNablaTau[ 1 ] + mU( 7, aIndex ) * mdNablaZetadZeta[ 1 ] + mV( 7, aIndex ) * mdNablaXidZeta[ 1 ] + mW( 7, aIndex ) * mdNablaTaudZeta[ 1 ] ;
            mFzeta( 2, 7 ) = mUzeta( 7, aIndex ) * mNablaZeta[ 2 ] + mVzeta( 7, aIndex ) * mNablaXi[ 2 ] + mWzeta( 7, aIndex ) * mNablaTau[ 2 ] + mU( 7, aIndex ) * mdNablaZetadZeta[ 2 ] + mV( 7, aIndex ) * mdNablaXidZeta[ 2 ] + mW( 7, aIndex ) * mdNablaTaudZeta[ 2 ] ;

            // tau = -zeta-xi
            mFzeta( 0, 8 ) = -mFzeta( 0, 6 ) - mFzeta( 0, 7 ) ;
            mFzeta( 1, 8 ) = -mFzeta( 1, 6 ) - mFzeta( 1, 7 ) ;
            mFzeta( 2, 8 ) = -mFzeta( 2, 6 ) - mFzeta( 2, 7 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFzeta( 0, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 0 ] + mVzeta( 9, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 9, aIndex ) * mNablaEta[ 0 ] + mU( 9, aIndex ) * mdNablaXidZeta[ 0 ] + mV( 9, aIndex ) * mdNablaZetadZeta[ 0 ] + mW( 9, aIndex ) * mdNablaEtadZeta[ 0 ] ;
            mFzeta( 1, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 1 ] + mVzeta( 9, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 9, aIndex ) * mNablaEta[ 1 ] + mU( 9, aIndex ) * mdNablaXidZeta[ 1 ] + mV( 9, aIndex ) * mdNablaZetadZeta[ 1 ] + mW( 9, aIndex ) * mdNablaEtadZeta[ 1 ] ;
            mFzeta( 2, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 2 ] + mVzeta( 9, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 9, aIndex ) * mNablaEta[ 2 ] + mU( 9, aIndex ) * mdNablaXidZeta[ 2 ] + mV( 9, aIndex ) * mdNablaZetadZeta[ 2 ] + mW( 9, aIndex ) * mdNablaEtadZeta[ 2 ] ;

            // zeta
            mFzeta( 0, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 0 ] + mVzeta( 10, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 10, aIndex ) * mNablaEta[ 0 ] + mU( 10, aIndex ) * mdNablaXidZeta[ 0 ] + mV( 10, aIndex ) * mdNablaZetadZeta[ 0 ] + mW( 10, aIndex ) * mdNablaEtadZeta[ 0 ] ;
            mFzeta( 1, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 1 ] + mVzeta( 10, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 10, aIndex ) * mNablaEta[ 1 ] + mU( 10, aIndex ) * mdNablaXidZeta[ 1 ] + mV( 10, aIndex ) * mdNablaZetadZeta[ 1 ] + mW( 10, aIndex ) * mdNablaEtadZeta[ 1 ] ;
            mFzeta( 2, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 2 ] + mVzeta( 10, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 10, aIndex ) * mNablaEta[ 2 ] + mU( 10, aIndex ) * mdNablaXidZeta[ 2 ] + mV( 10, aIndex ) * mdNablaZetadZeta[ 2 ] + mW( 10, aIndex ) * mdNablaEtadZeta[ 2 ] ;
            // eta= = -xi-zeta
            mFzeta( 0, 11 ) = -mFzeta( 0, 9 ) - mFzeta( 0, 10 ) ;
            mFzeta( 1, 11 ) = -mFzeta( 1, 9 ) - mFzeta( 1, 10 ) ;
            mFzeta( 2, 11 ) = -mFzeta( 2, 9 ) - mFzeta( 2, 10 ) ;
        }

//------------------------------------------------------------------------------

        void
        EF_TET10::compute_face_derivatives_straight( const uint aIndex )
        {
            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFxi( 0, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 0 ] + mVxi(  0, aIndex ) * mNablaEta[ 0 ] + mWxi( 0, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 1 ] + mVxi(  0, aIndex ) * mNablaEta[ 1 ] + mWxi( 0, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 0 ) = mUxi(  0, aIndex ) * mNablaXi[ 2 ] + mVxi(  0, aIndex ) * mNablaEta[ 2 ] + mWxi( 0, aIndex ) * mNablaTau[ 2 ] ;

            // eta
            mFxi( 0, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 0 ] + mVxi(  1, aIndex ) * mNablaEta[ 0 ] + mWxi( 1, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 1 ] + mVxi(  1, aIndex ) * mNablaEta[ 1 ] + mWxi( 1, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 1 ) = mUxi(  1, aIndex ) * mNablaXi[ 2 ] + mVxi(  1, aIndex ) * mNablaEta[ 2 ] + mWxi( 1, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -xi-eta
            mFxi( 0, 2 ) = -mFxi( 0, 0 ) - mFxi( 0, 1 ) ;
            mFxi( 1, 2 ) = -mFxi( 1, 0 ) - mFxi( 1, 1 ) ;
            mFxi( 2, 2 ) = -mFxi( 2, 0 ) - mFxi( 2, 1 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFxi( 0, 3 ) = mUxi( 3, aIndex ) * mNablaEta[ 0 ] + mVxi( 3, aIndex ) * mNablaZeta[ 0 ] + mWxi( 3, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 3 ) = mUxi( 3, aIndex ) * mNablaEta[ 1 ] + mVxi( 3, aIndex ) * mNablaZeta[ 1 ] + mWxi( 3, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 3 ) = mUxi( 3, aIndex ) * mNablaEta[ 2 ] + mVxi( 3, aIndex ) * mNablaZeta[ 2 ] + mWxi( 3, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFxi( 0, 4 ) = mUxi( 4, aIndex ) * mNablaEta[ 0 ] + mVxi( 4, aIndex ) * mNablaZeta[ 0 ] + mWxi( 4, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 4 ) = mUxi( 4, aIndex ) * mNablaEta[ 1 ] + mVxi( 4, aIndex ) * mNablaZeta[ 1 ] + mWxi( 4, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 4 ) = mUxi( 4, aIndex ) * mNablaEta[ 2 ] + mVxi( 4, aIndex ) * mNablaZeta[ 2 ] + mWxi( 4, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -eta-zeta
            mFxi( 0, 5 ) = -mFxi( 0, 3 ) - mFxi( 0, 4 ) ;
            mFxi( 1, 5 ) = -mFxi( 1, 3 ) - mFxi( 1, 4 ) ;
            mFxi( 2, 5 ) = -mFxi( 2, 3 ) - mFxi( 2, 4 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mFxi( 0, 6 ) = mUxi( 6, aIndex ) * mNablaZeta[ 0 ] + mVxi( 6, aIndex ) * mNablaXi[ 0 ] + mWxi( 6, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 6 ) = mUxi( 6, aIndex ) * mNablaZeta[ 1 ] + mVxi( 6, aIndex ) * mNablaXi[ 1 ] + mWxi( 6, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 6 ) = mUxi( 6, aIndex ) * mNablaZeta[ 2 ] + mVxi( 6, aIndex ) * mNablaXi[ 2 ] + mWxi( 6, aIndex ) * mNablaTau[ 2 ] ;

            // xi
            mFxi( 0, 7 ) = mUxi( 7, aIndex ) * mNablaZeta[ 0 ] + mVxi( 7, aIndex ) * mNablaXi[ 0 ] + mWxi( 7, aIndex ) * mNablaTau[ 0 ] ;
            mFxi( 1, 7 ) = mUxi( 7, aIndex ) * mNablaZeta[ 1 ] + mVxi( 7, aIndex ) * mNablaXi[ 1 ] + mWxi( 7, aIndex ) * mNablaTau[ 1 ] ;
            mFxi( 2, 7 ) = mUxi( 7, aIndex ) * mNablaZeta[ 2 ] + mVxi( 7, aIndex ) * mNablaXi[ 2 ] + mWxi( 7, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -zeta-xi
            mFxi( 0, 8 ) = -mFxi( 0, 6 ) - mFxi( 0, 7 ) ;
            mFxi( 1, 8 ) = -mFxi( 1, 6 ) - mFxi( 1, 7 ) ;
            mFxi( 2, 8 ) = -mFxi( 2, 6 ) - mFxi( 2, 7 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFxi( 0, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 0 ] + mVxi( 9, aIndex ) * mNablaZeta[ 0 ] + mWxi( 9, aIndex ) * mNablaEta[ 0 ] ;
            mFxi( 1, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 1 ] + mVxi( 9, aIndex ) * mNablaZeta[ 1 ] + mWxi( 9, aIndex ) * mNablaEta[ 1 ] ;
            mFxi( 2, 9 ) = mUxi( 9, aIndex ) * mNablaXi[ 2 ] + mVxi( 9, aIndex ) * mNablaZeta[ 2 ] + mWxi( 9, aIndex ) * mNablaEta[ 2 ] ;

            // zeta
            mFxi( 0, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 0 ] + mVxi( 10, aIndex ) * mNablaZeta[ 0 ] + mWxi( 10, aIndex ) * mNablaEta[ 0 ] ;
            mFxi( 1, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 1 ] + mVxi( 10, aIndex ) * mNablaZeta[ 1 ] + mWxi( 10, aIndex ) * mNablaEta[ 1 ] ;
            mFxi( 2, 10 ) = mUxi( 10, aIndex ) * mNablaXi[ 2 ] + mVxi( 10, aIndex ) * mNablaZeta[ 2 ] + mWxi( 10, aIndex ) * mNablaEta[ 2 ] ;

            // eta= -xi-zeta
            mFxi( 0, 11 ) = -mFxi( 0, 9 ) - mFxi( 0, 10 ) ;
            mFxi( 1, 11 ) = -mFxi( 1, 9 ) - mFxi( 1, 10 ) ;
            mFxi( 2, 11 ) = -mFxi( 2, 9 ) - mFxi( 2, 10 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFeta( 0, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 0 ] + mVeta(  0, aIndex ) * mNablaEta[ 0 ] + mWeta( 0, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 1 ] + mVeta(  0, aIndex ) * mNablaEta[ 1 ] + mWeta( 0, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 0 ) = mUeta(  0, aIndex ) * mNablaXi[ 2 ] + mVeta(  0, aIndex ) * mNablaEta[ 2 ] + mWeta( 0, aIndex ) * mNablaTau[ 2 ] ;

            // eta
            mFeta( 0, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 0 ] + mVeta(  1, aIndex ) * mNablaEta[ 0 ] + mWeta( 1, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 1 ] + mVeta(  1, aIndex ) * mNablaEta[ 1 ] + mWeta( 1, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 1 ) = mUeta(  1, aIndex ) * mNablaXi[ 2 ] + mVeta(  1, aIndex ) * mNablaEta[ 2 ] + mWeta( 1, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -xi-eta
            mFeta( 0, 2 ) = -mFeta( 0, 0 ) - mFeta( 0, 1 ) ;
            mFeta( 1, 2 ) = -mFeta( 1, 0 ) - mFeta( 1, 1 ) ;
            mFeta( 2, 2 ) = -mFeta( 2, 0 ) - mFeta( 2, 1 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFeta( 0, 3 ) = mUeta( 3, aIndex ) * mNablaEta[ 0 ] + mVeta( 3, aIndex ) * mNablaZeta[ 0 ] + mWeta( 3, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 3 ) = mUeta( 3, aIndex ) * mNablaEta[ 1 ] + mVeta( 3, aIndex ) * mNablaZeta[ 1 ] + mWeta( 3, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 3 ) = mUeta( 3, aIndex ) * mNablaEta[ 2 ] + mVeta( 3, aIndex ) * mNablaZeta[ 2 ] + mWeta( 3, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFeta( 0, 4 ) = mUeta( 4, aIndex ) * mNablaEta[ 0 ] + mVeta( 4, aIndex ) * mNablaZeta[ 0 ] + mWeta( 4, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 4 ) = mUeta( 4, aIndex ) * mNablaEta[ 1 ] + mVeta( 4, aIndex ) * mNablaZeta[ 1 ] + mWeta( 4, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 4 ) = mUeta( 4, aIndex ) * mNablaEta[ 2 ] + mVeta( 4, aIndex ) * mNablaZeta[ 2 ] + mWeta( 4, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -eta-zeta
            mFeta( 0, 5 ) = -mFeta( 0, 3 ) - mFeta( 0, 4 ) ;
            mFeta( 1, 5 ) = -mFeta( 1, 3 ) - mFeta( 1, 4 ) ;
            mFeta( 2, 5 ) = -mFeta( 2, 3 ) - mFeta( 2, 4 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mFeta( 0, 6 ) = mUeta( 6, aIndex ) * mNablaZeta[ 0 ] + mVeta( 6, aIndex ) * mNablaXi[ 0 ] + mWeta( 6, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 6 ) = mUeta( 6, aIndex ) * mNablaZeta[ 1 ] + mVeta( 6, aIndex ) * mNablaXi[ 1 ] + mWeta( 6, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 6 ) = mUeta( 6, aIndex ) * mNablaZeta[ 2 ] + mVeta( 6, aIndex ) * mNablaXi[ 2 ] + mWeta( 6, aIndex ) * mNablaTau[ 2 ] ;

            // xi
            mFeta( 0, 7 ) = mUeta( 7, aIndex ) * mNablaZeta[ 0 ] + mVeta( 7, aIndex ) * mNablaXi[ 0 ] + mWeta( 7, aIndex ) * mNablaTau[ 0 ] ;
            mFeta( 1, 7 ) = mUeta( 7, aIndex ) * mNablaZeta[ 1 ] + mVeta( 7, aIndex ) * mNablaXi[ 1 ] + mWeta( 7, aIndex ) * mNablaTau[ 1 ] ;
            mFeta( 2, 7 ) = mUeta( 7, aIndex ) * mNablaZeta[ 2 ] + mVeta( 7, aIndex ) * mNablaXi[ 2 ] + mWeta( 7, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -zeta-xi
            mFeta( 0, 8 ) = -mFeta( 0, 6 ) - mFeta( 0, 7 ) ;
            mFeta( 1, 8 ) = -mFeta( 1, 6 ) - mFeta( 1, 7 ) ;
            mFeta( 2, 8 ) = -mFeta( 2, 6 ) - mFeta( 2, 7 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFeta( 0, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 0 ] + mVeta( 9, aIndex ) * mNablaZeta[ 0 ] + mWeta( 9, aIndex ) * mNablaEta[ 0 ] ;
            mFeta( 1, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 1 ] + mVeta( 9, aIndex ) * mNablaZeta[ 1 ] + mWeta( 9, aIndex ) * mNablaEta[ 1 ] ;
            mFeta( 2, 9 ) = mUeta( 9, aIndex ) * mNablaXi[ 2 ] + mVeta( 9, aIndex ) * mNablaZeta[ 2 ] + mWeta( 9, aIndex ) * mNablaEta[ 2 ] ;

            // zeta
            mFeta( 0, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 0 ] + mVeta( 10, aIndex ) * mNablaZeta[ 0 ] + mWeta( 10, aIndex ) * mNablaEta[ 0 ] ;
            mFeta( 1, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 1 ] + mVeta( 10, aIndex ) * mNablaZeta[ 1 ] + mWeta( 10, aIndex ) * mNablaEta[ 1 ] ;
            mFeta( 2, 10 ) = mUeta( 10, aIndex ) * mNablaXi[ 2 ] + mVeta( 10, aIndex ) * mNablaZeta[ 2 ] + mWeta( 10, aIndex ) * mNablaEta[ 2 ] ;

            // eta= -xi-zeta
            mFeta( 0, 11 ) = -mFeta( 0, 9 ) - mFeta( 0, 10 ) ;
            mFeta( 1, 11 ) = -mFeta( 1, 9 ) - mFeta( 1, 10 ) ;
            mFeta( 2, 11 ) = -mFeta( 2, 9 ) - mFeta( 2, 10 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 0 xi->eta->tau
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFzeta( 0, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 0 ] + mVzeta(  0, aIndex ) * mNablaEta[ 0 ] + mWzeta( 0, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 1 ] + mVzeta(  0, aIndex ) * mNablaEta[ 1 ] + mWzeta( 0, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 0 ) = mUzeta(  0, aIndex ) * mNablaXi[ 2 ] + mVzeta(  0, aIndex ) * mNablaEta[ 2 ] + mWzeta( 0, aIndex ) * mNablaTau[ 2 ] ;

            // eta
            mFzeta( 0, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 0 ] + mVzeta(  1, aIndex ) * mNablaEta[ 0 ] + mWzeta( 1, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 1 ] + mVzeta(  1, aIndex ) * mNablaEta[ 1 ] + mWzeta( 1, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 1 ) = mUzeta(  1, aIndex ) * mNablaXi[ 2 ] + mVzeta(  1, aIndex ) * mNablaEta[ 2 ] + mWzeta( 1, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -xi-eta
            mFzeta( 0, 2 ) = -mFzeta( 0, 0 ) - mFzeta( 0, 1 ) ;
            mFzeta( 1, 2 ) = -mFzeta( 1, 0 ) - mFzeta( 1, 1 ) ;
            mFzeta( 2, 2 ) = -mFzeta( 2, 0 ) - mFzeta( 2, 1 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 1 eta->zeta->tau
            // - - - - - - - - - - - - - - - - - - -

            // eta
            mFzeta( 0, 3 ) = mUzeta( 3, aIndex ) * mNablaEta[ 0 ] + mVzeta( 3, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 3, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 3 ) = mUzeta( 3, aIndex ) * mNablaEta[ 1 ] + mVzeta( 3, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 3, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 3 ) = mUzeta( 3, aIndex ) * mNablaEta[ 2 ] + mVzeta( 3, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 3, aIndex ) * mNablaTau[ 2 ] ;

            // zeta
            mFzeta( 0, 4 ) = mUzeta( 4, aIndex ) * mNablaEta[ 0 ] + mVzeta( 4, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 4, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 4 ) = mUzeta( 4, aIndex ) * mNablaEta[ 1 ] + mVzeta( 4, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 4, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 4 ) = mUzeta( 4, aIndex ) * mNablaEta[ 2 ] + mVzeta( 4, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 4, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -eta-zeta
            mFzeta( 0, 5 ) = -mFzeta( 0, 3 ) - mFzeta( 0, 4 ) ;
            mFzeta( 1, 5 ) = -mFzeta( 1, 3 ) - mFzeta( 1, 4 ) ;
            mFzeta( 2, 5 ) = -mFzeta( 2, 3 ) - mFzeta( 2, 4 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 2 zeta->xi->tau
            // - - - - - - - - - - - - - - - - - - -

            // zeta
            mFzeta( 0, 6 ) = mUzeta( 6, aIndex ) * mNablaZeta[ 0 ] + mVzeta( 6, aIndex ) * mNablaXi[ 0 ] + mWzeta( 6, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 6 ) = mUzeta( 6, aIndex ) * mNablaZeta[ 1 ] + mVzeta( 6, aIndex ) * mNablaXi[ 1 ] + mWzeta( 6, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 6 ) = mUzeta( 6, aIndex ) * mNablaZeta[ 2 ] + mVzeta( 6, aIndex ) * mNablaXi[ 2 ] + mWzeta( 6, aIndex ) * mNablaTau[ 2 ] ;

            // xi
            mFzeta( 0, 7 ) = mUzeta( 7, aIndex ) * mNablaZeta[ 0 ] + mVzeta( 7, aIndex ) * mNablaXi[ 0 ] + mWzeta( 7, aIndex ) * mNablaTau[ 0 ] ;
            mFzeta( 1, 7 ) = mUzeta( 7, aIndex ) * mNablaZeta[ 1 ] + mVzeta( 7, aIndex ) * mNablaXi[ 1 ] + mWzeta( 7, aIndex ) * mNablaTau[ 1 ] ;
            mFzeta( 2, 7 ) = mUzeta( 7, aIndex ) * mNablaZeta[ 2 ] + mVzeta( 7, aIndex ) * mNablaXi[ 2 ] + mWzeta( 7, aIndex ) * mNablaTau[ 2 ] ;

            // tau = -zeta-xi
            mFzeta( 0, 8 ) = -mFzeta( 0, 6 ) - mFzeta( 0, 7 ) ;
            mFzeta( 1, 8 ) = -mFzeta( 1, 6 ) - mFzeta( 1, 7 ) ;
            mFzeta( 2, 8 ) = -mFzeta( 2, 6 ) - mFzeta( 2, 7 ) ;

            // - - - - - - - - - - - - - - - - - - -
            // FACE 3 xi->zeta->eta
            // - - - - - - - - - - - - - - - - - - -

            // xi
            mFzeta( 0, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 0 ] + mVzeta( 9, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 9, aIndex ) * mNablaEta[ 0 ] ;
            mFzeta( 1, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 1 ] + mVzeta( 9, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 9, aIndex ) * mNablaEta[ 1 ] ;
            mFzeta( 2, 9 ) = mUzeta( 9, aIndex ) * mNablaXi[ 2 ] + mVzeta( 9, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 9, aIndex ) * mNablaEta[ 2 ] ;

            // zeta
            mFzeta( 0, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 0 ] + mVzeta( 10, aIndex ) * mNablaZeta[ 0 ] + mWzeta( 10, aIndex ) * mNablaEta[ 0 ] ;
            mFzeta( 1, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 1 ] + mVzeta( 10, aIndex ) * mNablaZeta[ 1 ] + mWzeta( 10, aIndex ) * mNablaEta[ 1 ] ;
            mFzeta( 2, 10 ) = mUzeta( 10, aIndex ) * mNablaXi[ 2 ] + mVzeta( 10, aIndex ) * mNablaZeta[ 2 ] + mWzeta( 10, aIndex ) * mNablaEta[ 2 ] ;

            // eta= -xi-zeta
            mFzeta( 0, 11 ) = -mFzeta( 0, 9 ) - mFzeta( 0, 10 ) ;
            mFzeta( 1, 11 ) = -mFzeta( 1, 9 ) - mFzeta( 1, 10 ) ;
            mFzeta( 2, 11 ) = -mFzeta( 2, 9 ) - mFzeta( 2, 10 ) ;

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
                    case( 1 ) : // node 0 of face is identical
                    {
                        aE( 0, tCount ) = aF( 0, tOff );
                        aE( 1, tCount ) = aF( 1, tOff );
                        aE( 2, tCount ) = aF( 2, tOff );
                        ++tCount ;

                        aE( 0, tCount ) = aF( 0, tOff+2 );
                        aE( 1, tCount ) = aF( 1, tOff+2 );
                        aE( 2, tCount ) = aF( 2, tOff+2 );
                        ++tCount ;

                        break ;
                    }
                    case( 2 ) : // node 1 of face is identical
                    {
                        aE( 0, tCount ) = aF( 0, tOff+1 );
                        aE( 1, tCount ) = aF( 1, tOff+1 );
                        aE( 2, tCount ) = aF( 2, tOff+1 );
                        ++tCount ;

                        aE( 0, tCount ) = aF( 0, tOff );
                        aE( 1, tCount ) = aF( 1, tOff );
                        aE( 2, tCount ) = aF( 2, tOff );
                        ++tCount ;

                        break ;
                    }
                    case( 3 ) : // node 2 of face is identical
                    {
                        aE( 0, tCount ) = aF( 0, tOff+2 );
                        aE( 1, tCount ) = aF( 1, tOff+2 );
                        aE( 2, tCount ) = aF( 2, tOff+2 );
                        ++tCount ;

                        aE( 0, tCount ) = aF( 0, tOff+1 );
                        aE( 1, tCount ) = aF( 1, tOff+1 );
                        aE( 2, tCount ) = aF( 2, tOff+1 );
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
        EF_TET10::print_memory()
        {
            std::cout << " a      " << mM[ 27 ] << std::endl ;
            std::cout << " b      " << mM[ 28 ] << std::endl ;
            std::cout << " c      " << mM[ 29 ] << std::endl ;
            std::cout << " d      " << mM[ 30 ] << std::endl ;
            std::cout << " e      " << mM[ 31 ] << std::endl ;
            std::cout << " f      " << mM[ 32 ] << std::endl ;
            std::cout << " g      " << mM[ 33 ] << std::endl ;
            std::cout << " h      " << mM[ 34 ] << std::endl ;
            std::cout << " i      " << mM[ 35 ] << std::endl << std::endl ;


            std::cout << " a_xi   " << mM[  0 ] << std::endl ;
            std::cout << " b_xi   " << mM[  1 ] << std::endl ;
            std::cout << " c_xi   " << mM[  2 ] << std::endl ;
            std::cout << " d_xi   " << mM[  3 ] << std::endl ;
            std::cout << " e_xi   " << mM[  4 ] << std::endl ;
            std::cout << " f_xi   " << mM[  5 ] << std::endl ;
            std::cout << " g_xi   " << mM[  6 ] << std::endl ;
            std::cout << " h_xi   " << mM[  7 ] << std::endl ;
            std::cout << " i_xi   " << mM[  8 ] << std::endl << std::endl ;

            std::cout << " a_eta  " << mM[  9 ] << std::endl ;
            std::cout << " b_eta  " << mM[ 10 ] << std::endl ;
            std::cout << " c_eta  " << mM[ 11 ] << std::endl ;
            std::cout << " d_eta  " << mM[ 12 ] << std::endl ;
            std::cout << " e_eta  " << mM[ 13 ] << std::endl ;
            std::cout << " f_eta  " << mM[ 14 ] << std::endl ;
            std::cout << " g_eta  " << mM[ 15 ] << std::endl ;
            std::cout << " h_eta  " << mM[ 16 ] << std::endl ;
            std::cout << " i_eta  " << mM[ 17 ] << std::endl << std::endl ;

            std::cout << " a_zeta " << mM[ 18 ] << std::endl ;
            std::cout << " b_zeta " << mM[ 19 ] << std::endl ;
            std::cout << " c_zeta " << mM[ 20 ] << std::endl ;
            std::cout << " d_zeta " << mM[ 21 ] << std::endl ;
            std::cout << " e_zeta " << mM[ 22 ] << std::endl ;
            std::cout << " f_zeta " << mM[ 23 ] << std::endl ;
            std::cout << " g_zeta " << mM[ 24 ] << std::endl ;
            std::cout << " h_zeta " << mM[ 25 ] << std::endl ;
            std::cout << " i_zeta " << mM[ 26 ] << std::endl << std::endl ;


            std::cout << " t      " << mDetJ << std::endl ;
            std::cout << " 1/t^2  " << mM[ 36 ] << std::endl ;
            std::cout << " t_xi   " << mM[ 37 ] << std::endl ;
            std::cout << " t_eta  " << mM[ 38 ] << std::endl  ;
            std::cout << " t_zeta " << mM[ 39 ] << std::endl << std::endl ;


            std::cout << " l      " << mM[ 40 ] << std::endl ;
            std::cout << " m      " << mM[ 41 ] << std::endl ;
            std::cout << " n      " << mM[ 42 ] << std::endl ;
            std::cout << " p      " << mM[ 43 ] << std::endl ;
            std::cout << " q      " << mM[ 44 ] << std::endl ;
            std::cout << " r      " << mM[ 45 ] << std::endl ;
            std::cout << " u      " << mM[ 46 ] << std::endl ;
            std::cout << " v      " << mM[ 47 ] << std::endl ;
            std::cout << " w      " << mM[ 48 ] << std::endl << std::endl ;


            std::cout << " l_xi   " << mM[ 49 ] << std::endl ;
            std::cout << " m_xi   " << mM[ 50 ] << std::endl ;
            std::cout << " n_xi   " << mM[ 51 ] << std::endl ;
            std::cout << " p_xi   " << mM[ 52 ] << std::endl ;
            std::cout << " q_xi   " << mM[ 53 ] << std::endl ;
            std::cout << " r_xi   " << mM[ 54 ] << std::endl ;
            std::cout << " u_xi   " << mM[ 55 ] << std::endl ;
            std::cout << " v_xi   " << mM[ 56 ] << std::endl ;
            std::cout << " w_xi   " << mM[ 57 ] << std::endl << std::endl ;

            std::cout << " l_eta  " << mM[ 58 ] << std::endl ;
            std::cout << " m_eta  " << mM[ 59 ] << std::endl ;
            std::cout << " n_eta  " << mM[ 60 ] << std::endl ;
            std::cout << " p_eta  " << mM[ 61 ] << std::endl ;
            std::cout << " q_eta  " << mM[ 62 ] << std::endl ;
            std::cout << " r_eta  " << mM[ 63 ] << std::endl ;
            std::cout << " u_eta  " << mM[ 64 ] << std::endl ;
            std::cout << " v_eta  " << mM[ 65 ] << std::endl ;
            std::cout << " w_eta  " << mM[ 66 ] << std::endl << std::endl ;

            std::cout << " l_zeta " << mM[ 67 ] << std::endl ;
            std::cout << " m_zeta " << mM[ 68 ] << std::endl ;
            std::cout << " n_zeta " << mM[ 69 ] << std::endl ;
            std::cout << " p_zeta " << mM[ 70 ] << std::endl ;
            std::cout << " q_zeta " << mM[ 71 ] << std::endl ;
            std::cout << " r_zeta " << mM[ 72 ] << std::endl ;
            std::cout << " u_zeta " << mM[ 73 ] << std::endl ;
            std::cout << " v_zeta " << mM[ 74 ] << std::endl ;
            std::cout << " w_zeta " << mM[ 75 ] << std::endl << std::endl ;

            std::cout << "nabla" << std::endl ;
            std::cout << mNablaXi[ 0 ] << " " << mNablaEta[ 0 ] << " " << mNablaZeta[ 0 ] << " " << mNablaTau[ 0 ] << std::endl ;
            std::cout << mNablaXi[ 1 ] << " " << mNablaEta[ 1 ] << " " << mNablaZeta[ 1 ] << " " << mNablaTau[ 1 ] << std::endl ;
            std::cout << mNablaXi[ 2 ] << " " << mNablaEta[ 2 ] << " " << mNablaZeta[ 2 ] << " " << mNablaTau[ 2 ] << std::endl ;

            std::cout << "nabla_xi" << std::endl ;
            std::cout << mdNablaXidXi[ 0 ] << " " << mdNablaEtadXi[ 0 ] << " " << mdNablaZetadXi[ 0 ] << " " << mdNablaTaudXi[ 0 ] << std::endl ;
            std::cout << mdNablaXidXi[ 1 ] << " " << mdNablaEtadXi[ 1 ] << " " << mdNablaZetadXi[ 1 ] << " " << mdNablaTaudXi[ 1 ] << std::endl ;
            std::cout << mdNablaXidXi[ 2 ] << " " << mdNablaEtadXi[ 2 ] << " " << mdNablaZetadXi[ 2 ] << " " << mdNablaTaudXi[ 2 ] << std::endl ;

            std::cout << "nabla_eta" << std::endl ;
            std::cout << mdNablaXidEta[ 0 ] << " " << mdNablaEtadEta[ 0 ] << " " << mdNablaZetadEta[ 0 ] << " " << mdNablaTaudEta[ 0 ] << std::endl ;
            std::cout << mdNablaXidEta[ 1 ] << " " << mdNablaEtadEta[ 1 ] << " " << mdNablaZetadEta[ 1 ] << " " << mdNablaTaudEta[ 1 ] << std::endl ;
            std::cout << mdNablaXidEta[ 2 ] << " " << mdNablaEtadEta[ 2 ] << " " << mdNablaZetadEta[ 2 ] << " " << mdNablaTaudEta[ 2 ] << std::endl ;

            std::cout << "nabla_zeta" << std::endl ;
            std::cout << mdNablaXidZeta[ 0 ] << " " << mdNablaEtadZeta[ 0 ] << " " << mdNablaZetadZeta[ 0 ] << " " << mdNablaTaudZeta[ 0 ] << std::endl ;
            std::cout << mdNablaXidZeta[ 1 ] << " " << mdNablaEtadZeta[ 1 ] << " " << mdNablaZetadZeta[ 1 ] << " " << mdNablaTaudZeta[ 1 ] << std::endl ;
            std::cout << mdNablaXidZeta[ 2 ] << " " << mdNablaEtadZeta[ 2 ] << " " << mdNablaZetadZeta[ 2 ] << " " << mdNablaTaudZeta[ 2 ] << std::endl ;


        }

//------------------------------------------------------------------------------
    }
}