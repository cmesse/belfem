//
// Created by christian on 12/2/21.
//

#include "cl_EF_TRI6.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_dot.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_TRI6::EF_TRI6()
        {
            mNodeCoords.set_size( 6, 2 );

            mJ.set_size( 2, 2 );
            mInvJ.set_size( 2, 2 );
            mX.set_size( 6 );
            mY.set_size( 6 );

            mE.set_size( 2, 8 );

            mC.set_size( 1, 8 );

            mB.set_size( 2, 8 );
            mCurlA.set_size( 2, 6 );
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::precompute( const Matrix< real > & aXi )
        {
            // number of points
            uint tN = aXi.n_cols() ;

            mNxi.set_size( 6, tN );
            mNeta.set_size( 6, tN );

            mG.set_size( 8, tN );
            mH.set_size( 8, tN );
            mGxi.set_size( 8, tN );
            mHxi.set_size( 8, tN );
            mGeta.set_size( 8, tN );
            mHeta.set_size( 8, tN );

            for( uint k=0; k<tN; ++k )
            {
                real   xi = aXi( 0, k );
                real  eta = aXi( 1, k );
                real   xi2 = xi + xi ;
                real  eta2 = eta + eta ;

                real   xi3 = xi2 + xi ;
                real  eta3 = eta2 + eta ;

                real   xi4 = xi2 + xi2 ;
                real  eta4 = eta2 + eta2 ;

                real   xi8 = xi4 + xi4 ;
                real  eta8 = eta4 + eta4 ;

                // dN/dxi
                mNxi( 0, k ) = xi4 - 1.;
                mNxi( 1, k ) = 0.;
                mNxi( 2, k ) = xi4 + eta4 - 3.;
                mNxi( 3, k ) = eta4;
                mNxi( 4, k ) =  -eta4;
                mNxi( 5, k ) = 4. - xi8 - eta4;

                mNeta( 0, k ) = 0. ;
                mNeta( 1, k ) = eta4 - 1. ;
                mNeta( 2, k ) = xi4 + eta4 - 3.;
                mNeta( 3, k ) = xi4 ;
                mNeta( 4, k ) = 4.  - xi4 - eta8 ;
                mNeta( 5, k ) = - xi4 ;

                // G
                mG( 0, k ) = eta*(1.-xi4);
                mG( 1, k ) = eta2*( 1.-eta2 );
                mG( 2, k ) = eta2*( 1.-eta2 );
                mG( 3, k ) = eta4*( eta+xi )-eta3 ;
                mG( 4, k ) = eta*(eta4+xi4-6. )+2.-xi3 ;
                mG( 5, k ) = eta+xi3-eta*xi4-1.;
                mG( 6, k ) = eta4*( eta-xi-1.);
                mG( 7, k ) = eta4*( 2.-eta2-xi );
                //mG( 8, k ) = eta4*( eta+xi2-1.0);

                // dG/dxi
                mGxi( 0, k ) = -eta4;
                mGxi( 1, k ) = 0. ;
                mGxi( 2, k ) = 0. ;
                mGxi( 3, k ) = eta4;
                mGxi( 4, k ) = eta4-3.;
                mGxi( 5, k ) = 3.-eta4;
                mGxi( 6, k ) = -eta4;
                mGxi( 7, k ) = -eta4;
                //mGxi( 8, k ) =  eta8;

                // dG/deta
                mGeta( 0, k ) = 1.-xi4 ;
                mGeta( 1, k ) = 2.-eta8;
                mGeta( 2, k ) = 2.-eta8;
                mGeta( 3, k ) = eta8+xi4-3.;
                mGeta( 4, k ) = eta8+xi4-6.;
                mGeta( 5, k ) = 1.-xi4 ;
                mGeta( 6, k ) = eta8-xi4-4.;
                mGeta( 7, k ) = 8.-xi4-eta8-eta8 ;
                //mGeta( 8, k ) = eta8+xi8-4.;

                // H
                mH( 0, k ) = xi2*( xi2-1.) ;
                mH( 1, k ) = xi*( eta4-1.);
                mH( 2, k ) = eta*( xi4-3.)-xi+1.;
                mH( 3, k ) = xi*( 6.-xi4)-eta*( xi4-3.)-2.0;
                mH( 4, k ) = xi3-eta*xi4-xi*xi4;
                mH( 5, k ) = xi*( xi4-2.);
                mH( 6, k ) = xi4*( xi-eta-1.);
                mH( 7, k ) = xi4*( eta2+xi-1.);
                //mH( 8, k ) = xi4*( 2.-eta-xi2);

                // dH/dxi
                mHxi( 0, k ) = xi8-2.;
                mHxi( 1, k ) = eta4-1.;
                mHxi( 2, k ) = eta4-1.;
                mHxi( 3, k ) = 6.-eta4-xi8;
                mHxi( 4, k ) = 3.-eta4-xi8 ;
                mHxi( 5, k ) = xi8-2.;
                mHxi( 6, k ) = xi8-eta4-4.;
                mHxi( 7, k ) = eta8+xi8-4.;
                //mHxi( 8, k ) = 8.-eta4-xi8-xi8 ;

                // dH/deta
                mHeta( 0, k ) = 0. ;
                mHeta( 1, k ) =  xi4;
                mHeta( 2, k ) =  xi4-3.;
                mHeta( 3, k ) = 3.-xi4;
                mHeta( 4, k ) = -xi4;
                mHeta( 5, k ) = 0. ;
                mHeta( 6, k ) = -xi4;
                mHeta( 7, k ) =  xi8;
                //mHeta( 8, k ) = -xi4;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::link( Element  * aElement,
                       const bool aOperatorCurlH,
                       const bool aOperatorGradPhi,
                       const bool aOperatorCurlA )
        {
            // make sure that this is the correct element type
            BELFEM_ASSERT( aElement->element()->type() == ElementType::TRI6,
                          "Element %lu is not of type TRI6",
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

            // get y-coordinates
            mY( 0 ) = tElement->node( 0 )->y();
            mY( 1 ) = tElement->node( 1 )->y();
            mY( 2 ) = tElement->node( 2 )->y();
            mY( 3 ) = tElement->node( 3 )->y();
            mY( 4 ) = tElement->node( 4 )->y();
            mY( 5 ) = tElement->node( 5 )->y();

            // reset the index
            mLastJ          = BELFEM_UINT_MAX ;
            mLastNabla      = BELFEM_UINT_MAX ;
            mLastNablaDeriv = BELFEM_UINT_MAX ;

            // check if element is curved
            if( tElement->is_curved() )
            {
                if( aOperatorCurlH )
                {
                    // get the directions of the edges
                    aElement->edge_directions( mS );

                    // derivatives of jacobian matrix to xi
                    mW[ 0 ] = 4. * ( mX( 0 )+mX( 2 )-mX( 5 )-mX( 5 ) ); // a_xi
                    mW[ 1 ] = 4. * ( mY( 0 )+mY( 2 )-mY( 5 )-mY( 5 ) ); // b_xi
                    mW[ 2 ] = 4. * ( mX( 2 )+mX( 3 )-mX( 4 )-mX( 5 ) ); // c_xi
                    mW[ 3 ] = 4. * ( mY( 2 )+mY( 3 )-mY( 4 )-mY( 5 ) ); // d_xi

                    mW[ 4 ] = 4. * ( mX( 2 )+mX( 3 )-mX( 4 )-mX( 5 ) ); // a_eta
                    mW[ 5 ] = 4. * ( mY( 2 )+mY( 3 )-mY( 4 )-mY( 5 ) ); // b_eta
                    mW[ 6 ] = 4. * ( mX( 1 )+mX( 2 )-mX( 4 )-mX( 4 ) ); // c_eta
                    mW[ 7 ] = 4. * ( mY( 1 )+mY( 2 )-mY( 4 )-mY( 4 ) ); // d_eta

                    // link functions
                    mFunInterpolation = & EF_TRI6::E_curved ;
                    mFunCurl          = & EF_TRI6::C_curved ;
                }
                else
                {
                    mFunInterpolation = nullptr ;
                    mFunCurl          = nullptr ;
                    mC.fill( BELFEM_QUIET_NAN );
                }
                if( aOperatorGradPhi || aOperatorCurlA )
                {
                    mFunGrad = & EF_TRI6::B_curved ;
                }
                else
                {
                    mFunGrad = nullptr ;
                    mB.fill( BELFEM_QUIET_NAN );
                    mInvJ.fill( BELFEM_QUIET_NAN );
                }
            }
            else
            {
                if( aOperatorCurlH )
                {
                    // get the directions of the edges
                    aElement->edge_directions( mS );

                    // compute nabla ( it is constant for this element )
                    this->compute_nabla( 0 );

                    // link functions
                    mFunInterpolation = & EF_TRI6::compute_E ;
                    mFunCurl          = & EF_TRI6::C_straight ;
                }
                else
                {
                    mFunInterpolation = nullptr ;
                    mFunCurl          = nullptr ;
                    mC.fill( BELFEM_QUIET_NAN );
                }
                if( aOperatorGradPhi || aOperatorCurlA )
                {
                    if( ! aOperatorCurlH ) this->compute_J( 0 );
                    mFunGrad = & EF_TRI6::B_straight ;
                }
                else
                {
                    mFunGrad = nullptr ;
                    mB.fill( BELFEM_QUIET_NAN );
                    mCurlA.fill( BELFEM_QUIET_NAN );
                    mInvJ.fill( BELFEM_QUIET_NAN );
                    mJ.fill( BELFEM_QUIET_NAN );
                }
            }


        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::E_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );
            this->compute_E( aIndex );
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::C_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );
            this->compute_nabla_derivatives( aIndex );
            this->compute_E_xi_curved( aIndex );
            this->compute_E_eta_curved( aIndex );
            this->compute_C();
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::C_straight( const uint aIndex )
        {
            this->compute_E_xi_straight( aIndex );
            this->compute_E_eta_straight( aIndex );
            this->compute_C();
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_J( const uint aIndex )
        {
            if( aIndex != mLastJ )
            {
                // Jacobian
                mJ( 0, 0 ) = dot( mNxi.col( aIndex ), mX );
                mJ( 1, 0 ) = dot( mNeta.col( aIndex ), mX );
                mJ( 0, 1 ) = dot( mNxi.col( aIndex ), mY );
                mJ( 1, 1 ) = dot( mNeta.col( aIndex ), mY );

                // determinant
                mDetJ =   det( mJ );

                // absolute value
                mAbsDetJ = std::abs( mDetJ );

                // inverse
                mInvJ = inv( mJ );

                mLastJ = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_nabla( const uint aIndex )
        {
            if( aIndex != mLastNabla )
            {
                // J = [ a b; c d ];
                mW[  8 ] = dot( mNxi.col( aIndex ), mX );  // a
                mW[  9 ] = dot( mNxi.col( aIndex ), mY );  // b
                mW[ 10 ] = dot( mNeta.col( aIndex ), mX ); // c
                mW[ 11 ] = dot( mNeta.col( aIndex ), mY ); // d

                // determinant t = a*d - b*c
                mW[ 12 ] = mW[  8 ]*mW[ 11 ] - mW[  9 ]*mW[ 10 ];

                // copy data
                mDetJ    = mW[ 12 ];
                mAbsDetJ = std::abs( mDetJ );

                mJ( 0, 0 ) = mW[  8 ] ;
                mJ( 1, 0 ) = mW[ 10 ] ;
                mJ( 0, 1 ) = mW[  9 ] ;
                mJ( 1, 1 ) = mW[ 11 ] ;

                // temporary variable, will overwrite later
                mW[ 15 ] = 1. / mDetJ ;

                // Nablas
                mNablaXi[ 0 ]   =  mW[ 11 ] * mW[ 15 ] ;
                mNablaXi[ 1 ]   = -mW[ 10 ] * mW[ 15 ] ;
                mNablaEta[ 0 ]  = -mW[  9 ] * mW[ 15 ] ;
                mNablaEta[ 1 ]  =  mW[  8 ] * mW[ 15 ] ;

                // remember index
                mLastNabla = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_nabla_derivatives( const uint aIndex )
        {
            if( aIndex != mLastNablaDeriv )
            {
                // dt_dxi = a_xi*d + a*d_xi - b_xi*c - b*c_xi ;
                mW[ 13 ] = mW[  0 ]*mW[ 11 ] + mW[  8 ]*mW[  3 ] - mW[  1 ]*mW[ 10 ] - mW[  9 ]*mW[  2 ] ;

                // dt_deta = a_eta*d + a*d_eta - c_eta*b - c*b_eta ;
                mW[ 14 ] = mW[  4 ]*mW[ 11 ] + mW[  8 ]*mW[  7 ] - mW[  5 ]*mW[ 10 ] - mW[  9 ]*mW[  6 ] ;

                // 1/t^2;
                mW[ 15 ] = 1. / ( mW[ 12 ] * mW[ 12 ] );

                // iJ = [ p q; r s ] = [ d -b; -c; a ] / det(J)^2
                // dp_dxi = d_xi * t - d * t_xi
                mdNablaXidXi[ 0 ] = ( mW[  3 ]*mW[ 12 ] - mW[ 11 ]*mW[ 13 ] ) * mW[ 15 ] ;

                // dr_dxi = c * t_xi - t * c_xi
                mdNablaXidXi[ 1 ] = ( mW[ 10 ]*mW[ 13 ] - mW[  2 ]*mW[ 12 ] ) * mW[ 15 ] ;

                // dq_dxi = b * t_xi - t * b_xi
                mdNablaEtadXi[ 0 ] = ( mW[  9 ]*mW[ 13 ] - mW[  1 ]*mW[ 12 ] ) * mW[ 15 ] ;

                // ds_dxi = t * a_xi - a * t_xi
                mdNablaEtadXi[ 1 ] = ( mW[  0 ]*mW[ 12 ] - mW[  8 ]*mW[ 13 ] ) * mW[ 15 ] ;

                // dp_deta = t * d_eta - d * t_eta
                mdNablaXidEta[ 0 ] = ( mW[  7 ]*mW[ 12 ] - mW[ 11 ]*mW[ 14 ] ) * mW[ 15 ] ;

                // dr_deta = c * t_eta - t * c_eta ;
                mdNablaXidEta[ 1 ] = ( mW[  10 ]*mW[ 14 ] - mW[  6 ]*mW[ 12 ] ) * mW[ 15 ] ;

                // dq_deta = b * t_eta - t * b_eta
                mdNablaEtadEta[ 0 ] = ( mW[  9 ]*mW[ 14 ] - mW[  5 ]*mW[ 12 ] ) * mW[ 15 ] ;

                // ds_deta = t * a_eta - a * t_eta ;
                mdNablaEtadEta[ 1 ] = ( mW[  4 ]*mW[ 12 ] - mW[  8 ]*mW[ 14 ] ) * mW[ 15 ] ;

                mLastNablaDeriv = aIndex ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_E( const uint aIndex )
        {
            mE( 0, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaXi[ 0 ] + mH( 0, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 0 ) = mS[ 0 ] * ( mG( 0, aIndex ) * mNablaXi[ 1 ] + mH( 0, aIndex ) * mNablaEta[ 1 ] );
            mE( 0, 1 ) = mS[ 0 ] * ( mG( 1, aIndex ) * mNablaXi[ 0 ] + mH( 1, aIndex ) * mNablaEta[ 0 ] );
            mE( 1, 1 ) = mS[ 0 ] * ( mG( 1, aIndex ) * mNablaXi[ 1 ] + mH( 1, aIndex ) * mNablaEta[ 1 ] );

            mE( 0, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaXi[ 0 ] + mH( 2, aIndex )  * mNablaEta[ 0 ] );
            mE( 1, 2 ) = mS[ 1 ] * ( mG( 2, aIndex ) * mNablaXi[ 1 ] + mH( 2, aIndex )  * mNablaEta[ 1 ] );
            mE( 0, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaXi[ 0 ] + mH( 3, aIndex )  * mNablaEta[ 0 ] );
            mE( 1, 3 ) = mS[ 1 ] * ( mG( 3, aIndex ) * mNablaXi[ 1 ] + mH( 3, aIndex )  * mNablaEta[ 1 ] );

            mE( 0, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 0 ] + mH( 4, aIndex )  * mNablaEta[ 0 ] );
            mE( 1, 4 ) = mS[ 2 ] * ( mG( 4, aIndex ) * mNablaXi[ 1 ] + mH( 4, aIndex )  * mNablaEta[ 1 ] );
            mE( 0, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 0 ] + mH( 5, aIndex )  * mNablaEta[ 0 ] );
            mE( 1, 5 ) = mS[ 2 ] * ( mG( 5, aIndex ) * mNablaXi[ 1 ] + mH( 5, aIndex )  * mNablaEta[ 1 ] );

            mE( 0, 6 ) = mG( 6, aIndex ) * mNablaXi[ 0 ] + mH( 6, aIndex ) * mNablaEta[ 0 ];
            mE( 1, 6 ) = mG( 6, aIndex ) * mNablaXi[ 1 ] + mH( 6, aIndex ) * mNablaEta[ 1 ];

            mE( 0, 7 ) = mG( 7, aIndex ) * mNablaXi[ 0 ] + mH( 7, aIndex ) * mNablaEta[ 0 ];
            mE( 1, 7 ) = mG( 7, aIndex ) * mNablaXi[ 1 ] + mH( 7, aIndex ) * mNablaEta[ 1 ];
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_E_xi_curved( const uint aIndex )
        {
            uint tCount = 0 ;

            for( uint k=0; k<8; ++k )
            {
                mExi[  tCount++ ] =
                      mGxi( k, aIndex ) * mNablaXi[ 0 ]
                    + mG( k, aIndex ) * mdNablaXidXi[ 0 ]
                    + mHxi( k, aIndex ) * mNablaEta[ 0 ]
                    + mH( k, aIndex ) * mdNablaEtadXi[ 0 ] ;

                mExi[  tCount++ ] =
                      mGxi( k, aIndex ) * mNablaXi[ 1 ]
                    + mG( k, aIndex ) * mdNablaXidXi[ 1 ]
                    + mHxi( k, aIndex ) * mNablaEta[ 1 ]
                    + mH( k, aIndex ) * mdNablaEtadXi[ 1 ] ;
            }

        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_E_xi_straight( const uint aIndex )
        {
            uint tCount = 0 ;

            for( uint k=0; k<8; ++k )
            {
                mExi[  tCount++ ] =
                          mGxi( k, aIndex ) * mNablaXi[ 0 ]
                        + mHxi( k, aIndex ) * mNablaEta[ 0 ] ;

                mExi[  tCount++ ] =
                        mGxi( k, aIndex ) * mNablaXi[ 1 ]
                        + mHxi( k, aIndex ) * mNablaEta[ 1 ] ;
            }

        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_E_eta_curved( const uint aIndex )
        {
            uint tCount = 0 ;

            for( uint k=0; k<8; ++k )
            {
                mEeta[  tCount++ ] =
                          mGeta( k, aIndex ) * mNablaXi[ 0 ]
                        + mG( k, aIndex ) * mdNablaXidEta[ 0 ]
                        + mHeta( k, aIndex ) * mNablaEta[ 0 ]
                        + mH( k, aIndex ) * mdNablaEtadEta[ 0 ] ;

                mEeta[  tCount++ ] =
                          mGeta( k, aIndex ) * mNablaXi[ 1 ]
                        + mG( k, aIndex ) * mdNablaXidEta[ 1 ]
                        + mHeta( k, aIndex ) * mNablaEta[ 1 ]
                        + mH( k, aIndex ) * mdNablaEtadEta[ 1 ] ;
            }

        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_E_eta_straight( const uint aIndex )
        {
            uint tCount = 0 ;
            for( uint k=0; k<8; ++k )
            {
                mEeta[  tCount++ ] =
                          mGeta( k, aIndex ) * mNablaXi[ 0 ]
                        + mHeta( k, aIndex ) * mNablaEta[ 0 ] ;

                mEeta[  tCount++ ] =
                          mGeta( k, aIndex ) * mNablaXi[ 1 ]
                        + mHeta( k, aIndex ) * mNablaEta[ 1 ] ;
            }
        }


//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_C()
        {
            mC( 0, 0 ) = mS[ 0 ] * (
                      mNablaXi[ 0 ]  * mExi[  1 ]
                    - mNablaXi[ 1 ]  * mExi[  0 ]
                    + mNablaEta[ 0 ] * mEeta[  1 ]
                    - mNablaEta[ 1 ] * mEeta[  0 ] ) ;

            mC( 0, 1 ) = mS[ 0 ] * (
                      mNablaXi[ 0 ]  * mExi[  3 ]
                    - mNablaXi[ 1 ]  * mExi[  2 ]
                    + mNablaEta[ 0 ] * mEeta[  3 ]
                    - mNablaEta[ 1 ] * mEeta[  2 ] ) ;

            mC( 0, 2 ) = mS[ 1 ] * (
                      mNablaXi[ 0 ]  * mExi[  5 ]
                    - mNablaXi[ 1 ]  * mExi[  4 ]
                    + mNablaEta[ 0 ] * mEeta[  5 ]
                    - mNablaEta[ 1 ] * mEeta[  4 ] ) ;

            mC( 0, 3 ) = mS[ 1 ] * (
                      mNablaXi[ 0 ]  * mExi[  7 ]
                    - mNablaXi[ 1 ]  * mExi[  6 ]
                    + mNablaEta[ 0 ] * mEeta[  7 ]
                    - mNablaEta[ 1 ] * mEeta[  6 ] ) ;

            mC( 0, 4 ) = mS[ 2 ] * (
                      mNablaXi[ 0 ]  * mExi[  9 ]
                    - mNablaXi[ 1 ]  * mExi[  8 ]
                    + mNablaEta[ 0 ] * mEeta[  9 ]
                    - mNablaEta[ 1 ] * mEeta[  8 ] ) ;

            mC( 0, 5 ) = mS[ 2 ] * (
                      mNablaXi[ 0 ]  * mExi[ 11 ]
                    - mNablaXi[ 1 ]  * mExi[ 10 ]
                    + mNablaEta[ 0 ] * mEeta[ 11 ]
                    - mNablaEta[ 1 ] * mEeta[ 10 ] ) ;

            mC( 0, 6 ) =
                      mNablaXi[ 0 ]  * mExi[ 13 ]
                    - mNablaXi[ 1 ]  * mExi[ 12 ]
                    + mNablaEta[ 0 ] * mEeta[ 13 ]
                    - mNablaEta[ 1 ] * mEeta[ 12 ] ;

            mC( 0, 7 ) =
                      mNablaXi[ 0 ]  * mExi[ 15 ]
                    - mNablaXi[ 1 ]  * mExi[ 14 ]
                    + mNablaEta[ 0 ] * mEeta[ 15 ]
                    - mNablaEta[ 1 ] * mEeta[ 14 ] ;
        }

//------------------------------------------------------------------------------
    }
}