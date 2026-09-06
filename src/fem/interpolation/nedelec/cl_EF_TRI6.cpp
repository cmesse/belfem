//
// Created by christian on 12/2/21.
//

#include "nedelec/cl_EF_TRI6.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "fn_dot.hpp"
#include "fn_det.hpp"
#include "fn_inv2.hpp"

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
            mGrad.set_size( 4, 8 );

            mNumDofs = 8 ;
            mSumW = 0.5 ;
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
        EF_TRI6::link( Element  * aElement )
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

            // get the directions of the edges
            aElement->edge_directions( mS );

            // check if element is curved
            if( tElement->is_curved() )
            {
                // link functions
                mFunInterpolation = & EF_TRI6::E_curved ;
                mFunCurl          = & EF_TRI6::C_curved ;
                mFunGrad          = & EF_TRI6::G_curved ;

                // contracted second derivatives of the geometry map,
                // constant because the map is quadratic; table from
                // cl_IF_TRI6.hpp d2NdXi2, IF-Voigt rows ( xx, yy, xy )
                static const real tD2[ 3 ][ 6 ] = {
                        {  4.,  0.,  4.,  0.,  0., -8. },
                        {  0.,  4.,  4.,  0., -8.,  0. },
                        {  0.,  0.,  4.,  4., -4., -4. } };
                for ( uint v = 0; v < 3; ++v )
                {
                    mCurv[ v ][ 0 ] = 0.0 ;
                    mCurv[ v ][ 1 ] = 0.0 ;
                    for ( uint n = 0; n < 6; ++n )
                    {
                        mCurv[ v ][ 0 ] += tD2[ v ][ n ] * mX( n );
                        mCurv[ v ][ 1 ] += tD2[ v ][ n ] * mY( n );
                    }
                }
            }
            else
            {
                // compute nabla ( it is constant for this element )
                this->compute_nabla( 0 );

                // link functions
                mFunInterpolation = & EF_TRI6::compute_E ;
                mFunCurl          = & EF_TRI6::C_straight ;
                mFunGrad          = & EF_TRI6::G_straight ;

                mInvJ.fill( BELFEM_QUIET_NAN );
                mJ.fill( BELFEM_QUIET_NAN );

                // the hess channel must never be read on this branch
                // ( S = 0 on the affine map, and mJ above is NaN ) — NaN
                // rather than uninitialized so a missed skip screams
                for ( uint v = 0; v < 3; ++v )
                {
                    mCurv[ v ][ 0 ] = BELFEM_QUIET_NAN ;
                    mCurv[ v ][ 1 ] = BELFEM_QUIET_NAN ;
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
            this->compute_E_xi( aIndex );
            this->compute_E_eta( aIndex );
            this->compute_C();
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::C_straight( const uint aIndex )
        {
            this->compute_E_xi( aIndex );
            this->compute_E_eta( aIndex );
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

                // inverse, written in place
                inv2( mJ, mInvJ );

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
        EF_TRI6::compute_E_xi( const uint aIndex )
        {
            uint tCount = 0 ;

            for( uint k=0; k<8; ++k )
            {
                mExi[  tCount++ ] =
                      mGxi( k, aIndex ) * mNablaXi[ 0 ]
                    + mHxi( k, aIndex ) * mNablaEta[ 0 ];

                mExi[  tCount++ ] =
                      mGxi( k, aIndex ) * mNablaXi[ 1 ]
                    + mHxi( k, aIndex ) * mNablaEta[ 1 ];
            }

        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::compute_E_eta( const uint aIndex )
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

        void
        EF_TRI6::G_term1( const uint aIndex )
        {
            // reference partials with the nablas held fixed — the same
            // tables compute_C chains; they are UNSIGNED for this class,
            // so the edge signs are applied here ( interior dofs 6, 7
            // carry no sign )
            this->compute_E_xi( aIndex );
            this->compute_E_eta( aIndex );

            for ( uint k = 0; k < 8; ++k )
            {
                const real tSig = k < 6 ? mS[ k / 2 ] : 1.0 ;
                for ( uint j = 0; j < 2; ++j )        // field component
                {
                    for ( uint i = 0; i < 2; ++i )    // derivative direction
                    {
                        mGrad( i + 2 * j, k ) = tSig *
                              ( mExi[ 2 * k + j ]  * mNablaXi[ i ]
                              + mEeta[ 2 * k + j ] * mNablaEta[ i ] );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::G_straight( const uint aIndex )
        {
            // frozen nablas from link(); the affine map has S = 0, so the
            // hess channel does not exist here and MUST not be touched:
            // mCurv and mJ are NaN on this branch by design
            this->G_term1( aIndex );
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6::G_curved( const uint aIndex )
        {
            this->compute_nabla( aIndex );
            this->G_term1( aIndex );

            // E at this point for the Q-folded hess channel ( refreshes
            // mE — the documented curved-path contract )
            this->compute_E( aIndex );

            // the two inverse-map Hessians H^(d) = -A * Cm_d * A^T,
            // dof-independent at this point; A columns are the nabla
            // vectors — NEVER the mJ / mInvJ members
            const real * tA[ 2 ] = { mNablaXi, mNablaEta };
            real tH[ 2 ][ 4 ];
            for ( uint d = 0; d < 2; ++d )
            {
                // Cm entries from the IF-Voigt rows ( xx, yy, xy ),
                // contracted with nabla_d — the diagonal pairs are LIVE
                // for the quadratic map
                const real tCxx = mCurv[ 0 ][ 0 ] * tA[ d ][ 0 ]
                                + mCurv[ 0 ][ 1 ] * tA[ d ][ 1 ];
                const real tCyy = mCurv[ 1 ][ 0 ] * tA[ d ][ 0 ]
                                + mCurv[ 1 ][ 1 ] * tA[ d ][ 1 ];
                const real tCxy = mCurv[ 2 ][ 0 ] * tA[ d ][ 0 ]
                                + mCurv[ 2 ][ 1 ] * tA[ d ][ 1 ];

                for ( uint j = 0; j < 2; ++j )        // field component
                {
                    for ( uint i = 0; i < 2; ++i )    // derivative direction
                    {
                        tH[ d ][ i + 2 * j ] = -(
                              tA[ 0 ][ i ] * ( tCxx * tA[ 0 ][ j ]
                                             + tCxy * tA[ 1 ][ j ] )
                            + tA[ 1 ][ i ] * ( tCxy * tA[ 0 ][ j ]
                                             + tCyy * tA[ 1 ][ j ] ) );
                    }
                }
            }

            // Q-folded channel: Q^d = sum_i mW[ 8 + 2*d + i ] * e_i —
            // compute_nabla's own Jacobian dump, fresh for aIndex
            for ( uint k = 0; k < 8; ++k )
            {
                real tQ[ 2 ];
                for ( uint d = 0; d < 2; ++d )
                {
                    tQ[ d ] = mW[ 8 + 2 * d ]     * mE( 0, k )
                            + mW[ 8 + 2 * d + 1 ] * mE( 1, k );
                }
                for ( uint r = 0; r < 4; ++r )
                {
                    mGrad( r, k ) += tQ[ 0 ] * tH[ 0 ][ r ]
                                   + tQ[ 1 ] * tH[ 1 ][ r ];
                }
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TRI6::G( const uint aIndex )
        {
            ( this->*mFunGrad )( aIndex );
            return mGrad ;
        }

    }
}