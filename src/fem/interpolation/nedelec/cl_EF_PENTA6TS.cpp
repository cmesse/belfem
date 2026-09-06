//
// Created by christian on 8/20/25.
//
//
// Created by christian on 12/3/21.
//

#include "nedelec/cl_EF_PENTA6TS.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "fn_det.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_PENTA6TS::EF_PENTA6TS()
        {
            mJ.set_size( 2, 3 );
            mE.set_size( 3, 6 );
            mC.set_size( 3, 6 );
            mGrad.set_size( 9, 6 );

            mGram.set_size( 2, 2 );
            mInvGram.set_size( 2, 2 );
            mPseudoInvJ.set_size( 3, 2 );
            mInvJ.set_size( 3, 3 );
            mCoeffs.set_size( 3, 3 );

            mNumDofs = 6 ;
            mSumW = 1.0 ;
        }

//------------------------------------------------------------------------------

        void
        EF_PENTA6TS::link( Element * aElement )
        {
            // grab nodes
            const mesh::Node * tNode0 = aElement->facet()->node( 0 );
            const mesh::Node * tNode1 = aElement->facet()->node( 1 );
            const mesh::Node * tNode2 = aElement->facet()->node( 2 );

            // compute jacobian
            mJ( 0, 0 ) = tNode0->x() - tNode2->x();
            mJ( 1, 0 ) = tNode1->x() - tNode2->x();

            mJ( 0, 1 ) = tNode0->y() - tNode2->y();
            mJ( 1, 1 ) = tNode1->y() - tNode2->y();

            mJ( 0, 2 ) = tNode0->z() - tNode2->z();
            mJ( 1, 2 ) = tNode1->z() - tNode2->z();


            // compute the gram matrix
            mGram = mJ * trans( mJ );

            // invert the gram matrix in place
            inv2( mGram, mInvGram );

            // compute the pseudoinverse of J
            mPseudoInvJ = trans( mJ ) * mInvGram;

            // compute the determinant
            mDetGram = det( mGram );

            // compute the normal
            mN[ 0 ] = mJ( 0, 1 ) * mJ( 1, 2 ) - mJ( 0, 2 ) * mJ( 1, 1 );
            mN[ 1 ] = mJ( 0, 2 ) * mJ( 1, 0 ) - mJ( 0, 0 ) * mJ( 1, 2 );
            mN[ 2 ] = mJ( 0, 0 ) * mJ( 1, 1 ) - mJ( 0, 1 ) * mJ( 1, 0 );

            // doubled triangle surface
            real tNorm = std::sqrt( mN[ 0 ] * mN[ 0 ] + mN[ 1 ] * mN[ 1 ] + mN[ 2 ] * mN[ 2 ] );

            // normalize
            mN[ 0 ] /= tNorm ;
            mN[ 1 ] /= tNorm ;
            mN[ 2 ] /= tNorm ;

            // actual triangle surface
            mSurface = 0.5 * tNorm ;

            // grab edge signs and store them in container
            aElement->edge_directions( mS );

            // compute nabla values
            mNablaXi[ 0 ]   = mPseudoInvJ( 0, 0 ) ;
            mNablaXi[ 1 ]   = mPseudoInvJ( 1, 0 ) ;
            mNablaXi[ 2 ]   = mPseudoInvJ( 2, 0 ) ;

            mNablaEta[ 0 ]  = mPseudoInvJ( 0, 1 ) ;
            mNablaEta[ 1 ]  = mPseudoInvJ( 1, 1 ) ;
            mNablaEta[ 2 ]  = mPseudoInvJ( 2, 1 ) ;

            mNablaZeta[ 0 ] = - mNablaXi[ 0 ] - mNablaEta[ 0 ] ;
            mNablaZeta[ 1 ] = - mNablaXi[ 1 ] - mNablaEta[ 1 ] ;
            mNablaZeta[ 2 ] = - mNablaXi[ 2 ] - mNablaEta[ 2 ] ;

            // constant in-plane gradient tensor for G(), layout i + 3*j;
            // W = nXi (x) nEta - nEta (x) nXi, identical for all three edge
            // pairs because nZeta = -nXi - nEta. Deliberately placed here,
            // outside the block below where tJ aliases mCoeffs
            for ( uint tCol = 0; tCol < 3; ++tCol )        // field component
            {
                for ( uint tRow = 0; tRow < 3; ++tRow )    // derivative direction
                {
                    mGradW[ tRow + 3 * tCol ] =
                          mNablaXi[ tRow ] * mNablaEta[ tCol ]
                        - mNablaEta[ tRow ] * mNablaXi[ tCol ];
                }
            }

            // try to get the exact thickness from the mesh
            mThickness = aElement->parent()->parent()->mesh()->block( aElement->element()->block_id() )->thickness();

            if ( std::isnan( mThickness ) )
            {
               // if the exact thickness is not set, maybe because this is not a maxwell problem, we can compute
               // the thickness from the element on the mesh

                // center of first triangle
                real tX = aElement->element()->node( 0 )->x() + aElement->element()->node( 1 )->x() + aElement->element()->node( 2 )->x();
                real tY = aElement->element()->node( 0 )->y() + aElement->element()->node( 1 )->y() + aElement->element()->node( 2 )->y();
                real tZ = aElement->element()->node( 0 )->z() + aElement->element()->node( 1 )->z() + aElement->element()->node( 2 )->z();

                // minus center of second triangle
                tX -= aElement->element()->node( 3 )->x() + aElement->element()->node( 4 )->x() + aElement->element()->node( 5 )->x();
                tY -= aElement->element()->node( 3 )->y() + aElement->element()->node( 4 )->y() + aElement->element()->node( 5 )->y();
                tZ -= aElement->element()->node( 3 )->z() + aElement->element()->node( 4 )->z() + aElement->element()->node( 5 )->z();

                // thickness is the distance between both points
                mThickness = std::sqrt( tX*tX + tY*tY + tZ*tZ ) / 3.0 ;
            }



            // termporarily use mCoeffs as work matrix
            Matrix< real > & tJ = mCoeffs ;

            // the first two rows are identical to the rows in the reduced Jacobian
            tJ.set_row( 0, mJ.row( 0 ) );
            tJ.set_row( 1, mJ.row( 1 ) );

            // the final row is built using the thickness matrix (it refers to the tau-parameter)
            tJ( 2, 0 ) = 0.5 * mThickness * mN[ 0 ] ;
            tJ( 2, 1 ) = 0.5 * mThickness * mN[ 1 ] ;
            tJ( 2, 2 ) = 0.5 * mThickness * mN[ 2 ] ;

            // coompute the inverse in place
            inv3( tJ, mInvJ );

            // with the inverse, we can now compute the coefficients, which are fixed
            // for a linear element
            mCoeffs( 0, 0 ) = mNablaEta[1]*mInvJ(2,2)
                                                -mNablaEta[2]*mInvJ(1,2);

            mCoeffs( 1, 0 ) = mNablaEta[2]*mInvJ(0,2)
                                                -mNablaEta[0]*mInvJ(2,2);

            mCoeffs( 2, 0 ) = mNablaEta[0]*mInvJ(1,2)
                                                -mNablaEta[1]*mInvJ(0,2);


            mCoeffs( 0, 1 ) = mNablaXi[2]*mInvJ(1,2)
                                                -mNablaXi[1]*mInvJ(2,2);

            mCoeffs( 1, 1 ) = mNablaXi[0]*mInvJ(2,2)
                                                -mNablaXi[2]*mInvJ(0,2);

            mCoeffs( 2, 1 ) = mNablaXi[1]*mInvJ(0,2)
                                                -mNablaXi[0]*mInvJ(1,2);


            mCoeffs( 0, 2 ) = mNablaXi[2]*mInvJ(1,1)
                                                +mNablaEta[1]*mInvJ(2,0)
                                                -mNablaXi[1]*mInvJ(2,1)
                                                -mNablaEta[2]*mInvJ(1,0) ;

            mCoeffs( 1, 2 ) =  mNablaEta[2]*mInvJ(0,0)
                                                 +mNablaXi[0]*mInvJ(2,1)
                                                 -mNablaEta[0]*mInvJ(2,0)
                                                 -mNablaXi[2]*mInvJ(0,1);

            mCoeffs( 2, 2 ) = mNablaXi[1]*mInvJ(0,1)
                                                +mNablaEta[0]*mInvJ(1,0)
                                                -mNablaEta[1]*mInvJ(0,0)
                                                -mNablaXi[0]*mInvJ(1,1);

            // afterwards, we need to scale
            mCoeffs *= 0.5 ;

            // in the code, we use mDetJ for the volume increment
            // since the sum of all weights in a wedge 1, this is actually the element volume
            mDetJ = mThickness * mSurface ;

            // compute the absolute value
            mAbsDetJ = std::abs( mDetJ );
        }

//------------------------------------------------------------------------------

        void
        EF_PENTA6TS::precompute( const Matrix< real > & aXi )
        {
            // get number of integration points
            uint tN = aXi.n_cols();

            // compute factor for thickness position
            mF.set_size( 2, tN );
            for( uint k=0; k<tN; ++k )
            {
                mF( 0, k ) = 0.5 * ( 1.-aXi( 2, k ) );
                mF( 1, k ) = 0.5 * ( 1.+aXi( 2, k ) );
            }

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


            // parameter coordinates for
            mCurlPars.set_size( tN, {} );
            for ( uint k=0; k<tN; ++k )
            {
                real xi = aXi( 0, k );
                real eta = aXi( 1, k );
                real tau = aXi( 2, k );

                Matrix< real > & B = mCurlPars( k );
                B.set_size( 3, 6 );

                B( 0, 0 ) = xi;
                B( 1, 0 ) = eta;
                B( 2, 0 ) = tau-1.0;

                B( 0, 1 ) = xi-1.0;
                B( 1, 1 ) = eta;
                B( 2, 1 ) = tau-1.0;

                B( 0, 2 ) = xi;
                B( 1, 2 ) = eta-1.0;
                B( 2, 2 ) = tau-1.0;

                B( 0, 3 ) = -xi;
                B( 1, 3 ) = -eta;
                B( 2, 3 ) = -tau-1.0;

                B( 0, 4 ) = 1.0-xi;
                B( 1, 4 ) = -eta;
                B( 2, 4 ) = -tau-1.0;

                B( 0, 5 ) = -xi;
                B( 1, 5 ) = 1.0-eta;
                B( 2, 5 ) = -tau-1.0;
            }
        }

//------------------------------------------------------------------------------

        const Matrix <real> &
        EF_PENTA6TS::E( const uint aIndex )
        {
            // edge values without size and scaling

            // edge 0
            mE( 0, 0 ) = mG( 0, aIndex ) * mNablaEta[ 0 ] - mH( 0, aIndex ) * mNablaXi[ 0 ];
            mE( 1, 0 ) = mG( 0, aIndex ) * mNablaEta[ 1 ] - mH( 0, aIndex ) * mNablaXi[ 1 ];
            mE( 2, 0 ) = mG( 0, aIndex ) * mNablaEta[ 2 ] - mH( 0, aIndex ) * mNablaXi[ 2 ];

            // edge 1
            mE( 0, 1 ) = mG( 1, aIndex ) * mNablaZeta[ 0 ] - mH( 1, aIndex ) * mNablaEta[ 0 ];
            mE( 1, 1 ) = mG( 1, aIndex ) * mNablaZeta[ 1 ] - mH( 1, aIndex ) * mNablaEta[ 1 ];
            mE( 2, 1 ) = mG( 1, aIndex ) * mNablaZeta[ 2 ] - mH( 1, aIndex ) * mNablaEta[ 2 ];

            // edge 2
            mE( 0, 2 ) = mG( 2, aIndex ) * mNablaXi[ 0 ] - mH( 2, aIndex ) * mNablaZeta[ 0 ];
            mE( 1, 2 ) = mG( 2, aIndex ) * mNablaXi[ 1 ] - mH( 2, aIndex ) * mNablaZeta[ 1 ];
            mE( 2, 2 ) = mG( 2, aIndex ) * mNablaXi[ 2 ] - mH( 2, aIndex ) * mNablaZeta[ 2 ];

            // duplicate and finalize upper face
            mE( 0, 3 ) = mS[ 3 ] * mE( 0, 0 ) * mF( 1, aIndex );
            mE( 1, 3 ) = mS[ 3 ] * mE( 1, 0 ) * mF( 1, aIndex );
            mE( 2, 3 ) = mS[ 3 ] * mE( 2, 0 ) * mF( 1, aIndex );

            mE( 0, 4 ) = mS[ 4 ] * mE( 0, 1 ) * mF( 1, aIndex );
            mE( 1, 4 ) = mS[ 4 ] * mE( 1, 1 ) * mF( 1, aIndex );
            mE( 2, 4 ) = mS[ 4 ] * mE( 2, 1 ) * mF( 1, aIndex );

            mE( 0, 5 ) = mS[ 5 ] * mE( 0, 2 ) * mF( 1, aIndex );
            mE( 1, 5 ) = mS[ 5 ] * mE( 1, 2 ) * mF( 1, aIndex );
            mE( 2, 5 ) = mS[ 5 ] * mE( 2, 2 ) * mF( 1, aIndex );

            // finalize lower face
            mE( 0, 0 ) *= mS[ 0 ] * mF( 0, aIndex );
            mE( 1, 0 ) *= mS[ 0 ] * mF( 0, aIndex );
            mE( 2, 0 ) *= mS[ 0 ] * mF( 0, aIndex );

            mE( 0, 1 ) *= mS[ 1 ] * mF( 0, aIndex );
            mE( 1, 1 ) *= mS[ 1 ] * mF( 0, aIndex );
            mE( 2, 1 ) *= mS[ 1 ] * mF( 0, aIndex );

            mE( 0, 2 ) *= mS[ 2 ] * mF( 0, aIndex );
            mE( 1, 2 ) *= mS[ 2 ] * mF( 0, aIndex );
            mE( 2, 2 ) *= mS[ 2 ] * mF( 0, aIndex );

            return mE ;
        }

        const Matrix <real> &
        EF_PENTA6TS::C( const uint aIndex )
        {
            // compute the coefficients
            mC = mCoeffs * mCurlPars( aIndex );

            // adjust for edge orientations
            for ( uint k=0; k<6; ++k)
            {
                mC( 0, k ) *= mS[ k ] ;
                mC( 1, k ) *= mS[ k ] ;
                mC( 2, k ) *= mS[ k ] ;
            }

            return mC ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_PENTA6TS::G( const uint aIndex )
        {
            // unsigned w_m at this point, same pair table as E() —
            // recomputed here on purpose: calling E() would clobber mE, and
            // recovering w from mE fails on the faces where F vanishes
            real tW[ 3 ][ 3 ];
            for ( uint i = 0; i < 3; ++i )
            {
                tW[ 0 ][ i ] = mG( 0, aIndex ) * mNablaEta[ i ]
                             - mH( 0, aIndex ) * mNablaXi[ i ];
                tW[ 1 ][ i ] = mG( 1, aIndex ) * mNablaZeta[ i ]
                             - mH( 1, aIndex ) * mNablaEta[ i ];
                tW[ 2 ][ i ] = mG( 2, aIndex ) * mNablaXi[ i ]
                             - mH( 2, aIndex ) * mNablaZeta[ i ];
            }

            // thickness weights and their gradient contribution:
            // grad( e_m )     = s_m     ( F0 * W - 1/2 nTau (x) w_m )
            // grad( e_{m+3} ) = s_{m+3} ( F1 * W + 1/2 nTau (x) w_m )
            // with nTau_i = mInvJ( i, 2 ), the same column mCoeffs uses
            const real tF0 = mF( 0, aIndex );
            const real tF1 = mF( 1, aIndex );

            for ( uint tEdge = 0; tEdge < 3; ++tEdge )
            {
                for ( uint tCol = 0; tCol < 3; ++tCol )        // field component
                {
                    for ( uint tRow = 0; tRow < 3; ++tRow )    // derivative direction
                    {
                        const real tWij = mGradW[ tRow + 3 * tCol ];
                        const real tTau = 0.5 * mInvJ( tRow, 2 )
                                        * tW[ tEdge ][ tCol ];

                        mGrad( tRow + 3 * tCol, tEdge ) =
                            mS[ tEdge ] * ( tF0 * tWij - tTau );
                        mGrad( tRow + 3 * tCol, tEdge + 3 ) =
                            mS[ tEdge + 3 ] * ( tF1 * tWij + tTau );
                    }
                }
            }
            return mGrad ;
        }

    }
}