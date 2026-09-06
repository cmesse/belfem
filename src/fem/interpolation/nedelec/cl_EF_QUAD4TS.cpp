//
// Created by claude on 1/11/25.
//

#include "nedelec/cl_EF_QUAD4TS.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_QUAD4TS::EF_QUAD4TS()
        {
            mJ.set_size( 1, 2 );
            mE.set_size( 2, 2 );
            mC.set_size( 1, 2 );
            mGrad.set_size( 4, 2 );

            mGram.set_size( 1, 1 );
            mPseudoInvJ.set_size( 2, 1 );
            mCoeffs.set_size( 1, 1 );

            mNumDofs = 2 ;
            mSumW = 4.0 ;
        }

//------------------------------------------------------------------------------

        void
        EF_QUAD4TS::link( Element * aElement )
        {
            // grab nodes from the facet (LINE2 element)
            const mesh::Node * tNode0 = aElement->facet()->node( 0 );
            const mesh::Node * tNode1 = aElement->facet()->node( 1 );

            // compute jacobian (1x2 matrix for a line in 2D)
            mJ( 0, 0 ) = tNode1->x() - tNode0->x();
            mJ( 0, 1 ) = tNode1->y() - tNode0->y();

            // compute the gram matrix (1x1 for a line)
            mGram = mJ * trans( mJ );

            // using the pseudoinverse here. the gram matrix of a line is 1x1,
            // so J+ = J^T ( J J^T )^-1 collapses to a scalar division
            mPseudoInvJ = trans( mJ ) / mGram( 0, 0 );

            // compute the determinant
            mDetGram = det( mGram );

            // compute the normal (perpendicular to the line in 2D)
            // rotate by 90 degrees: (dx, dy) -> (-dy, dx)
            mN[ 0 ] = -mJ( 0, 1 );
            mN[ 1 ] =  mJ( 0, 0 );

            // tangent length
            real tNorm = std::sqrt( mN[ 0 ] * mN[ 0 ] + mN[ 1 ] * mN[ 1 ] );

            // normalize
            mN[ 0 ] /= tNorm ;
            mN[ 1 ] /= tNorm ;

            // actual line length
            mLength = tNorm ;

            // grab edge signs and store them in container
            aElement->edge_directions( mS );

            // mNablaXi = gradient along the line (Whitney edge function)
            mNablaXi[ 0 ]  = mPseudoInvJ( 0, 0 ) ;
            mNablaXi[ 1 ]  = mPseudoInvJ( 1, 0 ) ;

            // vector from the bottom curve center to the top curve center;
            // defines the through-thickness direction (the stacking side is
            // not guaranteed to be +mN)
            real tDx = 0.5 * ( aElement->element()->node( 2 )->x() + aElement->element()->node( 3 )->x()
                             - aElement->element()->node( 0 )->x() - aElement->element()->node( 1 )->x() );
            real tDy = 0.5 * ( aElement->element()->node( 2 )->y() + aElement->element()->node( 3 )->y()
                             - aElement->element()->node( 0 )->y() - aElement->element()->node( 1 )->y() );

            // geometric layer thickness
            real tD = std::sqrt( tDx*tDx + tDy*tDy );

            BELFEM_ASSERT( tD > 0.0, "Degenerate thin shell element %lu",
                           ( long unsigned int ) aElement->element()->id() );

            // try to get the exact thickness from the mesh
            mThickness = aElement->parent()->parent()->mesh()->block( aElement->element()->block_id() )->thickness();

            if ( std::isnan( mThickness ) )
            {
                // if the exact thickness is not set, maybe because this is not
                // a maxwell problem, we use the geometric thickness
                mThickness = tD ;
            }

            // mNablaEta = gradient of the thickness parameter eta in [-1,+1]
            real tThicknessScale = 2.0 / ( mThickness * tD );
            mNablaEta[ 0 ] = tDx * tThicknessScale ;
            mNablaEta[ 1 ] = tDy * tThicknessScale ;

            // coefficient for the curl: ( nabla eta x nabla xi )_z
            mCoeffs( 0, 0 ) = mNablaEta[ 0 ] * mNablaXi[ 1 ]
                            - mNablaEta[ 1 ] * mNablaXi[ 0 ];

            // gradient operator ( constant for this element, layout in
            // EdgeFunction::mGrad ): e_k = s_k F_k( eta ) w with w = nXi,
            // so grad( e_k ) = s_k F_k' nEta (x) nXi, F' = -+ 1/2.
            // mGrad is sized in the constructor, never resized here
            for ( uint tCol = 0; tCol < 2; ++tCol )        // field component
            {
                for ( uint tRow = 0; tRow < 2; ++tRow )    // derivative direction
                {
                    mGrad( tRow + 2 * tCol, 0 ) = -0.5 * mS[ 0 ]
                        * mNablaEta[ tRow ] * mNablaXi[ tCol ];
                    mGrad( tRow + 2 * tCol, 1 ) =  0.5 * mS[ 1 ]
                        * mNablaEta[ tRow ] * mNablaXi[ tCol ];
                }
            }

            // in the code, we use mDetJ for the volume increment;
            // the quadrature weights on the quad sum to 4, so this is
            // the element volume (area in 2D) divided by 4
            mDetJ = 0.25 * mThickness * mLength ;

            // compute the absolute value
            mAbsDetJ = std::abs( mDetJ );
        }

//------------------------------------------------------------------------------

        void
        EF_QUAD4TS::precompute( const Matrix< real > & aXi )
        {
            // get number of integration points
            uint tN = aXi.n_cols();

            // compute factor for thickness position (eta coordinate)
            mF.set_size( 2, tN );
            for( uint k=0; k<tN; ++k )
            {
                mF( 0, k ) = 0.5 * ( 1.-aXi( 1, k ) );
                mF( 1, k ) = 0.5 * ( 1.+aXi( 1, k ) );
            }

            // Curl parameters are derivatives of thickness weight functions
            mCurlPars.set_size( tN, {} );
            for ( uint k=0; k<tN; ++k )
            {
                Matrix< real > & B = mCurlPars( k );
                B.set_size( 1, 2 );

                // derivatives of the thickness weights
                // f0 = (1-eta)/2 and f1 = (1+eta)/2
                B( 0, 0 ) = -0.5 ;
                B( 0, 1 ) = +0.5 ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix <real> &
        EF_QUAD4TS::E( const uint aIndex )
        {
            // edge values without size and scaling

            // Compute base edge vector (tangent along the line)
            real tEx = mNablaXi[ 0 ];
            real tEy = mNablaXi[ 1 ];

            // Bottom edge (eta=-1)
            mE( 0, 0 ) = mS[ 0 ] * tEx * mF( 0, aIndex );
            mE( 1, 0 ) = mS[ 0 ] * tEy * mF( 0, aIndex );

            // Top edge (eta=+1); same tangent convention as the bottom edge,
            // so an interface edge shared with the next layer keeps the
            // tangential field continuous (matches the 3D thin shells)
            mE( 0, 1 ) = mS[ 1 ] * tEx * mF( 1, aIndex );
            mE( 1, 1 ) = mS[ 1 ] * tEy * mF( 1, aIndex );

            return mE ;
        }

        const Matrix <real> &
        EF_QUAD4TS::C( const uint aIndex )
        {
            // compute the coefficients
            mC = mCoeffs * mCurlPars( aIndex );

            // adjust for edge orientations
            for ( uint k=0; k<2; ++k)
            {
                mC( 0, k ) *= mS[ k ] ;
            }

            return mC ;
        }

//------------------------------------------------------------------------------
    }
}
