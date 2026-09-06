/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "nedelec/cl_EF_HEX8TB.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"

namespace belfem
{
    namespace fem
    {
        EF_HEX8TB::EF_HEX8TB()
        {
            mJ.set_size( 3, 3 );
            mInvJ.set_size( 3, 3 );

            mNumDofs = 4 ;

            // Reference-element volume of the cube [-1,1]^3 with standard
            // Gauss integration, consistent with the "reference-element
            // volume" convention documented in
            // src/fem/interpolation/doc/nedelec.md.
            mSumW = 8.0 ;

            mE.set_size( 3, 4 );
            mC.set_size( 3, 4 );
            mGrad.set_size( 9, 4 );

            mU.set_size( 3 );
            mV.set_size( 3 );
            mW.set_size( 3 );

            mP.set_size( 3 );
            mQ.set_size( 3 );

            for ( uint k=0; k<4; ++k )
            {
                mS[ k ] = 0.0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        EF_HEX8TB::link( Element * aElement )
        {
            // the two thin dimensions are exact physical data from the
            // blocks; the node positions of the wall are partly a drawing
            // proxy and must not feed the metric
            Mesh * tMesh = aElement->parent()->parent()->mesh() ;

            // the width sits on the connector block itself
            mWidth = tMesh->block(
                     aElement->element()->block_id() )->thickness() ;

            // the wall spans exactly one shell layer block; we recover its
            // thickness through the recovery facet, which always carries
            // the id of the wall element plus one and has that block's
            // element as its master. The physical tag is NOT used here:
            // the material machinery overwrites it
            mThickness = tMesh->block(
                     tMesh->facet( aElement->element()->id() + 1 )
                          ->master()->block_id() )->thickness() ;

            BELFEM_ASSERT( mWidth   > 0.0, "EF_HEX8TB: width not set correctly" );
            BELFEM_ASSERT( mThickness > 0.0, "EF_HEX8TB: thickness not set correctly" );

            aElement->edge_directions( mS );

            // midpoints of the two cross-section faces
            mU.fill( 0 );
            mV.fill( 0 );
            aElement->element()->node( 0 )->get_coords( mW ); mU += mW ;
            aElement->element()->node( 1 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 2 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 3 )->get_coords( mW ); mU += mW ;
            aElement->element()->node( 4 )->get_coords( mW ); mU += mW ;
            aElement->element()->node( 5 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 6 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 7 )->get_coords( mW ); mU += mW ;
            mU *= 0.25 ;
            mV *= 0.25 ;

            // unit tangent t along +xi, length from the midpoint chord
            mW = mV - mU ;
            mLength = norm( mW );
            BELFEM_ASSERT( mLength > 0.0, "EF_HEX8TB: element %lu has zero length",
                ( long unsigned int ) aElement->element()->id() );
            mU = mW ;
            mU /= mLength ;

            // raw layer direction from the midpoints of the bottom and top
            // faces: mid{4,5,6,7} - mid{0,1,2,3}
            mV.fill( 0 );
            aElement->element()->node( 0 )->get_coords( mW ); mV -= mW ;
            aElement->element()->node( 1 )->get_coords( mW ); mV -= mW ;
            aElement->element()->node( 2 )->get_coords( mW ); mV -= mW ;
            aElement->element()->node( 3 )->get_coords( mW ); mV -= mW ;
            aElement->element()->node( 4 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 5 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 6 )->get_coords( mW ); mV += mW ;
            aElement->element()->node( 7 )->get_coords( mW ); mV += mW ;
            mV *= 0.25 ;

            // unit normal n: orthogonalize against t, then normalize
            mV -= dot( mV, mU ) * mU ;
            real tNorm = norm( mV );
            BELFEM_ASSERT( tNorm > 0.0,
                "EF_HEX8TB: element %lu is degenerate, layer direction parallel to tangent",
                ( long unsigned int ) aElement->element()->id() );
            mV /= tNorm ;

            // unit binormal b = n x t, so that ( t, b, n ) is right-handed
            mW = cross( mV, mU );

            // geometry Jacobian, rows = d x / d ( xi, eta, zeta );
            // the frame is exact: t scales with the length, b with the
            // width, n with the thickness
            for ( uint i=0; i<3; ++i )
            {
                mJ( 0, i ) = 0.5 * mLength    * mU( i );
                mJ( 1, i ) = 0.5 * mWidth     * mW( i );
                mJ( 2, i ) = 0.5 * mThickness * mV( i );
            }

            // Nabla, columns = nabla xi, nabla eta, nabla zeta;
            // analytic inverse since the frame is orthonormal
            for ( uint i=0; i<3; ++i )
            {
                mInvJ( i, 0 ) = 2.0 * mU( i ) / mLength ;
                mInvJ( i, 1 ) = 2.0 * mW( i ) / mWidth ;
                mInvJ( i, 2 ) = 2.0 * mV( i ) / mThickness ;
            }

            mCrossSection = mWidth * mThickness ;
            mVolume = mCrossSection * mLength ;

            mDetJ = 0.125 * mVolume ;
            mAbsDetJ = mDetJ ; // positive by construction

            // constant curl factors:
            // curl( F_k nabla xi ) = dF_k/deta * ( nabla eta x nabla xi )
            //                      + dF_k/dzeta * ( nabla zeta x nabla xi )
            // with nabla eta x nabla xi = -4/(w*L) * n
            // and nabla zeta x nabla xi =  4/(d*L) * b
            mP = mV * ( -4.0 / ( mWidth * mLength ) );
            mQ = mW * (  4.0 / ( mThickness * mLength ) );
        }

//------------------------------------------------------------------------------

        void
        EF_HEX8TB::precompute( const Matrix< real > & aXi )
        {
            uint n = aXi.n_cols();

            mF.set_size( n, {{}} );
            mFxi.set_size( n, {{}} );

            for ( uint k=0; k<n; ++k )
            {
                real  eta = aXi( 1, k );
                real zeta = aXi( 2, k );

                Vector< real > & F = mF( k );

                F.set_size( mNumDofs );

                // Scale 1/8 = 1/4 (reference-QUAD Nédélec factor for unit edge
                // circulation on [-1,1]^2) × 1/2 (through-thickness ζ extrusion
                // over [-1,1]). Together they yield ∫ E_k · dℓ_j = s_k · δ_{jk}.
                F( 0 ) =   ( 1. - eta ) * ( 1. - zeta );
                F( 1 ) =   ( 1. + eta ) * ( 1. - zeta );
                F( 2 ) =   ( 1. - eta ) * ( 1. + zeta );
                F( 3 ) =   ( 1. + eta ) * ( 1. + zeta );

                F *= 0.125;

                Matrix< real > & dF = mFxi( k );
                dF.set_size( 3, mNumDofs, 0.0 );

                dF( 1, 0 ) = zeta - 1. ;
                dF( 2, 0 ) = eta - 1. ;

                dF( 1, 1 ) = 1. - zeta ;
                dF( 2, 1 ) = -1. - eta ;

                dF( 1, 2 ) = -1. - zeta ;
                dF( 2, 2 ) = 1. - eta ;

                dF( 1, 3 ) = 1. + zeta ;
                dF( 2, 3 ) = 1. + eta ;

                dF *= 0.125;
            }
        }

//------------------------------------------------------------------------------

        const Matrix <real> &
        EF_HEX8TB::E( const uint aIndex )
        {
            // all four dofs point along the curve:
            // E_k = s_k * F_k( eta, zeta ) * nabla xi
            const Vector< real > & F = mF( aIndex );

            mE( 0, 0 ) = mS[ 0 ] * F( 0 ) * mInvJ( 0, 0 );
            mE( 1, 0 ) = mS[ 0 ] * F( 0 ) * mInvJ( 1, 0 );
            mE( 2, 0 ) = mS[ 0 ] * F( 0 ) * mInvJ( 2, 0 );

            mE( 0, 1 ) = mS[ 1 ] * F( 1 ) * mInvJ( 0, 0 );
            mE( 1, 1 ) = mS[ 1 ] * F( 1 ) * mInvJ( 1, 0 );
            mE( 2, 1 ) = mS[ 1 ] * F( 1 ) * mInvJ( 2, 0 );

            mE( 0, 2 ) = mS[ 2 ] * F( 2 ) * mInvJ( 0, 0 );
            mE( 1, 2 ) = mS[ 2 ] * F( 2 ) * mInvJ( 1, 0 );
            mE( 2, 2 ) = mS[ 2 ] * F( 2 ) * mInvJ( 2, 0 );

            mE( 0, 3 ) = mS[ 3 ] * F( 3 ) * mInvJ( 0, 0 );
            mE( 1, 3 ) = mS[ 3 ] * F( 3 ) * mInvJ( 1, 0 );
            mE( 2, 3 ) = mS[ 3 ] * F( 3 ) * mInvJ( 2, 0 );

            return mE ;
        }

//------------------------------------------------------------------------------

        const Matrix <real> &
        EF_HEX8TB::C( const uint aIndex )
        {
            // curl( E_k ) = s_k * ( dF_k/deta * mP + dF_k/dzeta * mQ );
            // the current has no tangential component: the wall carries the
            // commutation current between the layers, transport current
            // along the tape lives in the thin shell
            const Matrix< real > & dF = mFxi( aIndex );

            mC( 0, 0 ) = mS[ 0 ] * ( dF( 1, 0 ) * mP( 0 ) + dF( 2, 0 ) * mQ( 0 ) );
            mC( 1, 0 ) = mS[ 0 ] * ( dF( 1, 0 ) * mP( 1 ) + dF( 2, 0 ) * mQ( 1 ) );
            mC( 2, 0 ) = mS[ 0 ] * ( dF( 1, 0 ) * mP( 2 ) + dF( 2, 0 ) * mQ( 2 ) );

            mC( 0, 1 ) = mS[ 1 ] * ( dF( 1, 1 ) * mP( 0 ) + dF( 2, 1 ) * mQ( 0 ) );
            mC( 1, 1 ) = mS[ 1 ] * ( dF( 1, 1 ) * mP( 1 ) + dF( 2, 1 ) * mQ( 1 ) );
            mC( 2, 1 ) = mS[ 1 ] * ( dF( 1, 1 ) * mP( 2 ) + dF( 2, 1 ) * mQ( 2 ) );

            mC( 0, 2 ) = mS[ 2 ] * ( dF( 1, 2 ) * mP( 0 ) + dF( 2, 2 ) * mQ( 0 ) );
            mC( 1, 2 ) = mS[ 2 ] * ( dF( 1, 2 ) * mP( 1 ) + dF( 2, 2 ) * mQ( 1 ) );
            mC( 2, 2 ) = mS[ 2 ] * ( dF( 1, 2 ) * mP( 2 ) + dF( 2, 2 ) * mQ( 2 ) );

            mC( 0, 3 ) = mS[ 3 ] * ( dF( 1, 3 ) * mP( 0 ) + dF( 2, 3 ) * mQ( 0 ) );
            mC( 1, 3 ) = mS[ 3 ] * ( dF( 1, 3 ) * mP( 1 ) + dF( 2, 3 ) * mQ( 1 ) );
            mC( 2, 3 ) = mS[ 3 ] * ( dF( 1, 3 ) * mP( 2 ) + dF( 2, 3 ) * mQ( 2 ) );

            return mC ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_HEX8TB::G( const uint aIndex )
        {
            // grad( e_k ) = s_k * ( grad F_k ) (x) ( nabla xi ) — the frame
            // is an element constant from link(); F has no xi dependence,
            // so grad F chains from the eta and zeta rows only
            const Matrix< real > & dF = mFxi( aIndex );

            for ( uint k = 0; k < 4; ++k )
            {
                for ( uint j = 0; j < 3; ++j )       // field component
                {
                    for ( uint i = 0; i < 3; ++i )   // derivative direction
                    {
                        mGrad( i + 3 * j, k ) = mS[ k ] *
                                ( mInvJ( i, 1 ) * dF( 1, k )
                                + mInvJ( i, 2 ) * dF( 2, k ) ) * mInvJ( j, 0 );
                    }
                }
            }

            return mGrad ;
        }

    }
}
