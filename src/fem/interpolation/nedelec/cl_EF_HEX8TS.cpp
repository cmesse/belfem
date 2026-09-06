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
#include "nedelec/cl_EF_HEX8TS.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "fn_norm.hpp"
#include "fn_inv2.hpp"
#include "fn_cross.hpp"
#include "fn_trans.hpp"

// Lagrange template specializations live in these headers (without the
// `inline` keyword), so include them only here, not in cl_EF_HEX8TS.hpp,
// to avoid multiple-definition errors when the header is included from
// other translation units.

#include "lagrange/cl_IF_QUAD4.hpp"

namespace belfem
{
    namespace fem
    {
        EF_HEX8TS::EF_HEX8TS()
        {
            mJ.set_size( 3, 3 );
            mInvJ.set_size( 3, 3 );

            mJ2.set_size( 2, 3 );
            mG.set_size( 2, 2  );
            mInvG.set_size( 2, 2 );
            mX.set_size( 4, 3 );
            mNumDofs = 8 ;

            // Reference-element volume of the cube [-1,1]^3 with standard
            // Gauss integration. Currently unused by any assembly site
            // (sum_w() is not consumed downstream), but kept consistent
            // with the "reference-element volume" convention documented
            // in src/fem/interpolation/doc/nedelec.md.
            mSumW = 8.0 ;

            mQuad4 = new InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >;

            mE.set_size( 3, 8 );
            mC.set_size( 3, 8 );
            mGrad.set_size( 9, 8 );

            mNx.set_size( 3, 8, 0. );

            mEx.set_size( 3, 8, 0. );
            mEy.set_size( 3, 8, 0. );
            mEz.set_size( 3, 8, 0. );

            mU.set_size( 3 );
            mV.set_size( 3 );
            mW.set_size( 3 );

            for ( uint k=0; k<8; ++k )
            {
                mS[ k ] = 0.0 ;
            }
        }

        EF_HEX8TS::~EF_HEX8TS()
        {
            if ( mQuad4 != nullptr ) delete mQuad4;
        }

        void
        EF_HEX8TS::link( Element * aElement )
        {
            // thin shell
            mThickness = aElement->parent()->parent()->mesh()->block (
                 aElement->element()->block_id() )->thickness() ;

            BELFEM_ASSERT( mThickness > 0.0, "EF_HEX8TS: thickness not set correctly" );

            // populate node coordinates
            for ( uint k=0; k<4; ++k )
            {
                mesh::Node * tNode = aElement->facet()->node( k );
                mX( k, 0 ) = tNode->x( 0 );
                mX( k, 1 ) = tNode->x( 1 );
                mX( k, 2 ) = tNode->x( 2 );
            }

            aElement->edge_directions( mS );

            // surface computation

            mU( 0 ) = mX( 1,0 )-mX( 0,0 );
            mU( 1 ) = mX( 1,1 )-mX( 0,1 );
            mU( 2 ) = mX( 1,2 )-mX( 0,2 );

            mV( 0 ) = mX( 3,0 )-mX( 0,0 );
            mV( 1 ) = mX( 3,1 )-mX( 0,1 );
            mV( 2 ) = mX( 3,2 )-mX( 0,2 );

            mSurface = norm( cross( mU, mV ) );

            mVolume = mSurface * mThickness;
            mLastIndex = BELFEM_UINT_MAX ;
        }

        void
        EF_HEX8TS::precompute( const Matrix< real > & aXi )
        {
            uint n = aXi.n_cols();

            mF.set_size( n, {{}} );
            mFxi.set_size( n, {{}} );

            mQuad4Nxi.set_size( n, {{}} );


            Vector< real > tXi( 3 );
            for ( uint k=0; k<n; ++k )
            {
                mQuad4->dNdXi( aXi.col( k ), mQuad4Nxi( k ) );

                real   xi = aXi( 0, k );
                real  eta = aXi( 1, k );
                real zeta = aXi( 2, k );

                Vector< real > & F = mF( k );

                F.set_size( 8 );

                // Scale 1/8 = 1/4 (reference-QUAD Nédélec factor for unit edge
                // circulation on [-1,1]^2) × 1/2 (through-thickness ζ extrusion
                // over [-1,1]). Together they yield ∫ E_k · dℓ_j = s_k · δ_{jk}.
                F( 0 ) =   ( 1.0 - eta ) * ( 1. - zeta );
                F( 1 ) =   ( 1.0 + xi )  * ( 1. - zeta );
                F( 2 ) =  -( 1.0 + eta ) * ( 1. - zeta );
                F( 3 ) =  -( 1.0 - xi )  * ( 1. - zeta );

                F( 4 ) =   ( 1.0 - eta ) * ( 1. + zeta );
                F( 5 ) =   ( 1.0 + xi )  * ( 1. + zeta );
                F( 6 ) =  -( 1.0 + eta ) * ( 1. + zeta );
                F( 7 ) =  -( 1.0 - xi )  * ( 1. + zeta );

                F *= 0.125;

                Matrix< real > & dF = mFxi( k );
                dF.set_size( 3, 8, 0.0 );

                dF( 1, 0 ) = zeta-1.0;
                dF( 2, 0 ) = eta-1.0;

                dF( 0, 1 ) = 1.0-zeta;
                dF( 2, 1 ) = -xi-1.0;

                dF( 1, 2 ) =  zeta-1.0;
                dF( 2, 2 ) =  eta+1.0;

                dF( 0, 3 ) = 1.0-zeta;
                dF( 2, 3 ) = 1.0-xi;

                dF( 1, 4 ) = -zeta-1.0;
                dF( 2, 4 ) = 1.0-eta;

                dF( 0, 5 ) = zeta+1.0;
                dF( 2, 5 ) = xi+1.0;

                dF( 1, 6 ) = -zeta-1.0;
                dF( 2, 6 ) = -eta-1.0;

                dF( 0, 7 ) = zeta+1.0;
                dF( 2, 7 ) = xi-1.0;

                dF *= 0.125;
            }
        }

        const Matrix <real> &
        EF_HEX8TS::E( const uint aIndex )
        {
            this->update_nabla( aIndex );

            Vector< real > & F = mF( aIndex );
            Matrix< real > & Nabla = mInvJ ;

            mE( 0, 0 ) = mS[ 0 ] * F ( 0 ) * Nabla( 0, 0 );
            mE( 1, 0 ) = mS[ 0 ] * F ( 0 ) * Nabla( 1, 0 );
            mE( 2, 0 ) = mS[ 0 ] * F ( 0 ) * Nabla( 2, 0 );

            mE( 0, 1 ) = mS[ 1 ] * F ( 1 ) * Nabla( 0, 1 );
            mE( 1, 1 ) = mS[ 1 ] * F ( 1 ) * Nabla( 1, 1 );
            mE( 2, 1 ) = mS[ 1 ] * F ( 1 ) * Nabla( 2, 1 );

            mE( 0, 2 ) = mS[ 2 ] * F ( 2 ) * Nabla( 0, 0 );
            mE( 1, 2 ) = mS[ 2 ] * F ( 2 ) * Nabla( 1, 0 );
            mE( 2, 2 ) = mS[ 2 ] * F ( 2 ) * Nabla( 2, 0 );

            mE( 0, 3 ) = mS[ 3 ] * F ( 3 ) * Nabla( 0, 1 );
            mE( 1, 3 ) = mS[ 3 ] * F ( 3 ) * Nabla( 1, 1 );
            mE( 2, 3 ) = mS[ 3 ] * F ( 3 ) * Nabla( 2, 1 );

            mE( 0, 4 ) = mS[ 4 ] * F ( 4 ) * Nabla( 0, 0 );
            mE( 1, 4 ) = mS[ 4 ] * F ( 4 ) * Nabla( 1, 0 );
            mE( 2, 4 ) = mS[ 4 ] * F ( 4 ) * Nabla( 2, 0 );

            mE( 0, 5 ) = mS[ 5 ] * F ( 5 ) * Nabla( 0, 1 );
            mE( 1, 5 ) = mS[ 5 ] * F ( 5 ) * Nabla( 1, 1 );
            mE( 2, 5 ) = mS[ 5 ] * F ( 5 ) * Nabla( 2, 1 );

            mE( 0, 6 ) = mS[ 6 ] * F ( 6 ) * Nabla( 0, 0 );
            mE( 1, 6 ) = mS[ 6 ] * F ( 6 ) * Nabla( 1, 0 );
            mE( 2, 6 ) = mS[ 6 ] * F ( 6 ) * Nabla( 2, 0 );

            mE( 0, 7 ) = mS[ 7 ] * F ( 7 ) * Nabla( 0, 1 );
            mE( 1, 7 ) = mS[ 7 ] * F ( 7 ) * Nabla( 1, 1 );
            mE( 2, 7 ) = mS[ 7 ] * F ( 7 ) * Nabla( 2, 1 );

            return mE ;
        }


        const Matrix<real> &
        EF_HEX8TS::C( const uint aIndex )
        {
            this->update_nabla( aIndex );

            Matrix<real> & Nabla = mInvJ;

            Matrix<real> & dF = mNx ;
            mNx = Nabla * mFxi( aIndex ) ;

            mEx( 0, 0 ) = mS[ 0 ] * dF( 0, 0 ) * Nabla( 0, 0 );
            mEx( 1, 0 ) = mS[ 0 ] * dF( 0, 0 ) * Nabla( 1, 0 );
            mEx( 2, 0 ) = mS[ 0 ] * dF( 0, 0 ) * Nabla( 2, 0 );

            mEx( 0, 1 ) = mS[ 1 ] * dF( 0, 1 ) * Nabla( 0, 1 );
            mEx( 1, 1 ) = mS[ 1 ] * dF( 0, 1 ) * Nabla( 1, 1 );
            mEx( 2, 1 ) = mS[ 1 ] * dF( 0, 1 ) * Nabla( 2, 1 );

            mEx( 0, 2 ) = mS[ 2 ] * dF( 0, 2 ) * Nabla( 0, 0 );
            mEx( 1, 2 ) = mS[ 2 ] * dF( 0, 2 ) * Nabla( 1, 0 );
            mEx( 2, 2 ) = mS[ 2 ] * dF( 0, 2 ) * Nabla( 2, 0 );

            mEx( 0, 3 ) = mS[ 3 ] * dF( 0, 3 ) * Nabla( 0, 1 );
            mEx( 1, 3 ) = mS[ 3 ] * dF( 0, 3 ) * Nabla( 1, 1 );
            mEx( 2, 3 ) = mS[ 3 ] * dF( 0, 3 ) * Nabla( 2, 1 );

            mEx( 0, 4 ) = mS[ 4 ] * dF( 0, 4 ) * Nabla( 0, 0 );
            mEx( 1, 4 ) = mS[ 4 ] * dF( 0, 4 ) * Nabla( 1, 0 );
            mEx( 2, 4 ) = mS[ 4 ] * dF( 0, 4 ) * Nabla( 2, 0 );

            mEx( 0, 5 ) = mS[ 5 ] * dF( 0, 5 ) * Nabla( 0, 1 );
            mEx( 1, 5 ) = mS[ 5 ] * dF( 0, 5 ) * Nabla( 1, 1 );
            mEx( 2, 5 ) = mS[ 5 ] * dF( 0, 5 ) * Nabla( 2, 1 );

            mEx( 0, 6 ) = mS[ 6 ] * dF( 0, 6 ) * Nabla( 0, 0 );
            mEx( 1, 6 ) = mS[ 6 ] * dF( 0, 6 ) * Nabla( 1, 0 );
            mEx( 2, 6 ) = mS[ 6 ] * dF( 0, 6 ) * Nabla( 2, 0 );

            mEx( 0, 7 ) = mS[ 7 ] * dF( 0, 7 ) * Nabla( 0, 1 );
            mEx( 1, 7 ) = mS[ 7 ] * dF( 0, 7 ) * Nabla( 1, 1 );
            mEx( 2, 7 ) = mS[ 7 ] * dF( 0, 7 ) * Nabla( 2, 1 );

            mEy( 0, 0 ) = mS[ 0 ] * dF( 1, 0 ) * Nabla( 0, 0 );
            mEy( 1, 0 ) = mS[ 0 ] * dF( 1, 0 ) * Nabla( 1, 0 );
            mEy( 2, 0 ) = mS[ 0 ] * dF( 1, 0 ) * Nabla( 2, 0 );

            mEy( 0, 1 ) = mS[ 1 ] * dF( 1, 1 ) * Nabla( 0, 1 );
            mEy( 1, 1 ) = mS[ 1 ] * dF( 1, 1 ) * Nabla( 1, 1 );
            mEy( 2, 1 ) = mS[ 1 ] * dF( 1, 1 ) * Nabla( 2, 1 );

            mEy( 0, 2 ) = mS[ 2 ] * dF( 1, 2 ) * Nabla( 0, 0 );
            mEy( 1, 2 ) = mS[ 2 ] * dF( 1, 2 ) * Nabla( 1, 0 );
            mEy( 2, 2 ) = mS[ 2 ] * dF( 1, 2 ) * Nabla( 2, 0 );

            mEy( 0, 3 ) = mS[ 3 ] * dF( 1, 3 ) * Nabla( 0, 1 );
            mEy( 1, 3 ) = mS[ 3 ] * dF( 1, 3 ) * Nabla( 1, 1 );
            mEy( 2, 3 ) = mS[ 3 ] * dF( 1, 3 ) * Nabla( 2, 1 );

            mEy( 0, 4 ) = mS[ 4 ] * dF( 1, 4 ) * Nabla( 0, 0 );
            mEy( 1, 4 ) = mS[ 4 ] * dF( 1, 4 ) * Nabla( 1, 0 );
            mEy( 2, 4 ) = mS[ 4 ] * dF( 1, 4 ) * Nabla( 2, 0 );

            mEy( 0, 5 ) = mS[ 5 ] * dF( 1, 5 ) * Nabla( 0, 1 );
            mEy( 1, 5 ) = mS[ 5 ] * dF( 1, 5 ) * Nabla( 1, 1 );
            mEy( 2, 5 ) = mS[ 5 ] * dF( 1, 5 ) * Nabla( 2, 1 );

            mEy( 0, 6 ) = mS[ 6 ] * dF( 1, 6 ) * Nabla( 0, 0 );
            mEy( 1, 6 ) = mS[ 6 ] * dF( 1, 6 ) * Nabla( 1, 0 );
            mEy( 2, 6 ) = mS[ 6 ] * dF( 1, 6 ) * Nabla( 2, 0 );

            mEy( 0, 7 ) = mS[ 7 ] * dF( 1, 7 ) * Nabla( 0, 1 );
            mEy( 1, 7 ) = mS[ 7 ] * dF( 1, 7 ) * Nabla( 1, 1 );
            mEy( 2, 7 ) = mS[ 7 ] * dF( 1, 7 ) * Nabla( 2, 1 );

            mEz( 0, 0 ) = mS[ 0 ] * dF( 2, 0 ) * Nabla( 0, 0 );
            mEz( 1, 0 ) = mS[ 0 ] * dF( 2, 0 ) * Nabla( 1, 0 );
            mEz( 2, 0 ) = mS[ 0 ] * dF( 2, 0 ) * Nabla( 2, 0 );

            mEz( 0, 1 ) = mS[ 1 ] * dF( 2, 1 ) * Nabla( 0, 1 );
            mEz( 1, 1 ) = mS[ 1 ] * dF( 2, 1 ) * Nabla( 1, 1 );
            mEz( 2, 1 ) = mS[ 1 ] * dF( 2, 1 ) * Nabla( 2, 1 );

            mEz( 0, 2 ) = mS[ 2 ] * dF( 2, 2 ) * Nabla( 0, 0 );
            mEz( 1, 2 ) = mS[ 2 ] * dF( 2, 2 ) * Nabla( 1, 0 );
            mEz( 2, 2 ) = mS[ 2 ] * dF( 2, 2 ) * Nabla( 2, 0 );

            mEz( 0, 3 ) = mS[ 3 ] * dF( 2, 3 ) * Nabla( 0, 1 );
            mEz( 1, 3 ) = mS[ 3 ] * dF( 2, 3 ) * Nabla( 1, 1 );
            mEz( 2, 3 ) = mS[ 3 ] * dF( 2, 3 ) * Nabla( 2, 1 );

            mEz( 0, 4 ) = mS[ 4 ] * dF( 2, 4 ) * Nabla( 0, 0 );
            mEz( 1, 4 ) = mS[ 4 ] * dF( 2, 4 ) * Nabla( 1, 0 );
            mEz( 2, 4 ) = mS[ 4 ] * dF( 2, 4 ) * Nabla( 2, 0 );

            mEz( 0, 5 ) = mS[ 5 ] * dF( 2, 5 ) * Nabla( 0, 1 );
            mEz( 1, 5 ) = mS[ 5 ] * dF( 2, 5 ) * Nabla( 1, 1 );
            mEz( 2, 5 ) = mS[ 5 ] * dF( 2, 5 ) * Nabla( 2, 1 );

            mEz( 0, 6 ) = mS[ 6 ] * dF( 2, 6 ) * Nabla( 0, 0 );
            mEz( 1, 6 ) = mS[ 6 ] * dF( 2, 6 ) * Nabla( 1, 0 );
            mEz( 2, 6 ) = mS[ 6 ] * dF( 2, 6 ) * Nabla( 2, 0 );

            mEz( 0, 7 ) = mS[ 7 ] * dF( 2, 7 ) * Nabla( 0, 1 );
            mEz( 1, 7 ) = mS[ 7 ] * dF( 2, 7 ) * Nabla( 1, 1 );
            mEz( 2, 7 ) = mS[ 7 ] * dF( 2, 7 ) * Nabla( 2, 1 );


            mC.set_row( 0, mEy.row( 2 ) - mEz.row( 1 ) );
            mC.set_row( 1, mEz.row( 0 ) - mEx.row( 2 ) );
            mC.set_row( 2, mEx.row( 1 ) - mEy.row( 0 ) );

            return mC;
        }

        void
        EF_HEX8TS::update_nabla( const uint aIndex )
        {
            if ( aIndex == mLastIndex ) return ;

            mLastIndex = aIndex;

            mJ2 = mQuad4Nxi( aIndex ) * mX;

            // [ dx/dxi dy/dxi dz/dxi ]
            mU( 0 ) = mJ2( 0, 0 ) ;
            mU( 1 ) = mJ2( 0, 1 ) ;
            mU( 2 ) = mJ2( 0, 2 ) ;

            // [ dx/deta dy/deta dz/deta ]
            mV( 0 ) = mJ2( 1, 0 ) ;
            mV( 1 ) = mJ2( 1, 1 ) ;
            mV( 2 ) = mJ2( 1, 2 ) ;

            mW = cross( mU, mV );


            // normal component  // [ dx/dzeta dy/dzeta dz/dzeta ]

            real tS = norm( mW );

            mDetJ = 0.5 * tS * mThickness ;

            mW /= mDetJ ;

            mG = mJ2 * trans( mJ2 ) ;

            inv2( mG, mInvG );

            mInvJ( 0, 0 ) = mJ2( 0, 0 ) * mInvG( 0, 0 ) + mJ2( 1, 0 ) * mInvG( 1, 0 );
            mInvJ( 1, 0 ) = mJ2( 0, 1 ) * mInvG( 0, 0 ) + mJ2( 1, 1 ) * mInvG( 1, 0 );
            mInvJ( 2, 0 ) = mJ2( 0, 2 ) * mInvG( 0, 0 ) + mJ2( 1, 2 ) * mInvG( 1, 0 );

            mInvJ( 0, 1 ) = mJ2( 0, 0 ) * mInvG( 0, 1 ) + mJ2( 1, 0 ) * mInvG( 1, 1 );
            mInvJ( 1, 1 ) = mJ2( 0, 1 ) * mInvG( 0, 1 ) + mJ2( 1, 1 ) * mInvG( 1, 1 );
            mInvJ( 2, 1 ) = mJ2( 0, 2 ) * mInvG( 0, 1 ) + mJ2( 1, 2 ) * mInvG( 1, 1 );


            mInvJ( 0, 2 ) = mW( 0 );
            mInvJ( 1, 2 ) = mW( 1 );
            mInvJ( 2, 2 ) = mW( 2 );

            mAbsDetJ = mDetJ ; // we know that it is positive
        }
//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_HEX8TS::G( const uint aIndex )
        {
            // grad( e_k ) = s_k * ( grad F_k ) (x) nabla_{a(k)} with the
            // per-point nablas held fixed — the same ingredients C() uses,
            // so antisym( G ) reproduces C on every geometry; scope of the
            // omitted Hessian term is documented at the declaration
            this->update_nabla( aIndex );

            mNx = mInvJ * mFxi( aIndex );

            // nabla column per dof, same alternation as E()
            const uint tAxis[ 8 ] = { 0, 1, 0, 1, 0, 1, 0, 1 };

            for ( uint k = 0; k < 8; ++k )
            {
                const uint a = tAxis[ k ];

                for ( uint j = 0; j < 3; ++j )       // field component
                {
                    for ( uint i = 0; i < 3; ++i )   // derivative direction
                    {
                        mGrad( i + 3 * j, k ) =
                                mS[ k ] * mNx( i, k ) * mInvJ( j, a );
                    }
                }
            }

            return mGrad ;
        }

    }
}