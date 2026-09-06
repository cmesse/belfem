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
#include "cl_EF_HEX8.hpp"
#include "cl_IF_InterpolationFunctionTemplate.hpp"

// Lagrange template specializations for HEX8 live in this header (without the
// `inline` keyword on the vtable-emitting functions). Include it here — not
// only in the factory — so this TU emits the vtable with the specialized
// N()/dNdXi(), not the generic throwing ones.
#include "lagrange/cl_IF_HEX8.hpp"
#include "fn_cross.hpp"
namespace belfem
{
    namespace fem
    {
        EF_HEX8::EF_HEX8()
        {
            mJ.set_size( 3, 3 );
            mInvJ.set_size( 3, 3 );
            mX.set_size( 8, 3 );
            mNumDofs = 12 ;

            // Reference-element volume of the cube [-1,1]^3 with standard
            // Gauss integration. Currently unused by any assembly site
            // (sum_w() is not consumed downstream), but kept consistent
            // with the "reference-element volume" convention documented
            // in src/fem/interpolation/doc/nedelec.md.
            mSumW = 8.0 ;

            mHex8 = new InterpolationFunctionTemplate<
                GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >;

            mE.set_size( 3, 12 );
            mC.set_size( 3, 12 );
            mGrad.set_size( 9, 12 );

            mNx.set_size( 3, 12, 0. );

            mEx.set_size( 3, 12, 0. );
            mEy.set_size( 3, 12, 0. );
            mEz.set_size( 3, 12, 0. );


            for ( uint k=0; k<12; ++k )
            {
                mS[ k ] = 0.0 ;
            }

            mU.set_size( 3 );
            mV.set_size( 3 );
        }

        EF_HEX8::~EF_HEX8()
        {
            if ( mHex8 != nullptr )  delete mHex8;
        }

        void
        EF_HEX8::precompute( const Matrix< real > & aXi )
        {
            uint n = aXi.n_cols();

            mF.set_size( n, {{}} );
            mFxi.set_size( n, {{}} );

            mHex8Nxi.set_size( n, {{}} );
            mHex8Nxi2.set_size( n, {{}} );

            Vector< real > tXi( 3 );
            for ( uint k=0; k<n; ++k )
            {
                mHex8->dNdXi( aXi.col( k ), mHex8Nxi( k ) );

                // Voigt second derivatives for G(); d2NdXi2 resizes its
                // work matrix, which is why this call lives in precompute
                // and never in G() itself
                mHex8->d2NdXi2( aXi.col( k ), mHex8Nxi2( k ) );

                real   xi = aXi( 0, k );
                real  eta = aXi( 1, k );
                real zeta = aXi( 2, k );

                Vector< real > & F = mF( k );

                F.set_size( 12 );

                // Scale 1/8 = 1/4 (reference-QUAD Nédélec factor for unit edge
                // circulation on [-1,1]^2) × 1/2 (through-thickness ζ extrusion
                // over [-1,1]). Together they yield ∫ E_k · dℓ_j = s_k · δ_{jk}.
                //
                // Signs follow the mesh::HEX8 edge orientations
                // (see cl_Element_HEX8.hpp):
                //   bottom loop (0..3): 0→1, 1→2, 2→3, 3→0 at zeta=-1
                //   top loop    (4..7): 4→5, 5→6, 6→7, 7→4 at zeta=+1
                //   verticals  (8..11): 0→4, 1→5, 2→6, 3→7
                F(0) =  (1.-eta)*(1.-zeta);
                F(1) =  (1.+xi) *(1.-zeta);
                F(2) = -(1.+eta)*(1.-zeta);
                F(3) = -(1.-xi) *(1.-zeta);

                F(4) = -(eta-1.)*(zeta+1.);
                F(5) =  (xi+1.)*(zeta+1.);
                F(6) = -(eta+1.)*(zeta+1.);
                F(7) =  (xi-1.)*(zeta+1.);


                F(8) = (eta-1.)*(xi-1.);
                F(9) = -(eta-1.)*(xi+1.);
                F(10) = (eta+1.)*(xi+1.);
                F(11) = -(eta+1.)*(xi-1.);

                F *= 0.125;

                Matrix< real > & dF = mFxi( k );
                dF.set_size( 3, 12 );

                dF(0,0) = 0.;
                dF(1,0) = zeta-1.;
                dF(2,0) = eta-1.;


                dF(0,1) = 1.-zeta;
                dF(1,1) = 0.;
                dF(2,1) = -xi-1.;

                dF(0,2) = 0.;
                dF(1,2) = zeta-1.;
                dF(2,2) = eta+1.;

                dF(0,3) = 1.-zeta;
                dF(1,3) = 0.;
                dF(2,3) = 1.-xi;

                dF(0,4) = 0.;
                dF(1,4) = -zeta-1.;
                dF(2,4) = -eta+1.;

                dF(0,5) = zeta+1.;
                dF(1,5) = 0.;
                dF(2,5) = xi+1 ;

                dF(0,6) = 0.;
                dF(1,6) = -zeta-1.;
                dF(2,6) = -eta-1.;

                dF(0,7) = zeta+1.;
                dF(1,7) = 0.;
                dF(2,7) = xi-1.;

                dF(0,8) = eta-1.;
                dF(1,8) = xi-1.;
                dF(2,8) = 0.;

                dF(0,9) = -eta+1.;
                dF(1,9) = -xi-1.;
                dF(2,9) = 0.;

                dF(0,10) = eta+1.;
                dF(1,10) = xi+1.;
                dF(2,10) = 0.;

                dF(0,11) = -eta-1.;
                dF(1,11) = -xi+1.;
                dF(2,11) = 0.;

                dF *= 0.125;
            }
        }

        void EF_HEX8::link( Element * aElement )
        {
            // populate node coordinates
            for ( uint k=0; k<8; ++k )
            {
                mesh::Node * tNode = aElement->element()->node( k );
                mX( k, 0 ) = tNode->x( 0 );
                mX( k, 1 ) = tNode->x( 1 );
                mX( k, 2 ) = tNode->x( 2 );
            }

            // grab directions from edge
            aElement->edge_directions( mS );

            this->estimate_volume();

            // reset nabla cache
            mLastIndex = BELFEM_UINT_MAX ;

        }

        void EF_HEX8::estimate_volume()
        {
            // test thickness in X-direction
            mU.fill( 0 );
            for ( uint k : mXLengthA )
            {
                mU( 0 ) += mX( k,0 );
                mU( 1 ) += mX( k,1 );
                mU( 2 ) += mX( k,2 );
            }
            mV.fill( 0 );
            for ( uint k : mXLengthB )
            {
                mV( 0 ) += mX( k,0 );
                mV( 1 ) += mX( k,1 );
                mV( 2 ) += mX( k,2 );
            }
            real h_xi = 0.25 * norm( mU - mV );

            // test thickness in Y-direction
            mU.fill( 0 );
            for ( uint k : mYLengthA )
            {
                mU( 0 ) += mX( k,0 );
                mU( 1 ) += mX( k,1 );
                mU( 2 ) += mX( k,2 );
            }
            mV.fill( 0 );
            for ( uint k : mYLengthB )
            {
                mV( 0 ) += mX( k,0 );
                mV( 1 ) += mX( k,1 );
                mV( 2 ) += mX( k,2 );
            }
            real h_eta = 0.25 * norm( mU - mV );

            // test thickness in Z-direction
            mU.fill( 0 );
            for ( uint k : mZLengthA )
            {
                mU( 0 ) += mX( k,0 );
                mU( 1 ) += mX( k,1 );
                mU( 2 ) += mX( k,2 );
            }
            mV.fill( 0 );
            for ( uint k : mZLengthB )
            {
                mV( 0 ) += mX( k,0 );
                mV( 1 ) += mX( k,1 );
                mV( 2 ) += mX( k,2 );
            }
            real h_zeta = 0.25 * norm( mU - mV );

            real h_min = std::min( h_xi, std::min( h_eta, h_zeta ) );
            real h_max = std::max( h_xi, std::max( h_eta, h_zeta ) );

            if ( h_min == h_xi )
            {

                mU( 0 ) = mX( 4,0 )-mX( 0,0 );
                mU( 1 ) = mX( 4,1 )-mX( 0,1 );
                mU( 2 ) = mX( 4,2 )-mX( 0,2 );

                mV( 0 ) = mX( 3,0 )-mX( 0,0 );
                mV( 1 ) = mX( 3,1 )-mX( 0,1 );
                mV( 2 ) = mX( 3,2 )-mX( 0,2 );
            }
            else if ( h_min == h_eta )
            {
                mU( 0 ) = mX( 1,0 )-mX( 0,0 );
                mU( 1 ) = mX( 1,1 )-mX( 0,1 );
                mU( 2 ) = mX( 1,2 )-mX( 0,2 );

                mV( 0 ) = mX( 4,0 )-mX( 0,0 );
                mV( 1 ) = mX( 4,1 )-mX( 0,1 );
                mV( 2 ) = mX( 4,2 )-mX( 0,2 );
            }
            else
            {
                mU( 0 ) = mX( 1,0 )-mX( 0,0 );
                mU( 1 ) = mX( 1,1 )-mX( 0,1 );
                mU( 2 ) = mX( 1,2 )-mX( 0,2 );

                mV( 0 ) = mX( 3,0 )-mX( 0,0 );
                mV( 1 ) = mX( 3,1 )-mX( 0,1 );
                mV( 2 ) = mX( 3,2 )-mX( 0,2 );
            }

            real tSurface = norm( cross( mU, mV ) );

            mApproxVolume = h_min * tSurface;

            // a thickness of less than 10 micrometers indicates that this element is very distorted
            // if so, we must use the approximate volume to get a reasonable jacobian
            mSmallVolume = h_min < 1e-5 || h_max / h_min > 100 ;
        }

        const Matrix <real> &
        EF_HEX8::E( const uint aIndex )
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

            mE( 0, 8 ) = mS[ 8 ] * F ( 8 ) * Nabla( 0, 2 );
            mE( 1, 8 ) = mS[ 8 ] * F ( 8 ) * Nabla( 1, 2 );
            mE( 2, 8 ) = mS[ 8 ] * F ( 8 ) * Nabla( 2, 2 );

            mE( 0, 9 ) = mS[ 9 ] * F ( 9 ) * Nabla( 0, 2 );
            mE( 1, 9 ) = mS[ 9 ] * F ( 9 ) * Nabla( 1, 2 );
            mE( 2, 9 ) = mS[ 9 ] * F ( 9 ) * Nabla( 2, 2 );

            mE( 0, 10 ) = mS[ 10 ] * F ( 10 ) * Nabla( 0, 2 );
            mE( 1, 10 ) = mS[ 10 ] * F ( 10 ) * Nabla( 1, 2 );
            mE( 2, 10 ) = mS[ 10 ] * F ( 10 ) * Nabla( 2, 2 );

            mE( 0, 11 ) = mS[ 11 ] * F ( 11 ) * Nabla( 0, 2 );
            mE( 1, 11 ) = mS[ 11 ] * F ( 11 ) * Nabla( 1, 2 );
            mE( 2, 11 ) = mS[ 11 ] * F ( 11 ) * Nabla( 2, 2 );


            return mE ;
        }

         const Matrix<real> &
         EF_HEX8::C( const uint aIndex )
        {
            this->update_nabla( aIndex );

            Matrix<real> & Nabla = mInvJ;

            Matrix<real> & dF = mNx ;
            mNx = mInvJ * mFxi( aIndex ) ;

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

            mEx( 0, 8 ) = mS[ 8 ] * dF( 0, 8 ) * Nabla( 0, 2 );
            mEx( 1, 8 ) = mS[ 8 ] * dF( 0, 8 ) * Nabla( 1, 2 );
            mEx( 2, 8 ) = mS[ 8 ] * dF( 0, 8 ) * Nabla( 2, 2 );

            mEx( 0, 9 ) = mS[ 9 ] * dF( 0, 9 ) * Nabla( 0, 2 );
            mEx( 1, 9 ) = mS[ 9 ] * dF( 0, 9 ) * Nabla( 1, 2 );
            mEx( 2, 9 ) = mS[ 9 ] * dF( 0, 9 ) * Nabla( 2, 2 );

            mEx( 0, 10 ) = mS[ 10 ] * dF( 0, 10 ) * Nabla( 0, 2 );
            mEx( 1, 10 ) = mS[ 10 ] * dF( 0, 10 ) * Nabla( 1, 2 );
            mEx( 2, 10 ) = mS[ 10 ] * dF( 0, 10 ) * Nabla( 2, 2 );

            mEx( 0, 11 ) = mS[ 11 ] * dF( 0, 11 ) * Nabla( 0, 2 );
            mEx( 1, 11 ) = mS[ 11 ] * dF( 0, 11 ) * Nabla( 1, 2 );
            mEx( 2, 11 ) = mS[ 11 ] * dF( 0, 11 ) * Nabla( 2, 2 );


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

            mEy( 0, 8 ) = mS[ 8 ] * dF( 1, 8 ) * Nabla( 0, 2 );
            mEy( 1, 8 ) = mS[ 8 ] * dF( 1, 8 ) * Nabla( 1, 2 );
            mEy( 2, 8 ) = mS[ 8 ] * dF( 1, 8 ) * Nabla( 2, 2 );

            mEy( 0, 9 ) = mS[ 9 ] * dF( 1, 9 ) * Nabla( 0, 2 );
            mEy( 1, 9 ) = mS[ 9 ] * dF( 1, 9 ) * Nabla( 1, 2 );
            mEy( 2, 9 ) = mS[ 9 ] * dF( 1, 9 ) * Nabla( 2, 2 );

            mEy( 0, 10 ) = mS[ 10 ] * dF( 1, 10 ) * Nabla( 0, 2 );
            mEy( 1, 10 ) = mS[ 10 ] * dF( 1, 10 ) * Nabla( 1, 2 );
            mEy( 2, 10 ) = mS[ 10 ] * dF( 1, 10 ) * Nabla( 2, 2 );

            mEy( 0, 11 ) = mS[ 11 ] * dF( 1, 11 ) * Nabla( 0, 2 );
            mEy( 1, 11 ) = mS[ 11 ] * dF( 1, 11 ) * Nabla( 1, 2 );
            mEy( 2, 11 ) = mS[ 11 ] * dF( 1, 11 ) * Nabla( 2, 2 );



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

            mEz( 0, 8 ) = mS[ 8 ] * dF( 2, 8 ) * Nabla( 0, 2 );
            mEz( 1, 8 ) = mS[ 8 ] * dF( 2, 8 ) * Nabla( 1, 2 );
            mEz( 2, 8 ) = mS[ 8 ] * dF( 2, 8 ) * Nabla( 2, 2 );

            mEz( 0, 9 ) = mS[ 9 ] * dF( 2, 9 ) * Nabla( 0, 2 );
            mEz( 1, 9 ) = mS[ 9 ] * dF( 2, 9 ) * Nabla( 1, 2 );
            mEz( 2, 9 ) = mS[ 9 ] * dF( 2, 9 ) * Nabla( 2, 2);

            mEz( 0, 10 ) = mS[ 10 ] * dF( 2, 10 ) * Nabla( 0, 2 );
            mEz( 1, 10 ) = mS[ 10 ] * dF( 2, 10 ) * Nabla( 1, 2 );
            mEz( 2, 10 ) = mS[ 10 ] * dF( 2, 10 ) * Nabla( 2, 2 );

            mEz( 0, 11 ) = mS[ 11 ] * dF( 2, 11 ) * Nabla( 0, 2 );
            mEz( 1, 11 ) = mS[ 11 ] * dF( 2, 11 ) * Nabla( 1, 2 );
            mEz( 2, 11 ) = mS[ 11 ] * dF( 2, 11 ) * Nabla( 2, 2 );

            // fix ( 2026-08-23 ): the rows were assembled with the
            // operands swapped, returning -curl. With mEx/mEy/mEz holding
            // d(e_i)/dx, d(e_i)/dy, d(e_i)/dz in their rows i, the curl is
            //   ( curl e )_x = d(e_z)/dy - d(e_y)/dz = mEy.row(2) - mEz.row(1)
            // and cyclic — the same assembly EF_HEX8TS::C() always had.
            // Sign pinned by the Stokes-consistency test in the G battery.
            mC.set_row( 0, mEy.row( 2 ) - mEz.row( 1 ) );
            mC.set_row( 1, mEz.row( 0 ) - mEx.row( 2 ) );
            mC.set_row( 2, mEx.row( 1 ) - mEy.row( 0 ) );

            return mC;
        }

        void
        EF_HEX8::update_nabla( const uint aIndex )
        {
            if ( aIndex == mLastIndex ) return ;

            mLastIndex = aIndex;

            mJ = mHex8Nxi( aIndex ) * mX;

            mInvJ( 0, 0 ) = mJ( 1, 1 ) * mJ( 2, 2 ) - mJ( 1, 2 ) * mJ( 2, 1 ) ;
            mInvJ( 1, 0 ) = mJ( 1, 2 ) * mJ( 2, 0 ) - mJ( 1, 0 ) * mJ( 2, 2 ) ;
            mInvJ( 2, 0 ) = mJ( 1, 0 ) * mJ( 2, 1 ) - mJ( 1, 1 ) * mJ( 2, 0 ) ;

            mInvJ( 0, 1 ) = mJ( 0, 2 ) * mJ( 2, 1 ) - mJ( 0, 1 ) * mJ( 2, 2 ) ;
            mInvJ( 1, 1 ) = mJ( 0, 0 ) * mJ( 2, 2 ) - mJ( 0, 2 ) * mJ( 2, 0 ) ;
            mInvJ( 2, 1 ) = mJ( 0, 1 ) * mJ( 2, 0 ) - mJ( 0, 0 ) * mJ( 2, 1 ) ;

            mInvJ( 0, 2 ) = mJ( 0, 1 ) * mJ( 1, 2 ) - mJ( 0, 2 ) * mJ( 1, 1 ) ;
            mInvJ( 1, 2 ) = mJ( 0, 2 ) * mJ( 1, 0 ) - mJ( 0, 0 ) * mJ( 1, 2 ) ;
            mInvJ( 2, 2 ) = mJ( 0, 0 ) * mJ( 1, 1 ) - mJ( 0, 1 ) * mJ( 1, 0 ) ;

            // fix for very small elements
            real tDet = det( mJ ) ;
            real tAbsDetJ = std::abs( tDet ) ;

            if ( mSmallVolume )
            {
                mDetJ = 0.125 * mApproxVolume ;
                mAbsDetJ = mDetJ ;
            }
            else
            {
                mDetJ = tDet;
                mAbsDetJ = tAbsDetJ;
            }

            mInvJ /= tAbsDetJ < BELFEM_EPSILON ? mDetJ : tDet ;
        }
//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_HEX8::G( const uint aIndex )
        {
            // grad( e_k ) = s_k [ grad( F_k ) (x) grad( xi_a )
            //                   + F_k * Hess( xi_a ) ]
            // The second term is the inverse-map Hessian the curl cancels
            // ( symmetric ) and the gradient keeps. Derivation + probe:
            // exchange thread coulomb_gauge_stepC_hex8.md.
            this->update_nabla( aIndex );

            Matrix< real > & A = mInvJ ;

            // physical gradient of the scalar factors, unsigned — the same
            // product C() builds; C() recomputes it unconditionally, so
            // writing mNx here is safe
            mNx = mInvJ * mFxi( aIndex ) ;

            // mixed second derivatives of the map:
            // tS[ pair ][ m ] = sum_n d2N( voigt-row, n ) * mX( n, m )
            // with pair 0 = ( xi, eta ) = Voigt row 5,
            //      pair 1 = ( xi, zeta ) = Voigt row 4,
            //      pair 2 = ( eta, zeta ) = Voigt row 3
            const Matrix< real > & tD2 = mHex8Nxi2( aIndex );
            const uint tVoigtRow[ 3 ] = { 5, 4, 3 };
            real tS[ 3 ][ 3 ];
            for ( uint p = 0; p < 3; ++p )
            {
                for ( uint m = 0; m < 3; ++m )
                {
                    real tSum = 0.0 ;
                    for ( uint n = 0; n < 8; ++n )
                    {
                        tSum += tD2( tVoigtRow[ p ], n ) * mX( n, m );
                    }
                    tS[ p ][ m ] = tSum ;
                }
            }

            // Hessians of the inverse map, one per parameter axis:
            //   Cm( d, e ) = S^(de) . A( :, a )   ( symmetric, zero diag )
            //   H^(a)      = - A * Cm * A^T
            // ONE form only — Cm is already symmetric, adding the pair
            // outer-product symmetrization on top would double-count H
            // ( and the curl tie could not see it )
            real tH[ 3 ][ 9 ];
            for ( uint a = 0; a < 3; ++a )
            {
                real tCm[ 3 ][ 3 ];
                real tC01 = 0.0 ;
                real tC02 = 0.0 ;
                real tC12 = 0.0 ;
                for ( uint m = 0; m < 3; ++m )
                {
                    tC01 += tS[ 0 ][ m ] * A( m, a );
                    tC02 += tS[ 1 ][ m ] * A( m, a );
                    tC12 += tS[ 2 ][ m ] * A( m, a );
                }
                tCm[ 0 ][ 0 ] = 0.0 ;  tCm[ 0 ][ 1 ] = tC01 ; tCm[ 0 ][ 2 ] = tC02 ;
                tCm[ 1 ][ 0 ] = tC01 ; tCm[ 1 ][ 1 ] = 0.0 ;  tCm[ 1 ][ 2 ] = tC12 ;
                tCm[ 2 ][ 0 ] = tC02 ; tCm[ 2 ][ 1 ] = tC12 ; tCm[ 2 ][ 2 ] = 0.0 ;

                for ( uint j = 0; j < 3; ++j )        // field component
                {
                    for ( uint i = 0; i < 3; ++i )    // derivative direction
                    {
                        real tSum = 0.0 ;
                        for ( uint d = 0; d < 3; ++d )
                        {
                            for ( uint e = 0; e < 3; ++e )
                            {
                                tSum += A( i, d ) * tCm[ d ][ e ] * A( j, e );
                            }
                        }
                        tH[ a ][ i + 3 * j ] = -tSum ;
                    }
                }
            }

            // assemble; axis table read from E(): edges 0,2,4,6 -> grad xi,
            // 1,3,5,7 -> grad eta, 8..11 -> grad zeta. mNx and mF are
            // unsigned; mS applied once on the whole bracket
            const uint tAxis[ 12 ] = { 0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2 };
            const Vector< real > & F = mF( aIndex );

            for ( uint k = 0; k < 12; ++k )
            {
                const uint a = tAxis[ k ];
                for ( uint j = 0; j < 3; ++j )        // field component
                {
                    for ( uint i = 0; i < 3; ++i )    // derivative direction
                    {
                        mGrad( i + 3 * j, k ) = mS[ k ] *
                              ( mNx( i, k ) * A( j, a )
                              + F( k ) * tH[ a ][ i + 3 * j ] );
                    }
                }
            }
            return mGrad ;
        }

    }
}