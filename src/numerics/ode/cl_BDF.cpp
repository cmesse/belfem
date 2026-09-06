/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_BDF.hpp"

#include "assert.hpp"

#include "../../linalg/lapack/fn_gesv.hpp"

namespace belfem
{
    namespace ode
    {
        BDF::BDF( ShiftRegister< real > & aH )
                : mH( aH )
        {
            mS.set_size( aH.capacity() + 1 );
            mPivot.set_size( aH.capacity() + 1 );
            BELFEM_ERROR(  mH.capacity() < 7, "BDF functions higher than 6 are not recommended" );
        }

        void
        BDF::compute_coefficients()
        {
            size_t n = mH.size() + 1 ;

            // compute the timesteps
            mS( 0 ) = 0.0 ;

            for ( size_t k=1; k<n; ++k )
            {
                mS( k ) = mS( k-1 ) - mH( k-1 );
            }
            // normalize timesteps
            mS /= mH( 0 );


            // populate the vandermonde matrix
            mVandermonde.set_size( n, n );
            for ( size_t j=0; j<n; ++j )
            {
                // first row is always 1
                mVandermonde( 0, j ) = 1.0;

                for ( size_t i=1; i<n; ++i )
                {
                    mVandermonde( i, j ) = mVandermonde( i-1, j ) * mS( j );
                }
            }

            // right-hand side: unit vector selecting the first-derivative row
            mCoefficients.set_size( n, 0.0 );
            mCoefficients( 1 ) = 1.0 ;

            // compute the alpha coefficients using LAPACK
            // note the mVandermonde is overwritten after this call
            gesv( mVandermonde, mCoefficients, mPivot );
        }

        real
        BDF::deval( const ShiftRegister< real > & aY, const bool aUpdateCoefficients )
        {
            BELFEM_ASSERT( aY.capacity() == mH.capacity()+1, "ShiftRegister capacity does not match BDF order"  );

            if ( aUpdateCoefficients ) this->compute_coefficients();

            real aF = 0.0 ;

            for ( size_t k=0; k<aY.size(); ++k )
            {
                aF += mCoefficients( k ) * aY( k );
            }

            aF /= mH( 0 );

            return aF ;
        }

        real
        BDF::eval( ShiftRegister< real > & aY , const real aF, const bool aUpdateCoefficients )
        {
            BELFEM_ASSERT( aY.capacity() == mH.capacity() + 1, "ShiftRegister capacity does not match BDF order"  );

            if ( aUpdateCoefficients ) this->compute_coefficients();

            real tA = mH( 0 ) * aF ;

            for ( size_t k=1; k<aY.size(); ++k )
            {
                tA -= mCoefficients( k ) * aY( k );
            }

            aY( 0 ) = tA / mCoefficients( 0 );

            return aY( 0 );
        }

        const Vector< real > &
        BDF::deval( const ShiftRegister< Vector< real > > & aY, const bool aUpdateCoefficients )
        {
            BELFEM_ASSERT( aY.capacity() == mH.capacity()+1, "ShiftRegister capacity does not match BDF order"  );

            if ( aUpdateCoefficients ) this->compute_coefficients();

            mF.set_size( aY( 0 ).length(), 0.0 );

            for ( size_t k=0; k<aY.size(); ++k )
            {
                mF += mCoefficients( k ) * aY( k );
            }

            mF /= mH( 0 );

            return mF ;
        }

        const Vector< real > &
        BDF::eval( ShiftRegister< Vector< real > > & aY , const Vector< real > & aF, const bool aUpdateCoefficients )
        {
            BELFEM_ASSERT( aY.capacity() == mH.capacity()+1, "ShiftRegister capacity does not match BDF order"  );

            if ( aUpdateCoefficients ) this->compute_coefficients();

            Vector< real > & tA = aY( 0 );

            tA = mH( 0 ) * aF ;

            for ( size_t k=1; k<aY.size(); ++k )
            {
                tA -= mCoefficients( k ) * aY( k );
            }

            aY( 0 )/= mCoefficients( 0 );

            return aY( 0 );
        }

    }
}