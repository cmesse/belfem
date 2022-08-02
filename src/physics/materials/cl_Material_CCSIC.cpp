//
// Created by Christian Messe on 14.04.20.
//
#include "cl_Material_CCSIC.hpp"
#include "fn_polyval.hpp"
#include "cl_Spline.hpp"
#include "fn_linspace.hpp"
#include "cl_SpMatrix.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        CCSIC::CCSIC() :
                OrthotropicMaterial( MaterialType::CCSIC )
        {
            this->create_heat_spline();

            // thermal conductivity
            // source: 10.18419/opus-9381
            mLambda1 = { -2.2075e-6, 4.9853e-3, 1.5640e1 };
            mLambda3 = { -1.0737e-03, 9.1618 };

            // optical emmittance
            // source: AST5-CT-2006-030729
            mEpsilon = { 5.030e-11, -4.3510e-07, 9.2540e-04, 0.33 };

            // for tangential expansion factor
            // source: 10.18419/opus-9381
            mMu1 = { 4.88398e-23, -5.55961e-19, 2.55950e-15, -6.08651e-12, 8.56421e-9, -3.66672e-6, 0.0e0 };
            mMu3 = { 4.04275e-16, -1.83263e-12, 4.37780e-9, 3.04970e-6, 0.0e0 };

        }
//----------------------------------------------------------------------------

        CCSIC::~CCSIC()
        {
            // delete the heatspline
            delete mHeat;
        }
//----------------------------------------------------------------------------

        void
        CCSIC::create_heat_spline()
        {
            // polynomial for T > 100 K
            // source: 10.18419/opus-9381
            Vector <real> tCp1 = {
                    1.14724476175820e-29,
                    -2.05114046370026e-25,
                    1.60007993165233e-21,
                    -7.13605193838964e-18,
                    2.00057453457363e-14,
                    -3.64192547014549e-11,
                    4.26478682127158e-08,
                    -2.99969103955274e-05,
                    9.64233891564731e-03,
                    1.54004632030134e00,
                    -8.80717557838845e01 };

            Vector <real> tCp0 = {
                    -3.44485322804155e-7,
                    1.36623742541241e-2,
                    0.0,
                    0.0
            };

            // correct interface
            real tCpA = polyval( tCp1, 100.0 );
            real tCpB = polyval( tCp0, 100.0 );
            tCp0( 3 ) = tCpA - tCpB;

            // steps for Temperatures
            uint tN = 81;
            Vector <real> tT = linspace( 0.0, 4000.0, tN );

            Vector <real> tCp( tN );

            // populate vector
            for ( uint k = 0; k < tN; ++k )
            {
                if ( tT( k ) < 100.0 )
                {
                    tCp( k ) = polyval( tCp0, tT( k ));
                }
                else
                {
                    tCp( k ) = polyval( tCp1, tT( k ));
                }
            }

            // create help matrix
            SpMatrix tHelpMatrix;

            // initialize help matrix
            spline::create_helpmatrix( tN, 50.0, tHelpMatrix );

            mHeat = new Spline( tT, tCp, tHelpMatrix );
        }

//----------------------------------------------------------------------------

        /**
          * Elasticity Matrix in 3D
          */
        void
        CCSIC::C( Matrix <real> & aC, const real aT ) const
        {
            this->compute_C3D( mE1,
                               mE2,
                               mE3,
                               mNu23,
                               mNu12,
                               mNu31,
                               mG23,
                               mG12,
                               mG13,
                               aC );
        }


//----------------------------------------------------------------------------

        /**
         * Elasticity matrix in plane stress
         */
        void
        CCSIC::C_ps( Matrix <real> & aC, const real aT ) const
        {
            this->compute_C_ps(
                    mE1,
                    mE3,
                    mNu13,
                    mG13,
                    aC );
        }

//----------------------------------------------------------------------------

        /**
         * Elasticity matrix in rotation symmetry
         */
        void
        CCSIC::C_rot( Matrix <real> & aC, const real aT ) const
        {
            this->compute_C_rot(
                    mE1,
                    mE2,
                    mE3,
                    mNu23,
                    mNu12,
                    mNu31,
                    mG12,
                    aC );
        }

//----------------------------------------------------------------------------

        real
        CCSIC::rho( const real aT ) const
        {
            return mRho;
        }

//----------------------------------------------------------------------------

        real
        CCSIC::c( const real aT ) const
        {
            return mHeat->eval( aT );
        }

//----------------------------------------------------------------------------

        void
        CCSIC::lambda_p( Matrix <real> & aLambda, const real aT ) const
        {
            aLambda( 0, 0 ) = polyval( mLambda1, aT );
            aLambda( 1, 0 ) = 0.0;
            aLambda( 0, 1 ) = 0.0;
            aLambda( 1, 1 ) = polyval( mLambda3, aT );
        }

//----------------------------------------------------------------------------

        void
        CCSIC::lambda_3d( Matrix <real> & aLambda, const real aT ) const
        {
            real tLambda1 = polyval( mLambda1, aT );

            aLambda( 0, 0 ) = tLambda1;
            aLambda( 1, 0 ) = 0.0;
            aLambda( 2, 0 ) = 0.0;

            aLambda( 0, 1 ) = 0.0;
            aLambda( 1, 1 ) = tLambda1;
            aLambda( 2, 1 ) = 0.0;

            aLambda( 0, 2 ) = 0.0;
            aLambda( 1, 2 ) = 0.0;
            aLambda( 2, 2 ) = polyval( mLambda3, aT );
        }

//----------------------------------------------------------------------------

        real
        CCSIC::epsilon( const real aT ) const
        {
            return polyval( mEpsilon, aT );
        }

//----------------------------------------------------------------------------

        void
        CCSIC::mu_ps( Matrix <real> & aMu, const real aT, const real aTref ) const
        {
            BELFEM_ASSERT( aMu.n_cols() == 2 && aMu.n_rows() == 2,
                          "Expansion matrix must be allocated as 2x2" );

            aMu( 0, 0 ) = std::exp( polyval( mMu1, aT ) - polyval( mMu1, aTref )) - 1.0;
            aMu( 1, 0 ) = 0.0;
            aMu( 0, 1 ) = 0.0;
            aMu( 1, 1 ) = std::exp( polyval( mMu3, aT ) - polyval( mMu3, aTref )) - 1.0;
        }

//----------------------------------------------------------------------------

        void
        CCSIC::mu( Matrix <real> & aMu, const real aT, const real aTref ) const
        {
            BELFEM_ASSERT( aMu.n_cols() == 3 && aMu.n_rows() == 3,
                          "Expansion matrix must be allocated as 3x3" );

            aMu( 0, 0 ) = std::exp( polyval( mMu1, aT ) - polyval( mMu1, aTref )) - 1.0;
            aMu( 1, 0 ) = 0.0;
            aMu( 2, 0 ) = 0.0;

            aMu( 0, 1 ) = 0.0;
            aMu( 1, 1 ) = aMu( 0, 0 );
            aMu( 2, 1 ) = 0.0;

            aMu( 0, 2 ) = 0.0;
            aMu( 1, 2 ) = 0.0;
            aMu( 2, 2 ) = std::exp( polyval( mMu3, aT ) - polyval( mMu3, aTref )) - 1.0;

        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */