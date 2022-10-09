//
// Created by christian on 4/22/22.
//

#include "cl_Material_Silver.hpp"
#include "nist_functions.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_linspace.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Silver::Silver() :
                IsotropicMaterial( MaterialType::Silver )
        {
            mTmax = 1235.0 ;

            mKohlerXmin = nist::extend_kohler( mKohlerA, mKohlerKmin, 4.0, mKohlerB ) ;

            mDensityPoly = { 0 };
            this->use_splines( true );

            this->create_density_poly();
            this->create_conductivity_polys() ;
            this->create_specific_heat_polys();
            this->create_specific_heat_spline() ;
            this->create_mech_polys() ;

            this->set_rrr( mRRR );

            mHasMechanical = true ;
            mHasResistivity = true ;
            mHasThermal = true ;

        }

//----------------------------------------------------------------------------

        Silver::~Silver()
        {
            delete mCpSpline ;
            delete mRhoSpline ;
        }

//----------------------------------------------------------------------------

        void
        Silver::set_rrr( const real aRRR )
        {
            mRRR = aRRR ;
            this->create_resistivity_polys() ;
            this->create_resistivity_spline();
            mRhoRef = this->rho_el0_nist( mTref );
        }

//----------------------------------------------------------------------------

        void
        Silver::create_density_poly()
        {
            // doi:10.1088/1742-6596/1677/1/012161
            mDensityPoly = {  -6.939999845777782e-05,  -0.604610761361504, 10687.00566698209 };
        }

//----------------------------------------------------------------------------

        void
        Silver::create_specific_heat_polys()
        {

            // polynomials based in Smith and Fickett, 10.6028/jres.100.012
            mSpecificHeatPoly1 = {  4.65503e2, -2.92254e2,   6.75689e1,
                                   -7.19040e0,  3.54982e-1, -5.12896e-3,
                                    2.55087e-5 };

            mSpecificHeatPoly3 = {   2.516291e5, -1.885135e4, 4.378403e2,
                                    -1.279849,  4.873041e-3, -9.720354e-6,
                                     8.049353e-9};


            // based on data from Touloukian
            mSpecificHeatPoly5 = { 2.093920e-8, -3.584665e-5, 7.431956e-2,
                                   2.147682e2 };

            // very low temperatue condition
            create_beam_poly( 0.0,
                              0.0,
                              0.0,
                              mSwitchCT0,
                              nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                              mSpecificHeatPoly0 );

            // glue conditions
            create_beam_poly( mSwitchCT1,
                              nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                              mSwitchCT2,
                              nist::cp_janaf( mSpecificHeatPoly3, mSwitchCT2 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT2 ),
                              mSpecificHeatPoly2 );

            create_beam_poly( mSwitchCT3,
                              nist::cp_janaf( mSpecificHeatPoly3, mSwitchCT3 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT3 ),
                              mSwitchCT4,
                              polyval( mSpecificHeatPoly5, mSwitchCT4 ),
                              dpolyval( mSpecificHeatPoly5, mSwitchCT4 ),
                              mSpecificHeatPoly4 );


        }

//--------------------------------------------------------------------------

        void
        Silver::create_specific_heat_spline()
        {
            uint tN = 701 ;
            real tDeltaT = 2.0 ;

            // populate data
            Vector< real > tT( tN );
            Vector< real > tC( tN );

            tT( 0 ) = 0.0 ;
            tC( 0 ) = 0.0 ;

            for( uint k=1; k<tN; ++k )
            {
                tT( k ) = tT(k-1) + tDeltaT ;
                tC( k ) = this->c_poly( tT( k ) );
            }

            SpMatrix tA;
            spline::create_helpmatrix( tN, tDeltaT, tA  );
            mCpSpline = new Spline( tT, tC, tA );

        }

//----------------------------------------------------------------------------

        void
        Silver::create_resistivity_polys()
        {
            real tScale = 1e12 ;
            uint tN = 101 ;

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 4 K <= T < 10 K
            // - - - - - - - - - - - - - - - - - - - -

            Vector< real > tT( tN );
            linspace( mSwitchRT0,mSwitchRT1 , tN, tT );
            tT( 0 ) = BELFEM_EPSILON ;
            Vector< real > tR( tN );
            for( uint i=0; i<tN; ++i )
            {
                tR( i ) = nist::rho0( mPrho, mRRR, tT( i ) );
            }
            tR *= tScale ;
            polyfit( tT, tR, 3, mResistivityPoly1 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 0 K <= T < 4 K
            // - - - - - - - - - - - - - - - - - - - -
            nist::extend_resistivity( mResistivityPoly1, mSwitchRT0, mResistivityPoly0 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 15 K <= T < 50 K
            // - - - - - - - - - - - - - - - - - - - -

            linspace( mSwitchRT2,mSwitchRT3 , tN, tT );
            for( uint i=0; i<tN; ++i )
            {
                tR( i ) = nist::rho0( mPrho, mRRR, tT( i ) );
            }
            tR *= tScale ;
            polyfit( tT, tR, 3, mResistivityPoly3 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 10 K <= T < 15 K
            // - - - - - - - - - - - - - - - - - - - -

            create_beam_poly( mSwitchRT4,
                              polyval( mResistivityPoly1, mSwitchRT4 ),
                              dpolyval( mResistivityPoly1, mSwitchRT4 ),
                              mSwitchRT2,
                              polyval( mResistivityPoly3, mSwitchRT2 ),
                              dpolyval( mResistivityPoly3, mSwitchRT2 ),
                              mResistivityPoly2 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 75 K <= T < T max K
            // - - - - - - - - - - - - - - - - - - - -

            linspace( mSwitchRT4,mTmax , tN, tT );
            for( uint i=0; i<tN; ++i )
            {
                tR( i ) = nist::rho0( mPrho, mRRR, tT( i ) );
            }
            tR *= tScale ;
            polyfit( tT, tR, 2, mResistivityPoly5 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 50 K <= T < 75 max K
            // - - - - - - - - - - - - - - - - - - - -

            create_beam_poly( mSwitchRT3,
                              polyval( mResistivityPoly3, mSwitchRT3 ),
                              dpolyval( mResistivityPoly3, mSwitchRT3 ),
                              mSwitchRT4,
                              polyval( mResistivityPoly5, mSwitchRT4 ),
                              dpolyval( mResistivityPoly5, mSwitchRT4 ),
                              mResistivityPoly4 );

            // backscale
            mResistivityPoly0 /= tScale ;
            mResistivityPoly1 /= tScale ;
            mResistivityPoly2 /= tScale ;
            mResistivityPoly3 /= tScale ;
            mResistivityPoly4 /= tScale ;
            mResistivityPoly5 /= tScale ;
        }

//--------------------------------------------------------------------------

        void
        Silver::create_resistivity_spline()
        {
            uint tN = 701 ;
            real tDeltaT = 2.0 ;

            // populate data
            Vector< real > tT( tN );

            tT( 0 ) = 0.0 ;
            for( uint k=1; k<tN; ++k )
            {
                tT( k ) = tT(k-1) + tDeltaT ;
            }

            SpMatrix tA;
            spline::create_helpmatrix( tN, tDeltaT, tA  );
            Vector< real > tR( tN );

            for( uint k=0; k<tN; ++k )
            {
                tR( k ) = this->rho_el0_nist( tT( k ) );
            }

            if( mRhoSpline != nullptr )
            {
                delete mRhoSpline ;
            }

            mRhoSpline = new Spline( tT, tR, tA );
        }

//----------------------------------------------------------------------------

        void
        Silver::create_mech_polys()
        {
            // Dataset from Smith and Fickett
            Vector< real > tT1 = { 1.176, 16.614, 31.42, 44.943, 57.827,
                                   68.136, 79.737, 91.338, 101.653, 113.897,
                                   125.498, 136.456, 148.057, 159.658,
                                   171.259, 182.86, 194.461, 205.42, 217.02,
                                   228.622, 240.223, 251.823, 263.424,
                                   275.025, 285.983, 297.584 } ;

            Vector< real > tE1 = { 91.2484, 91.1831, 90.9606, 90.6931, 90.3582,
                                   90.0678, 89.6879, 89.3079, 88.9278, 88.5479,
                                   88.1679, 87.7878, 87.4079, 87.0279, 86.6479,
                                   86.268, 85.888, 85.4855, 85.128, 84.7256,
                                   84.3456, 83.9881, 83.6081, 83.2281,
                                   82.8481, 82.4681};
            tE1 *= 1e9;

            polyfit( tT1, tE1, 4, mYoungPoly0 );

            // dataset from Blanke
            Vector< real > tT3 = { 97.44, 181.75, 262.42, 354.95, 438.31,
                                   579.35, 636.12, 748.69, 812.76, 885.05,
                                   958.24, 1027.73, 1093.54, 1159.31 };

            Vector< real > tE3 = { 89.733, 87.154, 84.832, 81.61, 78.583,
                                   73.047, 70.666, 65.391, 62.368, 58.831,
                                   54.974, 50.861, 46.75, 42.126};

            tE3 *= 1E9 ;
            polyfit( tT3, tE3, 4, mYoungPoly2 );

            create_beam_poly( mSwitchET0,
                              polyval( mYoungPoly0, mSwitchET0 ),
                              dpolyval( mYoungPoly0, mSwitchET0 ),
                              mSwitchET1,
                              polyval( mYoungPoly2, mSwitchET1 ),
                              dpolyval( mYoungPoly2, mSwitchET1 ),
                              mYoungPoly1 );

            // based on dataset from Smith and Fickett
            mNuPoly = { 2.6457e-5, 0.35969 };
        }

//----------------------------------------------------------------------------

        real
        Silver::c_poly( const real aT ) const
        {
            if( aT < mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if ( aT < mSwitchCT1 )
            {
                return nist::cp_janaf( mSpecificHeatPoly1, aT );
            }
            else if ( aT < mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if ( aT < mSwitchCT3 )
            {
                return nist::cp_janaf( mSpecificHeatPoly3, aT );
            }
            else if ( aT < mSwitchCT4 )
            {
                return polyval( mSpecificHeatPoly4, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly5, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Silver::E( const real aT ) const
        {
            if( aT < mSwitchCT0 )
            {
                return polyval( mYoungPoly0, aT );
            }
            else if ( aT < mSwitchCT1 )
            {
                return polyval( mYoungPoly1, aT );
            }
            else
            {
                return polyval( mYoungPoly2, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Silver::nu( const real aT ) const
        {
            return polyval( mNuPoly, aT );
        }

//----------------------------------------------------------------------------

        real
        Silver::G( const real aT ) const
        {
            return this->E( aT ) / ( 2.0 * ( 1 + this->nu( aT ) ) );
        }

//----------------------------------------------------------------------------

        real
        Silver::rho_el0_poly( const real aT ) const
        {
            if( aT < mSwitchRT0 )
            {
                return polyval( mResistivityPoly0, aT );
            }
            else if( aT < mSwitchRT1 )
            {
                return polyval( mResistivityPoly1, aT );
            }
            else if( aT < mSwitchRT2 )
            {
                return polyval( mResistivityPoly2, aT );
            }
            else if( aT < mSwitchRT3 )
            {
                return polyval( mResistivityPoly3, aT );
            }
            else if( aT < mSwitchRT4 )
            {
                return polyval( mResistivityPoly4, aT );
            }
            else
            {
                return polyval( mResistivityPoly5, aT );
            }
        }

//----------------------------------------------------------------------------

        void
        Silver::create_conductivity_polys()
        {

            // poly for w1, T< 50
            mThermalConductivityPoly1 = { -7.5168e-9, 1.0500e-6, -1.9953e-6, 0 } ;
            mThermalConductivityPoly3 = { 9.028e-7, 1.999e-3};

            // create connecting poly
            create_beam_poly(
                    mSwitchLT0,
                    polyval( mThermalConductivityPoly1, mSwitchLT0 ),
                    dpolyval( mThermalConductivityPoly1, mSwitchLT0 ),
                    mSwitchLT1,
                    polyval( mThermalConductivityPoly3, mSwitchLT1 ),
                    dpolyval( mThermalConductivityPoly3, mSwitchLT1 ),
                    mThermalConductivityPoly2 );
        }

//----------------------------------------------------------------------------

        // electric resistance
        real
        Silver::rho_el ( const real aJ, const real aT, const real aB, const real aAngle ) const
        {
            // resistivity for zero magnetic field
            real tRho0 = this->rho_el0( aT );

            return tRho0 * ( 1.0 + nist::kohler( mKohlerA, mKohlerB,
                                                 mKohlerXmin, mRhoRef / tRho0 * std::abs( aB ) ) );
        }

//----------------------------------------------------------------------------

        real
        Silver::lambda( const real aT ) const
        {

            real tW0 = 0.7323 / ( mRRR * aT );
            real tW1 = aT < mSwitchLT0 ?
                    polyval( mThermalConductivityPoly1, aT ) :
                       ( aT < mSwitchLT1 ?
                         polyval( mThermalConductivityPoly2, aT ) :
                         polyval( mThermalConductivityPoly3, aT ) ) ;

            return 1.0 / ( tW0 + tW1 );
        }

//----------------------------------------------------------------------------

        real
        Silver::lambda( const real aT, const real aB ) const
        {

            return this->lambda( aT )  * ( 1.0 + nist::kohler( mKohlerA, mKohlerB, mKohlerXmin,
                                                               mRhoRef / this->rho_el0_spline( aT )  * std::abs( aB ) ) );
        }

//----------------------------------------------------------------------------

        real
        Silver::rho_el0_nist( const real aT ) const
        {
            if( aT < mSwitchRT0 )
            {
                return polyval( mResistivityPoly0, aT );
            }
            else
            {
                return nist::rho0( mPrho, mRRR, aT );
            }
        }
//----------------------------------------------------------------------------

        void
        Silver::use_splines( const bool aSwitch )
        {
            if( aSwitch )
            {
                mCpFunction   = & Silver::c_spline ;
                mRho0Function = & Silver::rho_el0_spline ;
            }
            else
            {
                mCpFunction = & Silver::c_poly ;
                mRho0Function = & Silver::rho_el0_nist ;
            }
        }

//----------------------------------------------------------------------------
    }
}