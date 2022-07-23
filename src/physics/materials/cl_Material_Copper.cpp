//
// Created by Christian Messe on 29.07.20.
//

#include "cl_Material_Copper.hpp"

#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_linspace.hpp"
#include "cl_SpMatrix.hpp"
#include "nist_functions.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Copper::Copper() :
                IsotropicMaterial( MaterialType::Copper )
        {
            // set maximum temperature
            mTmax = 1357.15;
            mKohlerXmin = nist::extend_kohler( mKohlerA, mKohlerKmin, 8.0, mKohlerB ) ;

            this->create_density_poly();
            this->create_mech_polys();
            this->create_expansion_polys();
            this->create_specific_heat_polys();
            this->create_specific_heat_spline();

            this->set_rrr( mRRR );




            mHasThermal = true;
            mHasMechanical = true;
            mHasExpansion = true;
            mHasResistivity = true ;
        }

//----------------------------------------------------------------------------

        Copper::~Copper()
        {
            delete mCpSpline ;
            delete mRhoSpline ;

        }
//----------------------------------------------------------------------------

        void
        Copper::set_rrr( const real aRRR )
        {
            mRRR = aRRR ;

            //this->create_resistivity_polys() ;
            this->create_resistivity_spline() ;
            mRhoRef = this->rho_el0_spline( mTref );

            this->create_conductivity_polys() ;
        }


//----------------------------------------------------------------------------

        void
        Copper::create_specific_heat_polys()
        {
            // from CUDI for T < 9.441 K
            mSpecificHeatPoly0 = { -0.0308, 7.229 , - 2.1286, 101.89, 2.5631 } ;
            mSpecificHeatPoly0 /= 8920.0 ; // <-- density used in CUDI

           /*
            // dataset from MATPRO
            //  SPECIFIC HEAT OF COPPER
            // (SEE "BUBBLE CHAMBER GROUP-DATA HANDBOOK",VIII B.1-B.5)

            DATA TH  / 10., 12., 15., 16., 18., 20., 25., 30., 40., 50., 60.,&
                    70., 80., 90.,100.,120.,140.,160.,180.,200.,220.,240.,&
                    260.,280.,300.,500.,1000./
            DATA CAL / .86,1.42, 2.7, 3.4, 5. , 7.7, 16., 27., 60., 99.,137.,&
                    173.,205.,232.,254.,288.,313.,332.,346.,356.,364.,371.,&
                    376.,381.,386.,408.,408./
            */

            Vector< real > tT2  = { 10., 12., 15., 16., 18., 20., 25., 30., 40., 50., 60.  };
            Vector< real > tCp2 = { .86,1.42, 2.7, 3.4, 5. , 7.7, 16., 27., 60., 99., 137. };
            polyfit( tT2, tCp2, 4, mSpecificHeatPoly2 );

            // create connecting poly
            create_beam_poly(
                    mSwitchCT0,
                    polyval( mSpecificHeatPoly0, mSwitchCT0 ),
                    dpolyval( mSpecificHeatPoly0, mSwitchCT0 ),
                    mSwitchCT1,
                    polyval( mSpecificHeatPoly2, mSwitchCT1 ),
                    dpolyval( mSpecificHeatPoly2, mSwitchCT1 ),
                    mSpecificHeatPoly1 );

            Vector< real > tT4 = { 90., 100.,120.,140.,160.,180.,200.,220. };
            Vector< real > tCp4 = { 232., 254.,288.,313.,332.,346.,356.,364. };
            polyfit( tT4, tCp4, 3, mSpecificHeatPoly4 );

            // create connecting poly
            create_beam_poly(
                    mSwitchCT2,
                    polyval( mSpecificHeatPoly2, mSwitchCT2 ),
                    dpolyval( mSpecificHeatPoly2, mSwitchCT2 ),
                    mSwitchCT3,
                    polyval( mSpecificHeatPoly4, mSwitchCT3 ),
                    dpolyval( mSpecificHeatPoly4, mSwitchCT3 ),
                    mSpecificHeatPoly3 );

            Vector< real > tT6 = { 200.,220.,240., 260.,280.,300. };
            Vector< real > tC6 = { 356.,364.,371., 376.,381.,386., };

            polyfit( tT6, tC6, 3, mSpecificHeatPoly6 );

            // create connecting poly
            create_beam_poly(
                    mSwitchCT4,
                    polyval( mSpecificHeatPoly4, mSwitchCT4 ),
                    dpolyval( mSpecificHeatPoly4, mSwitchCT4 ),
                    mSwitchCT5,
                    polyval( mSpecificHeatPoly6, mSwitchCT5 ),
                    dpolyval( mSpecificHeatPoly6, mSwitchCT5 ),
                    mSpecificHeatPoly5 );

            // data from
            // Dinsdale, A T: SGTE data for pure elements. CALPHAD 15(1991)317/425.
            Vector <real> tT8 = { 295.05, 371.42, 474.53, 572.69, 672.81,
                                  772.92, 874.02, 971.15, 1071.26, 1169.39,
                                  1271.46, 1355.71 };

            Vector <real> tC8 = { 384.8, 397.0, 407.8, 419.0, 426.8, 434.0,
                                  440.8, 446.6, 452.8, 459.8, 463.6, 468.6 };

            polyfit( tT8, tC8, 3, mSpecificHeatPoly8 );

            // create connecting poly
            create_beam_poly(
                    mSwitchCT6,
                    polyval( mSpecificHeatPoly6, mSwitchCT6 ),
                    dpolyval( mSpecificHeatPoly6, mSwitchCT6 ),
                    mSwitchCT7,
                    polyval( mSpecificHeatPoly8, mSwitchCT7 ),
                    dpolyval( mSpecificHeatPoly8, mSwitchCT7 ),
                    mSpecificHeatPoly7 );

            mSwitchCT8 = -mSpecificHeatPoly8( 1 ) / ( 3. * mSpecificHeatPoly8( 0 ) );

            mSpecificHeatPoly9.set_size( 2 );
            mSpecificHeatPoly9( 0 ) = dpolyval( mSpecificHeatPoly8, mSwitchCT8 );
            mSpecificHeatPoly9( 1 ) = polyval( mSpecificHeatPoly8, mSwitchCT8 )
                    -  mSpecificHeatPoly9( 0 ) * mSwitchCT8;

        }
//----------------------------------------------------------------------------

        void
        Copper::create_specific_heat_spline()
        {
            uint tN = 701 ;
            real tDeltaT = 2.0 ;

            SpMatrix tA;
            spline::create_helpmatrix( tN, tDeltaT, tA );

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

            mCpSpline = new Spline( tT, tC, tA, 0, 0, 0 );
        }

//----------------------------------------------------------------------------

        void
        Copper::create_conductivity_polys()
        {
            mBeta = mRhoRef / ( mPlambda( 0 ) * mRRR );
            mPlambda( 7 ) = mPlambda( 8 ) / std::pow( mBeta / mPlambda( 9 ), mPlambda( 10 ) );



            /*
            // Touloukian
            Vector <real> tT6 = { 298.58, 375.36, 474.38, 576.4, 674.4, 773.42, 874.43,
                                  974.44, 1073.45, 1171.43, 1271.45, 1356.28 };
            Vector <real> tL6 = { 399.9, 394.9, 389.5, 382.0, 376.5, 370.6, 363.1,
                                  356.7, 350.7, 343.3, 337.3, 329.8 }; */


            // compute the temperatures for switching the lookup tables
            mSwitchLambdaT1 = this->compute_T_lambda_peak() ;
            mSwitchLambdaT0 = mSwitchLambdaT1 *  0.9 ;

            mSwitchLambdaT2 = mSwitchLambdaT1 *  1.1 ;
            mSwitchLambdaT3 = mSwitchLambdaT1 *  3.0 ;

            mSwitchLambdaT4 = mSwitchLambdaT1 *  5.0 ;
            mSwitchLambdaT5 = mSwitchLambdaT1 *  8.0 ;

            mSwitchLambdaT6 = mSwitchLambdaT1 * 10.0 ;
            mSwitchLambdaT7 = mSwitchLambdaT1 * 15.0 ;

            mSwitchLambdaT8 = mSwitchLambdaT1 * 20.0 ;

            // container for datasets
            uint tN = 100 ;
            Vector< real > tT( tN );
            Vector< real > tK( tN );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 0
            // - - - - - - - - - - - - - - - - - - - -

            // poplulate dataset
            linspace( BELFEM_EPSILON, mSwitchLambdaT0, tN, tT );
            for( uint k=0; k<tN; ++k )
            {
                tK( k ) = nist::lambda0( mPlambda, mBeta, tT( k ));
            }

            // create polynomial
            polyfit( tT, tK, 3, mThermalConductivityPoly0 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 3
            // - - - - - - - - - - - - - - - - - - - -

            // poplulate dataset
            linspace( mSwitchLambdaT2, mSwitchLambdaT3, tN, tT );
            for( uint k=0; k<tN; ++k )
            {
                tK( k ) = nist::lambda0( mPlambda, mBeta, tT( k ));
            }

            // create polynomial
            polyfit( tT, tK, 6, mThermalConductivityPoly3 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 5
            // - - - - - - - - - - - - - - - - - - - -

            // poplulate dataset
            linspace( mSwitchLambdaT4, mSwitchLambdaT5, tN, tT );
            for( uint k=0; k<tN; ++k )
            {
                tK( k )= nist::lambda0( mPlambda, mBeta, tT( k ));
            }

            // create polynomial
            polyfit( tT, tK, 4, mThermalConductivityPoly5 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 7
            // - - - - - - - - - - - - - - - - - - - -

            // poplulate dataset
            linspace( mSwitchLambdaT6, mSwitchLambdaT7, tN, tT );
            for( uint k=0; k<tN; ++k )
            {
                tK( k ) = nist::lambda0( mPlambda, mBeta, tT( k ));
            }

            // create polynomial
            polyfit( tT, tK, 4, mThermalConductivityPoly7 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 9
            // - - - - - - - - - - - - - - - - - - - -

            // poplulate dataset
            linspace( mSwitchLambdaT8, std::max( mTmax, 21. * mSwitchLambdaT1 ), tN, tT );
            for( uint k=0; k<tN; ++k )
            {
                tK( k ) = nist::lambda0( mPlambda, mBeta, tT( k ));
            }

            // create polynomial
            polyfit( tT, tK, 3, mThermalConductivityPoly9 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 1
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchLambdaT0,
                    polyval( mThermalConductivityPoly0, mSwitchLambdaT0 ),
                    dpolyval( mThermalConductivityPoly0, mSwitchLambdaT0 ),
                    mSwitchLambdaT1,
                    nist::lambda0( mPlambda, mBeta, mSwitchLambdaT1),
                    0.0,
                    mThermalConductivityPoly1 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 2
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchLambdaT1,
                    nist::lambda0( mPlambda, mBeta, mSwitchLambdaT1),
                    0.0,
                    mSwitchLambdaT2,
                    polyval( mThermalConductivityPoly3, mSwitchLambdaT2 ),
                    dpolyval( mThermalConductivityPoly3, mSwitchLambdaT2 ),
                    mThermalConductivityPoly2 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 4
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchLambdaT3,
                    polyval( mThermalConductivityPoly3, mSwitchLambdaT3 ),
                    dpolyval( mThermalConductivityPoly3, mSwitchLambdaT3 ),
                    mSwitchLambdaT4,
                    polyval( mThermalConductivityPoly5, mSwitchLambdaT4 ),
                    dpolyval( mThermalConductivityPoly5, mSwitchLambdaT4 ),
                    mThermalConductivityPoly4 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 6
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchLambdaT5,
                    polyval( mThermalConductivityPoly5, mSwitchLambdaT5 ),
                    dpolyval( mThermalConductivityPoly5, mSwitchLambdaT5 ),
                    mSwitchLambdaT6,
                    polyval( mThermalConductivityPoly7, mSwitchLambdaT6 ),
                    dpolyval( mThermalConductivityPoly7, mSwitchLambdaT6 ),
                    mThermalConductivityPoly6 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval 8
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchLambdaT7,
                    polyval( mThermalConductivityPoly7, mSwitchLambdaT7 ),
                    dpolyval( mThermalConductivityPoly7, mSwitchLambdaT7 ),
                    mSwitchLambdaT8,
                    polyval( mThermalConductivityPoly9, mSwitchLambdaT8 ),
                    dpolyval( mThermalConductivityPoly9, mSwitchLambdaT8 ),
                    mThermalConductivityPoly8 );
        }

//----------------------------------------------------------------------------

        real
        Copper::c_poly( const real aT ) const
        {
            if( aT < mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if( aT < mSwitchCT1 )
            {
                return polyval( mSpecificHeatPoly1, aT );
            }
            else if( aT < mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if( aT < mSwitchCT3 )
            {
                return polyval( mSpecificHeatPoly3, aT );
            }
            else if( aT < mSwitchCT4 )
            {
                return polyval( mSpecificHeatPoly4, aT );
            }
            else if( aT < mSwitchCT5 )
            {
                return polyval( mSpecificHeatPoly5, aT );
            }
            else if( aT < mSwitchCT6 )
            {
                return polyval( mSpecificHeatPoly6, aT );
            }
            else if( aT < mSwitchCT7 )
            {
                return polyval( mSpecificHeatPoly7, aT );
            }
            else if( aT < mSwitchCT8 )
            {
                return polyval( mSpecificHeatPoly8, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly9, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Copper::lambda( const real aT ) const
        {
            if( aT < mSwitchLambdaT0 )
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
            else if( aT < mSwitchLambdaT1 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else if( aT < mSwitchLambdaT2 )
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
            else if( aT < mSwitchLambdaT3 )
            {
                return polyval( mThermalConductivityPoly3, aT );
            }
            else if( aT < mSwitchLambdaT4 )
            {
                return polyval( mThermalConductivityPoly4, aT );
            }
            else if( aT < mSwitchLambdaT5 )
            {
                return polyval( mThermalConductivityPoly5, aT );
            }
            else if( aT < mSwitchLambdaT6 )
            {
                return polyval( mThermalConductivityPoly6, aT );
            }
            else if( aT < mSwitchLambdaT7 )
            {
                return polyval( mThermalConductivityPoly7, aT );
            }
            else if( aT < mSwitchLambdaT8 )
            {
                return polyval( mThermalConductivityPoly8, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly9, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Copper::lambda( const real aT, const real aB ) const
        {
            return this->lambda( aT )  * ( 1.0 + nist::kohler( mKohlerA, mKohlerB, mKohlerXmin,
                                                               mRhoRef / this->rho_el0_spline( aT )  * std::abs( aB ) ) );
        }

//----------------------------------------------------------------------------

        // electric resistance
        real
        Copper::rho_el ( const real aJ, const real aT, const real aB ) const
        {
            // resistivity for zero magnetic field
            real tRho0 = this->rho_el0_spline( aT );

            return tRho0 * ( 1.0 + nist::kohler( mKohlerA, mKohlerB, mKohlerXmin, mRhoRef / tRho0 * std::abs( aB ) ) );
        }

//----------------------------------------------------------------------------

        void
        Copper::create_density_poly()
        {
            // Henderson, J B; Hagemann, L and Blumm, J: Thermophysical properties of copper:
            // Netzsch Report 820.030/96 TPS No 4; Netzsch GmbH SeIb9 Germany, March 1996.

            Vector <real> tT = { 296.02, 374.94, 472.87, 574.78, 674.66, 775.56, 872.45,
                                 974.34, 1074.21, 1173.08, 1271.95, 1357.82 };
            Vector <real> tR = { 8929.96, 8890.37, 8848.47, 8804.24, 8738.98, 8687.74,
                                 8627.15, 8568.91, 8498.98, 8431.39, 8361.46, 8291.51 };

            polyfit( tT, tR, 2, mDensityPoly );
        }

//----------------------------------------------------------------------------

        void
        Copper::create_mech_polys()
        {

            // from doi : 10.1002/pssa.2210660209
            Vector <real> tT0 = { 5., 20., 40., 60., 80., 100., 120., 140., 160.,
                                  180., 200., 220., 240., 260., 280., 295. };

            Vector <real> tE0 = { 138.62, 138.57, 138.4, 137.95, 137.3, 136.55, 135.71,
                                  134.9, 134.03, 133.12, 132.33, 131.43, 130.61, 129.8,
                                  128.9, 128.17 };

            Vector <real> tG0 = { 51.72, 51.7, 51.63, 51.45, 51.19, 50.89, 50.55, 50.23,
                                  49.88, 49.52, 49.21, 48.86, 48.53, 48.21, 47.86, 47.57 };

            // fix units
            tE0 *= 1e9;
            tG0 *= 1e9;

            // create cryo poly
            polyfit( tT0, tE0, 4, mYoungPoly0 );
            polyfit( tT0, tG0, 4, mShearPoly0 );


            // based on doi: 10.1063/1.3253150 with extrapolation of nu from 10.1002/pssa.2210660209
            mYoungPoly2 = { -6.44239901E+07, 1.47580491E+11 };
            mShearPoly2 = { -2.47434399E+07, 5.43180029E+10 };

            create_beam_poly(
                    mSwitchMechT0,
                    polyval( mYoungPoly0, mSwitchMechT0 ),
                    dpolyval( mYoungPoly0, mSwitchMechT0 ),
                    mSwitchMechT1,
                    polyval( mYoungPoly2, mSwitchMechT1 ),
                    dpolyval( mYoungPoly2, mSwitchMechT1 ),
                    mYoungPoly1 );

            create_beam_poly(
                    mSwitchMechT0,
                    polyval( mShearPoly0, mSwitchMechT0 ),
                    dpolyval( mShearPoly0, mSwitchMechT0 ),
                    mSwitchMechT1,
                    polyval( mShearPoly2, mSwitchMechT1 ),
                    dpolyval( mShearPoly2, mSwitchMechT1 ),
                    mShearPoly1 );
        }
//----------------------------------------------------------------------------

        void
        Copper::create_expansion_polys()
        {
            // from  doi: 10.1063/1.1658614
            mThermalExpansionPoly = { 4.6314475E-14, -8.7810700E-12, 6.1242350E-08, 0.0, 0.0 };

        }

//----------------------------------------------------------------------------

        real
        Copper::E( const real aT ) const
        {
            if ( aT > mSwitchMechT1 )
            {
                return polyval( mYoungPoly2, aT );
            }
            else if ( aT > mSwitchMechT0 )
            {
                return polyval( mYoungPoly1, aT );
            }
            else
            {
                return polyval( mYoungPoly0, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Copper::G( const real aT ) const
        {
            if ( aT > mSwitchMechT1 )
            {
                return polyval( mShearPoly2, aT );
            }
            else if ( aT > mSwitchMechT0 )
            {
                return polyval( mShearPoly1, aT );
            }
            else
            {
                return polyval( mShearPoly0, aT );
            }
        }

//---------------------------------------------------------------------------

        void
        Copper::create_resistivity_polys()
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
            polyfit( tT, tR, 5, mResistivityPoly1 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 0 K <= T < 4 K
            // - - - - - - - - - - - - - - - - - - - -
            nist::extend_resistivity( mResistivityPoly1, mSwitchRT0, mResistivityPoly0 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 15 K <= T < 50 K
            // - - - - - - - - - - - - - - - - - - - -
            linspace( mSwitchRT2, mSwitchRT3, tN, tT );
            for( uint i=0; i<tN; ++i )
            {
                tR( i ) =  nist::rho0( mPrho, mRRR, tT( i ) );
            }
            tR *= tScale ;
            polyfit( tT, tR, 4, mResistivityPoly3 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 10 K <= T < 15 K
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchRT1,
                    polyval( mResistivityPoly1, mSwitchRT1 ),
                    dpolyval( mResistivityPoly1, mSwitchRT1 ),
                    mSwitchRT2,
                    polyval( mResistivityPoly3, mSwitchRT2 ),
                    dpolyval( mResistivityPoly3, mSwitchRT2 ),
                    mResistivityPoly2 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 75 K <= T < 200 K
            // - - - - - - - - - - - - - - - - - - - -
            linspace( mSwitchRT4, mSwitchRT5, tN, tT );
            tT( 0 ) = BELFEM_EPSILON ;
            for( uint i=0; i<tN; ++i )
            {
                tR( i ) = nist::rho0( mPrho, mRRR, tT( i ) );
            }
            tR *= tScale ;
            polyfit( tT, tR, 4, mResistivityPoly5 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 50 K <= T < 75 K
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchRT3,
                    polyval( mResistivityPoly3, mSwitchRT3 ),
                    dpolyval( mResistivityPoly3, mSwitchRT3 ),
                    mSwitchRT4,
                    polyval( mResistivityPoly5, mSwitchRT4 ),
                    dpolyval( mResistivityPoly5, mSwitchRT4 ),
                    mResistivityPoly4 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 250 K <= T < 2000 K
            // - - - - - - - - - - - - - - - - - - - -
            linspace( mSwitchRT6, mTmax, tN, tT );
            for( uint i=0; i<tN; ++i )
            {
                tR( i ) =  nist::rho0( mPrho, mRRR, tT( i ) );
            }
            tR *= tScale ;
            polyfit( tT, tR, 3, mResistivityPoly7 );

            // - - - - - - - - - - - - - - - - - - - -
            // create interval for 200 K <= T < 250 K
            // - - - - - - - - - - - - - - - - - - - -
            create_beam_poly(
                    mSwitchRT5,
                    polyval( mResistivityPoly5, mSwitchRT5 ),
                    dpolyval( mResistivityPoly5, mSwitchRT5 ),
                    mSwitchRT6,
                    polyval( mResistivityPoly7, mSwitchRT6 ),
                    dpolyval( mResistivityPoly7, mSwitchRT6 ),
                    mResistivityPoly6 );

            // - - - - - - - - - - - - - - - - - - - -
            // back scaling
            // - - - - - - - - - - - - - - - - - - - -
            mResistivityPoly0 /= tScale ;
            mResistivityPoly1 /= tScale ;
            mResistivityPoly2 /= tScale ;
            mResistivityPoly3 /= tScale ;
            mResistivityPoly4 /= tScale ;
            mResistivityPoly5 /= tScale ;
            mResistivityPoly6 /= tScale ;
            mResistivityPoly7 /= tScale ;
        }

//----------------------------------------------------------------------------

        void
        Copper::create_resistivity_spline()
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

            tR( 0 ) = 0.0 ;

            for( uint k=1; k<tN; ++k )
            {
                tR( k ) = nist::rho0( mPrho, mRRR, tT( k ) );
            }

            if( mRhoSpline != nullptr )
            {
                delete mRhoSpline ;
            }

            mRhoSpline = new Spline( tT, tR, tA, 0, 0, 0 );
        }

//----------------------------------------------------------------------------

        real
        Copper::rho_el0_poly( const real aT ) const
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
            else if( aT < mSwitchRT5 )
            {
                return polyval( mResistivityPoly5, aT );
            }
            else if( aT < mSwitchRT6 )
            {
                return polyval( mResistivityPoly6, aT );
            }
            else
            {
                return polyval( mResistivityPoly7, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Copper::compute_T_lambda_peak()
        {

            // we do a bisection to find the peak temperature.
            // this precomputed value should provide a first good estimate
            // now we compute that to numerical precision
            real tT0 = 1.0/polyval( mQ, std::log(mRRR ) );
            real tT1 = tT0 * 1.05 ;
            tT0 *= 0.95;
            real aT = BELFEM_QUIET_NAN ;
            real tF0 =  nist::dinvlambda0dT( mPlambda, mBeta, tT0 );
            real tF = nist::dinvlambda0dT( mPlambda, mBeta, tT1 );
            uint tCount = 0 ;

            while( std::abs( tT0 - tT1 ) > 100 * BELFEM_EPS )
            {
                aT = 0.5 * ( tT0 + tT1 );
                tF = nist::dinvlambda0dT( mPlambda, mBeta, aT );
                if( tF0*tF > 0.0 )
                {
                    tT0 = aT ;
                    tF0 = tF ;
                }
                else
                {
                    tT1 = aT ;
                }
                BELFEM_ERROR( tCount++ < 1000,
                              "failed to generate conductivity polynomial for Copper: too many iterations");
            }
           return aT ;
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */