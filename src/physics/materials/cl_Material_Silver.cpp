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

            this->create_specific_heat_polys();
            this->create_mech_polys() ;

            this->set_rrr( mRRR );

            mHasMechanical = true ;
            mHasResistivity = true ;

        }

//----------------------------------------------------------------------------

        void
        Silver::set_rrr( const real aRRR )
        {
            this->create_resistivity_polys() ;

        }

//----------------------------------------------------------------------------

        void
        Silver::create_specific_heat_polys()
        {

            // polynomials based in Smith and Fickett, 10.6028/jres.100.012
            mSpecificHeatPoly1 = {  1.685959e+3,
                                   -7.397296e+2,
                                    1.292855e+2,
                                   -1.126873e+1,
                                    4.915730e-1,
                                   -7.351090e-3,
                                    3.940887e-5 };

            mSpecificHeatPoly3 = {   2.173720e+5,
                                    -1.728826e+4,
                                     4.101784e+2,
                                    -1.035062,
                                     3.722819e-3,
                                    -6.981773e-6,
                                     5.452852e-9 };

            // very low temperatue condition
            create_beam_poly( 0.0,
                              0.0,
                              0.0,
                              mSwitchCT0,
                              nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                              mSpecificHeatPoly0 );

            // glue condition
            create_beam_poly( mSwitchCT1,
                              nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                              mSwitchCT2,
                              nist::cp_janaf( mSpecificHeatPoly3, mSwitchCT2 ),
                              nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT2 ),
                              mSpecificHeatPoly2 );

            // linear function for high temperatures
            mSpecificHeatPoly4.set_size( 2 );
            mSpecificHeatPoly4( 0 )
                = nist::dcpdT_janaf( mSpecificHeatPoly3, mSwitchCT3 ) ;
            mSpecificHeatPoly4( 1 ) = nist::cp_janaf( mSpecificHeatPoly3, mSwitchCT3 )
                                            -mSpecificHeatPoly4( 0 ) * mSwitchCT3 ;

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
                              polyval( mResistivityPoly1, mSwitchRT4 ),
                              mSwitchRT2,
                              polyval( mResistivityPoly3, mSwitchRT2 ),
                              polyval( mResistivityPoly3, mSwitchRT2 ),
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
                              polyval( mResistivityPoly3, mSwitchRT3 ),
                              mSwitchRT4,
                              polyval( mResistivityPoly5, mSwitchRT4 ),
                              polyval( mResistivityPoly5, mSwitchRT4 ),
                              mResistivityPoly4 );

            // backscale
            mResistivityPoly0 /= tScale ;
            mResistivityPoly1 /= tScale ;
            mResistivityPoly2 /= tScale ;
            mResistivityPoly3 /= tScale ;
            mResistivityPoly4 /= tScale ;
            mResistivityPoly5 /= tScale ;
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
        Silver::c( const real aT ) const
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
            else
            {
                return polyval( mSpecificHeatPoly4, aT );
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

        }

//----------------------------------------------------------------------------

        // electric resistance
        real
        Silver::rho_el ( const real aJ, const real aT, const real aB ) const
        {
            // resistivity for zero magnetic field
            real tRho0 = this->rho_el0_poly( aT );

            return tRho0 * ( 1.0 + nist::kohler( mKohlerA, mKohlerB, mKohlerXmin, mRhoRef / tRho0 * std::abs( aB ) ) );
        }


    }
}