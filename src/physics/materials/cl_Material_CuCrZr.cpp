//
// Created by Christian Messe on 13.12.20.
//

#include "cl_Material_CuCrZr.hpp"

#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        CuCrZr::CuCrZr() :
                IsotropicMaterial( MaterialType::CuCrZr )
        {
            // set maximum temperature
            mTmax = 1357.15;

            this->create_specific_heat_polys();
            this->create_conductivity_polys();

            mDensityPoly.set_size( 1, 8910.0 );

            mHasThermal = true;
            mHasMechanical = false;
            mHasExpansion = false;


        }

//----------------------------------------------------------------------------

        void
        CuCrZr::create_specific_heat_polys()
        {
            // from datasheet
            Vector< real > tT1  = { 293.15, 373.15, 473.15, 573.15 };
            Vector< real > tCp1 = { 370.,450.,480.,500. };

            // scaled datapoints, based on NIST data
            // seriously, this is an educated guess
            Vector< real > tT3  = { 600., 700., 800., 900., 1000., 1100., 1200., 1300. };
            Vector< real > tCp3 = { 505.988249, 516.1399393, 525.9220264, 536.9163789,
                                    549.209915, 565.4277502, 585.7943597, 616.3137886};

            polyfit( tT1, tCp1, 2, mSpecificHeatPoly1 );
            polyfit( tT3, tCp3, 3, mSpecificHeatPoly3 );

            create_beam_poly( 0.0,
                              0.0,
                              0.0,
                              mSwitchCT0,
                              polyval( mSpecificHeatPoly1, mSwitchCT0 ),
                              dpolyval( mSpecificHeatPoly1, mSwitchCT0 ),
                              mSpecificHeatPoly0 ) ;

            create_beam_poly( mSwitchCT1,
                              polyval( mSpecificHeatPoly1, mSwitchCT1 ),
                              dpolyval( mSpecificHeatPoly1, mSwitchCT1 ),
                              mSwitchCT2,
                              polyval( mSpecificHeatPoly3, mSwitchCT2 ),
                              dpolyval( mSpecificHeatPoly3, mSwitchCT2 ),
                              mSpecificHeatPoly2 ) ;
        }

//----------------------------------------------------------------------------

        void
        CuCrZr::create_conductivity_polys()
        {
            // from DOI 10.1007/s12289-008-0-0201-2
            // values for CuCrZ-X
            Vector< real > tT1 = { 4.710872684, 5.493947427, 7.49485782,
                                   9.496372497, 11.48299017, 15.01758601,
                                   18.99948826, 23.00878032, 26.95516112,
                                   33.49074226, 42.49867926, 52.48564091 };
            Vector< real > tL1 = { 25.66203344, 31.1157454, 41.04163214,
                                   52.48090385, 61.55508358, 79.94170575,
                                   99.10247202, 120.9652524, 137.5477206,
                                   159.5530278, 190.4852666, 219.9829798 };


            Vector< real > tT3 = { 67.41066014, 87.49830253, 109.8665383,
                                   130.0758725, 150.1062525, 170.1147154,
                                   189.9059987, 229.6338645, 269.8329798,
                                   309.9787065 };

            Vector< real > tL3 = { 239.8309732, 252.9243739, 262.6306573,
                                   272.1040387, 276.9687905, 286.326543,
                                   291.4455687, 303.9723136, 314.9384765,
                                   319.8589631};
            // from Datasheet Deutsches Kupferinstitut,
            // except last point, which is from 0.1007/s10853-010-4542-0
            //Vector< real > tT5 = { 373.15, 473.15, 573.15, 673.15, 1072.627 };

            //Vector< real > tL5 = { 315., 324., 333., 336., 350.462 };
            // manually drawn polynomial

            mThermalConductivityPoly5 = {  9.126151E-08,
                                          -2.354669E-04,
                                           2.366939E-01,
                                           2.590615E+02 };


            polyfit( tT1, tL1, 2, mThermalConductivityPoly1 );
            polyfit( tT3, tL3, 4, mThermalConductivityPoly3 );
            //polyfit( tT5, tL5, 2, mThermalConductivityPoly5 );

            create_beam_poly( 0.0,
                              0.0,
                              0.0,
                              mSwitchLambdaT0,
                              polyval( mSpecificHeatPoly1, mSwitchLambdaT0 ),
                              dpolyval( mSpecificHeatPoly1, mSwitchLambdaT0 ),
                              mThermalConductivityPoly0 ) ;

            create_beam_poly( mSwitchLambdaT1,
                              polyval( mThermalConductivityPoly1, mSwitchLambdaT1 ),
                              dpolyval( mThermalConductivityPoly1, mSwitchLambdaT1 ),
                              mSwitchLambdaT2,
                              polyval( mThermalConductivityPoly3, mSwitchLambdaT2 ),
                              dpolyval( mThermalConductivityPoly3, mSwitchLambdaT2 ),
                              mThermalConductivityPoly2 ) ;

            create_beam_poly( mSwitchLambdaT3,
                              polyval( mThermalConductivityPoly3, mSwitchLambdaT3 ),
                              dpolyval( mThermalConductivityPoly3, mSwitchLambdaT3 ),
                              mSwitchLambdaT4,
                              polyval( mThermalConductivityPoly5, mSwitchLambdaT4 ),
                              dpolyval( mThermalConductivityPoly5, mSwitchLambdaT4 ),
                              mThermalConductivityPoly4 ) ;

        }

//----------------------------------------------------------------------------

        real
        CuCrZr::c( const real aT ) const
        {
            if( aT > mTmax )
            {
                return polyval( mSpecificHeatPoly3, mTmax ) ;
            }
            else if ( aT > mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly3, aT );
            }
            else if ( aT > mSwitchCT1 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if ( aT > mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly1, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        CuCrZr::lambda( const real aT ) const
        {
            if( aT > mTmax )
            {
                return polyval( mThermalConductivityPoly5, mTmax );
            }
            else if ( aT > mSwitchLambdaT4 )
            {
                return polyval( mThermalConductivityPoly5, aT );
            }
            else if ( aT > mSwitchLambdaT3 )
            {
                return polyval( mThermalConductivityPoly4, aT );
            }
            else if ( aT > mSwitchLambdaT2 )
            {
                return polyval( mThermalConductivityPoly3, aT );
            }
            else if ( aT > mSwitchLambdaT1 )
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
            else if ( aT > mSwitchLambdaT0 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
        }

//----------------------------------------------------------------------------

    } /* end namespace material */
}  /* end namespace belfem */