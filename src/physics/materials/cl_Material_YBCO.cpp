//
// Created by christian on 7/27/22.
//


#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_Material_YBCO.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        YBCO::YBCO( ) :
                IsotropicMaterial( MaterialType::YBCO )
        {
            mTmax = 1273.15 ;

            this->create_conductivity_polys() ;
            this->create_specific_heat_polys();
            mDensityPoly.set_size( 1, 6300.0 );
            mHasThermal = true ;
        }

//----------------------------------------------------------------------------

        void
        YBCO::create_specific_heat_polys()
        {
            // based on data from 10.1109/TASC.2014.2341180
            mSpecificHeatPoly0 = {  6.905665e-07,
                                   -2.340218e-04,
                                    2.733142e-02,
                                    1.194272e-01,
                                    0.0 };

            mSpecificHeatPoly2 = {  2.895108e-06,
                                   -4.506863e-03,
                                    2.042874,
                                    -4.258366e1 };

            // create connecting poly
            create_beam_poly(
                    mSwitchCT0,
                    polyval( mSpecificHeatPoly0, mSwitchCT0 ),
                    dpolyval( mSpecificHeatPoly0, mSwitchCT0 ),
                    mSwitchCT1,
                    polyval( mSpecificHeatPoly2, mSwitchCT1 ),
                    dpolyval( mSpecificHeatPoly2, mSwitchCT1 ),
                    mSpecificHeatPoly1 );

            // create extrapolation poly
            mSpecificHeatPoly3.set_size( 2 );
            mSpecificHeatPoly3( 0 ) = dpolyval( mSpecificHeatPoly2, mSwitchCT2 );
            mSpecificHeatPoly3( 1 ) = polyval( mSpecificHeatPoly2, mSwitchCT2 )
                    - mSpecificHeatPoly3( 0 ) * mSwitchCT2 ;

        }

//----------------------------------------------------------------------------

        void
        YBCO::create_conductivity_polys()
        {
            // based on data from 10.1109/TASC.2014.2341180

            // from 60 - 100
            mThermalConductivityPoly2 = {  -2.475512e-5, 7.031881e-3,  -6.925542e-1, 2.883919e1 } ;

            // from 150 - 200
            mThermalConductivityPoly4 = {  -9.600000E-04, 4.798333E+00 } ;

            // peak value
            real tKpeak = 8.296353 ;

            create_beam_poly(
                    0,
                    0,
                    0,
                    mSwitchLambdaT0,
                    tKpeak,
                    0,
                    mThermalConductivityPoly0 );

            create_beam_poly(
                    mSwitchLambdaT0,
                    tKpeak,
                    0,
                    mSwitchLambdaT1,
                    polyval( mThermalConductivityPoly2,  mSwitchLambdaT1 ),
                    dpolyval( mThermalConductivityPoly2,  mSwitchLambdaT1 ),
                    mThermalConductivityPoly1 );

            create_beam_poly(
                    mSwitchLambdaT2,
                    polyval( mThermalConductivityPoly2,  mSwitchLambdaT2 ),
                    dpolyval( mThermalConductivityPoly2,  mSwitchLambdaT2 ),
                    mSwitchLambdaT3,
                    polyval( mThermalConductivityPoly4,  mSwitchLambdaT3),
                    dpolyval( mThermalConductivityPoly4,  mSwitchLambdaT3 ),
                    mThermalConductivityPoly3 );


        }

//----------------------------------------------------------------------------

        real
        YBCO::c( const real aT ) const
        {
            if ( aT < mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if ( aT < mSwitchCT1 )
            {
                return polyval( mSpecificHeatPoly1, aT );
            }
            else if ( aT < mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else
            {
                return polyval( mSpecificHeatPoly3, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        YBCO::lambda( const real aT ) const
        {
            if ( aT < mSwitchLambdaT0 )
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
            else if ( aT < mSwitchLambdaT1 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else if ( aT < mSwitchLambdaT2 )
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
            else if ( aT < mSwitchLambdaT3 )
            {
                return polyval( mThermalConductivityPoly3, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly4, aT );
            }
        }

//----------------------------------------------------------------------------
    }
}