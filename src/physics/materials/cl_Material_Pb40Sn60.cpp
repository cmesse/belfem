//
// Created by christian on 7/27/22.
//

#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_Material_Pb40Sn60.hpp"
#include "nist_functions.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Pb40Sn60::Pb40Sn60() :
                IsotropicMaterial( MaterialType::Pb40Sn60 )
        {
            mTmax = 456.15;
            mDensityPoly.set_size( 1, 8520.0 );

            this->create_specific_heat_polys();
            this->create_conductivity_polys();
            this->create_resistivity_polys();
            mHasThermal = true;
            mHasResistivity = true ;

        }

//----------------------------------------------------------------------------

        void
        Pb40Sn60::create_specific_heat_polys()
        {
            // nist poly for 4 - 20 K
            mSpecificHeatPoly1 = {
                     1.147156e2,
                    -1.427373e2,
                     6.842287e1,
                    -1.579399e1,
                     1.765134e0,
                    -7.233934e-2,
                     1.104775e-3 } ;

            // normal poly for 30 - 75 K
            mSpecificHeatPoly3 = {
                     2.853872e-4,
                    -6.335783e-2,
                     4.983330,
                    -2.091478e1 };

            // poly for T > 133 K
            mSpecificHeatPoly5 = { 6.118826e-2, 1.153249e2 };

            // create connecting poly
            create_beam_poly(
                    0,
                    0,
                    0,
                    mSwitchCT0,
                    nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                    nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT0 ),
                    mSpecificHeatPoly0 );

            create_beam_poly(
                    mSwitchCT1,
                    nist::cp_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                    nist::dcpdT_janaf( mSpecificHeatPoly1, mSwitchCT1 ),
                    mSwitchCT2,
                    polyval( mSpecificHeatPoly3, mSwitchCT2 ),
                    dpolyval( mSpecificHeatPoly3, mSwitchCT2 ),
                    mSpecificHeatPoly2 );

            create_beam_poly(
                    mSwitchCT3,
                    polyval( mSpecificHeatPoly3, mSwitchCT3 ),
                    dpolyval( mSpecificHeatPoly3, mSwitchCT3 ),
                    mSwitchCT4,
                    polyval( mSpecificHeatPoly5, mSwitchCT4 ),
                    dpolyval( mSpecificHeatPoly5, mSwitchCT4 ),
                    mSpecificHeatPoly4 );
        }

//----------------------------------------------------------------------------

        void
        Pb40Sn60::create_conductivity_polys()
        {
            // 4 < T < 8
            mThermalConductivityPoly1 = {
                      1.505199e-1,
                     -5.233981,
                      6.104678e1,
                     -1.281518e2 };

            // 20 < T < 60
            mThermalConductivityPoly4 = {
                    -3.416565e-4,
                     6.637951e-2,
                    -4.401279,
                     1.601836e2 };

            // 80 < T < 300
            mThermalConductivityPoly6 = {
                    -4.578020e-7,
                     3.743639e-4,
                    -1.022176e-1,
                     6.530614e1 };

            create_beam_poly(
                    0,
                    0,
                    0,
                    mSwitchLT0,
                    polyval( mThermalConductivityPoly1, mSwitchLT0 ),
                    dpolyval( mThermalConductivityPoly1, mSwitchLT0 ),
                    mThermalConductivityPoly0 );

            create_beam_poly(
                    mSwitchLT1,
                    polyval( mThermalConductivityPoly1, mSwitchLT1 ),
                    dpolyval( mThermalConductivityPoly1, mSwitchLT1 ),
                    mSwitchLT2,
                    mLambdaMax,
                    0,
                    mThermalConductivityPoly2 );

            create_beam_poly(
                    mSwitchLT2,
                    mLambdaMax,
                    0,
                    mSwitchLT3,
                    polyval( mThermalConductivityPoly4, mSwitchLT3 ),
                    dpolyval( mThermalConductivityPoly4, mSwitchLT3 ),
                    mThermalConductivityPoly3 );

            create_beam_poly(
                    mSwitchLT4,
                    polyval( mThermalConductivityPoly4, mSwitchLT4 ),
                    dpolyval( mThermalConductivityPoly4, mSwitchLT4 ),
                    mSwitchLT5,
                    polyval( mThermalConductivityPoly6, mSwitchLT5 ),
                    dpolyval( mThermalConductivityPoly6, mSwitchLT5 ),
                    mThermalConductivityPoly5 );

            mThermalConductivityPoly7.set_size( 2 );
            mThermalConductivityPoly7( 0 )
                = dpolyval( mThermalConductivityPoly6, mSwitchLT6 );
            mThermalConductivityPoly7( 1 )
                    = polyval( mThermalConductivityPoly6, mSwitchLT6 )
                            - mThermalConductivityPoly7( 0 ) * mSwitchLT6 ;

        }


//----------------------------------------------------------------------------

        void
        Pb40Sn60::create_resistivity_polys()
        {
           mResistivityPoly0 = { 1.160240e-14, 8.763652e-13, 1.963404E-10, 0.0 };
           mResistivityPoly2 = { 5.640467E-13, 3.344170E-10, - 2.768285E-10 };


           create_beam_poly(
                   mSwitchRT0,
                   polyval( mResistivityPoly0, mSwitchRT0 ),
                   dpolyval( mResistivityPoly0, mSwitchRT0 ),
                   mSwitchRT1,
                   polyval( mResistivityPoly2, mSwitchRT1 ),
                   dpolyval( mResistivityPoly2, mSwitchRT1 ),
                   mResistivityPoly1 );

        }
 //----------------------------------------------------------------------------

        real
        Pb40Sn60::c( const real aT ) const
        {
            if ( aT < mSwitchCT0 )
            {
                return polyval( mSpecificHeatPoly0, aT );
            }
            else if ( aT < mSwitchCT1 )
            {
                return nist::cp_janaf(  mSpecificHeatPoly1, aT );
            }
            else if ( aT < mSwitchCT2 )
            {
                return polyval( mSpecificHeatPoly2, aT );
            }
            else if ( aT < mSwitchCT3 )
            {
                return polyval( mSpecificHeatPoly3, aT );
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
        Pb40Sn60::lambda( const real aT ) const
        {
            if ( aT < mSwitchLT0 )
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
            else if ( aT < mSwitchLT1 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else if ( aT < mSwitchLT2 )
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
            else if ( aT < mSwitchLT3 )
            {
                return polyval( mThermalConductivityPoly3, aT );
            }
            else if ( aT < mSwitchLT4 )
            {
                return polyval( mThermalConductivityPoly4, aT );
            }
            else if ( aT < mSwitchLT5 )
            {
                return polyval( mThermalConductivityPoly5, aT );
            }
            else if ( aT < mSwitchLT6 )
            {
                return polyval( mThermalConductivityPoly6, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly7, aT );
            }
        }

//----------------------------------------------------------------------------

        real
        Pb40Sn60::rho_el ( const real aJ, const real aT, const real aB, const real aAngle ) const
        {
            if( aT < mSwitchRT0 )
            {
                return polyval( mResistivityPoly0, aT ) ;
            }
            else if ( aT < mSwitchRT1 )
            {
                return polyval( mResistivityPoly1, aT ) ; ;
            }
            else
            {
                return polyval( mResistivityPoly2, aT ) ;
            }
        }

//---------------------------------------------------------------------------
    }
}