//
// Created by Christian Messe on 14.12.20.
//

#include "cl_Material_Zirconia.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Zirconia::Zirconia() :
                IsotropicMaterial( MaterialType::Zirconia )
        {
            // set maximum temperature
            mTmax = 3000.0;


            this->create_specific_heat_polys();
            this->create_conductivity_polys();
            this->create_mech_polys() ;
            this->create_density_poly( 5890 );

            mHasMechanical = true;
            mHasThermal = true;
            mHasExpansion = true;
        }

//----------------------------------------------------------------------------

        real
        Zirconia::c( const real aT ) const
        {
            if ( aT > mTmax )
            {
                return polyval( mSpecificHeatPoly1, mTmax );
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

        /**
         * Poisson Number
         */
        real
        Zirconia::nu( const real aT ) const
        {
            return mPoisson;
        }

//----------------------------------------------------------------------------

        /**
         * Shear Modulus in Pa
         */
        real
        Zirconia::G( const real aT ) const
        {
            return this->E( aT ) / ( 2.0 * mPoisson + 2.0 );
        }

//----------------------------------------------------------------------------

        void
        Zirconia::create_specific_heat_polys()
        {
            // dataset based on JANAF tables and
            // Touloukian, Y. S. (Hrsg.): 1st. Edition. Bd. Volume 1 - Elements: Thermophysical
            //properties of high temperature solid materials. New York : Purdue
            //University,Thermophysical Properties Research Center, 1967

            mSpecificHeatPoly1 = { 7.36354e-14, -6.26045e-10, 2.08140e-6, -3.14632e-3, 2.23168e0 };

            // create low temperature poly
            create_beam_poly( 0.0,
                              0.0,
                              0.0,
                              mSwitchCT0,
                              polyval( mSpecificHeatPoly1, mSwitchCT0 ),
                              dpolyval( mSpecificHeatPoly1, mSwitchCT0 ),
                              mSpecificHeatPoly0 );
        }

//----------------------------------------------------------------------------

        void
        Zirconia::create_conductivity_polys()
        {
            // Based on data from Touloukian 1968, Blanke 1989 and Schlichting 2001

            mThermalConductivityPoly = { 5.13591e-10,
                                         -1.505743e-6,
                                         1.67486e-3,
                                         1.25613 };

        }

//----------------------------------------------------------------------------

        void
        Zirconia::create_mech_polys()
        {
            mYoungPoly = {
                    -3.09616e1,
                    1.18243e5,
                    -1.63052e8,
                    2.61414e11
                    };

            mPoisson = 0.296391 ;

            mThermalExpansionPoly.set_size( 2, 0.0 );
            mThermalExpansionPoly( 0 ) = 1.02762e-5 ;
        }

//----------------------------------------------------------------------------
    }
}