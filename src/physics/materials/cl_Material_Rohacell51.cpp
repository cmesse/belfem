//
// Created by Christian Messe on 17.07.20.
//

#include "cl_Material_AltraMat80.hpp"
#include "fn_polyval.hpp"
#include "fn_polyfit.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"

#include "cl_Material_Rohacell51.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Rohacell51::Rohacell51()  :
                IsotropicMaterial( MaterialType::Rohacell51 )
        {
            mTmax = 480.0;

            // (Verstraete, et al., 2010) ?
            Vector <real> tT0 = { 30.415, 39.4942, 50.3891, 64.9157, 77.1725,
                                  88.0674, 102.594, 119.39, 132.101, 142.542,
                                  157.069, 167.51, 180.22, 187.938, 200.649,
                                  207.458, 216.083, 223.8, 228.34, 235.603,
                                  243.32, 251.492, 260.571, 267.834, 276.913,
                                  284.63, 291.44, 297.795, 302.335, 308.236,
                                  313.684, 319.131, 324.125 };

            Vector <real> tL0 = { 0.00493421, 0.00578947, 0.00697368, 0.00861842,
                                  0.00973684, 0.0106579, 0.0119079, 0.0134868,
                                  0.0147368, 0.0153947, 0.0167105, 0.0175, 0.0184868,
                                  0.0192105, 0.0205263, 0.0211184, 0.0221053, 0.0227632,
                                  0.0234211, 0.0243421, 0.0248026, 0.0258553, 0.0271711,
                                  0.0282237, 0.0296053, 0.0305263, 0.0318421, 0.0325658,
                                  0.0330921, 0.0338816, 0.0348684, 0.0357895, 0.0365132 };

            // table as given in doi.org/10.2514/6.2012-3010
            Vector <real> tT2 = { 95.83, 134.94, 154.72, 174.44,
                                  228.89, 288.33, 367.22, 479.44 };

            Vector <real> tL2 = { 0.0127, 0.0171, 0.019, 0.0208,
                                  0.0272, 0.036, 0.0567, 0.0848 };

            // lower polynomial
            polyfit( tT0, tL0, 3, mThermalConductivityPoly0 );

            // upper polynomial
            polyfit( tT2, tL2, 2, mThermalConductivityPoly2 );

            // connecting polynomial
            create_beam_poly(
                    mSwitchLambdaT0, // T = 200 K
                    polyval( mThermalConductivityPoly0, mSwitchLambdaT0 ), // lambda @T=200 K
                    dpolyval( mThermalConductivityPoly0, mSwitchLambdaT0 ), // dlambdadT @T=200 K,
                    mSwitchLambdaT1, // T = 400 K
                    polyval( mThermalConductivityPoly2, mSwitchLambdaT1 ), // lambda @T=400 K
                    dpolyval( mThermalConductivityPoly2, mSwitchLambdaT1 ), // dlambdadT @T=400 K
                    mThermalConductivityPoly1 );


            // from datasheet

            mDensityPoly.set_size( 1, 52.1 );

            mYoungPoly.set_size( 1, 70.0e6 );
            mShearPoly.set_size( 1, 19.0e6 );

            // For temperatures < 300 K: data for Polymethacrylimide
            // tanken from doi.org/10.1063/1.555671
            // for temperatures > 300 K: data for 200 WF version, taken from
            // JACOB LANGER, Investigation and Simulation of the Heating Effects
            // in Sandwich Cores During Vibrational Loading

            // assumption:
            mSpecificHeatPoly = { 2.157253E-06, -5.645322E-03, 6.184665, 0.0 };

            mHasThermal = true;
            mHasMechanical = true;
        }

//----------------------------------------------------------------------------

        real
        Rohacell51::lambda( const real aT ) const
        {
            if ( aT < mSwitchLambdaT0 )
            {
                return polyval( mThermalConductivityPoly0, aT );
            }
            else if ( aT < mSwitchLambdaT1 )
            {
                return polyval( mThermalConductivityPoly1, aT );
            }
            else
            {
                return polyval( mThermalConductivityPoly2, aT );
            }
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */