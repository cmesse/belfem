//
// Created by Christian Messe on 23.01.21.
//

#include "cl_Material_Inconel750X.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        Inconel750X::Inconel750X()  :
                IsotropicMaterial( MaterialType::Inconel750X )
        {
            mSpecificHeatPoly = { 9.856161E-07, -1.858865E-03, 1.304636E+00, 1.806522E+02 } ;
            mThermalConductivityPoly = { 1.485812E-06, 1.185818E-02, 8.281965E+00 };
            mHasThermal  = true ;

            mYoungPoly = { -9.285698E+04, 2.333451E+07 , 2.141477E+11 };
            mShearPoly = { -2.369588E+04, -1.359715E+07, 8.751591E+10 };
            mDensityPoly.set_size( 1, 8276.29 ) ;

            mHasMechanical = true ;

            mThermalExpansionPoly =  { -6.2754359E-17, 2.8714439E-12,
                                       2.3017839E-10, 1.1435641E-05,
                                       -3.5056053E-03};

            mHasExpansion = true ;

            mTmax = 1393.0 ;
        }

//----------------------------------------------------------------------------
    }
}