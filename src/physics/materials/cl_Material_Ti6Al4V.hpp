//
// Created by Christian Messe on 02.08.20.
//

#ifndef BELFEM_CL_TI6AL4V_HPP
#define BELFEM_CL_TI6AL4V_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        class Ti6Al4V : public IsotropicMaterial
        {
            const real mSwitchCT0 = 75.0;
            const real mSwitchCT1 = 200.0;
            const real mSwitchCT2 = 298.0;
            const real mSwitchCT3 = 1253.0;
            const real mSwitchCT4 = 1283.0;

            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;
            Vector <real> mSpecificHeatPoly5;

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;
            Vector <real> mThermalConductivityPoly3;
            Vector <real> mThermalConductivityPoly4;
            Vector <real> mThermalConductivityPoly5;
            Vector <real> mThermalConductivityPoly6;

            const real mSwitchLambdaT0 = 50.0;
            const real mSwitchLambdaT1 = 75.0;
            const real mSwitchLambdaT2 = 200.0;
            const real mSwitchLambdaT3 = 298.0;
            const real mSwitchLambdaT4 = 1253.0;
            const real mSwitchLambdaT5 = 1283.0;

            Vector <real> mYoungPoly0;
            Vector <real> mYoungPoly1;
            Vector <real> mYoungPoly2;

            Vector <real> mShearPoly0;
            Vector <real> mShearPoly1;
            Vector <real> mShearPoly2;

            const real mSwitchMechT0 = 1253.0;
            const real mSwitchMechT1 = 1283.0;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Ti6Al4V();

            ~Ti6Al4V() = default;

//----------------------------------------------------------------------------

            real
            E( const real aT ) const;

//----------------------------------------------------------------------------

            real
            G( const real aT ) const;

//----------------------------------------------------------------------------

            real
            c( const real aT ) const;

//----------------------------------------------------------------------------

            real
            lambda( const real aT ) const;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_density_poly();

//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------

            void
            create_mech_polys();

//----------------------------------------------------------------------------

        };
//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_TI6AL4V_HPP
