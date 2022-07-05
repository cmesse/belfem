//
// Created by Christian Messe on 27.07.20.
//

#ifndef BELFEM_CL_INCONEL718_HPP
#define BELFEM_CL_INCONEL718_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        class Inconel718 : public IsotropicMaterial
        {
            const real mSwitchLambdaT0 = 70.0;
            const real mSwitchLambdaT1 = 300.0;

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;

            const real mSwitchCT0 = 50.0;
            const real mSwitchCT1 = 100.0;
            const real mSwitchCT2 = 250.0;
            const real mSwitchCT3 = 400.0;
            const real mSwitchCT4 = 950.0;
            const real mSwitchCT5 = 1500.0;

            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;
            Vector <real> mSpecificHeatPoly5;
            Vector <real> mSpecificHeatPoly6;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Inconel718();

//----------------------------------------------------------------------------

            ~Inconel718() = default;

//----------------------------------------------------------------------------

            real
            lambda( const real aT ) const;

//----------------------------------------------------------------------------

            real
            c( const real aT ) const;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_mech_polys();

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_INCONEL718_HPP
