//
// Created by christian on 7/28/22.
//

#ifndef BELFEM_CL_MATERIAL_SS314_HPP
#define BELFEM_CL_MATERIAL_SS314_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
        class SAE314 : public IsotropicMaterial
        {
            const real mSwitchCT0 = 10.0;
            const real mSwitchCT1 = 40.0;
            const real mSwitchCT2 = 60.0;
            const real mSwitchCT3 = 221.34;

            Vector< real > mSpecificHeatPoly0;
            Vector< real > mSpecificHeatPoly1;
            Vector< real > mSpecificHeatPoly2;
            Vector< real > mSpecificHeatPoly3;
            Vector< real > mSpecificHeatPoly4;

            const real mSwitchLT0 = 20.0;
            const real mSwitchLT1 = 40.0;
            const real mSwitchLT2 = 80.0;
            const real mSwitchLT3 = 100.0;
            const real mSwitchLT4 = 150.0;
            const real mSwitchLT5 = 250.0;

            Vector< real > mThermalConductivityPoly0;
            Vector< real > mThermalConductivityPoly1;
            Vector< real > mThermalConductivityPoly2;
            Vector< real > mThermalConductivityPoly3;
            Vector< real > mThermalConductivityPoly4;
            Vector< real > mThermalConductivityPoly5;
            Vector< real > mThermalConductivityPoly6;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            SAE314();

//----------------------------------------------------------------------------

            ~SAE314() = default;

//----------------------------------------------------------------------------

            /**
             * specific heat capacity in J /( kg*K)
             */
            real
            c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * thermal conductivity in W/(m*K)
             */
            real
            lambda( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_expansion_poly();

//----------------------------------------------------------------------------

            void
            create_mech_polys();

//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MATERIAL_SS314_HPP
