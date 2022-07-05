//
// Created by Christian Messe on 13.12.20.
//

#ifndef BELFEM_CL_MATERIAL_CUCRZR_HPP
#define BELFEM_CL_MATERIAL_CUCRZR_HPP

#include "cl_Material_Copper.hpp"

#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_linspace.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class CuCrZr: public IsotropicMaterial
        {

            const real mSwitchCT0 = 100.0;
            const real mSwitchCT1 = 400.0;
            const real mSwitchCT2 = 700.0;

            Vector< real > mSpecificHeatPoly0;
            Vector< real > mSpecificHeatPoly1;
            Vector< real > mSpecificHeatPoly2;
            Vector< real > mSpecificHeatPoly3;

            Vector< real > mThermalConductivityPoly0;
            Vector< real > mThermalConductivityPoly1;
            Vector< real > mThermalConductivityPoly2;
            Vector< real > mThermalConductivityPoly3;
            Vector< real > mThermalConductivityPoly4;
            Vector< real > mThermalConductivityPoly5;

            const real mSwitchLambdaT0 = 4.0;
            const real mSwitchLambdaT1 = 40.0;
            const real mSwitchLambdaT2 = 70.0;
            const real mSwitchLambdaT3 = 150.0;
            const real mSwitchLambdaT4 = 300.0;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            CuCrZr();

//----------------------------------------------------------------------------

            ~CuCrZr() = default;

//----------------------------------------------------------------------------

            /**
             * specific heat capacity in J/(kg*K)
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
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */


#endif //BELFEM_CL_MATERIAL_CUCRZR_HPP
