//
// Created by christian on 7/27/22.
//

#ifndef BELFEM_CL_MATERIAL_YBCO_HPP
#define BELFEM_CL_MATERIAL_YBCO_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
        class YBCO : public IsotropicMaterial
        {
            const real mSwitchCT0 = 80.0;
            const real mSwitchCT1 = 100.0;
            const real mSwitchCT2 = 300.0;

            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;

            const real mSwitchLambdaT0 = 40.4417 ;
            const real mSwitchLambdaT1 = 60.0 ;
            const real mSwitchLambdaT2 = 100.0 ;
            const real mSwitchLambdaT3 = 150.0 ;

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;
            Vector <real> mThermalConductivityPoly3;
            Vector <real> mThermalConductivityPoly4;

 //----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            YBCO() ;

//----------------------------------------------------------------------------

            ~YBCO() = default ;

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
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MATERIAL_YBCO_HPP
