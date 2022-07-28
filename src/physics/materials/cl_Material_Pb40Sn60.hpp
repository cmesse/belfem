//
// Created by christian on 7/27/22.
//

#ifndef BELFEM_CL_MATERIAL_PB40SN60_HPP
#define BELFEM_CL_MATERIAL_PB40SN60_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
        class Pb40Sn60 : public IsotropicMaterial
        {
            const real mSwitchCT0 = 4.0;
            const real mSwitchCT1 = 20.0;
            const real mSwitchCT2 = 30.0;
            const real mSwitchCT3 = 75.0;
            const real mSwitchCT4 = 133.0;

            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;
            Vector <real> mSpecificHeatPoly5;

            const real mSwitchLT0 = 4.0 ;
            const real mSwitchLT1 = 8.0 ;
            const real mSwitchLT2 = 12.0 ;
            const real mSwitchLT3 = 20.0 ;
            const real mSwitchLT4 = 60.0 ;
            const real mSwitchLT5 = 80.0 ;
            const real mSwitchLT6 = 300.0 ;
            const real mLambdaMax = 110.809108 ;

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;
            Vector <real> mThermalConductivityPoly3;
            Vector <real> mThermalConductivityPoly4;
            Vector <real> mThermalConductivityPoly5;
            Vector <real> mThermalConductivityPoly6;
            Vector <real> mThermalConductivityPoly7;

            const real mSwitchRT0 = 50.0 ;
            const real mSwitchRT1 = 200.0 ;

            Vector <real> mResistivityPoly0;
            Vector <real> mResistivityPoly1;
            Vector <real> mResistivityPoly2;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Pb40Sn60() ;

//----------------------------------------------------------------------------

            ~Pb40Sn60() = default ;

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

            real
            rho_el ( const real aJ, const real aT, const real aB ) const ;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------

            void
            create_resistivity_polys();

//----------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_MATERIAL_PB40SN60_HPP
