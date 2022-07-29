//
// Created by Christian Messe on 27.07.20.
//

#ifndef BELFEM_CL_ALUMINUM_HPP
#define BELFEM_CL_ALUMINUM_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class Aluminum : public IsotropicMaterial
        {

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;

            const real mSwitchLambdaT0 = 140.0;
            const real mSwitchLambdaT1 = 300.0;

            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;

            const real mSwitchCT0 = 50.0;
            const real mSwitchCT1 = 75.0;
            const real mSwitchCT2 = 200.0;
            const real mSwitchCT3 = 298.0;

            // nist data for electric resistivity by means of Kohler's rule
            const Vector < real > mA = { 0.01827, -0.1839, 0.6229, 0.3168, -2.662 } ;

            const Vector< real > mPrho = { 1.467e-8, 1.18777e-15, 3.56474, 3.18825e10, 1.09259, 49.7615, 3.01754, -0.557109	};

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Aluminum();

//----------------------------------------------------------------------------

            ~Aluminum() = default;

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
//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_ALUMINUM_HPP
