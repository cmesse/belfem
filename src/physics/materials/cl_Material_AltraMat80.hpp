//
// Created by Christian Messe on 02.07.20.
//

#ifndef BELFEM_CL_ALTRAMAT80_HPP
#define BELFEM_CL_ALTRAMAT80_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        class AltraMat80 : public IsotropicMaterial
        {
            // polynomial data extrapolated from Janaf
            // 20% mass Silicon Oxide, Cristobalite, High (SiO2)
            // 80% mass Aluminum Oxide, Gamma (Al2O3)
            // polynomial from 0 K to 200 K ;
            Vector <real> mSpecificHeatPoly0;

            // polynomial from 200 K to 800 K ;
            Vector <real> mSpecificHeatPoly1;

            // polynomial from 800 K to 1200 K ;
            Vector <real> mSpecificHeatPoly2;

            // polynomial from 1200 K to 3000 K ;
            Vector <real> mSpecificHeatPoly3;

            // switch temperatures
            const real mSwitchCpT0 = 200;
            const real mSwitchCpT1 = 800;
            const real mSwitchCpT2 = 1200;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AltraMat80();

//----------------------------------------------------------------------------

            ~AltraMat80() = default;

//----------------------------------------------------------------------------

            real
            c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        };
//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_ALTRAMAT80_HPP
