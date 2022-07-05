//
// Created by Christian Messe on 17.07.20.
//

#ifndef BELFEM_CL_ROHACELL51_HPP
#define BELFEM_CL_ROHACELL51_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        class Rohacell51 : public IsotropicMaterial
        {

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;

            const real mSwitchLambdaT0 = 200;
            const real mSwitchLambdaT1 = 400;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Rohacell51();

//----------------------------------------------------------------------------

            ~Rohacell51() = default;

//----------------------------------------------------------------------------

            /**
             * thermal conductivity in W/(m*K)
             */
            real
            lambda( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_ROHACELL51_HPP
