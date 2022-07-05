//
// Created by Christian Messe on 02.12.19.
//

#ifndef BELFEM_CL_ATMOSPHEREMODEL_ISA1976_HPP
#define BELFEM_CL_ATMOSPHEREMODEL_ISA1976_HPP

#include "typedefs.hpp"
#include "cl_AtmosphereModel.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace atmosphere
    {
        /**
         * Up tp 100 km Standard ISA 1976 atmosphere
         * higher atmosphere, see Sforza, Manned Spacecraft Design Principles,
         * Chapter 2
         */
        class AtmosphereModel_ISA1976 : public AtmosphereModel
        {
            // reference temperature in K
            const real mT0 = 288.15;

            // reference pressure in Pa
            const real mP0 = 1.01325e5;

            // Gas Specific constant for Air in J/(kg*K)
            const real mR = 287.05019205;

            // reference geopotential altitudes
            const Vector < real > mH;

            // reference lapse rate in K/m
            const Vector < real > mL;

            const real mMaxAltitude;

            // temperature at 100 km
            real mT100;

            // pressure at 100 km
            real mP100;

            // constant for hight altitude temperature equation
            real mTheta100;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

             AtmosphereModel_ISA1976();

//------------------------------------------------------------------------------

            ~AtmosphereModel_ISA1976() = default;

//------------------------------------------------------------------------------

            real
            max_altitude() const;

//------------------------------------------------------------------------------

            void
            compute_T_and_p(
                    const real & aAltitude,
                          real & aTemperature,
                          real & aPressure ) const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_ATMOSPHEREMODEL_ISA1976_HPP
