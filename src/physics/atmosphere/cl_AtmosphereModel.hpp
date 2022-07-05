//
// Created by Christian Messe on 02.12.19.
//

#ifndef BELFEM_CL_ATMOSPHEREMODEL_HPP
#define BELFEM_CL_ATMOSPHEREMODEL_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace atmosphere
    {
//------------------------------------------------------------------------------

        class AtmosphereModel
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            AtmosphereModel()        = default;

//------------------------------------------------------------------------------

            virtual ~AtmosphereModel() = default;

//------------------------------------------------------------------------------

            virtual real
            max_altitude() const = 0;

//------------------------------------------------------------------------------

            virtual void
            compute_T_and_p(
                    const real & aAltitude,
                          real & aTemperature,
                          real & aPressure ) const = 0;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_ATMOSPHEREMODEL_HPP
