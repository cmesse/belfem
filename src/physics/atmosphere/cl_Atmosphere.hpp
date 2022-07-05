//
// Created by Christian Messe on 02.12.19.
//

#ifndef BELFEM_CL_ATMOSPHERE_HPP
#define BELFEM_CL_ATMOSPHERE_HPP

#include "typedefs.hpp"
#include "en_AtmosphereType.hpp"
#include "cl_AtmosphereModel.hpp"
namespace belfem
{
    class Atmosphere
    {
        const AtmosphereType          mType;
        atmosphere::AtmosphereModel * mModel;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Atmosphere( const AtmosphereType aType );

//------------------------------------------------------------------------------

        ~Atmosphere();

//------------------------------------------------------------------------------

        void
        compute_T_and_P( const real & aAltitude,
                               real & aTemperature,
                               real & aPressure ) const;

//------------------------------------------------------------------------------

        // return the type
        inline AtmosphereType
        type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------
    };
}

#endif //BELFEM_CL_ATMOSPHERE_HPP
