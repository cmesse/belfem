//
// Created by Christian Messe on 23.01.21.
//

#ifndef BELFEM_CL_MATERIAL_INCONEL750X_HPP
#define BELFEM_CL_MATERIAL_INCONEL750X_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class Inconel750X : public IsotropicMaterial
        {
            // polynomial for specific heat
            Vector< real > mSpecificHeatPoly;

            // polynomial for thermal conductivity
            Vector< real > mThermalConductivityPoly;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Inconel750X() ;

//----------------------------------------------------------------------------

            ~Inconel750X() = default ;

//----------------------------------------------------------------------------

            real
            lambda( const real aT ) const;

//----------------------------------------------------------------------------

            real
            c( const real aT ) const;

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------


        inline real
        Inconel750X::lambda( const real aT ) const
        {
            return polyval( mThermalConductivityPoly, aT );
        }

//----------------------------------------------------------------------------

        inline real
        Inconel750X::c( const real aT ) const
        {
            return polyval( mSpecificHeatPoly, aT );
        }

//----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_MATERIAL_INCONEL750X_HPP
