//
// Created by Christian Messe on 14.12.20.
//

#ifndef BELFEM_CL_MATERIAL_ZIRCONIA_HPP
#define BELFEM_CL_MATERIAL_ZIRCONIA_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class Zirconia : public IsotropicMaterial
        {
            const real mSwitchCT0 = 100.0;

            Vector< real > mSpecificHeatPoly0;
            Vector< real > mSpecificHeatPoly1;

            real mPoisson ;

//----------------------------------------------------------------------------

        public:
//----------------------------------------------------------------------------

            Zirconia();

//----------------------------------------------------------------------------

            ~Zirconia() = default ;

//----------------------------------------------------------------------------

            /**
             * specific heat capacity in J/(kg*K)
             */
            real
            c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * Poisson Number
             */
            real
            nu( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * Shear Modulus in Pa
             */
            real
            G( const real aT=BELFEM_TREF ) const;

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
            create_mech_polys();

//----------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_MATERIAL_ZIRCONIA_HPP
