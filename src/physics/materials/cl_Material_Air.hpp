//
// Created by Christian Messe on 03.08.20.
//

#ifndef BELFEM_CL_MATERIAL_AIR_HPP
#define BELFEM_CL_MATERIAL_AIR_HPP
#include "constants.hpp"
#include "cl_Spline.hpp"
#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        /**
         * a fake material that wraps an air object
         */
        class Air  : public IsotropicMaterial
        {


            // spline for specific heat capacity
            Spline mHeatSpline ;

            // spline for thermal conductivity
            Spline mConductivitySpline ;

            // gas constant
            real mR ;

            // reference pressure
            real mP = BELFEM_PREF ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Air() ;

            ~Air() = default ;

//----------------------------------------------------------------------------
            /**
             * sets the pressure for the air material.
             * this funcition is only available for air, and must be accessed
             * via reinterpret_cast< >
             *
             * @param aP pressure in Pa for density computation
             */
             void
             set_reference_pressure( const real & aP );

//----------------------------------------------------------------------------

            inline real
            rho( const real aT = BELFEM_TREF ) const
            {
                 return mP / ( mR * aT );
            }

//----------------------------------------------------------------------------

            inline real
            c( const real aT = BELFEM_TREF ) const
            {
                return mHeatSpline.deval( aT );
            }

//----------------------------------------------------------------------------

            inline real
            lambda( const real aT = BELFEM_TREF ) const
            {
                return mConductivitySpline.eval( aT );
            }

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_splines();

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_MATERIAL_AIR_HPP
