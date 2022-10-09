//
// Created by Christian Messe on 18.04.22.
//

#ifndef BELFEM_CL_MATERIAL_HastelloyC276C276_HPP
#define BELFEM_CL_MATERIAL_HastelloyC276C276_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        class HastelloyC276 : public IsotropicMaterial
        {

            const real mSwitchCT1 = 20.0;
            const real mSwitchCT2 = 40.0;
            const real mSwitchCT3 = 90.0;
            const real mSwitchCT4 = 110.0;
            const real mSwitchCT5 = 298.15;


            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;
            Vector <real> mSpecificHeatPoly5;
            Vector <real> mSpecificHeatPoly6;

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;

            const real mSwitchLambdaT0 = 50.0;
            const real mSwitchLambdaT1 = 100.0;

            const real mSwitchRhoT0 = 25.0 ;
            const real mSwitchRhoT1 = 60.0 ;

            Vector <real> mElectricResistivityPoly0;
            Vector <real> mElectricResistivityPoly1;
            Vector <real> mElectricResistivityPoly2;

            Vector< real > mIntAlphaPoly ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            HastelloyC276();

//----------------------------------------------------------------------------

            ~HastelloyC276() = default;

//----------------------------------------------------------------------------

            /**
             * specific heat capacity in J/(kg*K)
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

            /**
             * Young's Modulus in Pa
             */
            real
            E( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * Shear Modulus in Pa
             */
            real
            G( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            // electric resistance
            real
            rho_el ( const real aJ=0, const real aT=0, const real aB=0, const real aAngle=0 ) const ;

            /**
             * thermal strain parameter
             *
             * mu = exp( int( a, T, aTref, aT ) ) - 1
             */
            real
            mu( const real aT=BELFEM_TREF, const real aTref=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_resistivity_polys();

//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------

            // void
            // create_density_poly();

//----------------------------------------------------------------------------

            // void
            // create_mech_polys();

//----------------------------------------------------------------------------

            void
            create_expansion_poly();

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------

        inline real
        HastelloyC276::E( const real aT ) const
        {
            return 205e9 ;
        }

//----------------------------------------------------------------------------

        inline real
        HastelloyC276::G( const real aT ) const
        {
            return 79e9 ;
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */

#endif //BELFEM_CL_MATERIAL_HastelloyC276C276_HPP
