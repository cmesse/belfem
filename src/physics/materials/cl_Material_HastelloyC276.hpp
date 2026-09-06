/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_MATERIAL_HASTELLOYC276_HPP
#define BELFEM_CL_MATERIAL_HASTELLOYC276_HPP

#include "cl_Material_Metal.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Hastelloy C-276 nickel-based superalloy implementation
         *
         * Hastelloy C-276 is a nickel-molybdenum-chromium alloy with excellent
         * corrosion resistance and good low-temperature properties.
         *
         * Composition (from 10.1016/j.cryogenics.2023.103776):
         * - Ni: balance (~57%)
         * - Mo: 14.18%
         * - Cr: 14.01%
         * - Fe: 4.76%
         * - W: 2.93%
         * - Co: 3.71%
         * - Mn: 0.53%
         * - Al: 0.26%
         *
         * Data sources:
         * - Mechanical: 10.1016/j.cryogenics.2006.01.014, MATWEB
         * - Electrical: 10.1063/1.2899058
         * - Thermal: 10.1063/1.2899058
         * - Magnetic susceptibility: χ(T) fit prepared in create_mu() but NOT
         *   routed — the material currently reports mu = mu0 ( see create_mu() )
         *
         * MaterialType: LookupAlloy (properties depend only on temperature, not field)
         */
        class HastelloyC276 : public Metal
        {
            // Mechanical properties
            Vector< real > mYoungPoly ;      // Young's modulus polynomial
            Vector< real > mPoissonPoly ;    // Poisson's ratio polynomial

            // Electrical resistivity (piecewise: polynomial and log-log)
            real mTRhoSwitch ;
            Cell< Vector< real > > mRhoPolys ;

            // Thermal expansion (polynomial with constant tail). NOT the split
            // temperature of the Grueneisen branch - this is the temperature ABOVE
            // which alpha is held at the polynomial's maximum. See
            // Material::alpha_switch_temperature() for the unrelated base class one.
            real mTAlphaPlateau ;
            Cell< Vector< real > > mAlphaPolys ;
            Vector< real > mThermalExpansionCryo ;

            // Magnetic susceptibility χ (piecewise in 1/T)
            Vector< real > mTMuSwitch ;
            Cell< Vector< real > > mMuPolys ;

            // Specific heat (complex piecewise with Debye tail)
            Vector< real > mTCpSwitch ;
            Cell< Vector< real > > mCpPolys ;

            // Thermal conductivity (phonon + electronic contributions)
            real mTLambdaSwitch ;
            Cell< Vector< real > > mLambdaPolys ;

        public:

            HastelloyC276();

            ~HastelloyC276() override;

        protected:

            real
            E_custom( const real T ) const override ;

            real
            nu_custom( const real T ) const override ;

            real
            rho_custom( const real T ) const override ;

            real
            alpha_custom(const real T) const override;

            real
            mu_custom( const real H, const real T ) const override;

            real
            cp_custom( const real T ) const override;

            real
            lambda_custom(const real T) const override;

        private:

            void
            set_constants();

            void
            create_mech();

            void
            create_rho();

            void
            create_alpha();

            void
            create_alpha_cryo();

            void
            create_mu();

            void
            create_cp();

            void
            create_lambda();
        };

        inline real
        HastelloyC276::E_custom( const real T ) const
        {
            return polyval( mYoungPoly, T );
        }

        inline real
        HastelloyC276::nu_custom( const real T ) const
        {
            return polyval( mPoissonPoly, T );
        }

        inline real
        HastelloyC276::alpha_custom( const real T ) const
        {
            if ( T < this->alpha_switch_temperature() )
            {
                // Grueneisen branch: alpha = C( T ) cp( T )
                return std::exp( polyval( mThermalExpansionCryo, T ) ) * this->cp( T );
            }
            else if ( T < mTAlphaPlateau )
            {
                return polyval( mAlphaPolys( 0 ), T );
            }
            else
            {
                return mAlphaPolys( 1 )( 0 );
            }
        }

    }
}

#endif //BELFEM_CL_MATERIAL_HASTELLOYC276_HPP