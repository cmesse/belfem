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

#ifndef BELFEM_CL_MATERIAL_YBCO_HPP
#define BELFEM_CL_MATERIAL_YBCO_HPP

#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief YBCO (YBa₂Cu₃O₇) high-temperature superconductor implementation
         *
         * YBa₂Cu₃O₇ is a cuprate superconductor with Tc ≈ 92.5 K (typical).
         * This implementation includes both normal-state and superconducting
         * properties.
         *
         * Composition: YBa₂Cu₃O₇ (13 atoms per formula unit)
         * Critical temperature: ~92.5 K (typical, can vary with oxygen content)
         * Density: ~6300 kg/m³
         *
         * Key features:
         * - Mechanical: Based on Lei & Ledbetter 1991 (NISTIR 3980)
         * - Thermal expansion: Salomons et al. 1987 (10.1016/0378-4363(87)90092-1)
         * - Specific heat: Baak 1989, Lang et al. 1988
         * - Electrical: Bloch-Grüneisen model with RRR ≈ 50
         * - Thermal conductivity: Callaway model with e-ph and umklapp scattering
         *   including d-wave superconducting gap effects
         *
         * IMPORTANT: For HTS materials, use rho_powerlaw() for resistivity
         * calculations in the superconducting state, not rho(T).
         *
         * MaterialType: HTS (ρ, λ depend on j, T, n×B and n·B)
         *
         * References:
         * - Sommerfeld et al. 2003 (10.1103/PhysRevB.67.174520) - fitted parameters
         * - Raman data: 10.1103/PhysRevB.80.064505
         */
        class YBCO : public Metal
        {
            // Reference values for scaling normalized curves
            Vector< real > mReferenceValues ;
            Vector< real > mCallawayModel ;

            // Young's modulus (Bezier for low T, polynomial for high T)
            Bezier * mYoungBezier = nullptr ;
            Vector< real > mYoungPoly ;

            // Poisson's ratio (piecewise polynomials)
            Cell< Vector< real > > mPoissionPolys ;

            // Thermal expansion (polynomial with constant tail)
            //! temperature ABOVE which alpha is held at the polynomial's maximum.
            //! Unrelated to Material::alpha_switch_temperature(), which is the split
            //! temperature of the Grueneisen branch and runs the other way.
            real mTAlphaPlateau ;
            Vector< real > mAlphaPoly ;
            Vector< real > mThermalExpansionCryo ;

            // Specific heat (4 piecewise regions in log-log space)
            Vector< real > mTCpSwitch ;
            Cell< Vector< real > > mCpPolys ;

        public:

            /**
             * @brief Constructor - initializes YBCO with typical properties
             *
             * Sets up:
             * - Tc = 92.5 K (typical)
             * - RRR = 50 (typical for thin films)
             * - Debye temperature θD = 275 K (fitted)
             * - Layer thickness = 1 μm (default)
             * - Young's modulus = 210 GPa @ RT (typical for thin films)
             * - Poisson's ratio = 0.25 @ RT (assumed reasonable value)
             */
            YBCO();

            /**
             * @brief Destructor - deletes Bezier curve for Young's modulus
             */
            ~YBCO() override;

            /**
             * @brief Test function for Callaway model parameter fitting
             * @param T Temperature [K]
             * @param aParams Model parameters to test
             * @return Thermal conductivity [W/(m·K)]
             *
             * Used for calibrating the Callaway thermal conductivity model.
             */
            real
            test_callaway( const real T , const Vector< real > & aParams );

        protected:

            /**
             * @brief Young's modulus (custom evaluation)
             * @param T Temperature [K]
             * @return Young's modulus [Pa]
             */
            real
            E_custom( const real T ) const override ;

            /**
             * @brief Poisson's ratio (custom evaluation)
             * @param T Temperature [K]
             * @return Poisson's ratio [-]
             */
            real
            nu_custom( const real T ) const override ;

            /**
             * @brief Thermal expansion coefficient (custom evaluation)
             * @param T Temperature [K]
             * @return Thermal expansion coefficient [1/K]
             */
            real
            alpha_custom(const real T) const override ;

            /**
             * @brief Specific heat capacity (custom evaluation)
             * @param T Temperature [K]
             * @return Specific heat capacity [J/(kg·K)]
             */
            real
            cp_custom(const real T) const override ;

            /**
             * @brief Thermal conductivity using Callaway model
             * @param T Temperature [K]
             * @return Thermal conductivity [W/(m·K)]
             *
             * Implements λ = λ_e + λ_ph where:
             * - λ_e: Electronic contribution (Wiedemann-Franz law)
             * - λ_ph: Phonon contribution (Callaway model with umklapp,
             *         boundary, electron-phonon scattering, and d-wave gap)
             *
             * Parameters fitted to Sommerfeld et al. 2003 data.
             */
            real
            lambda_custom(const real T) const override;

        private:

            /**
             * @brief Set Young's modulus reference value
             * @param aEref Reference Young's modulus [Pa]
             * @param aTref Reference temperature [K]
             */
            void
            set_young( const real aEref, const real aTref );

            /**
             * @brief Set Poisson's ratio reference value
             * @param aNuref Reference Poisson's ratio [-]
             * @param aTref Reference temperature [K]
             */
            void
            set_poisson( const real aNuref, const real aTref );

            /**
             * @brief Create mechanical property curves (E and ν)
             *
             * Based on Lei & Ledbetter 1991 (NISTIR 3980, Fig 4.3)
             */
            void
            create_mech();

            /**
             * @brief Create thermal expansion curve
             *
             * Based on Salomons et al. 1987 (10.1016/0378-4363(87)90092-1)
             */
            void
            create_expansion();

            void
            create_expansion_cryo();

            /**
             * @brief Create specific heat curve
             *
             * Based on Baak 1989, Lang et al. 1988, with 4 piecewise regions
             */
            void
            create_cp();
        };

    }
}

#endif //BELFEM_CL_MATERIAL_YBCO_HPP