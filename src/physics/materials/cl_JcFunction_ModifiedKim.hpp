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

#ifndef BELFEM_CL_JCFUNCTIONMODIFIEDKIM_HPP
#define BELFEM_CL_JCFUNCTIONMODIFIEDKIM_HPP

#include "cl_JcFunction.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Modified Kim model for critical current density
         *
         * Implements an analytical model for Jc(B, angle) based on the
         * modified Kim-Anderson relation with angular dependence.
         *
         * Model equation:
         * @code
         * Jc(B,θ) = Jc0 / [1 + √(k²·sin²(θ) + cos²(θ)) · B/B0]^α
         * @endcode
         *
         * where:
         * - Jc0: Critical current density at zero field [A/m²]
         * - B0: Characteristic field scale [T]
         * - k: Anisotropy parameter (k < 1 means field parallel to the
         *      c-axis / surface normal reduces Jc more than in-plane field;
         *      k = 1 is isotropic)
         * - α: Exponent controlling field dependence (typically 0.5-1.0)
         * - θ: Angle between surface normal and magnetic field [rad]
         *
         * Dependencies: normB, angleNxB
         *
         * Reference:
         * - 10.1088/0953-2048/24/6/065024 (Ekin scaling law with angular correction)
         *
         * Usage:
         * @code
         * JcFunction* jc = new JcFunctionModifiedKim(
         *     1e9,   // Jc0 = 1 GA/m² at zero field
         *     5.0,   // B0 = 5 T characteristic field
         *     0.5,   // k² = 0.5 anisotropy
         *     2.0    // α = 2.0 exponent
         * );
         * material->set_jc_function(jc);  // Material takes ownership
         * @endcode
         */
        class JcFunctionModifiedKim : public JcFunction
        {
            const real mJc0 ;    //!< Critical current at zero field [A/m²]
            const real mB0 ;     //!< Characteristic field scale [T]
            const real mK2 ;     //!< Anisotropy parameter k² [-]
            const real mAlpha ;  //!< Field dependence exponent α [-]

        public:

            /**
             * @brief Constructor for modified Kim model
             * @param jc0 Critical current density at zero field [A/m²]
             * @param B0 Characteristic field scale [T]
             * @param k2 Anisotropy parameter k² (dimensionless, typically 0.1-10)
             * @param alpha Field dependence exponent α (typically 0.5-2.0)
             *
             * Creates a modified Kim function with angular dependence.
             * The model is valid for normB > 0 and all angles.
             *
             * Reference: 10.1088/0953-2048/24/6/065024
             */
            JcFunctionModifiedKim(
                const real jc0,
                const real B0,
                const real k2,
                const real alpha ) :
            mJc0( jc0 ), mB0( B0 ), mK2( k2 ), mAlpha( alpha )
            {
                this->set_dependency( JcParameter::normB );
                this->set_dependency( JcParameter::angleNxB );
            }

            /**
             * @brief Destructor
             */
            ~JcFunctionModifiedKim() override = default;

            /**
             * @brief Evaluate Jc at given field and angle
             * @param normB Magnetic field magnitude [T]
             * @param angle Angle between surface normal and field [rad]
             * @return Critical current density Jc [A/m²]
             *
             * Evaluates the modified Kim model:
             * Jc = Jc0 / [1 + √(k²·sin²(θ) + cos²(θ)) · B/B0]^α
             *
             * Even in θ by construction (cos², sin² only), so the unfolded
             * angle θ ∈ [0, π] delivered since 2026-08-16 evaluates
             * identically to the historical folded [0, π/2] input — this law
             * needs no fold of its own.
             */
            real
            eval( const real normB, const real angle  ) const override
            {
                real c = std::cos( angle );
                real c2 = c * c;
                real s2 = 1.0 - c2;

                return mJc0 /
                    std::pow( 1. +
                        std::sqrt( ( mK2 * s2 + c2 ) ) * normB / mB0, mAlpha );

            }
        };
    }
}
#endif //BELFEM_CL_JCFUNCTIONMODIFIEDKIM_HPP