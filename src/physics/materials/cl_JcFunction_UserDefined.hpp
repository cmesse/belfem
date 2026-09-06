/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_JCFUNCTION_USERDEFINED_HPP
#define BELFEM_CL_JCFUNCTION_USERDEFINED_HPP

#include "cl_JcFunction.hpp"
#include "cl_Material.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief User-defined critical current density (Jc) or power law exponent (n) function
         *
         * This class provides a wrapper for user-defined functions that compute Jc or n
         * for superconducting materials. It supports both 2-parameter and 3-parameter
         * function signatures:
         *
         * - MatFunc2: f(Material*, normB, angle) - Field and angle dependent
         * - MatFunc3: f(Material*, normB, angle, T) - Field, angle, and temperature dependent
         *
         * ANGLE CONTRACT (changed 2026-08-16): the angle is the
         * UNFOLDED field-to-tape-normal angle θ ∈ [0, π] — before that change it
         * arrived folded into [0, π/2]. A user function that is even in θ
         * (cos², sin² only) is unaffected; one that used cos(θ) linearly or
         * assumed θ ≤ π/2 must fold internally now. Plugins compiled against
         * the old contract get the new angle silently — recompile and review.
         *
         * The user provides a function pointer that will be called to evaluate the property.
         * This allows complete customization while maintaining compatibility with the
         * material framework.
         *
         * IMPORTANT: This class is typically created and managed by UserDefinedMaterial.
         * The material owns the JcFunctionUserDefined instance and will delete it.
         *
         * Example usage (within a user-defined material library):
         * @code
         * // User's custom Jc function
         * real my_jc_function(const Material* mat, real B, real angle, real T) {
         *     // Custom implementation
         *     return 1e9 / (1.0 + B/5.0) * exp(-T/90.0);
         * }
         *
         * // In the material initialization function
         * extern "C" void MyMaterial_init(Material* mat) {
         *     mat->set_user_defined_function(
         *         MaterialProperty::jc,
         *         MaterialDependency::normB,
         *         MaterialDependency::angleNxB,
         *         MaterialDependency::T,
         *         &my_jc_function
         *     );
         * }
         * @endcode
         */
        class JcFunctionUserDefined : public JcFunction
        {
            Material * mMaterial ;  //!< Pointer to the parent material (not owned)

            MatFunc2  * mUserJcFunction2 = nullptr ;  //!< User function pointer (2 parameters)
            MatFunc3  * mUserJcFunction3 = nullptr ;  //!< User function pointer (3 parameters)

        public:

            /**
             * @brief Constructor
             * @param Mat Pointer to the parent material (not owned by this class)
             */
            JcFunctionUserDefined( Material * Mat ) : mMaterial( Mat ) {};

            JcFunctionUserDefined( const JcFunctionUserDefined & ) = delete;
            JcFunctionUserDefined & operator=( const JcFunctionUserDefined & ) = delete;

            ~JcFunctionUserDefined() override = default;

            /**
             * @brief Evaluate the function (2-parameter version)
             * @param normB Magnetic flux density magnitude [T]
             * @param angle Angle between field and normal direction [rad]
             * @return Critical current density [A/m²] or power law exponent [dimensionless]
             */
            real
            eval( const real normB, const real angle ) const override ;

            /**
             * @brief Evaluate the function (3-parameter version)
             * @param normB Magnetic flux density magnitude [T]
             * @param angle Angle between field and normal direction [rad]
             * @param T Temperature [K]
             * @return Critical current density [A/m²] or power law exponent [dimensionless]
             */
            real
            eval( const real normB, const real angle, const real T ) const override ;

            /**
             * @brief Set the user-defined function pointer (2-parameter version)
             * @param Function Pointer to user function with signature: real(Material*, real, real)
             */
            void
            set_jc_function( MatFunc2 * Function )
            {
                mUserJcFunction2 = Function ;
                this->set_dependency( JcParameter::normB );
                this->set_dependency( JcParameter::angleNxB );
            }

            /**
             * @brief Set the user-defined function pointer (3-parameter version)
             * @param Function Pointer to user function with signature: real(Material*, real, real, real)
             */
            void
            set_jc_function( MatFunc3 * Function )
            {
                mUserJcFunction3 = Function ;
                this->set_dependency( JcParameter::normB );
                this->set_dependency( JcParameter::angleNxB );
                this->set_dependency( JcParameter::T );
            }

        private:

            /**
             * @brief Internal wrapper for 2-parameter user function
             */
            real
            jc_user_2( const real normB, const real angle ) const
            {
                return mUserJcFunction2( mMaterial, normB, angle );
            }

            /**
             * @brief Internal wrapper for 3-parameter user function
             */
            real
            jc_user_3( const real normB, const real angle, const real T ) const
            {
                return mUserJcFunction3( mMaterial, normB, angle, T );
            }
        };

        inline real
        JcFunctionUserDefined::eval( const real normB, const real angle ) const
        {
            BELFEM_ERROR( mUserJcFunction2 != nullptr,
                "2-parameter user function is not set" );

            return jc_user_2( normB, angle );
        }

        inline real
        JcFunctionUserDefined::eval( const real normB, const real angle, const real T ) const
        {
            BELFEM_ERROR( mUserJcFunction3 != nullptr,
                "3-parameter user function is not set" );

            return jc_user_3( normB, angle, T );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_JCFUNCTION_USERDEFINED_HPP