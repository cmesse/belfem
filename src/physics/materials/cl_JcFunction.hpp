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

#ifndef BELFEM_CL_JCFUNCTION_HPP
#define BELFEM_CL_JCFUNCTION_HPP

/**
 * @file cl_JcFunction.hpp
 * @brief Base class hierarchy for critical current density (Jc) and n-value functions
 *
 * This file defines the abstract base class JcFunction and its parameter system
 * for representing field-, angle-, and temperature-dependent properties in
 * high-temperature superconductors (HTS).
 *
 * ## Usage
 *
 * JcFunction objects are created by MaterialFactory and assigned to Material
 * objects. The Material takes ownership and is responsible for deletion:
 *
 * @code
 * // Create analytical Jc function (modified Kim model)
 * JcFunction* jc = factory.create_jc_function(
 *     1e9,   // Jc0 [A/m²]
 *     5.0,   // B0 [T]
 *     0.5,   // k² anisotropy
 *     2.0    // α exponent
 * );
 * material->set_jc_function(jc);  // Material takes ownership
 *
 * // Or create from database file ( same creator for jc and n tables )
 * JcFunction* n_value = factory.create_jc_function(
 *     "path/to/n_value_data.h5", "n"
 * );
 * material->set_n_function(n_value);  // Material takes ownership
 * @endcode
 *
 * ## Derived Classes
 *
 * - JcFunctionModifiedKim: Analytical model with field and angle dependence
 * - JcFunctionDatabase: Lookup table with field, angle, and temperature dependence
 *
 * ## Important Notes
 *
 * 1. OWNERSHIP: Material takes ownership - do NOT delete manually
 * 2. UNITS: Jc in [A/m²], B in [T], angle in [rad], T in [K]
 * 3. Same class used for both Jc and n-value (identical dependencies)
 * 4. Constant values stored directly in Material, not as JcFunction objects
 *
 * Dependencies:
 * - normB: Magnitude of magnetic field [T]
 * - angleNxB: Angle between surface normal and field [rad]
 * - T: Temperature [K]
 */

#include "typedefs.hpp"
#include "cl_Bitset.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Parameters that Jc (or n) can depend on
         *
         * - normB: Magnitude of magnetic field [T]
         * - angleNxB: Angle between surface normal and field [rad]
         * - T: Temperature [K]
         */
        enum class JcParameter
        {
            normB          = 1,
            angleNxB       = 2,
            T              = 3,
            UNDEFINED      = 4
        };

        /**
         * @brief Bitset for tracking Jc function dependencies
         */
        typedef Bitset< static_cast< uint >( JcParameter::UNDEFINED ) > JcDependency;

//------------------------------------------------------------------------------

        /**
         * @brief Base class for critical current density (Jc) and n-value functions
         *
         * This abstract base class provides the interface for computing:
         * - Jc(B, angle, T): Critical current density [A/m²]
         * - n(B, angle, T): Power law exponent [-]
         *
         * The same class is used for both properties since they have identical
         * functional dependencies and evaluation signatures.
         *
         * ## Derived Classes
         *
         * - JcFunctionModifiedKim: Analytical modified Kim model
         *   * Dependencies: normB, angleNxB
         *   * Overrides: eval(B, angle)
         *
         * - JcFunctionDatabase: Lookup table from HDF5 file
         *   * Dependencies: normB, angleNxB, T
         *   * Overrides: eval(B, angle, T)
         *
         * ## Virtual Function Behavior
         *
         * Both eval() methods have default implementations that throw errors.
         * Derived classes override the appropriate method(s) based on their
         * dependency set:
         * - If depends on T: Override eval(B, angle, T)
         * - If no T dependence: Override eval(B, angle)
         *
         * The Material class calls the appropriate eval() based on the
         * dependency flags queried via depends_on().
         *
         * ## Ownership and Lifetime
         *
         * Material takes ownership of assigned JcFunction objects:
         * @code
         * JcFunction* jc = factory.create_jc_function(1e9, 5.0, 0.5, 2.0);
         * hts->set_jc_function(jc);  // Material deletes jc in destructor
         * // Do NOT call delete jc!
         * @endcode
         *
         * @ingroup grp_physics_materials
         * @see @ref physics_materials_materials_usage_guide
         */
        class JcFunction
        {
        protected:

            JcDependency mDependency ;

        public:

            /**
             * @brief Default constructor
             */
            JcFunction() = default;

            /**
             * @brief Virtual destructor
             */
            virtual
            ~JcFunction() = default;

            /**
             * @brief Check if function depends on a parameter
             * @param aParameter Parameter to check (normB, angleNxB, or T)
             * @return True if function depends on this parameter
             */
            bool
            depends_on( const JcParameter aParameter ) const;

            /**
             * @brief Smallest value this function can return, when known
             * @return min over the table for database-driven functions;
             *         NaN when no cheap bound exists ( analytic laws )
             *
             * Used by Material::set_n_function for a load-time sanity
             * warning on measured n tables that soften toward n = 1.
             */
            virtual real
            min_value() const
            {
                return BELFEM_QUIET_NAN ;
            }

            /**
             * @brief Check if function is constant (no dependencies)
             * @return True if function has no dependencies
             */
            bool
            is_constant() const ;

            /**
             * @brief Evaluate with field and angle dependence
             * @param normB Magnetic field magnitude [T]
             * @param angle Angle between surface normal and field [rad]
             * @return Jc or n value [A/m²] or [-]
             */
            virtual real
            eval(
                const real normB,
                const real angle
                ) const;

            /**
             * @brief Evaluate with full temperature, field, and angle dependence
             * @param normB Magnetic field magnitude [T]
             * @param angle Angle between surface normal and field [rad]
             * @param T Temperature [K]
             * @return Jc or n value [A/m²] or [-]
             */
            virtual real
            eval(
                const real normB,
                const real angle,
                const real T
                ) const;

            /**
             * @brief Derivative of eval with respect to the field magnitude
             * @param normB Magnetic field magnitude [T]
             * @param angle Angle between surface normal and field [rad]
             * @param T Temperature [K]
             * @return d(jc)/d|B| [A/m²/T] or d(n)/d|B| [1/T]
             *
             * The base implementation returns 0, which is EXACT for constant
             * functions and the conservative fallback for derived classes
             * that have not implemented an analytic form ( ModifiedKim,
             * UserDefined ): a zero here reproduces the pre-lookup-derivative tangent
             * for that function rather than inventing a wrong one. Unlike
             * eval(), it does not abort — the Newton consumer
             * ( add_rho_field_tangent ) treats 0 as "no field channel" and
             * early-outs.
             */
            virtual real
            deval_dB(
                const real normB,
                const real angle,
                const real T
                ) const;

            /**
             * @brief Derivative of eval with respect to the field angle
             * @return d(jc)/dθ [A/m²/rad] or d(n)/dθ [1/rad]
             *
             * Same zero-default contract as deval_dB. NOTE: for HTS the
             * angle is bn_angle ( field to tape normal, UNFOLDED [ 0, pi ]
             * since 2026-08-16 ), and the existing Newton kernel differentiates
             * bj_angle ( field to current, the metal Kohler variable ) — so
             * this hook must NOT be wired into add_rho_field_tangent's beta
             * channel until a signed bn_angle d(theta)/dq chain exists there
             * ( 2026-08-13 audit, both voices; known residual ).
             */
            virtual real
            deval_dbeta(
                const real normB,
                const real angle,
                const real T
                ) const;

            /**
             * @brief Derivative of eval with respect to temperature
             * @return d(jc)/dT [A/m²/K] or d(n)/dT [1/K]
             *
             * Same zero-default contract as deval_dB.
             */
            virtual real
            deval_dT(
                const real normB,
                const real angle,
                const real T
                ) const;

        protected:

            /**
             * @brief Mark function as depending on a parameter
             * @param aParameter Parameter to add to dependency set
             *
             * Used by derived class constructors to declare which parameters
             * their implementation depends on.
             */
            void
            set_dependency( const JcParameter aParameter );

        };

        inline bool
        JcFunction::depends_on( const JcParameter aParameter ) const
        {
            return mDependency.test( static_cast< uint >( aParameter ) );
        }

        inline bool
        JcFunction::is_constant() const
        {
            return mDependency.count() == 0 ;
        }

    }
}
#endif //BELFEM_CL_JCFUNCTION_HPP