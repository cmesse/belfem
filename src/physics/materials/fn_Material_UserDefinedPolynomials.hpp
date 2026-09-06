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

#ifndef BELFEM_FN_MATERIAL_USERDEFINEDPOLYNOMIALS_HPP
#define BELFEM_FN_MATERIAL_USERDEFINEDPOLYNOMIALS_HPP

#include "cl_Material.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @file fn_Material_UserDefinedPolynomials.hpp
         * @brief Wrapper functions for polynomial-based material property evaluation
         *
         * These inline functions provide the interface between user-defined materials
         * and the polynomial evaluation mechanism. They are used internally by
         * UserDefinedMaterial::set_user_defined_polynomial().
         *
         * Each function has the signature required by MatFunc1 or MatFunc2, allowing
         * them to be assigned as user-defined functions while delegating to the
         * material's polynomial evaluation method.
         */

//------------------------------------------------------------------------------

        /**
         * @brief Evaluate Young's modulus using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return Young's modulus [Pa]
         */
        inline real
        polyval_E( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::E, T );
        }

        /**
         * @brief Evaluate Poisson's ratio using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return Poisson's ratio [dimensionless]
         */
        inline real
        polyval_nu( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::nu, T );
        }

        /**
         * @brief Evaluate specific heat using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return Specific heat capacity [J/(kg·K)]
         */
        inline real
        polyval_cp( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::cp, T );
        }

        /**
         * @brief Evaluate thermal conductivity using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return Thermal conductivity [W/(m·K)]
         */
        inline real
        polyval_lambda( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::lambda, T );
        }

        /**
         * @brief Evaluate magnetic permeability using polynomial (temperature-dependent only)
         *
         * Note: The magnetic field parameter H is ignored. This function evaluates
         * mu(T) only, suitable for non-ferromagnetic materials where mu does not
         * depend on applied field.
         *
         * @param Mat Material pointer
         * @param H Magnetic field strength [A/m] (ignored)
         * @param T Temperature [K]
         * @return Magnetic permeability μ = ∂B/∂H [H/m] (absolute, not relative)
         */
        inline real
        polyval_mu( const Material * Mat, const real H, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::mu, T );
        }

        /**
         * @brief Evaluate electrical resistivity using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return Electrical resistivity [Ω·m]
         */
        inline real
        polyval_rho( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::rho, T );
        }

        /**
         * @brief Evaluate thermal expansion coefficient using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return Linear thermal expansion coefficient [1/K]
         */
        inline real
        polyval_alpha( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::alpha, T );
        }

        /**
         * @brief Evaluate 0.2% proof stress using polynomial
         * @param Mat Material pointer
         * @param T Temperature [K]
         * @return 0.2% proof stress [Pa]
         */
        inline real
        polyval_Rp02( const Material * Mat, const real T )
        {
            return Mat->evaluate_polynomial( MaterialProperty::Rp02, T );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_MATERIAL_USERDEFINEDPOLYNOMIALS_HPP