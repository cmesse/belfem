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

#ifndef BELFEM_CL_MATERIAL_USERDEFINED_HPP
#define BELFEM_CL_MATERIAL_USERDEFINED_HPP

#include "cl_Material.hpp"

// note: NEVER include cl_Vector.hpp or cl_Matrix.hpp here, or any class that uses it.
//       Doing so would break the API for the user defined materials.

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        template < typename T >
        T
        polyval( const std::vector< T > & coeffs, const T x )
        {
            const index_t n = coeffs.size() ;
            T y = coeffs[ 0 ];

            for( index_t k=1; k<n; ++k )
            {
                y *= x;
                y += coeffs[ k ];
            }

            return y;
        }

//------------------------------------------------------------------------------

        template < typename T >
        T
        dpolyval( const std::vector< T > & coeffs, const T x )
        {
            const index_t n = coeffs.size() - 1 ;
            T p = ( T ) n ;
            T dydx = p * coeffs[ 0 ];

            for( index_t k=1; k<n; ++k )
            {
                dydx *= x;
                p -= 1.0 ;
                dydx += p * coeffs[ k ];
            }

            return dydx;
        }

//------------------------------------------------------------------------------

        /**
         * @brief User-defined material loaded from external shared library
         *
         * This class enables users to define custom materials by implementing property
         * functions in a separate shared library (.so on Linux, .dylib on macOS).
         * The library is dynamically loaded at runtime using dlopen.
         *
         * LOADING MECHANISM:
         * The library must contain an initialization function with the signature:
         * @code
         * extern "C" void \<label\>_init(Material* mat);
         * @endcode
         * where \<label\> is the material label provided to the constructor.
         *
         * PROPERTY DEFINITION:
         * Within the init function, users can define material properties using:
         * - set_user_defined_function() - For custom property functions
         * - set_user_defined_polynomial() - For polynomial temperature dependencies
         * - set_constant() - For constant property values
         * - set_jc_function() / set_n_function() - For superconductor properties
         *
         * SUPPORTED PROPERTIES:
         * - Mechanical: E(T), nu(T), Rp02(T), alpha(T)
         * - Thermal: cp(T), lambda(T) or lambda(T,B,angle), rho(T) or rho(T,B,angle)
         * - Magnetic: mu(H,T)
         * - Superconducting: jc(B,angle) or jc(B,angle,T), n(B,angle) or n(B,angle,T)
         *
         * Example user library (mymat.cpp):
         * @code
         * #include "cl_Material.hpp"
         * using namespace belfem;
         *
         * // Custom resistivity function
         * real my_rho(const Material* mat, real T) {
         *     return 1.7e-8 * (1.0 + 0.004 * (T - 293.0));
         * }
         *
         * // Initialization function (must be extern "C")
         * extern "C" void MyAlloy_init(Material* mat) {
         *     // Set constant properties
         *     mat->set_constant(MaterialProperty::E, 200e9);
         *     mat->set_constant(MaterialProperty::nu, 0.3);
         *
         *     // Set custom resistivity function
         *     mat->set_user_defined_function(MaterialProperty::rho,
         *                                     MaterialDependency::T,
         *                                     &my_rho);
         *
         *     // Set polynomial for specific heat: cp = 385 + 0.12*T  (DESCENDING order: T¹, T⁰)
         *     std::vector<real> cp_coeffs = {0.12, 385.0};
         *     mat->set_user_defined_polynomial(MaterialProperty::cp, cp_coeffs);
         * }
         * @endcode
         *
         * Compile as shared library. An installed BELFEM keeps the module
         * layout under include/belfem, and the headers include each other by
         * bare name, so three directories are needed:
         * @code
         * g++ -shared -fPIC mymat.cpp -o libmymat.so \
         *     -I<prefix>/include/belfem/core \
         *     -I<prefix>/include/belfem/containers \
         *     -I<prefix>/include/belfem/physics/materials
         * @endcode
         * The shipped share/belfem/templates/UserMaterialTemplate.cmake does
         * this for both an installed prefix and a source tree.
         *
         * Usage in BELFEM:
         * @code
         * MaterialFactory factory;
         * Material* mat = factory.create_material("libmymat.so", "MyAlloy");
         * @endcode
         */
        class UserDefinedMaterial : public Material
        {
            void *  mHandle ;  //!< Handle to dynamically loaded library (from dlopen)

            Cell< MatFunc1 * > mUserFunctions ;  //!< Array of single-parameter user functions

            Cell< Cell< real > > mPolynomials ;  //!< Polynomial coefficients for each property

            real ( UserDefinedMaterial::*mMuFunction ) ( const real H, const real T ) const = nullptr ;  //!< Function pointer for mu evaluation

            MatFunc2 * mUserMuFunction = nullptr ;  //!< User-defined mu(H,T) function
            MatFunc3 * mUserRhoFunction = nullptr ;  //!< User-defined rho(T,B,angle) function
            MatFunc3 * mUserLambdaFunction = nullptr ;  //!< User-defined lambda(T,B,angle) function


            public :

                /**
                 * @brief Constructor - loads and initializes user-defined material
                 *
                 * Dynamically loads the shared library at aLibraryPath and calls the
                 * initialization function "\<aLabel\>_init" to configure material properties.
                 *
                 * @param aLibraryPath Path to shared library (.so, .dylib, .dll)
                 * @param aLabel Material label (must match init function name)
                 * @throws Error if library cannot be loaded or init function not found
                 */
                UserDefinedMaterial( const string & aLibraryPath, const string & aLabel ) ;

                /**
                 * @brief Destructor - closes the dynamically loaded library
                 */
                ~UserDefinedMaterial() override;

            /**
             * @brief Set a user-defined function with one dependency (typically temperature)
             *
             * Assigns a custom function for evaluating a material property. The function
             * receives the material pointer and one parameter (typically temperature).
             *
             * Supported properties: E, nu, cp, lambda, rho, alpha, Rp02
             *
             * @param Property The material property to define
             * @param Dependency The dependency type (must be MaterialDependency::T)
             * @param Function Pointer to user function: real(const Material*, real)
             * @throws Error if Dependency is not T, or if Property is jc or n
             *         (those take the field-dependent overloads below)
             *
             * Example:
             * @code
             * real my_cp(const Material* mat, real T) {
             *     return 385.0 + 0.12*T;
             * }
             * mat->set_user_defined_function(MaterialProperty::cp,
             *                                 MaterialDependency::T,
             *                                 &my_cp);
             * @endcode
             */
            void
            set_user_defined_function(
                const MaterialProperty   Property,
                const MaterialDependency Dependency,
                      MatFunc1 * Function ) override ;

            /**
             * @brief Set a user-defined function with two dependencies
             *
             * Assigns a custom function for evaluating a material property with two
             * parameters. Currently supports:
             * - mu(H, T): Magnetic permeability as function of field and temperature
             * - jc(B, angle): Critical current density as function of field and angle
             * - n(B, angle): Power law exponent as function of field and angle
             *
             * @param Property The material property to define (mu, jc, or n)
             * @param Dependency1 First dependency (normH for mu, normB for jc/n)
             * @param Dependency2 Second dependency (T for mu, angleNxB for jc/n)
             * @param Function Pointer to user function: real(const Material*, real, real)
             * @throws Error if property or dependency combination is not supported
             *
             * Example:
             * @code
             * real my_jc(const Material* mat, real B, real angle) {
             *     return 1e9 / (1.0 + B/5.0);
             * }
             * mat->set_user_defined_function(MaterialProperty::jc,
             *                                 MaterialDependency::normB,
             *                                 MaterialDependency::angleNxB,
             *                                 &my_jc);
             * @endcode
             */
            void
            set_user_defined_function(
                const MaterialProperty   Property,
                const MaterialDependency Dependency1,
                const MaterialDependency Dependency2,
                      MatFunc2 * Function ) override ;

            /**
             * @brief Set a user-defined function with three dependencies
             *
             * Assigns a custom function for evaluating a material property with three
             * parameters. Currently supports:
             * - lambda(T, B, angle): Thermal conductivity with field dependence
             * - rho(T, B, angle): Electrical resistivity with field dependence
             * - jc(B, angle, T): Critical current density with temperature
             * - n(B, angle, T): Power law exponent with temperature
             *
             * Argument order differs by property: lambda/rho take
             * (T, normB, angleBxJ); jc/n take (normB, angleNxB, T).
             *
             * @param Property The material property to define
             * @param Dependency1 First dependency (T for lambda/rho, normB for jc/n)
             * @param Dependency2 Second dependency (normB for lambda/rho, angleNxB for jc/n)
             * @param Dependency3 Third dependency (angleBxJ for lambda/rho, T for jc/n)
             * @param Function Pointer to user function: real(const Material*, real, real, real)
             * @throws Error if property or dependency combination is not supported
             *
             * Example:
             * @code
             * real my_rho(const Material* mat, real T, real B, real angle) {
             *     real rho0 = 1.7e-8 * (1.0 + 0.004*(T-293.0));
             *     real delta_rho = rho0 * 0.01 * B*B;  // Magnetoresistance
             *     return rho0 + delta_rho;
             * }
             * mat->set_user_defined_function(MaterialProperty::rho,
             *                                 MaterialDependency::T,
             *                                 MaterialDependency::normB,
             *                                 MaterialDependency::angleBxJ,
             *                                 &my_rho);
             * @endcode
             */
            void
            set_user_defined_function(
                const MaterialProperty   Property,
                const MaterialDependency Dependency1,
                const MaterialDependency Dependency2,
                const MaterialDependency Dependency3,
                      MatFunc3 * Function ) override ;

            /**
             * @brief Set a polynomial function for a material property
             *
             * Defines a property as a polynomial in temperature: f(T) = c₀T^n + c₁T^(n-1) + ... + c_n
             * This is a convenience function that internally creates a user function
             * using polyval().
             *
             * IMPORTANT: Coefficients are in DESCENDING order (MATLAB style), highest degree first.
             *
             * Supported properties: E, nu, cp, lambda, mu, rho, alpha, Rp02
             *
             * @param Property The material property to define
             * @param Coefficients Polynomial coefficients [c₀, c₁, c₂, ...] in DESCENDING order (highest degree first)
             * @throws Error if property does not support polynomial definition
             *
             * Example:
             * @code
             * // cp(T) = -1e-5*T² + 0.12*T + 385
             * std::vector<real> cp_coeffs = {-1e-5, 0.12, 385.0};  // Descending order: T², T¹, T⁰
             * mat->set_user_defined_polynomial(MaterialProperty::cp, cp_coeffs);
             * @endcode
             */
            void
            set_user_defined_polynomial(
                const MaterialProperty Property,
                const Cell< real > & Coefficients ) override;

            void
            set_user_defined_polynomial(
               const MaterialProperty Property,
               const std::vector< real > & Coefficients ) override ;

            /**
             * @brief Evaluate a polynomial function for a given property and temperature
             *
             * Internal method used by polynomial-based property evaluations.
             * Uses Horner's method for efficient polynomial evaluation.
             *
             * @param Property The material property
             * @param T Temperature [K]
             * @return Property value at temperature T
             * @throws Error if polynomial for this property is not defined
             */
            real
            evaluate_polynomial(
                const MaterialProperty  Property,
                const real T ) const override ;

            real
            evaluate_derivative_of_polynomial(
                const MaterialProperty  Property,
                const real T ) const override ;

        protected:

            //! @brief Custom Young's modulus evaluation - calls user function
            real
            E_custom(const real T) const override
            {
                return mUserFunctions(  static_cast<size_t>(MaterialProperty::E) )( this, T );
            }

            //! @brief Custom Poisson's ratio evaluation - calls user function
            real
            nu_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::nu))( this, T);
            }

            //! @brief Custom specific heat evaluation - calls user function
            real
            cp_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::cp))( this, T);
            }

            //! @brief Custom thermal conductivity evaluation (temperature only) - calls user function
            real
            lambda_custom( const real T ) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::lambda))( this, T);
            }

            //! @brief Custom thermal conductivity evaluation (with field) - calls user function
            real
            lambda_custom(const real T, const real normB, const real angle) const override
            {
                return ( * mUserLambdaFunction ) (this, T, normB, angle );
            }

            //! @brief Custom electrical resistivity evaluation (temperature only) - calls user function
            real
            rho_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::rho))( this, T );
            }

            //! @brief Custom electrical resistivity evaluation (with field) - calls user function
            real
            rho_kohler(const real T, const real normB, const real angle) const override
            {
                return ( * mUserRhoFunction ) (this, T, normB, angle );
            }

            //! @brief Custom jc evaluation (temperature only) - calls user function
            real
            jc_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::jc))( this, T );
            }

            //! @brief Custom n evaluation (temperature only) - calls user function
            real
            n_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::n))( this, T );
            }

            //! @brief Custom thermal expansion coefficient evaluation - calls user function
            real
            alpha_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::alpha))( this, T);
            }

            //! @brief Custom 0.2% proof stress evaluation - calls user function
            real
            Rp02_custom(const real T) const override
            {
                return mUserFunctions( static_cast<size_t>(MaterialProperty::Rp02))( this, T);
            }

            //! @brief Custom magnetic permeability evaluation - dispatches to polynomial or user function
            real
            mu_custom( const real H, const real T ) const override
            {
                return ( this->*mMuFunction )( H, T );
            }

        private:

            //! @brief Evaluate mu using polynomial (ignores H)
            real
            mu_polynomial( const real H, const real T ) const
            {
                return evaluate_polynomial( MaterialProperty::mu, T );
            }

            //! @brief Evaluate mu using user-provided function
            real
            mu_user( const real H, const real T ) const
            {
                return mUserMuFunction( this, H, T );
            }

            void
            set_user_defined_function( const MaterialProperty Property );

        };

//------------------------------------------------------------------------------

        inline real
        UserDefinedMaterial::evaluate_polynomial( const MaterialProperty  Property, const real T ) const
        {
            BELFEM_ASSERT( mPolynomials( static_cast<size_t>(Property) ).size() > 0,
                "return polynomial for %s of %s is not defined",
                to_string( Property ).c_str(),
                this->label().c_str()
            );

            return polyval( mPolynomials( static_cast<size_t>(Property) ).vector_data(), T );
        }

        inline real
        UserDefinedMaterial::evaluate_derivative_of_polynomial( const MaterialProperty  Property, const real T ) const
        {
            BELFEM_ASSERT( mPolynomials( static_cast<size_t>(Property) ).size() > 0,
                "return polynomial for %s of %s is not defined",
                to_string( Property ).c_str(),
                this->label().c_str()
            );

            return dpolyval( mPolynomials( static_cast<size_t>(Property) ).vector_data(), T );
        }

    }
}
#endif //BELFEM_CL_MATERIAL_USERDEFINED_HPP