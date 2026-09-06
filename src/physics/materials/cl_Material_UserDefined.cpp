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

#include <dlfcn.h>
#include "cl_Vector.hpp"
#include "cl_Material_UserDefined.hpp"
#include "fn_Material_UserDefinedPolynomials.hpp"
#include "cl_JcFunction_UserDefined.hpp"
#include "cl_JcFunction_ModifiedKim.hpp"
#include "cl_JcFunction_Database.hpp"

namespace belfem
{
    namespace  material
    {
        UserDefinedMaterial::UserDefinedMaterial(
            const string & aLibraryPath,
            const string & aLabel ) :
            Material( MaterialType::UserDefined )
        {
            this->set_label( aLabel );
            mUserFunctions.set_size( gNumMaterialProperties, nullptr );
            mPolynomials.set_size( gNumMaterialProperties, {} );

            mHandle = dlopen( aLibraryPath.c_str(), RTLD_LAZY );

            BELFEM_ERROR( mHandle != nullptr, "Could not load library %s ( error: %s )",
                aLibraryPath.c_str(), dlerror() );

            string tFunctionName = sprint( "%s_init", aLabel.c_str() );

            void (* tInitFunc )( Material * ) =
                reinterpret_cast< void(*)( Material * ) >( dlsym( mHandle, tFunctionName.c_str() ) );

            BELFEM_ERROR( tInitFunc != nullptr, "Could not find function %s in library %s ( error: %s )",
                 tFunctionName.c_str(), aLibraryPath.c_str(), dlerror() );

            tInitFunc( this );
        }

        UserDefinedMaterial::~UserDefinedMaterial()
        {
            dlclose( mHandle );
        }

        void
        UserDefinedMaterial::set_user_defined_function(
                const MaterialProperty   Property,
                const MaterialDependency Dependency,
                      MatFunc1 * Function )
        {
            BELFEM_ERROR( Dependency == MaterialDependency::T,
                "parameter of function %s for %s must be temperature in K" ,
                to_string( Property ).c_str(), this->label().c_str() );

            // the power law evaluates jc and n through a JcFunction; a T-only
            // callback registered here would sit in a slot nothing reads
            BELFEM_ERROR( Property != MaterialProperty::jc && Property != MaterialProperty::n,
                "%s of %s must depend on the field: register it with the "
                "( normB, angleNxB ) or ( normB, angleNxB, T ) overload",
                to_string( Property ).c_str(), this->label().c_str() );

            mUserFunctions( static_cast<size_t>(Property) ) = *Function;
            this->set_custom( Property );
        }


        void
        UserDefinedMaterial::set_user_defined_function(
                        const MaterialProperty   Property,
                        const MaterialDependency Dependency1,
                        const MaterialDependency Dependency2,
                              MatFunc2 * Function )
        {


            switch ( Property )
            {
                case MaterialProperty::mu:
                {
                    BELFEM_ERROR( Dependency1 == MaterialDependency::normH,
                        "first parameter of function %s for %s must be a magnetic field in T" ,
                        to_string( Property ).c_str(), this->label().c_str() );

                    BELFEM_ERROR( Dependency2 == MaterialDependency::T,
                        "second parameter of function %s for %s must be temperature in K" ,
                        to_string( Property ).c_str(), this->label().c_str() );
                    mUserMuFunction = Function;
                    mMuFunction = & UserDefinedMaterial::mu_custom ;
                    this->set_custom( Property );
                    break ;
                }
                case MaterialProperty::jc:
                {
                    BELFEM_ERROR( mJcFunction == nullptr,
                        "jc function of material %s is already set", this->label().c_str() );

                    JcFunctionUserDefined * tJc = new JcFunctionUserDefined( this );
                    tJc->set_jc_function( Function );
                    mJcFunction = tJc;
                    this->set_have( MaterialProperty::jc ) ;
                    this->set_dependency( MaterialProperty::jc, MaterialDependency::normB ) ;
                    this->set_dependency( MaterialProperty::jc, MaterialDependency::angleNxB ) ;
                    break ;
                }
                case MaterialProperty::n:
                {
                    BELFEM_ERROR( mNFunction == nullptr,
                        "n function of material %s is already set", this->label().c_str() );

                    JcFunctionUserDefined * tN = new JcFunctionUserDefined( this );
                    tN->set_jc_function( Function );
                    mNFunction = tN;
                    this->set_have( MaterialProperty::n ) ;
                    this->set_dependency( MaterialProperty::n, MaterialDependency::normB ) ;
                    this->set_dependency( MaterialProperty::n, MaterialDependency::angleNxB ) ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false,
                       "two-parameter model for property %s for %s is not supported",
                       to_string( Property ).c_str(), this->label().c_str() );
                }
            }
        }

        void
        UserDefinedMaterial::set_user_defined_function(
                       const MaterialProperty   Property,
                       const MaterialDependency Dependency1,
                       const MaterialDependency Dependency2,
                       const MaterialDependency Dependency3,
                             MatFunc3 * Function )
        {
            switch ( Property )
            {
                case MaterialProperty::lambda:
                case MaterialProperty::rho:
                {
                    BELFEM_ERROR( Dependency1 == MaterialDependency::T,
                        "first parameter of function %s for %s must be a temperature in K" ,
                        to_string( Property ).c_str(), this->label().c_str() );

                    BELFEM_ERROR( Dependency2 == MaterialDependency::normB,
                        "second parameter of function %s for %s must be a magnetic field in T" ,
                        to_string( Property ).c_str(), this->label().c_str() );

                    BELFEM_ERROR( Dependency3 == MaterialDependency::angleBxJ,
                        "third parameter of function %s for %s must be an angle in rad" ,
                        to_string( Property ).c_str(), this->label().c_str() );

                    this->reset_dependencies( Property );
                    this->set_dependency( Property, MaterialDependency::T );
                    this->set_dependency( Property, MaterialDependency::normB );
                    this->set_dependency( Property, MaterialDependency::angleBxJ );

                    if( Property == MaterialProperty::lambda )
                    {
                        mUserLambdaFunction = Function;
                    }
                    else
                    {
                        mUserRhoFunction = Function;
                    }
                    this->set_custom( Property );

                    break ;
                }
                case MaterialProperty::jc:
                case MaterialProperty::n:
                {
                    BELFEM_ERROR( Dependency1 == MaterialDependency::normB,
                       "first parameter of function %s for %s must be magnetic field in T" ,
                       to_string( Property ).c_str(), this->label().c_str() );

                    BELFEM_ERROR( Dependency2 == MaterialDependency::angleNxB,
                        "second parameter of function %s for %s must be an angle in rad" ,
                        to_string( Property ).c_str(), this->label().c_str() );

                    BELFEM_ERROR( Dependency3 == MaterialDependency::T,
                       "second parameter of function %s for %s must be a temperature in K" ,
                       to_string( Property ).c_str(), this->label().c_str() );

                    this->reset_dependencies( Property );
                    this->set_dependency( Property, MaterialDependency::normB );
                    this->set_dependency( Property, MaterialDependency::angleNxB );
                    this->set_dependency( Property, MaterialDependency::T );

                    if ( Property == MaterialProperty::jc )
                    {
                        BELFEM_ERROR( mJcFunction == nullptr,
                            "jc function of material %s is already set", this->label().c_str() );

                        this->set_have( MaterialProperty::jc ) ;
                        JcFunctionUserDefined * tJc = new JcFunctionUserDefined( this );
                        tJc->set_jc_function( Function );
                        mJcFunction = tJc;
                    }
                    else
                    {
                        BELFEM_ERROR( mNFunction == nullptr,
                            "n function of material %s is already set", this->label().c_str() );

                        this->set_have( MaterialProperty::n ) ;
                        JcFunctionUserDefined * tN = new JcFunctionUserDefined( this );
                        tN->set_jc_function( Function );
                        mNFunction = tN;
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false,
                    "three-parameter model for property %s for %s is not supported",
                    to_string( Property ).c_str(), this->label().c_str() );
                }
            }
        }

        void
        UserDefinedMaterial::set_user_defined_polynomial(
            const MaterialProperty Property,
            const Cell< real > & Coefficients )
        {
            mPolynomials( static_cast<size_t>(Property) ) = Coefficients ;

            this->set_user_defined_function( Property );
        }

        void
        UserDefinedMaterial::set_user_defined_polynomial(
            const MaterialProperty Property,
            const std::vector< real > & Coefficients )
        {

            mPolynomials( static_cast<size_t>(Property) ).vector_data() = Coefficients;

            this->set_user_defined_function( Property );
        }


        void
        UserDefinedMaterial::set_user_defined_function( const MaterialProperty Property )
        {
              switch ( Property )
            {
                case MaterialProperty::E :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_E );
                    break ;
                }
                case MaterialProperty::nu :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_nu );
                    break ;
                }
                case MaterialProperty::cp :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_cp );
                    break ;
                }
                case MaterialProperty::lambda :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_lambda );
                    break ;
                }
                case MaterialProperty::mu :
                {
                    this->set_user_defined_function( Property,
                        MaterialDependency::normH,
                        MaterialDependency::T,
                        polyval_mu );
                    mMuFunction = & UserDefinedMaterial::mu_polynomial ;
                    break ;
                }
                case MaterialProperty::rho :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_rho );
                    break ;
                }
                case MaterialProperty::alpha :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_alpha );
                    break ;
                }
                case MaterialProperty::Rp02 :
                {
                    this->set_user_defined_function( Property, MaterialDependency::T, polyval_Rp02 );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "polynomial model for property %s for %s is not supported",
                        to_string( Property ).c_str(), this->label().c_str() );
                }
            }
        }
    }
}