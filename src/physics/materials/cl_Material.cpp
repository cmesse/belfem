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

#include <iostream>
#include <dlfcn.h>

#include "assert.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "constants.hpp"

#include "fn_material_data_path.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"

#include <algorithm>
#include <cmath>
#include "cl_Material.hpp"
#include "cl_Vector.hpp"

#include "Spline_Enums.hpp"

namespace belfem
{
    string
    to_string( const MaterialProperty aProperty )
    {
        switch ( aProperty )
        {
            case MaterialProperty::density : return "density" ;
            case MaterialProperty::cp : return "cp" ;
            case MaterialProperty::lambda : return "lambda" ;
            case MaterialProperty::mu : return "mu" ;
            case MaterialProperty::rho : return "rho" ;
            case MaterialProperty::jc : return "jc" ;
            case MaterialProperty::n : return "n" ;
            case MaterialProperty::E : return "E" ;
            case MaterialProperty::nu : return "nu" ;
            case MaterialProperty::alpha : return "alpha" ;
            case MaterialProperty::Rp02 : return "Rp02" ;
            case MaterialProperty::debye : return "debye" ;
            case MaterialProperty::M : return "M" ;
            default : return "unknown" ;
        }
    }

    unit
    get_unit( const MaterialProperty aProperty )
    {
        // base units L, M, T, I, theta, N, J
        switch ( aProperty )
        {
           // density
           case MaterialProperty::density :
           {
             return  { -3, 1, 0, 0, 0, 0, 0 } ;
           }

           // specific heat capacity
           case MaterialProperty::cp :
           {
               return { 2, 0,  -2,   0,  -1,   0,   0 } ;
           }

           // thermal conductivity
           case MaterialProperty::lambda :
           {
               return { 1, 1, -3, 0 ,-1 , 0 , 0 } ;
           }

           // magnetic permeability, defined as ∂B/∂H
           case MaterialProperty::mu :
           {
               return { 1, 1, -2, -2, 0, 0, 0 } ;
           }

           // electric resistivity
           case MaterialProperty::rho :
           {
               return { 3, 1, -3, -2, 0, 0, 0 } ;
           }

           // critical current
           case MaterialProperty::jc :
           {
               return { -2, 0, 0, 1, 0, 0, 0 } ;
           }

           // exponent for power law
           case MaterialProperty::n :
           {
               return { 0, 0, 0, 0, 0, 0, 0 } ;
           }

           // youngs modulus
           case MaterialProperty::E :
           {
               return { -1, 1, -2, 0, 0, 0, 0 } ;
           }

           // Poissoin's ratio
           case MaterialProperty::nu :
           {
               return { 0, 0, 0, 0, 0, 0, 0 } ;
           }

           // thermal expansion coefficient, defined as (1/l)*∂l/∂T, not (1/l)*Δl/ΔT !
           case MaterialProperty::alpha :
           {
               return { 0, 0, 0, 0, -1, 0, 0 } ;
           }

           // yield stress
           case MaterialProperty::Rp02 :
           {
               return { -1, 1, -2, 0, 0, 0, 0 } ;
           }

           // debye temperature
           case MaterialProperty::debye :
           {
               return { 0, 0, 0, 0, 1, 0, 0 } ;
           }

           // molar mass
           case MaterialProperty::M:
           {
               return { 0, 1, 0, 0, 0, -1, 0 } ;
           }
           default :
           {
               BELFEM_ERROR( false, "Unknown material property ");
               return { 0, 0, 0, 0, 0, 0, 0 } ;
           }
       }
   }

//------------------------------------------------------------------------------

    Material::Material( const MaterialType aType, const bool aIsIsotropic ) :
      mCommRank( comm_rank() ),
      mType( aType ),
      mIsIsotropic( aIsIsotropic )
    {
        mPropertyDependencies.set_size( gNumMaterialProperties, nullptr );


        for ( size_t i = 0; i < gNumMaterialProperties; ++i )
        {
            mPropertyDependencies( i ) = new MaterialDependencyBitset ;
        }

        mConstantProperties.set_size( gNumMaterialProperties, BELFEM_QUIET_NAN );

        // assuming mu0 per default
        this->set_constant( MaterialProperty::mu, constant::mu0 );

        // assuming one atom per molecule
        this->set_constant( MaterialProperty::q, 1.0 );

        // default density scaling factor
        this->set_constant( MaterialProperty::density_correction, 1.0 );
    }

    Material::~Material()
    {

       for ( MaterialDependencyBitset * tDep : mPropertyDependencies )
       {
           delete tDep ;
       }

       if ( mJcFunction != nullptr ) delete mJcFunction ;
       if ( mNFunction != nullptr ) delete mNFunction ;

       // last: mDefectFunction points into this library, and the resistivity
       // laws call it on every residual evaluation. dlclose( nullptr ) is not
       // a POSIX no-op, and the handle is null for every material without a
       // defect, so the guard is required rather than tidy
       if ( mDefectHandle != nullptr )
       {
           dlclose( mDefectHandle );
           mDefectHandle = nullptr ;
       }

       // same for the heating plugin
       if ( mHeatHandle != nullptr )
       {
           dlclose( mHeatHandle );
           mHeatHandle = nullptr ;
       }
    }

    // set label
    void
    Material::set_label( const string & aLabel )
    {
       mLabel = aLabel ;
    }

    void
    Material::set_number( const string & aNumber )
    {
       mNumber = aNumber ;
    }

    void
    Material::set_RRR( const real RRR )
    {
        BELFEM_ERROR( false, "Material::set_RRR not implemented for this material" );
    }

    // sets the property flag
    void
    Material::set_have( const MaterialProperty aProperty, const bool aHave )
    {
        if ( aHave )
        {
            mHaveProperty.set(  static_cast< size_t >( aProperty ) );
        }
        else
        {
            mHaveProperty.reset(  static_cast< size_t >( aProperty ) );
        }
    }

    void
    Material::reset_dependencies( const MaterialProperty aProperty )
    {
       index_t tIndex = static_cast< size_t >( aProperty ) ;
       mConstantProperties( tIndex ) = BELFEM_QUIET_NAN ;
       mPropertyDependencies( tIndex )->reset() ;
    }

    // sets the dependency flag
    void
    Material::set_dependency( const MaterialProperty aProperty, const MaterialDependency aDependency )
    {
       index_t tIndex = static_cast< size_t >( aProperty ) ;
       mConstantProperties( tIndex ) = BELFEM_QUIET_NAN ;
       mPropertyDependencies( tIndex )->set( static_cast< size_t >( aDependency ) );
    }


//------------------------------------------------------------------------------

    void
    Material::set_constant( const MaterialProperty aProperty, const real aValue )
    {
        index_t tIndex = static_cast< size_t >( aProperty ) ;
        mConstantProperties( tIndex ) = aValue ;
        mHaveProperty.set( tIndex );

        // riva already selected -- gate a late constant n ( same rule as
        // set_n_function; keeps the check order-independent )
        if ( aProperty == MaterialProperty::n &&
             mResistivityLaw == ResistivityLaw::Riva )
        {
            this->check_riva_n_source() ;
        }


        bool tComputeTheta0 = false ;

        switch ( aProperty )
        {
            case MaterialProperty::M :
            {
                size_t k = static_cast< size_t >( MaterialProperty::R ) ;
                mConstantProperties( k ) = constant::Rm / aValue ;
                mHaveProperty.set( k );
                tComputeTheta0 = this->have( MaterialProperty::beta ) ;
                break ;
            }
            case MaterialProperty::R :
            {
                size_t k = static_cast< size_t >( MaterialProperty::M ) ;
                mConstantProperties( k ) = aValue * constant::Rm ;
                mHaveProperty.set( k );
                tComputeTheta0 = this->have( MaterialProperty::beta ) ;
                break ;
            }
            case MaterialProperty::beta :
            {
                tComputeTheta0 = this->have( MaterialProperty::R );
                break ;
            }
            case MaterialProperty::q :
            {
                // a compound sets its atom count after beta ( YBCO, Magnesia )
                tComputeTheta0 = this->have( MaterialProperty::R ) && this->have( MaterialProperty::beta );
                break ;
            }
            default:
            {
                // pass
            }
        }

        if ( tComputeTheta0 )
        {
            // Debye T^3 law per atom: beta = 12 pi^4 q R / ( 5 theta^3 ), so the
            // atoms per formula unit q enter for a compound ( MgO: 2, YBCO: 13 )
            real R    = this->constant_property( MaterialProperty::R ) ;
            real q    = this->constant_property( MaterialProperty::q ) ;
            real beta = this->constant_property( MaterialProperty::beta ) ;
            real pi4  = constant::pi * constant::pi * constant::pi * constant::pi ;

            real theta = std::pow( 2.4 * pi4 * q * R / beta, 1.0 / 3.0 ) ;
            this->set_constant( MaterialProperty::debye0K, theta ) ;
        }

        if ( tIndex >= gNumNonConstantMaterialProperties ) return ;

        mPropertyDependencies( tIndex )->reset() ;

        switch ( aProperty )
        {

            case MaterialProperty::cp :
            {
                mFunctionCp      = & Material::cp_const ;
                mFunctiondCpdT   = & Material::return_zero ;
                mFunctiond2CpdT2 = & Material::return_zero ;
                break ;
            }
            case MaterialProperty::lambda :
            {
                mFunctionLambda    = & Material::lambda_const ;
                mFunctiondLambdadT = & Material::return_zero ;
                break ;
            }
            case MaterialProperty::rho :
            {
                mFunctionRho   = & Material::rho_const ;
                mFunctiondRhodT = & Material::return_zero ;
                break;
            }
            case MaterialProperty::mu :
            {
                mFunctionH     = & Material::H_const ;
                mFunctionMu    = & Material::mu_const ;
                mFunctionDMuDH = & Material::dmudH_const ;
                break ;
            }
            case MaterialProperty::E :
            {
                mFunctionE = & Material::E_const ;
                break ;
            }
            case MaterialProperty::nu :
            {
               mFunctionNu = & Material::nu_const ;
               break;
            }
            case MaterialProperty::alpha :
            {
               mFunctionAlpha = & Material::alpha_const ;
               break;
            }
            case MaterialProperty::Rp02 :
            {
               mFunctionRp02 = & Material::Rp02_const ;
               break;
            }
            case MaterialProperty::debye :
            {
               mFunctionDebye = & Material::debye_const ;
               break;
            }
            default:
            {
               BELFEM_ERROR( false, "Unknown material property ");
            }
        }
    }
//------------------------------------------------------------------------------

    void
    Material::set_custom( const MaterialProperty aProperty )
    {
        // jc and n are evaluated through a JcFunction ( jc_eval / n_eval
        // route on its own dependency flags ), so a T-only custom callback
        // is never reached on the assembly path, which calls the full
        // ( normJ, T, normB, angleNxB ) overloads. The two former cases
        // wrote into the rho_i and debye dispatch slots. Checked before the
        // have-flag is raised so a rejected call leaves the material as it was
        BELFEM_ERROR( aProperty != MaterialProperty::jc && aProperty != MaterialProperty::n,
            "Material::set_custom must not be called for %s of %s: "
            "use set_jc_function() / set_n_function()",
            to_string( aProperty ).c_str(), mLabel.c_str() );

        index_t tIndex = static_cast< size_t >( aProperty ) ;
        mConstantProperties( tIndex ) = BELFEM_QUIET_NAN ;
        mHaveProperty.set( tIndex );
        mPropertyDependencies( tIndex )->reset() ;

        this->set_dependency( aProperty, MaterialDependency::T );

        switch ( aProperty )
        {
           case MaterialProperty::density :
           {
               mFunctionDensity = & Material::density_custom ;
               break ;
           }
           case MaterialProperty::M:
           {
               BELFEM_ERROR( false, "Material::set_custom must not be called for molar-property");
               break ;
           }
           case MaterialProperty::cp :
           {
               mFunctionCp    = & Material::cp_custom ;
               mFunctiondCpdT = & Material::dcpdT_custom ;
               mFunctiond2CpdT2 = & Material::d2cpdT2_custom ;
               break ;
           }
           case MaterialProperty::lambda :
           {
                mFunctionLambda    = & Material::lambda_custom ;
                mFunctiondLambdadT = & Material::dlambdadT_custom ;
                break ;
           }
           case MaterialProperty::rho :
           {
                mFunctionRho    = & Material::rho_custom ;
                mFunctiondRhodT = & Material::drhodT_custom ;
                break;
           }
           case MaterialProperty::mu :
           {
                this->set_dependency( aProperty, MaterialDependency::normH );
                mFunctionMu = & Material::mu_custom ;
                break ;
           }
           case MaterialProperty::E :
           {
                mFunctionE = & Material::E_custom ;
                break ;
           }
           case MaterialProperty::nu :
           {
               mFunctionNu = & Material::nu_custom ;
               break;
           }
           case MaterialProperty::alpha :
           {
               mFunctionAlpha = & Material::alpha_custom ;
               break;
           }
           case MaterialProperty::Rp02 :
           {
               mFunctionRp02 = & Material::Rp02_custom ;
               break;
           }
           case MaterialProperty::rho_i :
           {
               mFunctionRhoI = & Material::rho_i_custom ;
               break;
           }
           case MaterialProperty::debye :
           {
               mFunctionDebye = & Material::debye_custom ;
               break;
           }
           default:
           {
               BELFEM_ERROR( false, "Unsupported material property ");
           }
        }
    }
//------------------------------------------------------------------------------

    real
    Material::cp_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::cp_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::dcpdT_custom( const real T ) const
    {
        return this->dcpdT_finite_difference( T );
    }


    real
    Material::d2cpdT2_custom( const real T ) const
    {
        return this->d2cpdT2_finite_difference( T );
    }

    real
    Material::lambda_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::lambda_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::dlambdadT_custom( const real T ) const
    {
        return this->dlambdadT_finite_difference( T );
    }

    real
    Material::E_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::E_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::dEdT_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::dEdT_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::nu_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::nu_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::mu_custom( const real H, const real T ) const
    {
        BELFEM_ERROR( false, "Material::mu_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::alpha_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::alpha_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::Rp02_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::Rp02_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::debye_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::debye_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    inline real
    Material::rho_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::rho_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    inline real
    Material::jc_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::jc_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    inline real
    Material::n_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::n_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real Material::rho_i_custom( const real T ) const
    {
        BELFEM_ERROR( false, "Material::rho_i_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::rho_kohler( const real T, const real B, const real beta  ) const
    {
        BELFEM_ERROR( false, "Material::rho_kohler not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::drhodT_kohler( const real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::drhodT_kohler not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::drhodB_kohler( const real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::drhodB_kohler not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::drhodbeta_kohler( const real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::drhodbeta_kohler not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::rho_table( real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::rho_table not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::drhodT_table( real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::drhodT_table not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::drhodB_table( real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::drhodB_table not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::drhodbeta_table( real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::drhodbeta_table not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::lambda_custom( const real T, const real B, const real beta  ) const
    {
        BELFEM_ERROR( false, "Material::lambda_custom not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::lambda_table( real T, const real B, const real beta ) const
    {
        BELFEM_ERROR( false, "Material::lambda_table not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    void
    Material::set_bh_curve( const material::BhCurve * aCurve )
    {
        BELFEM_ERROR( false, "not implemented for material %s", mLabel.c_str() );
    }

    void
    Material::set_jc_function( const material::JcFunction * aFunction )
    {
        BELFEM_ERROR( mJcFunction == nullptr, "Material %s already has a jc function", mLabel.c_str() );

        mJcFunction = aFunction ;

        if ( aFunction->depends_on( material::JcParameter::normB ) )
        {
            this->set_dependency( MaterialProperty::jc, MaterialDependency::normB );
        }
        if ( aFunction->depends_on( material::JcParameter::angleNxB ) )
        {
            this->set_dependency( MaterialProperty::jc, MaterialDependency::angleNxB );
        }
        if ( aFunction->depends_on( material::JcParameter::T ) )
        {
            this->set_dependency( MaterialProperty::jc, MaterialDependency::T );
        }

        this->set_have( MaterialProperty::jc );

        // set default value of ecrit if not specified
        if ( ! this->have( MaterialProperty::ec ) )
        {
            this->set_constant( MaterialProperty::ec, 1e-4 );
        }
    }

    void
    Material::set_n_function( const material::JcFunction * aFunction )
    {
        BELFEM_ERROR( mNFunction == nullptr, "Material %s already has an n function", mLabel.c_str() );

        mNFunction = aFunction ;

        // load-time sanity check on measured n tables: near T_crit they
        // soften toward the ohmic limit. Values at or below 1 are floored by
        // n_eval, and the piecewise Bezier blend degenerates below n ~ 4 —
        // both are handled, but the user should know their table gets there.
        // The threshold is 1.01, not higher: the shipped tables deliberately
        // ramp to a 1.02 nodal floor at T_crit (with the interpolant held
        // above 1.01 by their Bernstein control net), so a stricter test
        // would fire on every correctly built table. min_value() is the
        // NODAL minimum, which does not bound the interpolant — below 1.01
        // the flooring may genuinely engage.
        real tNmin = aFunction->min_value() ;
        if ( comm_rank() == 0 && tNmin < 1.01 )  // NaN ( no bound known ) fails the test
        {
            message( InfoLevel::Default,
                "    Warning: the n table of material %s reaches down to n = %g.\n"
                "             Values <= 1 are floored to the ohmic limit ( n_eval );\n"
                "             the piecewise law's flux-flow blend degenerates below n = 4.",
                mLabel.c_str(), ( double ) tNmin );
        }

        this->set_have( MaterialProperty::n );

        // riva selected before the n source arrived ( plugin _init call
        // order ) -- gate the source now instead
        if ( mResistivityLaw == ResistivityLaw::Riva )
        {
            this->check_riva_n_source() ;
        }

        if ( aFunction->depends_on( material::JcParameter::normB ) )
        {
            this->set_dependency( MaterialProperty::n, MaterialDependency::normB );
        }
        if ( aFunction->depends_on( material::JcParameter::angleNxB ) )
        {
            this->set_dependency( MaterialProperty::n, MaterialDependency::angleNxB );
        }
        if ( aFunction->depends_on( material::JcParameter::T ) )
        {
            this->set_dependency( MaterialProperty::n, MaterialDependency::T );
        }
    }

    void
    Material::set_user_defined_function(
            const MaterialProperty   Property,
            const MaterialDependency Dependency,
                  MatFunc1 * Function )
    {
        BELFEM_ERROR( false, "Material::set_user_defined_function() is not implemented for %s", mLabel.c_str() );
    }

    void
    Material::set_user_defined_function(
        const MaterialProperty   Property,
        const MaterialDependency Dependency1,
        const MaterialDependency Dependency2,
              MatFunc2 * Function )
    {
        BELFEM_ERROR( false, "Material::set_user_defined_function() is not implemented for %s", mLabel.c_str() );
    }

    void
    Material::set_user_defined_function(
        const MaterialProperty   Property,
        const MaterialDependency Dependency1,
        const MaterialDependency Dependency2,
        const MaterialDependency Dependency3,
              MatFunc3 * Function )
    {
        BELFEM_ERROR( false, "Material::set_user_defined_function() is not implemented for %s", mLabel.c_str() );
    }

    void
    Material::set_user_defined_defect(
              DefectFunc * Function )
    {
        BELFEM_ERROR( this->have( MaterialProperty::jc ), "Material::set_user_defined_defect() requires the material %s to have Jc", mLabel.c_str() );
        mDefectFunction = Function ;
    }

    void
    Material::set_user_defined_heating(
              HeatFunc * Function )
    {
        mHeatFunction = Function ;
    }

    void
    Material::set_piecewise( const bool aUsePiecewise )
    {
        mResistivityLaw = aUsePiecewise ?
            ResistivityLaw::Piecewise : ResistivityLaw::PowerLaw ;
    }

    void
    Material::set_resistivity_law( const ResistivityLaw aLaw )
    {
        if ( aLaw == ResistivityLaw::Riva )
        {
            this->check_riva_n_source() ;
        }

        mResistivityLaw = aLaw ;
    }

    void
    Material::check_riva_n_source() const
    {
        // The riva law is total over its inputs: at runtime n <= 1 folds
        // into the ohmic limit ec/jc and the run continues ( powerlaws.hpp,
        // "The riva law" ). That totality is deliberate, but it also means
        // a mis-built n source would run silently as a plain resistor --
        // refuse at setup, where the failure is attributable. Called from
        // set_resistivity_law, set_n_function and set_constant( n ), so
        // the outcome does not depend on installation order.
        //
        // The table test is the stored nodal minimum over the WHOLE table,
        // including T > T_crit where riva never reads n. Deliberate: the
        // build pipeline floors the entire table at 1.02, so any stored
        // value at or below 1 marks a table that escaped it. Two accepted
        // limits: the nodal minimum does not bound the Lagrange
        // interpolant ( undershoot at eval time is what the runtime ohmic
        // branch remains for ), and min_value() is NaN for analytic and
        // plugin fits ( no cheap bound ) -- NaN fails the <= comparison
        // and passes unchecked.
        if ( mNFunction != nullptr )
        {
            const real tNmin = mNFunction->min_value() ;

            BELFEM_ERROR( ! ( tNmin <= 1.0 ),
                "The n table of material %s reaches down to n = %g. "
                "Under the riva law n <= 1 is the plain ohmic resistor ec/jc, "
                "so this table would silently stop superconducting there. "
                "The databases are built to stay above n = 1; this one is broken.",
                mLabel.c_str(), ( double ) tNmin ) ;
        }
        else if ( this->have( MaterialProperty::n ) )
        {
            // raw slot, not constant_property(): the accessor debug-asserts
            // on NaN with an unrelated message, and NaN is one of the
            // states this gate exists to name
            const real tN = mConstantProperties(
                    static_cast< index_t >( MaterialProperty::n ) ) ;

            BELFEM_ERROR( std::isfinite( tN ) && tN > 1.0,
                "A constant n for the riva law must be finite and greater "
                "than 1 ( material %s has n = %g ).",
                mLabel.c_str(), ( double ) tN ) ;
        }
    }

    void
    Material::read_defect( const string & aLibraryPath,
                    const string & aLabel )
    {
        // one defect plugin per material: a second call would strand the
        // first mapping while mDefectFunction silently moved to the new one
        BELFEM_ERROR( mDefectHandle == nullptr,
            "Material already loaded a defect plugin; cannot load %s on top of it",
            aLibraryPath.c_str() );

        const string tPath = material::data_file( aLibraryPath );

        mDefectHandle = dlopen( tPath.c_str(), RTLD_LAZY );

        BELFEM_ERROR( mDefectHandle != nullptr, "Could not load library %s ( error: %s )",
            tPath.c_str(), dlerror() );

        string tFunctionName = sprint( "%s_init", aLabel.c_str() );

        void (* tInitFunc )( Material * ) =
            reinterpret_cast< void(*)( Material * ) >( dlsym( mDefectHandle, tFunctionName.c_str() ) );

        // the resolved path, not the name as written: naming a file the
        // loader never opened is what made this error hard to act on
        BELFEM_ERROR( tInitFunc != nullptr, "Could not find function %s in library %s ( error: %s )",
             tFunctionName.c_str(), tPath.c_str(), dlerror() );

        tInitFunc( this );
    }

    void
    Material::read_heating( const string & aLibraryPath,
                const string & aLabel )
    {
        // one heating plugin per material: a second call would strand the
        // first mapping while mHeatFunction silently moved to the new one
        BELFEM_ERROR( mHeatHandle == nullptr,
            "Material already loaded a heating plugin; cannot load %s on top of it",
            aLibraryPath.c_str() );

        const string tPath = material::data_file( aLibraryPath );

        mHeatHandle = dlopen( tPath.c_str(), RTLD_LAZY );

        BELFEM_ERROR( mHeatHandle != nullptr, "Could not load library %s ( error: %s )",
            tPath.c_str(), dlerror() );

        string tFunctionName = sprint( "%s_init", aLabel.c_str() );

        void (* tInitFunc )( Material * ) =
            reinterpret_cast< void(*)( Material * ) >( dlsym( mHeatHandle, tFunctionName.c_str() ) );

        // the resolved path, not the name as written: naming a file the
        // loader never opened is what made this error hard to act on
        BELFEM_ERROR( tInitFunc != nullptr, "Could not find function %s in library %s ( error: %s )",
             tFunctionName.c_str(), tPath.c_str(), dlerror() );

        tInitFunc( this );
    }

    real
    Material::evaluate_polynomial( const MaterialProperty  Property,
      const real T ) const
    {
        BELFEM_ERROR( false, "Material::evaluate_polynomial() is not implemented for %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::evaluate_derivative_of_polynomial(
    const MaterialProperty  Property,
    const real T ) const
    {
        BELFEM_ERROR( false, "Material::evaluate_derivative_of_polynomial() is not implemented for %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

    void
    Material::set_user_defined_polynomial(
        const MaterialProperty Property,
        const Cell< real > & Coefficients )
    {
        BELFEM_ERROR( false, "Material::set_user_defined_polynomial() is not implemented for %s", mLabel.c_str() );
    }

    void
    Material::set_user_defined_polynomial(
        const MaterialProperty Property,
        const std::vector< real > & Coefficients )
    {
        BELFEM_ERROR( false, "Material::set_user_defined_polynomial() is not implemented for %s", mLabel.c_str() );
    }

    void
    Material::load_bh_curve( const material::BhCurve * aCurve )
    {
        this->set_bh_curve( aCurve );

        // reset constant value
        this->set_have( MaterialProperty::mu );
        this->reset_dependencies( MaterialProperty::mu );
        this->set_dependency( MaterialProperty::mu, MaterialDependency::normH );

        mFunctionH     = & Material::H_bhcurve ;
        mFunctionMu    = & Material::mu_bhcurve ;
        mFunctionDMuDH = & Material::dmudH_bhcurve ;
    }

    real
    Material::H_bhcurve( const real B ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::mu_bhcurve( real H, const real T ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
        return BELFEM_QUIET_NAN ;
    }


    void
    Material::dmudH_bhcurve( const real H, real & mu, real & dmudH ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
    }

    void
    Material::set_table_flags( const bool aFlag )
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
    }

    real
    Material::spline_property( const MaterialProperty aProperty, const real aX ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::dspline_property( const MaterialProperty aProperty, const real aX ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
        return BELFEM_QUIET_NAN ;
    }

    real
    Material::ddspline_property( const MaterialProperty aProperty, const real aX ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
        return BELFEM_QUIET_NAN ;
    }

    void
    Material::reset_spline( const MaterialProperty aProperty )
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
    }

    void
    Material::create_spline( real (Material::*aFunction)(const real aT) const,
        const MaterialProperty aProperty,
        const uint aStartBC,
        const uint aEndBC,
        const real adYdX0,
        const real adYdX1 )
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
    }

    real
    Material::l( const real T ) const
    {
        BELFEM_ERROR( false, "Not implemented for main material class" );
        return BELFEM_QUIET_NAN ;
    }

    void
    Material::create_spline(  const MaterialProperty aProperty, const real adYdX0, const real adYdX1 )
    {
        BELFEM_ASSERT( this->have( aProperty ), "Material %s does not have property %s", mLabel.c_str(), to_string(aProperty).c_str() );
        BELFEM_ASSERT( this->have( MaterialProperty::T_max ), "Material %s does not have a maximum temperature set (needed for %s).", mLabel.c_str(), to_string(aProperty).c_str() );

        real (Material::*tFunction)(const real aT) const = nullptr;

        spline::SplineBC tStartBC = std::isnan( adYdX0 ) ? spline::SplineBC::Parabolic : spline::SplineBC::Tangent ;
        spline::SplineBC tEndBC =  std::isnan( adYdX1 ) ?  spline::SplineBC::Parabolic : spline::SplineBC::Tangent ;

        this->reset_spline( aProperty );

        switch ( aProperty )
        {
            case MaterialProperty::density :
            {
                tFunction = & Material::density_custom ;
                break ;
            }
            case MaterialProperty::cp :
            {
                tFunction = & Material::cp_custom ;
                break ;
            }
            case MaterialProperty::lambda :
            {
                tFunction = & Material::lambda_custom ;
                break ;
            }
            case MaterialProperty::rho :
            {
                tFunction = & Material::rho_custom ;
                break ;
            }
            case MaterialProperty::mu :
            {
                BELFEM_ERROR( false, "Material::create_spline must not be called for mu-property");
                break;
            }
            case MaterialProperty::E :
            {
                tStartBC = spline::SplineBC::Tangent ;
                tFunction = & Material::E_custom ;
                break ;
            }
            case MaterialProperty::nu :
            {
                tStartBC = spline::SplineBC::Tangent ;
                tFunction = & Material::nu_custom ;
                break ;
            }
            case MaterialProperty::alpha :
            {
                tFunction = & Material::alpha_custom ;
                break;
            }
            case MaterialProperty::Rp02 :
            {
                tFunction = & Material::Rp02_custom ;
                break;
            }
            case MaterialProperty::rho_i :
            {
                tFunction = & Material::rho_i_custom ;
                break;
            }
            case MaterialProperty::debye :
            {
                tFunction = & Material::debye_custom ;
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Wrong material property ");
            }
        }

        this->create_spline( tFunction,
            aProperty,
            static_cast< uint > ( tStartBC ),
            static_cast< uint > ( tEndBC  ),
            adYdX0,
            adYdX1 );
    }

//------------------------------------------------------------------------------
}
