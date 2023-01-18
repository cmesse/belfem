//
// Created by Christian Messe on 13.11.19.
//

#include "assert.hpp"
#include "cl_Material.hpp"

namespace belfem
{
//----------------------------------------------------------------------------


    Material::Material( const MaterialType aType ) :
        mType( aType )
    {
        mLabel = to_string( aType );
    }

//----------------------------------------------------------------------------

    const string &
    Material::label() const
    {
        return mLabel;
    }

//----------------------------------------------------------------------------

    const string &
    Material::number() const
    {
        return mNumber ;
    }

//----------------------------------------------------------------------------
//      ELASTIC PROPERTIES
//----------------------------------------------------------------------------

    real
    Material::E( const real aT ) const
    {
        BELFEM_ERROR( false, "E is not implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    real
    Material::nu( const real aT ) const
    {
        BELFEM_ERROR( false, "nu is not implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    real
    Material::G( const real aT ) const
    {
        BELFEM_ERROR( false, "G is not implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    void
    Material::set_rrr( const real aRRR )
    {
        BELFEM_ERROR( false, "set_rrr is not implemented for material %s",
                      mLabel.c_str() );
    }

//----------------------------------------------------------------------------

    void
    Material::use_splines( const bool aSwitch )
    {
        BELFEM_ERROR( false, "use_splines is not implemented for material %s",
                      mLabel.c_str() );
    }

//----------------------------------------------------------------------------


    void
    Material::C(  Matrix< real > & aC, const real aT ) const
    {
        BELFEM_ERROR( false, "C is not implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------

    void
    Material::C_ps( Matrix< real > & aC, const real aT ) const
    {
        BELFEM_ERROR( false, "C_ps is not implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------

    void
    Material::C_rot( Matrix< real > & aC, const real aT  ) const
    {
        BELFEM_ERROR( false, "C_rot is not implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------

    real
    Material::rho( const real aT ) const
    {
        BELFEM_ERROR( false, "rho is implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------
//   THERMAL PROPERTIES
//----------------------------------------------------------------------------

    /**
     * thermal conductivity in W/(m*K)
     */
    real
    Material::lambda( const real aT ) const
    {
        BELFEM_ERROR( false, "lambda is implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    void
    Material::lambda_p( Matrix< real > & aLambda, const real aT ) const
    {
        BELFEM_ERROR( false, "lambda_p is implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------

    void
    Material::lambda_3d( Matrix< real > & aLambda, const real aT ) const
    {
        BELFEM_ERROR( false, "lambda is implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------


    real
    Material::c( const real aT ) const
    {
        BELFEM_ERROR( false, "c is not implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------
// Thermal Expansion
//----------------------------------------------------------------------------

    real
    Material::alpha( const real aT ) const
    {
        BELFEM_ERROR( false, "alpha is not implemented for material %s",
                 mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    real
    Material::mu( const real aT, const real aTref ) const
    {
        BELFEM_ERROR( false, "mu is not implemented for material %s",
                          mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    void
    Material::mu_ps( Matrix< real > & aMu, const real aT, const real aTref ) const
    {
        BELFEM_ERROR( false, "mu_ps is not implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------

    void
    Material::mu( Matrix< real > & aMu, const real aT, const real aTref  ) const
    {
        BELFEM_ERROR( false, "mu is not implemented for material %s",
                     mLabel.c_str() );
    }

//----------------------------------------------------------------------------
// Optical Properties
//----------------------------------------------------------------------------

    real
    Material::epsilon( const real aT ) const
    {
        BELFEM_ERROR( false, "epsilon is not implemented for material %s",
                     mLabel.c_str() );
        return 0.0;
    }

//----------------------------------------------------------------------------

    bool
    Material::has_thermal() const
    {
        BELFEM_ERROR( false, "has_thermal() not implemented for material %s", mLabel.c_str() );
        return false ;
    }

//----------------------------------------------------------------------------

    bool
    Material::has_mechanical() const
    {
        BELFEM_ERROR( false, "has_mechanical() not implemented for material %s", mLabel.c_str() );
        return false ;
    }

//----------------------------------------------------------------------------

    bool
    Material::has_thermal_expansion() const
    {
        BELFEM_ERROR( false, "has_thermal_expansion() not implemented for material %s", mLabel.c_str() );
        return false ;
    }

//----------------------------------------------------------------------------

    bool
    Material::has_electric_resistivity() const
    {
        BELFEM_ERROR( false, "has_thermal_expansion() not implemented for material %s", mLabel.c_str() );
        return false ;
    }


//----------------------------------------------------------------------------

    real
    Material::mu_r ( const real aH, const real aT ) const
    {
        BELFEM_ERROR( false, "mu_r() not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    // electric resistance
    real
    Material::j_crit( const real aT, const real aB, const real aAngle ) const
    {
        return  BELFEM_REAL_MAX ;
    }

//----------------------------------------------------------------------------

    void
    Material::compute_jcrit_and_rho(
                           real & aJcrit,
                           real & aRho,
                           const real aJ,
                           const real aT,
                           const real aB,
                           const real aAngle ) const
    {
        aJcrit = BELFEM_REAL_MAX ;
        aRho = this->rho_el( aJ, aT, aB, aAngle );
    }

//----------------------------------------------------------------------------

    real
    Material::creep_expinent_minus_1() const
    {
        return  0 ;
    }

//----------------------------------------------------------------------------

    // electric resistance
    real
    Material::rho_el ( const real aJ, const real aT, const real aB, const real aAngle ) const
    {
        BELFEM_ERROR( false, "rho_el() not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    Material::nu_s( const real aB, const real aT ) const
    {
        BELFEM_ERROR( false, "nu_s() not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    Material::mu_s( const real aB, const real aT ) const
    {
        BELFEM_ERROR( false, "mu_s() not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//-----------------------------------------------------------------------
}