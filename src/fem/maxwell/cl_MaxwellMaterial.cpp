//
// Created by christian on 9/21/21.
//
#include "stringtools.hpp"
#include "cl_MaterialFactory.hpp"
#include "cl_MaxwellMaterial.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    MaxwellMaterial::MaxwellMaterial( const string aLabel ) :
        IsotropicMaterial( MaterialType::Maxwell )
    {
        // set name
        mLabel = aLabel ;

        // temporary cell with label
        Cell< string > tWords = string_to_words( aLabel );

        if( tWords.size() == 1 )
        {
            // if there is only one word, this points to the one group with that name
            mDomainLabels.push( tWords( 0 ) );
        }
        else
        {
            // otherwise, all but the first words are targets
            for( index_t k=1; k<tWords.size(); ++k )
            {
                mDomainLabels.push( tWords( k ) );
            }
        }
        // set default values
        this->set_rho_el_const( BELFEM_REAL_MAX );

        this->set_mu_r( 1.0 );
    }

//----------------------------------------------------------------------------

    MaxwellMaterial::~MaxwellMaterial()
    {
        if( mNuSpline != nullptr )
        {
            delete mNuSpline ;
        }
        if( mMuSpline != nullptr )
        {
            delete mMuSpline ;
        }
        if( mRhoSpline != nullptr )
        {
            delete mRhoSpline ;
        }
        if( mThermalMaterial != nullptr )
        {
            delete mThermalMaterial ;
        }
        if( mJcData != nullptr )
        {
            delete mJcData ;
        }
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_rho_el_const( const real aRhoEl )
    {
        mResistivityLaw = ResistivityLaw::Constant ;

        mFunRho = & MaxwellMaterial::rho_el_const ;
        this->reset_parameters();
        mRhoc = aRhoEl ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_rho_el_ej( const real aEc, const real aJc, const real aN )
    {
        mResistivityLaw = ResistivityLaw::DependJ ;

        this->reset_parameters();

        mFunRho   = & MaxwellMaterial::rho_el_powerlaw ;
        mFunJcrit = & MaxwellMaterial::j_crit_const ;
        mFunMm1   = & MaxwellMaterial::nm1_const ;

        mRhoc  = aEc / aJc ;
        mEc    = aEc ;
        mJc    = aJc ;
        mNm1   = aN - 1.;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_rho_el_ejb(
            const real aEc,
            const real aJc0,
            const real aN0,
            const real aN1,
            const real aB0 )
    {
        mResistivityLaw = ResistivityLaw::DependJB ;

        mFunRho   = & MaxwellMaterial::rho_el_powerlaw ;
        mFunJcrit = & MaxwellMaterial::j_crit_b ;
        mFunMm1   = & MaxwellMaterial::nm1_b ;

        this->reset_parameters();

        mEc   = aEc ;
        mJc0  = aJc0 ;
        mN0m1 = aN0 - 1.;
        mN1m1 = aN1 - 1.;
        mB0   = aB0 ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_rho_el_ejt(
            const real aEc,
            const real aJc0,
            const real aN0,
            const real aT0,
            const real aTc )
    {
        mResistivityLaw = ResistivityLaw::DependJT ;
        mFunRho   = & MaxwellMaterial::rho_el_powerlaw ;
        mFunJcrit = & MaxwellMaterial::j_crit_t ;
        mFunMm1   = & MaxwellMaterial::nm1_t ;

        this->reset_parameters();

        mEc   = aEc ;
        mJc0  = aJc0 ;
        mN0m1 = aN0 - 1.;
        mT0   = aT0 ;
        mTc   = aTc ;

        mBeta =  aJc0 / ( 1. + aT0 / aTc );
        mAlpha  = -mBeta / aTc ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_rho_el_ejbt(
            const real aEc,
            const real aJc0,
            const real aN0,
            const real aN1,
            const real aB0,
            const real aT0,
            const real aTc )
    {
        mResistivityLaw = ResistivityLaw::DependJBT ;
        mFunRho   = & MaxwellMaterial::rho_el_powerlaw ;
        mFunJcrit = & MaxwellMaterial::j_crit_bt ;
        mFunMm1   = & MaxwellMaterial::nm1_bt ;

        this->reset_parameters();

        mEc   = aEc ;
        mJc0  = aJc0 ;
        mN0m1 = aN0 - 1.;
        mN1m1 = aN1 - 1.;
        mB0   = aB0 ;
        mT0   = aT0 ;
        mTc   = aTc ;

        mBeta =  aJc0 / ( 1. + aT0 / aTc );
        mAlpha  = -mBeta / aTc ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_j_crit( const string & aDatabase, const string & aDataset )
    {
        if( mJcData != nullptr )
        {
            delete mJcData ;
        }
        mJcData = new Database( aDatabase, aDataset );
        mFunJcrit = & MaxwellMaterial::j_crit_database ;
        mResistivityLaw = ResistivityLaw::DependJBT ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_mu_r( const real aMuR )
    {
        mPermeabilityLaw = PermeabilityLaw::Constant ;
        mFunMur = & MaxwellMaterial::mu_r_const ;
        mFunNus = & MaxwellMaterial::nu_s_const ;
        mMuR = aMuR ;
        mNuS = constant::nu0 / mMuR ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_nu_s( Spline * aMuSpline, Spline * aNuSpline )
    {
        mMuSpline = aMuSpline ;
        mNuSpline = aNuSpline ;

        mPermeabilityLaw = PermeabilityLaw::Spline ;
        mFunNus = & MaxwellMaterial::nu_s_spline ;
        mFunMus = & MaxwellMaterial::mu_s_spline ;
        mBmax = aMuSpline->x_max() ;
        mHmax = aMuSpline->eval( mBmax ) * mBmax ;

        // we don't want the mu_r to return a reasonable value
        // don't worry, mu_r should never be called in this mode
        // but if it does,
        mFunMur = nullptr ;
        mMuR = BELFEM_QUIET_NAN ;
        mNuS = BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_rho_el( Spline * aSpline )
    {
        mRhoSpline = aSpline ;
        mResistivityLaw = ResistivityLaw::DependT ;
        mFunRho = & MaxwellMaterial::rho_el_spline_t ;

        mFunJcrit = & MaxwellMaterial::j_crit_const ;
        mJc = 0 ;

        mPermeabilityLaw = PermeabilityLaw::Constant ;
        mFunMur = & MaxwellMaterial::mu_r_const ;
        mMuR = std::isnan( mMuR ) ? 1.0 : mMuR ;
        mNuS = constant::nu0 / mMuR ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_nu_s( const real aNu )
    {
        mPermeabilityLaw = PermeabilityLaw::Constant ;
        mFunNus = & MaxwellMaterial::nu_s_const ;
        mFunMur = & MaxwellMaterial::mu_r_const ;
        mNuS = aNu ;
        mMuR = constant::nu0 / aNu ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::reset_parameters()
    {
        // critical resistance in Ohm * m
        mRhoc  = BELFEM_QUIET_NAN ;

        // critical electric field in V/m
        mEc    = BELFEM_QUIET_NAN ;

        // critical current density in A/m^2
        mJc    = BELFEM_QUIET_NAN ;

        // critical temperature in K
        mTc    = BELFEM_QUIET_NAN ;

        // operating temperature in K
        mT0    = BELFEM_QUIET_NAN ;

        // critical current density at operating temperature A/m^2
        mJc0   = BELFEM_QUIET_NAN ;

        // help parameter for mJc = alpha * T + beta
        mAlpha = BELFEM_QUIET_NAN ;

        // help parameter for mJc = alpha * T + beta
        mBeta  = BELFEM_QUIET_NAN ;

        // flux creep exponent minus one
        mNm1   = BELFEM_QUIET_NAN ;

        // flux creep exponent minus one at operating temperature A/m^2 with b=0
        mN0m1  = BELFEM_QUIET_NAN ;

        // flux creep exponent minus one at operating temperature A/m^2 with b>>0
        mN1m1  = BELFEM_QUIET_NAN ;

        // constant used for magnetic field dependencies
        mB0    = BELFEM_QUIET_NAN ;

        // constant reduced magnetic permeability
        mMuR   = 1.0 ;

        // secant magnetic resistance
        mNuS   = BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_thermal_material( const MaterialType aType )
    {
        MaterialFactory tFactory ;
        if( mThermalMaterial != nullptr )
        {
            delete mThermalMaterial;
        }
        mThermalMaterial = tFactory.create_material( aType );
    }
//----------------------------------------------------------------------------

    void
    MaxwellMaterial::set_thermal_material( const string aLabel )
    {
        MaterialFactory tFactory ;
        if( mThermalMaterial != nullptr )
        {
            delete mThermalMaterial;
        }
        mThermalMaterial = tFactory.create_material( aLabel );
    }

//----------------------------------------------------------------------------
}