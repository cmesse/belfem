//
// Created by christian on 9/21/21.
//

#ifndef BELFEM_CL_MAXWELLMATERIAL_HPP
#define BELFEM_CL_MAXWELLMATERIAL_HPP

#include "cl_IsotropicMaterial.hpp"
#include "en_FEM_DomainType.hpp"
#include "cl_Spline.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    class MaxwellMaterial : public IsotropicMaterial
    {

        // spline for the permeability
        Spline * mNuSpline = nullptr ;

        // spline if resisivity is only temperature dependent
        Spline * mRhoSpline = nullptr ;

        // cell with labels of target groups
        Cell< string > mDomainLabels ;

        // constant value for mu_r
        real mMuR = 1.0 ;

        // constant value for nu_s
        real mNuS = constant::nu0 ;

        // critical resistance in Ohm * m
        real mRhoc  = BELFEM_QUIET_NAN ;

        // critical electric field in V/m
        real mEc    = BELFEM_QUIET_NAN ;

        // critical current density in A/m^2
        real mJc    = BELFEM_QUIET_NAN ;

        // critical temperature in K
        real mTc    = BELFEM_QUIET_NAN ;

        // operating temperature in K
        real mT0    = BELFEM_QUIET_NAN ;

        // critical current density at operating temperature A/m^2
        real mJc0   = BELFEM_QUIET_NAN ;

        // help parameter for mJc = alpha * T + beta
        real mAlpha = BELFEM_QUIET_NAN ;

        // help parameter for mJc = alpha * T + beta
        real mBeta  = BELFEM_QUIET_NAN ;

        // flux creep exponent minus one
        real mNm1   = BELFEM_QUIET_NAN ;

        // flux creep exponent minus one at operating temperature A/m^2 with b=0
        real mN0m1  = BELFEM_QUIET_NAN ;

        // flux creep exponent minus one at operating temperature A/m^2 with b>>0
        real mN1m1  = BELFEM_QUIET_NAN ;

        // constant used for magnetic field dependencies
        real mB0    = BELFEM_QUIET_NAN ;

        // constant for magnetizatoin offset in BH curve
        real mM0 = BELFEM_QUIET_NAN ;

        real
        ( MaxwellMaterial::*mFunMur )
        ( const real aH, const real aT  ) const;

        real
        ( MaxwellMaterial::*mFunNus )
        ( const real aB, const real aT  ) const;

        real
        ( MaxwellMaterial::*mFunRho )
        ( const real aJ, const real aT, const real aB ) const;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        MaxwellMaterial( const string aLabel );

        ~MaxwellMaterial();

        // magnetic permeability
        real
        mu_r ( const real aH=0, const real aT=BELFEM_TREF ) const ;

        // magnetic resistance
        real
        nu_s ( const real aB=0, const real aT=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        // electric resistance
        real
        rho_el ( const real aJ=0, const real aT=BELFEM_TREF, const real aB=0 ) const ;

//----------------------------------------------------------------------------

        void
        set_rho_el_const( const real aRhoEl );

//----------------------------------------------------------------------------

        void
        set_rho_el_ej( const real aEc, const real aJc, const real aN );

//----------------------------------------------------------------------------

        void
        set_rho_el_ejb( const real aEc, const real aJc0, const real aN0, const real aN1, const real aB0 );

//----------------------------------------------------------------------------

        void
        set_rho_el_ejt( const real aEc, const real aJc0, const real aN0, const real aT0, const real aTc );

//----------------------------------------------------------------------------

        void
        set_rho_el_ejbt( const real aEc, const real aJc0, const real aN0, const real aN1, const real aB0, const real aT0, const real aTc  );



//----------------------------------------------------------------------------

        void
        set_mu_r( const real aMuR );

//----------------------------------------------------------------------------

        void
        set_nu_s( Spline * aSpline );

//----------------------------------------------------------------------------

        void
        set_nu_s( const real aNu );

//----------------------------------------------------------------------------

        void
        set_rho_el( Spline * aSpline );

//----------------------------------------------------------------------------

        void
        set_m0( const real & aM0 );

        real
        m0() const ;

//----------------------------------------------------------------------------

        /**
         * list of names this material is connected to
         */
        const Cell< string > &
        domain_labels() const ;

//----------------------------------------------------------------------------

        /**
         * return the critical current density
         */
         real
         jc() const;

//----------------------------------------------------------------------------
    private:
//----------------------------------------------------------------------------

        void
        reset_parameters();

//----------------------------------------------------------------------------

        real
        mu_r_const ( const real aH=0, const real aT=BELFEM_TREF ) const ;

        real
        nu_s_const ( const real aB=0, const real aT=BELFEM_TREF ) const ;

        real
        nu_s_spline ( const real aB=0, const real aT=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        real
        rho_el_const ( const real aJ, const real aT=BELFEM_TREF, const real aB=0 ) const ;

        real
        rho_el_powerlaw_ej( const real aJ, const real aT=BELFEM_TREF, const real aB=0  ) const ;

        real
        rho_el_powerlaw_ejt(const real aJ, const real aT=BELFEM_TREF, const real aB=0 ) const ;

        real
        rho_el_powerlaw_ejb( const real aJ, const real aT=BELFEM_TREF, const real aB=0  ) const ;

        real
        rho_el_powerlaw_ejbt( const real aJ, const real aT=BELFEM_TREF, const real aB=0  ) const ;

        real
        rho_el_spline_t( const real aJ, const real aT=BELFEM_TREF, const real aB=0  ) const ;

//---------------------------------------------------------------------------
    };

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el( const real aJ, const real aT, const real aB ) const
        {
            return ( this->*mFunRho )( aJ, aT, aB );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::mu_r( const real aH, const real aT ) const
        {
            return ( this->*mFunMur )( aH, aT );
        }
//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::nu_s ( const real aB, const real aT ) const
        {
            return ( this->*mFunNus )( aB, aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_const( const real aJ, const real aT, const real aB ) const
        {
            return mRhoc ;
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_powerlaw_ej( const real aJ, const real aT, const real aB ) const
        {
            return mRhoc * std::pow( aJ / ( mJc + BELFEM_EPS ), mNm1 );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_powerlaw_ejt( const real aJ, const real aT, const real aB ) const
        {
            // compute critical current
            real tJc = mAlpha * aT + mBeta ;

            // compute flux creep exponent
            return ( mEc / ( tJc + BELFEM_EPS ) ) * std::pow( aJ / tJc, mN0m1 * mT0 / aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_powerlaw_ejb( const real aJ, const real aT, const real aB  ) const
        {
            // compute critical current
            real tJc = mJc0 / ( 1. + aB / mB0 );

            return ( mEc / ( tJc + BELFEM_EPS ) ) * std::pow( aJ / tJc,
                                             mN1m1 + ( mN0m1-mN1m1 ) / ( 1. +  aB / mB0 ) );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_powerlaw_ejbt( const real aJ, const real aT, const real aB  ) const
        {
            // compute critical current
            real tJc = ( mAlpha * aT + mBeta ) / ( 1. + aB / mB0 );

            return ( mEc / ( tJc + BELFEM_EPS ) ) * std::pow( aJ / tJc,
                                             ( mN1m1 + ( mN0m1-mN1m1 ) / ( 1. +  aB / mB0 ) )
                                             * ( mT0 / aT ) );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_spline_t( const real aJ, const real aT, const real aB ) const
        {
            return mRhoSpline->eval( aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::mu_r_const( const real aH, const real aT ) const
        {
            BELFEM_ASSERT( mPermeabilityLaw == PermeabilityLaw::Constant,
                "Nope. You are doing it wrong! Never all mu_r() this way if it is not constant.");

            return mMuR ;
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::nu_s_const ( const real aB, const real aT ) const
        {
            return mNuS ;
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::nu_s_spline ( const real aB, const real aT ) const
        {
            return std::exp( mNuSpline->eval( aB ) );
        }

//----------------------------------------------------------------------------

        inline const Cell< string > &
        MaxwellMaterial::domain_labels() const
        {
            return mDomainLabels ;
        }

//----------------------------------------------------------------------------

        inline void
        MaxwellMaterial::set_m0( const real & aM0 )
        {
            mM0 = aM0 ;
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::m0() const
        {
            return mM0 ;
        }

//----------------------------------------------------------------------------

        /**
         * return the critical current density
         */
        inline real
        MaxwellMaterial::jc() const
        {
            BELFEM_ASSERT( mResistivityLaw == ResistivityLaw::Constant
                || mResistivityLaw == ResistivityLaw::DependJ,
                "material law bust be either constant or e-j power to return jc");

            return mJc ;

        }

//----------------------------------------------------------------------------
}
#endif //BELFEM_CL_FEM_MAXWELLMATERIAL_HPP
