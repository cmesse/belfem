//
// Created by christian on 9/21/21.
//

#ifndef BELFEM_CL_MAXWELLMATERIAL_HPP
#define BELFEM_CL_MAXWELLMATERIAL_HPP

#include "cl_IsotropicMaterial.hpp"
#include "en_FEM_DomainType.hpp"
#include "cl_Spline.hpp"
#include "cl_Database.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    class MaxwellMaterial : public IsotropicMaterial
    {
        // spline for the permeability
        Spline * mNuSpline = nullptr ;

        // spline if resisivity is only temperature dependent
        Spline * mRhoSpline = nullptr ;

        // for jc database
        Database * mJcData = nullptr ;

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
        real mBmax  = BELFEM_QUIET_NAN ;
        real mHmax  = BELFEM_QUIET_NAN ;

        // constant for magnetizatoin offset in BH curve
        real mM0 = BELFEM_QUIET_NAN ;

        // maximum resistance for powerlaw
        real mRhoMax = 1e-10 ;

        // help material for thermal data
        Material * mThermalMaterial = nullptr ;

        real
        ( MaxwellMaterial::*mFunMur )
        ( const real aH, const real aT  ) const;

        real
        ( MaxwellMaterial::*mFunNus )
        ( const real aB, const real aT  ) const;

        real
        ( MaxwellMaterial::*mFunMus )
                ( const real aB, const real aT  ) const;

        real
        ( MaxwellMaterial::*mFunRho )
        ( const real aJ, const real aT, const real aB, const real aAngle ) const;

        real
        ( MaxwellMaterial::*mFunRhoJcrit )
                ( const real aJ, const real aJc, const real aT, const real aB, const real aAngle ) const;

        real
        ( MaxwellMaterial::*mFunJcrit )
                (  const real aT, const real aB, const real aAngle ) const;

        real
        ( MaxwellMaterial::*mFunMm1 )
                (  const real aT, const real aB  ) const;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        MaxwellMaterial( const string aLabel );

        ~MaxwellMaterial();

        // magnetic permeability
        real
        mu_r ( const real aH=0, const real aT=BELFEM_TREF ) const ;

        // magnetic permeability
        real
        mu_s ( const real aH=0, const real aT=BELFEM_TREF ) const ;

        // magnetic resistance
        real
        nu_s ( const real aB=0, const real aT=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        // electric resistance
        real
        rho_el ( const real aJ=0, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0 ) const ;

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
        set_j_crit( const string & aDatabase, const string & aDataset );

//----------------------------------------------------------------------------

        void
        set_mu_r( const real aMuR );

//----------------------------------------------------------------------------

        void
        set_nu_s( Spline * aNuSpline );

//----------------------------------------------------------------------------

        void
        set_nu_s( const real aNu );

//----------------------------------------------------------------------------

        void
        set_rho_el( Spline * aSpline );

//----------------------------------------------------------------------------

        void
        set_m0( const real & aM0 );

//----------------------------------------------------------------------------

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
         j_crit( const real aT=0, const real aB=0, const real aAlpha=0 ) const;

//----------------------------------------------------------------------------

        /**
          * Density in kg/m^3
          */
        real
        rho( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * specific heat capacity in J/(kg*K)
         */
        real
        c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * thermal conductivity in W/(m*K)
         */
        real
        lambda( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * thermal conductivity in W/(m*K)
         */
        real
        lambda( const real aT, const real aB ) const;

//----------------------------------------------------------------------------

        void
        set_thermal_material( const MaterialType aType );

//----------------------------------------------------------------------------

        void
        set_thermal_material( const string aLabel );

//----------------------------------------------------------------------------

        real
        creep_expinent_minus_1() const ;

//----------------------------------------------------------------------------

        void
        compute_jcrit_and_rho( real & aJcrit,
                                                real & aRho,
                                                const real aJ=0,
                                                const real aT=0,
                                                const real aB=0,
                                                const real aAngle=0 ) const ;

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

        real
        mu_s_spline ( const real aH=0, const real aT=BELFEM_TREF ) const ;

        real
        j_crit_const( const real aT, const real aB, const real aAngle ) const ;

        real
        j_crit_t( const real aT, const real aB, const real aAngle ) const ;

        real
        j_crit_b( const real aT, const real aB, const real aAngle ) const ;

        real
        j_crit_bt( const real aT, const real aB, const real aAngle ) const ;

        real
        j_crit_database( const real aT, const real aB, const real aAngle ) const ;

        real
        nm1_const( const real aT, const real aB ) const ;

        real
        nm1_t( const real aT, const real aB ) const ;

        real
        nm1_b( const real aT, const real aB ) const ;

        real
        nm1_bt( const real aT, const real aB ) const ;

//----------------------------------------------------------------------------

        real
        rho_el_const ( const real aJ, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0 ) const ;

        real
        rho_el_powerlaw( const real aJ, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0   ) const ;

        real
        rho_el_spline_t( const real aJ, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0   ) const ;
//----------------------------------------------------------------------------

        real
        rho_el_jc_const ( const real aJ, const real aJc, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0 ) const ;

        real
        rho_el_jc_powerlaw( const real aJ, const real aJc, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0) const ;

        real
        rho_el_jc_spline_t( const real aJ, const real aJc, const real aT=BELFEM_TREF, const real aB=0, const real aAngle=0 ) const ;

//---------------------------------------------------------------------------
    };

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el( const real aJ, const real aT, const real aB, const real aAngle ) const
        {
            return ( this->*mFunRho )( std::abs( aJ ), aT, std::abs( aB ), aAngle );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::mu_r( const real aH, const real aT ) const
        {
            return ( this->*mFunMur )( std::abs( aH ), aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::nu_s ( const real aB, const real aT ) const
        {
            return ( this->*mFunNus )( std::abs( aB ), aT );
        }
//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::mu_s ( const real aH, const real aT ) const
        {
            return ( this->*mFunMus )( std::abs( aH ), aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_const( const real aJ, const real aT, const real aB, const real aAngle ) const
        {
            return mRhoc ;
        }

//----------------------------------------------------------------------------

        inline void
        MaxwellMaterial::compute_jcrit_and_rho( real & aJcrit,
                               real & aRho,
                               const real aJ,
                               const real aT,
                               const real aB,
                               const real aAngle ) const
        {
            aJcrit =  ( this->*mFunJcrit )( aT, aB, aAngle );
            aRho = ( this->*mFunRhoJcrit )( aJ, aJcrit, aT, aB, aAngle );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_powerlaw( const real aJ, const real aT, const real aB, const real aAngle ) const
        {
            real tJc = ( this->*mFunJcrit )( aT, aB, aAngle );

            return  std::min( ( mEc / tJc ) * std::pow( aJ / ( tJc + BELFEM_EPSILON ), ( this->*mFunMm1 )( aT, aB ) ) , mRhoMax );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_jc_const ( const real aJ, const real aJc, const real aT, const real aB, const real aAngle ) const
        {
            return mRhoc ;
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_jc_powerlaw( const real aJ, const real aJc, const real aT, const real aB, const real aAngle ) const
        {
            return  std::min( ( mEc / aJc ) * std::pow( aJ / ( aJc + BELFEM_EPSILON ), ( this->*mFunMm1 )( aT, aB ) ) , mRhoMax );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_jc_spline_t( const real aJ, const real aJc, const real aT, const real aB, const real aAngle ) const
        {
            return  mRhoSpline->eval( aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho_el_spline_t( const real aJ, const real aT, const real aB, const real aAngle ) const
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
            if( aB < mBmax )
            {
                return std::exp( mNuSpline->eval( aB ));
            }
            else
            {
                return ( mHmax + constant::nu0 * ( aB - mBmax ) ) / aB ;
            }
        }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::mu_s_spline ( const real aH, const real aT ) const
    {
        if( aH < BELFEM_EPSILON )
        {
            return  0 ;
        }
        else if( aH < mHmax )
        {

            real tF = 1 ;
            uint tCount = 0 ;


            real tB ;
            real tB0 = constant::mu0 * aH ;
            real tF0 = this->nu_s( tB0 ) * tB0 - aH ;
            real tB1 = mBmax ;


            while( std::abs( tF ) > aH *1e-12 )
            {
                tB = 0.5 * ( tB1 + tB0 );

                tF = this->nu_s( tB ) * tB - aH ;

                if ( tF0 * tF  > 0 )
                {
                    tB0 = tB;
                    tF0 = tF;
                }
                else
                {
                    tB1 = tB;
                }

                BELFEM_ERROR( tCount++ < 200, "Too many iterations.");
            }
            return tB / aH ;
        }
        else
        {
            return ( mBmax + constant::mu0 * ( aH - mHmax ) ) / aH ;
        }
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
        MaxwellMaterial::j_crit( const real aT, const real aB, const real aAngle ) const
        {
            return ( this->*mFunJcrit )( aT, aB, aAngle );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::rho( const real aT ) const
        {
            BELFEM_ASSERT( mThermalMaterial != nullptr,
                "thermal material not set for %s", mLabel.c_str() );

            return mThermalMaterial->rho( aT );
        }

//---------------------------------------------------------------------------

        inline real
        MaxwellMaterial::c( const real aT ) const
        {
            BELFEM_ASSERT( mThermalMaterial != nullptr,
                           "thermal material not set for %s", mLabel.c_str() );

            return mThermalMaterial->c( aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::lambda( const real aT ) const
        {
            BELFEM_ASSERT( mThermalMaterial != nullptr,
                           "thermal material not set for %s", mLabel.c_str() );

            return mThermalMaterial->lambda( aT );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::lambda( const real aT, const real aB ) const
        {
            BELFEM_ASSERT( mThermalMaterial != nullptr,
                           "thermal material not set for %s", mLabel.c_str() );

            return mThermalMaterial->lambda( aT, aB );
        }

//----------------------------------------------------------------------------

        inline real
        MaxwellMaterial::creep_expinent_minus_1() const
        {
            return mNm1 ;
        }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::j_crit_const( const real aT, const real aB, const real aAngle ) const
    {
        return mJc ;
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::j_crit_t( const real aT, const real aB, const real aAngle ) const
    {
        return mAlpha * aT + mBeta ;
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::j_crit_b( const real aT, const real aB, const real aAngle ) const
    {
        return mJc0 / ( 1. + aB / mB0 );
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::j_crit_bt( const real aT, const real aB, const real aAngle ) const
    {
        return ( mAlpha * aT + mBeta ) / ( 1. + aB / mB0 );
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::j_crit_database( const real aT, const real aB, const real aAngle ) const
    {
        return mJcData->compute( aAngle, aT, aB );
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::nm1_const( const real aT, const real aB ) const
    {
        return mNm1 ;
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::nm1_t( const real aT, const real aB ) const
    {
        return mN0m1 * mT0 / aT ;
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::nm1_b( const real aT, const real aB ) const
    {
        return mN1m1 + ( mN0m1-mN1m1 ) / ( 1. +  aB / mB0 ) ;
    }

//----------------------------------------------------------------------------

    inline real
    MaxwellMaterial::nm1_bt( const real aT, const real aB ) const
    {
        return ( mN1m1 + ( mN0m1-mN1m1 ) / ( 1. +  aB / mB0 ) ) * ( mT0 / aT ) ;
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_CL_FEM_MAXWELLMATERIAL_HPP
