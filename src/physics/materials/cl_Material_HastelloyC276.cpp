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

#include "cl_Material_HastelloyC276.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "cl_Material_Abundance.hpp"
#include "debye.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * Hastelloy C-276 material property implementation
         *
         * This alloy exhibits complex temperature-dependent behavior requiring
         * multiple piecewise polynomial regions for accurate representation.
         *
         * Key features:
         * - Electrical resistivity: Low-T polynomial + log-log fit (10.1063/1.2899058)
         * - Thermal expansion: Polynomial with constant tail at high T
         * - Magnetic susceptibility: χ(T) modeled as χ(1/T) with smooth transitions,
         *   prepared in create_mu() but NOT routed — mu stays mu0 ( see the note there )
         * - Specific heat: 6 piecewise regions + Debye tail for high T
         * - Thermal conductivity: λ = λ_phonon + λ_electron (Wiedemann-Franz)
         *
         * Molar mass and impurity parameter computed from composition.
         */
        HastelloyC276::HastelloyC276() :
            Metal( "HastelloyC276", MaterialType::LookupAlloy )
        {
            this->set_constants();         // Composition, density, melting point

            // create_cp() samples cp_from_debye(), whose dilation term reads
            // alpha( T ), E( T ) and nu( T ): the polynomial alpha, the moduli
            // and rho come first, cp after them, and only then the cryogenic
            // alpha branch, which in turn is fitted against cp
            this->create_alpha();          // Thermal expansion coefficient ( polynomial )
            this->create_rho();            // Electrical resistivity
            this->create_mech();           // Young's modulus, Poisson's ratio
            this->create_mu();             // Magnetic susceptibility fit ( NOT routed, mu = mu0 )

            this->create_cp();             // Specific heat capacity
            this->create_alpha_cryo();     // Grueneisen branch of alpha below the split

            this->create_lambda();         // Thermal conductivity
        }

        HastelloyC276::~HastelloyC276()
        {

        }

        void HastelloyC276::set_constants()
        {
            this->set_constant( MaterialProperty::T_max, 1643.15 );
            this->set_constant( MaterialProperty::ref_density, 8890.0);
            this->set_constant( MaterialProperty::T_ref_density, gTroom );


            Abundance tAbundance ;
            real M ;
            real Gamma ;

            // composition see 10.1016/j.cryogenics.2023.103776
            tAbundance.compute_molar_mass_and_impurity_from_masses(
                { "Ni", "Mo", "Cr", "Fe", "W", "Co", "Mn", "Al" },
                { 0, 0.1418, 0.1401, 0.0476, 0.0293, 0.0371, 0.0053,  0.0026 },
                M,
                Gamma );

            this->set_constant( MaterialProperty::M, M );
            this->set_constant( MaterialProperty::Gamma, Gamma );
        }

        void
        HastelloyC276::create_mech()
        {
            // based on data from 10.1016/j.cryogenics.2006.01.014
            // completed with data from MATWEB
            mYoungPoly = {-3.90589006630972E04, 0.0, 2.04065547724810E+11 };

            // we know that nu @ RT = 0.307, at melting temperature, nu->0.5.
            // We back off slightly from 0.5 because K = E/(3·(1-2ν)) is sampled
            // by the cp spline up to T_max via cp_from_debye, and ν=0.5 makes
            // K diverge, poisoning the spline with NaN.
            real T0 = gTroom ;
            real T1 = this->constant_property( MaterialProperty::T_max );

            real nu0 = 0.307 ;
            real nu1 = 0.49 ;

            // moreover, assuming that dν/dT = 0 at 0K, we can obtain
            // ν(T) = a·T² + c with ν(T0)=ν0, ν(T1)=ν1, hence
            //   a = (ν1-ν0) / (T1²-T0²)
            //   c = (ν0·T1² - ν1·T0²) / (T1²-T0²)
            real a = ( nu1-nu0 ) / (T1*T1-T0*T0) ;
            real c = (nu0*T1*T1-nu1*T0*T0) / (T1*T1-T0*T0) ;
            mPoissonPoly = { a, 0.0, c };

            this->set_have( MaterialProperty::E );
            this->set_have( MaterialProperty::nu );

            this->set_custom( MaterialProperty::E );
            this->set_custom( MaterialProperty::nu );
        }

        void
        HastelloyC276::create_rho()
        {
            // fitted against curve from 10.1063/1.2899058
            mRhoPolys.set_size( 2 );

            // for data > 10.73 K, log-log curve
            mRhoPolys( 1 ) = {  1.93174404699903E-04, - 2.21920482428318E-03,+ 1.07884809281212E-02, - 2.40374094654130E-02,- 1.35906507676544E+01};

            // find temperature at which the curve turns
            real x = 2.373 ;

            real f = dpolyval( mRhoPolys( 1 ), x );
            real df = ddpolyval( mRhoPolys( 1 ), x );

            while ( abs( f ) > BELFEM_EPSILON )
            {
                x -= f/df;
                f = dpolyval( mRhoPolys( 1 ), x );
                df = ddpolyval( mRhoPolys( 1 ), x );
            }


            // now we need to build the curve for the low temperature
            real & T0 = mTRhoSwitch ;
            T0 = std::exp( x );
            f = std::exp( polyval( mRhoPolys( 1 ), x ));

            real a = 2.42719471085423E-13 ;
            real b = 0.0 ;
            real c = 3.69608432119349E-11 ;
            real d = 0.0 ;

            Vector< real > & p = mRhoPolys( 0 );
            p.set_size( 5 );

            p( 0 ) = a ;
            p( 1 ) = b-4.*T0*a;
            p( 2 ) = T0*(6.*T0*a-3.*b) + c ;
            p( 3 ) = T0 * ( T0*(3.*b- 4.*T0*a) - 2.*c ) + d ;
            p( 4 ) = T0*(T0*(T0*(a*T0-b)+c)-d) + f;

            this->set_have( MaterialProperty::rho );
            this->create_spline( MaterialProperty::rho );
        }

        void
        HastelloyC276::create_alpha()
        {
            mAlphaPolys.set_size( 2 );

            // based on data from 10.1063/1.2899058
            mAlphaPolys( 0 ) = {
                -7.26491387242953E-13,
                - 2.74765591990654E-12,
                + 1.12899447795587E-07,
                0.0
            };

            // the polynomial has a maximum, from which onwards we use a constant value
            real a = 3. * mAlphaPolys( 0 )( 0 ) ;
            real b = 2. * mAlphaPolys( 0 )( 1 ) ;
            real c =      mAlphaPolys( 0 )( 2 ) ;

            mTAlphaPlateau = ( -b - std::sqrt( b*b - 4.0*a*c ) ) / ( 2.0*a );

            real f = polyval( mAlphaPolys( 0 ), mTAlphaPlateau );
            mAlphaPolys( 1 ).set_size( 1, f );

            this->set_have( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::alpha );
            this->set_have( MaterialProperty::density );
            this->set_custom( MaterialProperty::density );

            // provisional spline on the plain polynomial: create_cp() needs
            // density( T ) and alpha( T ); create_alpha_cryo() replaces it
            this->set_alpha_switch_temperature( 0.0 );
            this->create_spline( MaterialProperty::alpha );
            this->spline( MaterialProperty::alpha )->create_integral( this->constant_property( MaterialProperty::T_ref_density ), 0.0 );
        }

        void
        HastelloyC276::create_alpha_cryo()
        {
            // The polynomial is linear in T near 0 K, whereas alpha must follow
            // cp ( T^3 plus a small electronic term ). Below the split the
            // Grueneisen branch takes over, anchored on the polynomial's VALUE at
            // the split; the derivatives there come from K( T ), because the
            // polynomial's own slope is not trustworthy ( dln(alpha/cp)/dT < 0
            // at every temperature ). The plateau lies above the split.
            real T = this->set_alpha_switch_temperature();
            this->create_cryo_expansion_anchored( polyval( mAlphaPolys( 0 ), T ), mThermalExpansionCryo );
        }

        real
        HastelloyC276::rho_custom( const real T ) const
        {
            if ( T < mTRhoSwitch )
            {
                return polyval( mRhoPolys( 0 ), T );
            }
            else
            {
                return std::exp(polyval( mRhoPolys( 1 ), std::log(T ) ) );
            }
        }

        void
        HastelloyC276::create_mu()
        {
            mMuPolys.set_size( 3 );



            mMuPolys( 0 ) = { 1.08056149E+00, 5.03655878E-03 };
            mMuPolys( 2 ).set_size( 1, 0.34698162 );

            real x0 = 0.3 ;
            real x1 = 0.375 ;


            create_beam_poly(
                x0,
                polyval( mMuPolys( 0 ), x0 ),
                dpolyval( mMuPolys( 0 ), x0 ),
                x1,
                mMuPolys( 2 )(0),
                0.0,
                mMuPolys( 1 ));

            mTMuSwitch = { 1./x1, 1./x0 };

            // NOT routed on purpose. This routine once registered the curve
            // under MaterialProperty::nu ( a slip that left mu at mu0 ), and
            // routing it under mu is not a one-line fix: mu_custom() evaluates
            // mMuPolys( 2 ) - the low-T constant - on the high-T branch, so it
            // returns chi = 0.347 at 300 K instead of the Curie-Weiss value;
            // set_custom( mu ) leaves dmudH on the constant path, which asserts
            // once mu is no longer a constant; and the Maxwell calculator would
            // switch this material to its non-constant-mu path. Until those
            // three are settled, the substrate stays at mu = mu0.
        }

        real
        HastelloyC276::mu_custom( const real H, const real T ) const
        {
            real chi ;

            if ( T < mTMuSwitch( 0  ) )
            {
                chi = mMuPolys( 2 )(0);
            }
            else if (  T < mTMuSwitch( 1  ) )
            {
                chi = polyval( mMuPolys( 1 ), 1./T );
            }
            else
            {
                chi = polyval( mMuPolys( 2 ), 1./T );
            }

            return constant::mu0 * ( chi + 1.0 );
        }

        void
        HastelloyC276::create_cp()
        {
            mTCpSwitch.set_size( 6 );

            real & T0 = mTCpSwitch( 0 ) ;
            real & T1 = mTCpSwitch( 1 ) ;
            real & T2 = mTCpSwitch( 2 ) ;
            real & T3 = mTCpSwitch( 3 ) ;
            real & T4 = mTCpSwitch( 4 ) ;
            real & T5 = mTCpSwitch( 5 ) ;

            mCpPolys.set_size( 6);

            Vector< real > & p0 = mCpPolys( 0 );
            Vector< real > & p1 = mCpPolys( 1 );
            Vector< real > & p2 = mCpPolys( 2 );
            Vector< real > & p3 = mCpPolys( 3 );
            Vector< real > & p4 = mCpPolys( 4 );
            Vector< real > & p5 = mCpPolys( 5 );

            p0 = {
                -2.76893713E-07,
                 6.52481421E-05,
                - 5.23091044E-03,
                3.90556958E-01,
            };

            T0 = std::sqrt(-p0(1)/(3.*p0(0)));

            // for 15 < T < 40 from 10.1063/1.2899058
            T1 = 15.0 ;

            p2 = { 5.99E-04, 0, 1.33E-01, 0 };

            // cv = γ * T + β * T³
            this->set_constant( MaterialProperty::beta, p2( 0 ) );
            this->set_constant( MaterialProperty::gamma, p2( 2 ) );

            // from 40 < T < 120
            p4 = { -1.16634965E-04, + 1.96837543E-02, + 2.10987137E+00, - 6.34829196E+01 };

            T2 = 19.0 ;

            real y0 = polyval( p0, T0*T0 );
            real dy0 = dpolyval( p0, T0*T0 ) * 2.0 * T0 * T0 + y0 ;
            y0 *= T0 ;

            create_beam_poly( T0, y0, dy0,
                T1, polyval( p2, T1 ), dpolyval( p2, T1 ), p1 );

            // use turning point to connect
            T3 = -p4(1)/(3.*p4(0));

            real x0 = 18 ;
            create_beam_poly(
              x0, polyval( p2, x0 ), dpolyval( p2, x0 ),
              T3, polyval( p4, T3 ), dpolyval( p4, T3 ), p3 );

            real f0 = T3 +p3(1)/(3.*p3(0));

            real x1 = 19 ;
            create_beam_poly(
                          x1, polyval( p2, x1 ), dpolyval( p2, x1 ),
                          T3, polyval( p4, T3 ), dpolyval( p4, T3 ), p3 );
            real f1 = T3 + p3(1)/(3.*p3(0));

            real x = x0 ;
            real f = f0 ;

            uint tCount = 0 ;

            while ( std::abs( f ) > 1e-12 && tCount++ < 100 )
            {
                x -= 0.95 * f0 *( x1-x0)/( f1-f0 );
                if ( x < x0 || x > x1 )
                {
                    x = 0.5 * ( x0 + x1 );
                }

                create_beam_poly(
                          x, polyval( p2, x ), dpolyval( p2, x ),
                          T3, polyval( p4, T3 ), dpolyval( p4, T3 ), p3 );
                f = T3 + p3(1)/(3.*p3(0));
                if ( f0 * f > 0.0 )
                {
                    x0 = x ;
                    f0 = f ;
                }
                else
                {
                    x1 = x ;
                    f1 = f ;
                }
            }
            T2 = x ;

            BELFEM_ERROR( tCount < 100, "Couldn't detect cp interval" );


            // we know that the literature value of cp at RT is 427
            // the debye temperature in this region is between 250 and 300 k
            real T = 293.15 ;
            real cp = 427 ;

            // The Debye temperature here is a fitting parameter that makes the
            // Debye cv model reproduce the literature cp at room temperature; it
            // is not the same as the Debye temperature derived from rho or from
            // the low-T T³ heat capacity, so do not interpret it physically.
            // The bracket must contain the root: with γ ≈ 0.133 J/(kg·K²) and
            // R/M ≈ 134 J/(kg·K), cp(293.15, θ) is decreasing in θ and only
            // crosses 427 J/(kg·K) above ≈ 350 K, so [300, 700] is safe.
            x0 = 300.0 ;
            f0 = this->cp_from_debye( T, x0 ) - cp ;
            x1 = 700.0 ;
            f1 = this->cp_from_debye( T, x1 ) - cp ;

            x = x0 ;
            f = -f0 ;

            tCount = 0 ;
            while ( std::abs( f ) > 1e-12 && tCount++ < 100 )
            {
                x -= 0.95 * f0 *( x1-x0)/( f1-f0 );
                if ( x < x0 || x > x1 )
                {
                    x = 0.5 * ( x0 + x1 );
                }

                f = this->cp_from_debye( T, x ) - cp ;
                if ( f0 * f > 0.0 )
                {
                    x0 = x ;
                    f0 = f ;
                }
                else
                {
                    x1 = x ;
                    f1 = f ;
                }
            }
            this->set_constant( MaterialProperty::debye, x );
            T4 = 100.0 ;
            T5 = 250.0 ;

            real dT = 1e-3 ;

            cp = this->cp_from_debye( T5, x ) ;
            real dcpdT = (this->cp_from_debye( T5+dT, x )
                - this->cp_from_debye( T5-dT, x )) / (2.0*dT);

            create_beam_poly(
               T4, polyval( p4, T4 ), dpolyval( p4, T4 ),
                T5, cp, dcpdT , p5 );

            this->set_have( MaterialProperty::cp );
            this->create_spline( MaterialProperty::cp, p0(p0.length()-1) );
        }

        void
        HastelloyC276::create_lambda()
        {
            mLambdaPolys.set_size( 2 );

            Vector< real > & p0 = mLambdaPolys( 0 );
            Vector< real > & p1 = mLambdaPolys( 1 );

            p0 = { -1.26684292E-05, - 6.33602091E-04, + 1.81777951E-01, 0.0 };

            real a = 3 * p0(0);
            real b = 2 * p0(1);
            real c = p0(2);

            mTLambdaSwitch = ( -b - std::sqrt( b*b - 4.0*a*c ) ) / ( 2.0*a );
            p1 = { polyval( p0, mTLambdaSwitch ) };

            this->set_have( MaterialProperty::lambda );
            this->create_spline( MaterialProperty::lambda );
        }

        real
        HastelloyC276::cp_custom( const real T ) const
        {
            if ( T < mTCpSwitch( 0 ) ) return polyval( mCpPolys( 0 ), T*T )*T;
            if ( T < mTCpSwitch( 1 ) ) return polyval( mCpPolys( 1 ), T );
            if ( T < mTCpSwitch( 2 ) ) return polyval( mCpPolys( 2 ), T );
            if ( T < mTCpSwitch( 3 ) ) return polyval( mCpPolys( 3 ), T );
            if ( T < mTCpSwitch( 4 ) ) return polyval( mCpPolys( 4 ), T );
            if ( T < mTCpSwitch( 5 ) ) return polyval( mCpPolys( 5 ), T );

            return this->cp_from_debye( T, this->constant_property( MaterialProperty::debye ) );
        }

        real
        HastelloyC276::lambda_custom( const real T ) const
        {
            if ( T == 0.0 ) return 0.0;

            real lambda_ph = T < mTLambdaSwitch ?
                polyval( mLambdaPolys( 0 ), T ) : mLambdaPolys( 1 )( 0 ) ;

            real lambda_e = constant::L0 * T / this->rho_custom( T );
            return lambda_ph + lambda_e;
        }

    }
}