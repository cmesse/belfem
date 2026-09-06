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

#include "cl_Material_YBCO.hpp"

#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "cl_Material_Abundance.hpp"
#include "debye.hpp"


namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * YBCO (YBa₂Cu₃O₇) material property implementation
         *
         * This is the most complex material in the database, combining:
         * - Normal-state metallic properties (Bloch-Grüneisen resistivity)
         * - Superconducting properties (d-wave gap, power law resistivity)
         * - Advanced thermal conductivity (Callaway model with multiple scattering mechanisms)
         *
         * Key implementation features:
         * - Mechanical: Normalized curves from Lei & Ledbetter scaled to user values
         * - Thermal expansion: Polynomial with constant tail at high T
         * - Specific heat: 4 piecewise regions (polynomial, log-log fits)
         * - Electrical: Bloch-Grüneisen + RRR model (fitted to Sommerfeld 2003)
         * - Thermal conductivity: Callaway model including:
         *   * Phonon-phonon umklapp scattering
         *   * Boundary scattering (layer thickness dependent)
         *   * Electron-phonon scattering
         *   * Optical phonon scattering (Raman active mode at 501 cm⁻¹)
         *   * D-wave superconducting gap effects
         *
         * Data sources:
         * - Lei & Ledbetter 1991 (NISTIR 3980) - mechanical
         * - Salomons et al. 1987 (10.1016/0378-4363(87)90092-1) - thermal expansion
         * - Baak 1989, Lang et al. 1988 - specific heat
         * - Sommerfeld et al. 2003 (10.1103/PhysRevB.67.174520) - fitted parameters
         */
        YBCO::YBCO() :
            Metal( "YBCO", MaterialType::HTS )
        {
            mReferenceValues.set_size( 12, 1.0 );

            this->set_constant( MaterialProperty::T_max, 400.0 );

            this->create_mech();        // Young's modulus, Poisson's ratio
            this->create_expansion();   // Thermal expansion coefficient ( polynomial )
            this->create_cp();          // Specific heat capacity

            Abundance tAbundance ;

            real M ;
            real Gamma ;
            tAbundance.compute_molar_mass_and_impurity_from_volumes(
                { "Y", "Ba", "Cu", "O" },
                {  1.0, 2.0, 3.0, 7.0 },
                M,
                Gamma );


            this->set_constant( MaterialProperty::M, M );

            // Gamma from the isotope composition ( 1.1e-5 ) is far below the
            // point-defect scattering of a coated conductor. The value used for
            // the Callaway kernel is an EFFECTIVE, fitted disorder parameter -
            // see the note in lambda_custom() on what the fit can and cannot
            // determine.
            ( void ) Gamma ;
            this->set_constant( MaterialProperty::Gamma, 1.0 );

            //this->set_constant( MaterialProperty::M, M );

            this->set_constant( MaterialProperty::ref_density, 6300 );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );



            // typical value
            this->set_constant( MaterialProperty::T_crit, 92.5 );

            // number of atoms per molecule
            this->set_constant( MaterialProperty::q, 13 ); // assuming YBa2Cu3O7

            // assuming a reasonable value for rho_i_ref and theta [fitted]

            // fitted against measurements
            // Sommerfeld et. al. Phys Rev 2003
            // 10.1103/PhysRevB.67.174520
            this->set_constant( MaterialProperty::debye, 275 );
            this->set_rho_i_ref( 273.15, 60e-8 );

            this->set_constant( MaterialProperty::layer_thickness, 1.0e-6 );
            this->set_custom( MaterialProperty::density );

            // we assume reasonable values for thin films
            this->set_young( 210e9, gTroom );
            this->set_poisson( 0.25, gTroom );

            // the cryogenic alpha branch needs cp and K( T ) - hence after cp,
            // E and nu - and must precede the lambda spline below, whose Callaway
            // input grueneisen( T ) reads alpha
            this->create_expansion_cryo();

            this->set_RRR( 50.0 );

            this->set_have( MaterialProperty::lambda );
            // start tangent of the electronic term L0 * T / ( rho_i + rho_0 )
            // at T -> 0; the Callaway phonon term starts flat
            this->create_spline( MaterialProperty::lambda,
                constant::L0 / this->constant_property( MaterialProperty::rho_0 ) );

        }

        YBCO::~YBCO()
        {
            if ( mYoungBezier != nullptr ) delete mYoungBezier ;
        }


        /**
         * Create mechanical properties using normalized curves
         *
         * Uses Lei & Ledbetter 1991 (NISTIR 3980, Fig 4.3) which provides
         * normalized E(T)/E(RT) and ν(T)/ν(RT) curves. These are then scaled
         * by user-specified reference values via set_young() and set_poisson().
         *
         * Young's modulus: Bezier curve for T < 100K, linear for T > 100K
         * Poisson's ratio: Cubic polynomial for T < 100K, linear for T > 100K
         */
        void YBCO::create_mech()
        {
            // This is the normalized curve from Lei, Ledbetter 1991, NISTIR 3980, Fig 4.3

            // Reference temperature
            real xref = 295. ;
            real yref = 1.0 ;

            // end temperature for this curve
            real x2 = 100.0 ;

            // value at 0 K
            real x0 = 0.0 ;
            real y0 = 1.06523215 ;

            // value at x1
            real y2 = 1.05097719 ;

            // center point
            real x1 = 0.5 * ( x2 - x0 );

            // linear extrapolation for control value
            real a = ( yref  - y2 ) / ( xref - x2 );
            real b = yref - a * xref;
            real y1 = a * x1 + b ;

            mYoungBezier = new Bezier();

            mYoungBezier->basis_x() = { x0, x1, x1, x2 };
            mYoungBezier->basis_y() = { y0, y0, y1, y2 };
            mYoungPoly = { a, b };

            // linear approximation for curve for nu, also from Fig 4.3
            a = 4.7519e-5 ;
            b = yref - a * xref ;
            y0 = 0.9859 ;
            y2 = a * x2 + b ;

            mPoissionPolys.set_size( 2, {} );
            create_beam_poly( x0,  y0, 0.0, x2, y2, a, mPoissionPolys( 0 ) );
            mPoissionPolys( 1 ) = { a, b };
        }

        void
        YBCO::set_young( const belfem::real aEref, const belfem::real aTref )
        {
            uint tIndex = static_cast< uint >( MaterialProperty::E );
            mReferenceValues( tIndex ) = 1.0 ;

            real tScale =  this->E_custom( aTref );
            mReferenceValues( tIndex ) = aEref / tScale ;
            this->set_have( MaterialProperty::E );
            this->create_spline( MaterialProperty::E, 0.0 );
        }

        void
        YBCO::set_poisson( const belfem::real aNuref, const belfem::real aTref )
        {
            uint tIndex = static_cast< uint >( MaterialProperty::nu );
            mReferenceValues( tIndex ) = 1.0 ;

            real tScale =  this->nu_custom( aTref );
            mReferenceValues( tIndex ) = aNuref / tScale ;

            BELFEM_ERROR( this->nu_custom( this->constant_property( MaterialProperty::T_max ) ) < 0.49,
                "Chosen function exceeds nu = 0.49 at maximum temperature.");

            this->set_have( MaterialProperty::nu );
            this->create_spline( MaterialProperty:: nu, 0.0 );
        }

        real
        YBCO::E_custom( const real T ) const
        {
            uint tIndex = static_cast< uint >( MaterialProperty::E );

            if ( T < mYoungBezier->basis_x()(3) )
            {
                return mReferenceValues( tIndex ) * mYoungBezier->y( T );
            }
            else
            {
                return mReferenceValues( tIndex ) * polyval( mYoungPoly, T );
            }
        }

        real
        YBCO::nu_custom( const real T ) const
        {
            uint tIndex = static_cast< uint >( MaterialProperty::nu );

            if ( T < mYoungBezier->basis_x()(3) )
            {
                return mReferenceValues( tIndex ) * polyval( mPoissionPolys( 0 ), T );
            }
            else
            {
                return mReferenceValues( tIndex ) * polyval( mPoissionPolys( 1 ), T );
            }
        }

        void
        YBCO::create_expansion()
        {
            // extracted from dataset of Salomons et. al, 1987, 10.1016/0378-4363(87)90092-1

            Vector< real > & p = mAlphaPoly ;

            p = { 2.00236453619024E-16, - 7.06531471374544E-14, - 1.23172613411674E-10, + 7.79099007152825E-08, 0.0 };

            // this polynomial has a maximum at around 350 K. We just chop the dataset at this point

            real & x = mTAlphaPlateau ;

            x = 350 ;

            real f = dpolyval( p, x );
            real df ;

            while ( abs( f ) > BELFEM_EPSILON )
            {
                f = dpolyval( p, x );
                df = ddpolyval( p, x );
                x -= f / df ;
            }

            mReferenceValues( static_cast< size_t >( MaterialProperty::alpha ) ) = polyval( p, x );

            this->set_have( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::alpha );

            // provisional spline on the plain polynomial, replaced by
            // create_expansion_cryo() once cp exists
            this->set_alpha_switch_temperature( 0.0 );
            this->create_spline( MaterialProperty::alpha );
            this->spline( MaterialProperty::alpha )->create_integral( 0.0 );
        }

        void
        YBCO::create_expansion_cryo()
        {
            // the quartic is linear in T near 0 K; below the split the Grueneisen
            // branch alpha = C( T ) cp( T ) takes over, anchored on the quartic's
            // VALUE at the split with derivatives from K( T ) - the quartic's own
            // slope is not trustworthy ( dln(alpha/cp)/dT < 0 everywhere ). The
            // plateau at 350 K lies above the split.
            const Vector< real > & p = mAlphaPoly ;
            real T = this->set_alpha_switch_temperature();
            this->create_cryo_expansion_anchored( polyval( p, T ), mThermalExpansionCryo );
        }

        real YBCO::alpha_custom( const real T ) const
        {
            if ( T < this->alpha_switch_temperature() )
            {
                // Grueneisen branch: alpha = C( T ) cp( T )
                return std::exp( polyval( mThermalExpansionCryo, T ) ) * this->cp( T );
            }
            else if ( T < mTAlphaPlateau )
            {
                return polyval( mAlphaPoly, T );
            }
            else
            {
                return mReferenceValues( static_cast< size_t >( MaterialProperty::alpha ) );
            }
        }

        void
        YBCO::create_cp()
        {


            mCpPolys.set_size( 4, {} );

            Vector< real > & p0 = mCpPolys( 0 );
            Vector< real > & p1 = mCpPolys( 1 );
            Vector< real > & p2 = mCpPolys( 2 );
            Vector< real > & p3 = mCpPolys( 3 );

            mTCpSwitch.set_size( 3 );
            real & T0 = mTCpSwitch( 0 );
            real & T1 = mTCpSwitch( 1 );
            real & T2 = mTCpSwitch( 2 );

            // cv = γ * T + β * T³
            // based on curve 1, extracted from Baak 1989, 10.1016/0921-4534(89)91125-8
            // valid for T < 20 K
            p0 = { 4.00111872035290E-04, 0., 7.31512890202218E-03, 0. };

            this->set_constant(  MaterialProperty::M, 0.6661937 );
            this->set_constant(  MaterialProperty::beta, p0( 0 ) );
            this->set_constant(  MaterialProperty::gamma, p0( 2 ) );

            // fitted for 50 < T < 350, using dataset of Lang et al, 1988, 10.1007/BF01312506
            p2 = {  1.25814185256104E-01, - 2.18941935692446E+00, + 1.31963053875634E+01, - 2.14028044611340E+01 };

            // this polynomial has a turning point at x ≈ 5.8

            T0 = 20.0 ;
            T1 = 50.0 ;

            real x0   = std::log( T0 );
            real x1   = std::log( T1 );
            real g0   = polyval( p0, T0 );
            real dg0  = dpolyval( p0, T0 );

            real f0 = std::log( g0 );
            real df0 = dg0 / g0 * T0 ;

            real f1 = polyval( p2, x1 );
            real df1 = dpolyval( p2, x1 );

            create_beam_poly( x0, f0, df0, x1, f1, df1, p1 );

            x0 = 5.79 ;
            f0 = ddpolyval( p2, x0 );

            x1 = 5.81 ;
            f1 = ddpolyval( p2, x1 );

            real f = f1 ;
            real x = 0.5 * ( x0 + x1 );

            while ( abs( f ) > BELFEM_EPSILON )
            {
                f = ddpolyval( p2, x );
                if ( f0 * f < 0 )
                {
                    x1 = x ;
                }
                else
                {
                    x0 = x ;
                    f0 = f ;
                }
                x = 0.5 * ( x0 + x1 );
            }

            T2 = std::exp( x );

            real a = dpolyval( p2, x );
            real b = polyval( p2, x ) - a * x ;

            p3 = { a, b };

            this->set_have( MaterialProperty::cp );
            this->create_spline( MaterialProperty::cp, 0.0 );
        }

        real YBCO::cp_custom( const real T ) const
        {
            if ( T <= 0 ) return 0.0;
            if ( T < mTCpSwitch( 0 ) )  return polyval( mCpPolys( 0 ), T );
            if ( T < mTCpSwitch( 1 ) )  return std::exp(polyval( mCpPolys( 1 ), std::log(T)));
            if ( T < mTCpSwitch( 2 ) )  return std::exp(polyval( mCpPolys( 2 ), std::log(T)));

            return std::exp(polyval( mCpPolys( 3 ), std::log(T)));
        }

        real
        YBCO::test_callaway( const real T,  const Vector< real > & aParams )
        {


            Vector< double > params( 20 );

            params( 0 ) = T ;
            params( 1 ) = this->debye( T );
            params( 2 ) = this->group_velocity( T );
            params( 3 ) = this->density( T );
            params( 4 ) = this->G( T );
            params( 5 ) = this->grueneisen( T  );

            params( 6 ) = this->constant_property( MaterialProperty::M );
            params( 7 ) = this->constant_property( MaterialProperty::Gamma );

            params( 8 ) =  this->constant_property( MaterialProperty::T_crit );

            real theta = 275.0 ;
            this->set_constant( MaterialProperty::debye, theta );
            //this->set_constant( MaterialProperty::debye0K, theta );

            this->set_rho_i_ref( 273.15, 60e-8, theta);


            params( 9 )  = this->constant_property( MaterialProperty::layer_thickness ); // same L0 as lambda_custom
            params( 10 ) = 4.0 ; //
            params( 11 ) = 26.14; // b parameter for Umklapp scattering
            params( 12 ) = 5.0; // d parameter for Umklapp scattering
            params( 13 ) = 10  ;
            params( 14 ) = 5 ; // factor for effective electron mass, typical, 10.1103/PhysRevLett.62.2317
            params( 15 ) = 2e27 ; // conduction electrons concentration, typical, 10.1103/PhysRevLett.62.2317
            params( 16 ) = 2.1; // Delta 0 factor for d-wave coupling [fitted]
            params( 17 ) = 0.0 ;
            params( 18 ) = aParams( 1 ); // lambda_opt [fitted]
            params( 19 ) = 501.0 ;  // Raman shift oxygen band peak 10.1103/PhysRevB.80.064505

            real RRR = 50.0;


            real k_ph ;
            int status = 0 ;
            callaway_conductivity(  params.data(), &k_ph, &status );

            // real Tc = params(9);
            real rho_0 = rho_i_custom(273.15) / RRR ;


            real rho_i = rho_i_custom(T);                    // ideal e-ph resistivity
            real k_e = constant::L0 * T / ( rho_i + rho_0 );            // normal state κ_e

            return k_e + k_ph  ;
        }

        /**
         * Thermal conductivity using Callaway model
         *
         * Implements λ = λ_e + λ_ph where:
         * - λ_e: Electronic contribution via Wiedemann-Franz law
         * - λ_ph: Phonon contribution via Callaway model
         *
         * The Callaway model includes:
         * 1. Umklapp scattering (phonon-phonon): τ_U ∝ exp(θD/(b·T))
         * 2. Boundary scattering: τ_B = v·L (layer thickness dependent)
         * 3. Electron-phonon scattering: temperature and carrier dependent
         * 4. Optical phonon scattering: Raman mode at 501 cm⁻¹
         * 5. D-wave superconducting gap: Δ₀/(kB·Tc) factor
         *
         * All parameters fitted to Sommerfeld et al. 2003 (10.1103/PhysRevB.67.174520)
         */
        real
        YBCO::lambda_custom(const real T) const
        {
            if ( T <= 0 ) return 0.0;

            double params[ 20 ];

            // Fitted against measurements from Sommerfeld et al. 2003
            // (10.1103/PhysRevB.67.174520).
            //
            // 2026-08-25: the optical channel's unit defect in debye.f90 was
            // fixed, and alpha( T ) below the split now follows cp, which made
            // grueneisen( T ) a sane ~2.1 instead of ~800 at 2 K. Both changes
            // opened the phonon channel wide, so b, d and Gamma were refitted.
            // The fit is DEGENERATE: the electronic term k_e = L0 T / rho( T )
            // below alone reproduces the data to 10 % rms, and every fit
            // reaches the same floor by suppressing k_ph ( ~ 1 W/(m K) ). The
            // phonon parameters are therefore interim values, not measurements;
            // the open modelling question is the electronic term, whose
            // resistivity anchor sits well below the literature and dominates.

            params[  0 ] = T ;                        // Temperature in K

            // temperature dependent properties
            params[  1 ] = this->debye( T );          // debye temperature in K
            params[  2 ] = this->group_velocity( T ); // vg in m/s
            params[  3 ] = this->density( T );        // density in kg/m³
            params[  4 ] = this->G( T );              // Shear Modulus in Pa
            params[  5 ] = this->grueneisen( T );

            // material constants
            params[  6 ] = this->constant_property( MaterialProperty::M ); // molar mass
            params[  7 ] = this->constant_property( MaterialProperty::Gamma ); // impurity parameter
            params[  8 ] = this->constant_property( MaterialProperty::T_crit ); // critical temperature

            // model specific constants
            params[  9 ] = this->constant_property( MaterialProperty::layer_thickness ); // thickness of probe
            params[ 10 ] = 4.0 ; //  exponent for the function, usually n=4

            // parameters for Umklapp scattering
            // C++ index | Fortran index | meaning
            //    11     | params(12)    | b parameter [fitted]
            //    12     | params(13)    | d parameter [fitted]
            params[ 11 ] = 26.14 ;
            params[ 12 ] = 5.0 ;

            // parameters for electron-phonon scattering
            params[ 13 ] = 10.0 ; // deformation potential in eV, typical value
            params[ 14 ] = 5. ; // factor for effective electron mass, typical, 10.1103/PhysRevLett.62.2317
            params[ 15 ] = 2e27 ; // conduction electrons concentration, typical, 10.1103/PhysRevLett.62.2317
            params[ 16 ] = 2.1 ; // Δ0/(kB*T),  parameter [ fitted ]
            params[ 17 ] = 0.0 ;   // λ_ph, value for acoustic electron-phonon coupling (not used)
            params[ 18 ] = 1.06 ;  // λ_opt, parameter for optical scattering [fitted]
            params[ 19 ] = 501.0 ; // omega_opt, Raman shift oxygen band peak [cm^-1], 10.1103/PhysRevB.80.064505


            real k_ph ;
            int status = 0 ;
            callaway_conductivity(  params, &k_ph, &status );

            BELFEM_ASSERT( status == 0, "Error in callaway_conductivity" );

            real rho_0 = this->constant_property( MaterialProperty::rho_0 );

            real rho_i = rho_i_custom(T);                    // ideal e-ph resistivity
            real k_e = constant::L0 * T / ( rho_i + rho_0 );  // normal state κ_e

            return k_e + k_ph  ;
        }
    }
}
