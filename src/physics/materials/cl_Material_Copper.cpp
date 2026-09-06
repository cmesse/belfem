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


#include "fn_gesv.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_create_fifth_order_beam_poly.hpp"
#include "cl_Material_Abundance.hpp"

#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "cl_Material_Copper.hpp"
namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * Copper material property implementation
         *
         * Data sources:
         * - Mechanical: copper.org, Blanke (Thermophysikalische Stoffgrößen, Springer 1989)
         * - Thermal expansion: Touloukian dataset
         * - Specific heat: Touloukian dataset
         * - Resistivity: Matula, Bloch-Grüneisen model with reference at 273.15K
         * - Thermal conductivity: Hust 1984 (10.6028/nbs.ir.84-3007)
         * - Magnetoresistance (Kohler): De Launey 1959, Benz 1969, Arentz 1982
         */
        Copper::Copper( const real RRR, const bool aBuildTables ) :
            Metal( "Copper", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();
            this->create_cp();             // Specific heat capacity
            this->create_alpha();          // Thermal expansion coefficiens


            // Reference resistivity from Matula, results in theta = 320 @ 25 C
            this->set_rho_i_ref( 273.15, 1.541e-8, 321.352 );

            this->create_debye();          // Debye temperature

            this->create_rho();

            // Thermal conductivity coefficients from Hust 1984, 10.6028/nbs.ir.84-3007
            this->set_lambda_coefficients( {0.1661, 1.754e-8, 2.763, 1102, -0.165, 70, 1.765, 0.838} );

            this->create_kohler() ;  // Magnetoresistance curves

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Copper.org and Wolfram Cloud
            this->create_mech( 134.447, 0.11492, 922.14, 293.15, 0.344 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }

        }

        Copper::~Copper()
        {
            if ( mKohlerLongBezier != nullptr ) delete mKohlerLongBezier ;

            if ( mThermalExpansion != nullptr )  delete mThermalExpansion ;
        }

        void
        Copper::set_constants()
        {
            // set the melting temperature
            this->set_constant( MaterialProperty::T_max, 1358.);

            // reference density at room temperature
            this->set_constant( MaterialProperty::ref_density, 8933. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );


            Abundance tAbundance ;
            this->set_constant( MaterialProperty::M, tAbundance.compute_molar_mass( "Cu" ) );

            // Sommerfeld and Debye coefficients of cp ≈ cv = γ * T + β * T³,
            // extracted from the Touloukian data. Must be set after M, because
            // writing beta derives debye0K from the molar mass.
            this->set_constant( MaterialProperty::gamma, 0.010969801407435 );
            this->set_constant( MaterialProperty::beta,  0.000753034913427 );

            // TPPM vol. 1
            this->set_constant( MaterialProperty::rho0_pure, 0.000851e-8 );

        }


        /**
         * Create thermal expansion coefficient α(T)
         *
         * A cubic Bezier for dL/L( T ) referred to 293.15 K ( Touloukian ),
         * alpha = L'/L above the split temperature; below it the Grueneisen
         * branch alpha = C( T )·cp( T ) from create_cryo_expansion(). The
         * alpha spline and its integral ( for l( T ) ) are built by
         * finish_cryo_expansion().
         */
        void
        Copper::create_alpha()
        {
            // mapped against Touloukian dataset
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 79.424524719, 293.15, 1358.0 };
            mThermalExpansion->basis_y() = { -0.340793558792, -0.340793558792, -0.197203414821, 2.173881258550 };
            mThermalExpansion->basis_y() *= 0.01 ;

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }

        void
        Copper::create_cp()
        {
            // fitted against Touloukian dataset
            Metal::create_cp( { 2.61467567494, 3.35435402912, 3.74800814310, 4.79681400831 },
                             {0.782267668755, 3.02229820851,5.02776431176, 5.66825740852 },
                             {4.79681400831, 5.53151720444, 6.59829615326, 7.21303166000 },
                             { 5.66825740852,  6.11693179231, 5.90224773827, 6.25860246932 } );

            // final routing for cp: the analytic Beziers, not the sampled spline.
            // The cryogenic alpha branch is fitted against cp and must see the
            // same curve that is evaluated at run time, so the routing has to be
            // settled here and not in create_debye()
            this->set_custom( MaterialProperty::cp );
        }

        void
        Copper::create_debye()
        {
            mTDebyeSwitch.set_size( 3 );

            real & T1 = mTDebyeSwitch( 0 );
            real & T2 = mTDebyeSwitch( 1 );
            real & T3 = mTDebyeSwitch( 2 );


            mDebyePolys.set_size( 4, {} );

            Vector< real > & p0 = mDebyePolys( 0 );
            Vector< real > & p1 = mDebyePolys( 1 );
            Vector< real > & p2 = mDebyePolys( 2 );
            Vector< real > & p3 = mDebyePolys( 3 );

            real y0 = this->constant_property( MaterialProperty::debye0K );

            // chosen value to get well behaved curves
            T1 = 14.0 ;
            real dT = 1e-6 ;
            real y1 =  this->compute_debye_from_cv( T1 );
            real dy1 = ( this->compute_debye_from_cv( T1+dT ) - this->compute_debye_from_cv( T1-dT ) ) / ( 2.0 * dT ) ;

            T2 = 50.0 ;

            // optimized to fit Matula dataset for 70 K < T < 298.15 K ;

            real Tm = 106.322773 ;
            real ym = 327.900623 ;

            // fitted against data from Matula for values of T  > 298.15 K ;
            Vector< real > x = { 273.15 , 298.15, 1100.0 };
            Vector< real > y = { 321.352, 320.0, 276.99 };

            T3 = x(0) ;
            real y3 = y(0) ;
            polyfit( x,y,2,p3 );

            // middle poly
            real x2 = Tm * 0.001 ;
            real x3 = T3 * 0.001 ;
            Matrix< real > V = {
                { x2*x2*x2*x2, x2*x2*x2, x2*x2, x2, 1. },
                          { 4.*x2*x2*x2, 3*x2*x2, 2.*x2, 1., 0 },
                          { x3*x3*x3*x3, x3*x3*x3, x3*x3, x3, 1. },
                          { 4.*x3*x3*x3, 3*x3*x3, 2.*x3, 1., 0 },
                          { 12.*x3*x3, 6.*x3, 2., 0., 0. }
            };

            p2 = { ym, 0., y3, dpolyval( p3, T3 ) * 1.e3, ddpolyval( p3, T3 ) * 1.e6 };
            Vector< int_t > tPivot( 5, 0. );
            gesv(  V, p2, tPivot );

            p2( 0 ) *= 1e-12 ;
            p2( 1 ) *= 1e-9 ;
            p2( 2 ) *= 1e-6 ;
            p2( 3 ) *= 1e-3 ;

            // ultra low poly
            create_beam_poly(
               0.0,
               y0,
               0.0,
               T1,
               y1,
               dy1,
               p0 );

            // low poly
            create_beam_poly(
               T1,
               y1,
               dy1,
               T2,
               polyval(p2,T2),
               dpolyval( p2,T2 ),
               p1 );

            // to obtain the derivatives, we use the analytic forms instead of the
            // splines. E is not touched here: create_mech() owns its routing
            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::cp );
            this->set_have( MaterialProperty::debye );

            this->create_spline( MaterialProperty::debye, 0.0 );
        }



        real
        Copper::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        real
        Copper::debye_custom( const real T ) const
        {
            if ( T <= 0 ) return this->constant_property( MaterialProperty::debye0K );

            if ( T < mTDebyeSwitch( 0 ) ) return polyval( mDebyePolys(0), T );
            if ( T < mTDebyeSwitch( 1 ) ) return polyval( mDebyePolys(1), T );
            if ( T < mTDebyeSwitch( 2 ) ) return polyval( mDebyePolys(2), T );

            return polyval( mDebyePolys( 3 ), T );
        }

        void
        Copper::create_kohler()
        {
            // fit to data from
            // De Launey 1959, 10.1016/0022-3697(59)90038-1
            // Benz 1969, 10.1063/1.1657896
            // Arentz. 1982, 10.1103/PhysRevB.26.2727
            mKohlerTransPolys.set_size( 3 );
            Vector< real > & q1 = mKohlerTransPolys( 1 );
            q1 = {  -2.3638899907E-02,1.3389925106, - 7.2085920334 };

            // fit to data from
            // De Launey 1959, 10.1016/0022-3697(59)90038-1
            // Clausecker 1969, 10.1007/BF02422537
            // Strom-Olsen 1967, https://www.jstor.org/stable/2415889
            // Arentz. 1982, 10.1103/PhysRevB.26.2727
            mKohlerLongBezier = new Bezier(
            { 2.3, 5.7581541250, 5.3025840143, 10.0 },
            {-6.7358206535,-1.7358206535, 0.2114121663,   0.2114121663 } );

            // polynomial for low field area
            mBSKohlerSwitch.set_size( 2 );

            real w = 2.3; // w = ln(x) w' = 1/x w'' = -1/x^2
            real x = std::exp( w );

            // transition points
            mBSKohlerSwitch( 0 ) = x ;
            mBSKohlerSwitch( 1 ) = std::exp( 10.0 ); // field saturates here


            mKohlerLongPolys.set_size( 2 );


            x = std::exp( w );

            // compute longitudinal polynomial for low field
            // let f = exp(g) and w = ln(x)
            real g      = mKohlerLongBezier->y( w );
            real dgdw   = mKohlerLongBezier->dydx( w ) ;
            real d2dgw2 = mKohlerLongBezier->d2ydx2( w ) ;
            real dgdx = dgdw / x ;
            real d2gdx2 = ( d2dgw2 - dgdw)/(x*x);

            real f = std::exp( g );
            real dfdx = f * dgdw / x ;
            real d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;

            Vector< real > & p0 = mKohlerLongPolys( 0 );
            p0.set_size( 4 );
            p0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            p0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            p0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            p0( 3 ) = 0.0 ;


            Vector< real > & p1 = mKohlerLongPolys( 1 );
            p1.set_size( 1, std::exp(mKohlerLongBezier->y( 10.0 )) );


            // compute transversal polynomial for low field
            g = polyval( q1, w );
            dgdw = dpolyval( q1, w );
            d2dgw2 = ddpolyval( q1, w );
            dgdx = dgdw / x ;
            d2gdx2 = ( d2dgw2 - dgdw)/(x*x);

            f = std::exp( g );
            dfdx = f * dgdw / x ;
            d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;
            Vector< real > & q0 = mKohlerTransPolys( 0 );
            q0.set_size( 4 );
            q0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            q0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            q0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            q0( 3 ) = 0.0 ;

            // compute transversal polynomial for high field
            Vector< real > & q2 = mKohlerTransPolys( 2 );
            q2.set_size( 2 );
            q2(0) = dpolyval( q1, 10.0 ) ;
            q2(1) = polyval( q1, 10.0 ) - q2(0)*10.0 ;

            this->set_kohler_dependencies();
        }



        real Copper::kohler( const real B, const real S, const real beta ) const
        {
            if ( B < BELFEM_EPSILON ) return 0.0 ;

            real BS = B * S ;

            real Along ;
            real Atrans ;
            if ( BS < mBSKohlerSwitch( 0 ) )
            {
                Along  = polyval( mKohlerLongPolys(0), BS );
                Atrans = polyval( mKohlerTransPolys(0), BS );
            }
            else if ( BS < mBSKohlerSwitch( 1 ) )
            {
                real w = std::log( BS ) ;
                Along  = std::exp( mKohlerLongBezier->y( w ) );
                Atrans =  std::exp( polyval( mKohlerTransPolys(1), w ) );
            }
            else
            {
                Along  =  mKohlerLongPolys(1)(0);
                Atrans = std::exp( polyval( mKohlerTransPolys(2), std::log( BS ) ) );
            }

            real c = std::cos( beta );
            real c2 = c*c ;
            real s2 = 1.0 - c2 ;

            // Pippard's angular interpolation formula
            return Along * c2 + Atrans * s2 ;
        }

    }
}