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

#include "fn_create_beam_poly.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"

#include "fn_polyfit.hpp"
#include "petsctools.hpp"

#include "cl_Material_Abundance.hpp"
#include "cl_Material_Silver.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * Silver material property implementation
         *
         * Data sources:
         * - Constants: Smith and Fickett (10.6028/jres.100.012)
         * - Mechanical: Fitted against experimental data
         * - Thermal expansion: Touloukian dataset
         * - Specific heat: Touloukian dataset
         * - Resistivity: Matula, Bloch-Grüneisen model with reference at 273.15K
         * - Thermal conductivity: Fitted to Mendelssohn 1952, Gerritsen 1952, Li 2015
         *   NOTE: Silver conductivity data with RRR is surprisingly sparse
         * - Magnetoresistance (Kohler): Lüthi 1960, Strom-Olsen 1967, Stout 1939, Iwasa 1993
         */
        Silver::Silver( const real RRR, const bool aBuildTables ) :
            Metal( "Silver", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();

            this->create_cp();             // Specific heat capacity

            this->create_alpha();          // Thermal expansion coefficient
            this->set_custom( MaterialProperty::density );

            // Reference resistivity from Matula
            this->set_rho_i_ref(  273.15, 1.466e-8 , 221.991 );

            this->create_debye();


            this->create_rho() ;

            /* Thermal conductivity coefficients fitted to sparse experimental data.
             * While Copper is very well documented, reliable conductivity data for
             * silver in conjunction with the RRR value are surprisingly sparse.
             * Main data sources:
             * - Mendelssohn 1952, sample Agm RRR 30.9, 10.1088/0370-1298/65/6/301
             * - Gerritsen 1952, sample AG2T, RRR 407, 10.1016/S0031-8914(56)90035-0
             * - Li 2015, 10.1088/1757-899X/102/1/012027
             */
            this->set_lambda_coefficients(
                {
                    2.11674345e-01,
                    3.12572182e-08,
                    2.94024132e+00,
                    9.90947611e+02,
                    -1.35325295e-01,
                    1.04277725e+01,
                    1.21386795e+00,
                    8.07441604e-01 } );

            this->create_kohler();  // Magnetoresistance curves

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Wolfram Cloud
            this->create_mech( 86.647, 0.06298, 587.35, 293.15, 0.37 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }

        }

        Silver::~Silver()
        {
            delete mThermalExpansion ;
        }

        void
        Silver::set_constants()
        {
            // constants taken from from Smith and Fickett, 0.6028/jres.100.012

            // set the melting temperature
            this->set_constant( MaterialProperty::T_max, 1235.08 );

            // reference density at room temperature
            this->set_constant( MaterialProperty::ref_density, 10492. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );

            Abundance tAbundance ;
            this->set_constant( MaterialProperty::M, tAbundance.compute_molar_mass( "Ag" ) );

            this->set_custom( MaterialProperty::density );

            // fitted against dataset from Touloukian
            this->set_constant( MaterialProperty::gamma, 5.87443823888e-3 );
            this->set_constant( MaterialProperty::beta, 1.53669606066e-3 );

            real R = this->constant_property( MaterialProperty::R ) ;
            real theta = this->constant_property( MaterialProperty::debye0K );
            real pi = constant::pi ;
            real beta = ( 2.4 * pi * pi * pi * pi * R ) / ( theta * theta * theta ) ;
            this->set_constant( MaterialProperty::beta, beta );

            // TPPM vol. 1
            this->set_constant( MaterialProperty::rho0_pure, 0.000620e-8   );

        }

        void
        Silver::create_alpha()
        {
            // mapped against Touloukian dataset
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 49.022633809, 293.15, 1235.08 };
            mThermalExpansion->basis_y() = {-0.434387147544, -0.434387147544, -0.183012211724, 2.126971848186 };
            mThermalExpansion->basis_y() *= 0.01 ;

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );
            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }

        void
        Silver::create_cp()
        {
            // fitted against Touloukian dataset
            Metal::create_cp( { 2.78202115300, 2.96259724117, 3.53602466952, 4.41755060119 },
                            { 2.07290202941, 2.66333881571, 4.64553042143, 5.15013111673 },
                            { 4.41755060119, 5.18492206955, 7.07069952285, 7.10918477700 },
                            { 5.15013111673, 5.58938782216, 5.41551370364, 5.77804911703 } );

            // final routing for cp: the analytic Beziers, not the sampled spline.
            // The cryogenic alpha branch is fitted against cp and must see the
            // same curve that is evaluated at run time, so the routing has to be
            // settled here and not in create_debye()
            this->set_custom( MaterialProperty::cp );
        }


        void Silver::create_debye()
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
            T1 = 12.0 ;
            real dT = 1e-6 ;
            real y1 =  this->compute_debye_from_cv( T1 );
            real dy1 = ( this->compute_debye_from_cv( T1+dT ) - this->compute_debye_from_cv( T1-dT ) ) / ( 2.0 * dT ) ;


            // optimized to fit Matula dataset for 70 K < T < 298.15 K ;

            real Tm = 30.921895 ;
            real ym = 226.686867 ;
            T2 = Tm ;

            // fitted against data from Matula for values of T  > 298.15 K ;
            Vector< real > x = { 273.15 , 298.15, 800.0 };
            Vector< real > y = { 221.991, 221.0, 198.937 };

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

        void
        Silver::create_kohler()
        {
            mBSKohlerSwitch.set_size( 2 );
            mKohlerLongPolys.set_size( 3 );
            mKohlerTransPolys.set_size( 3 );
            Vector< real > & p0 = mKohlerLongPolys( 0 );
            Vector< real > & p1 = mKohlerLongPolys( 1 );
            Vector< real > & p2 = mKohlerLongPolys( 2 );
            Vector< real > & q0 = mKohlerTransPolys( 0 );
            Vector< real > & q1 = mKohlerTransPolys( 1 );
            Vector< real > & q2 = mKohlerTransPolys( 2 );

            // longitudinal field, fit to data from
            // Lüthi 1960, PhD thesis ETH Zürich, Switzerland
            // Strom-Olsen 1967, https://www.jstor.org/stable/2415889
            p1 = { -1.46572872E-01, + 2.56099423E+00, - 1.09120064E+01 };

            // transversal field, fit to data from
            // Stout 1939, 10.1021/ja01871a006
            // Strom-Olsen 1967, https://www.jstor.org/stable/2415889
            // Iwasa 1993, 10.1016/0011-2275(93)90199-X
            // Li 2015, 10.1088/1757-899X/102/1/012027
            q1 = {2.48185096E-03, - 1.37777811E-01, + 2.44232600E+00, - 1.02282979E+01  };

            // find point where the distance between both curves is extreme
            real a = 3 * q1(0);
            real b = 2 * ( q1(1) - p1(0));
            real c = q1(2) - p1(1);
            real w0 = ( -b - std::sqrt( b*b - 4*a*c ) ) / (2*a);

            // low field longitudinal
            real g      = polyval( p1, w0 );
            real dgdw   = dpolyval( p1, w0 ) ;
            real d2dgw2 = ddpolyval( p1, w0 ) ;
            real x = std::exp( w0 );
            real dgdx = dgdw / x ;
            real d2gdx2 = ( d2dgw2 - dgdw)/(x*x);


            real f = std::exp( g );
            real dfdx = f * dgdw / x ;
            real d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;

            q0.set_size( 4 );
            q0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            q0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            q0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            q0( 3 ) = 0.0 ;

            // low field transversal
            g = polyval( q1, w0 );
            dgdw = dpolyval( q1, w0 ) ;
            d2dgw2 = ddpolyval( q1, w0 ) ;
            f = std::exp( g );
            dfdx = f * dgdw / x ;
            d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;
            q0.set_size( 4 );
            q0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            q0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            q0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            q0( 3 ) = 0.0 ;

            // low field parallel
            g = polyval( p1, w0 );
            dgdw = dpolyval( p1, w0 ) ;
            d2dgw2 = ddpolyval( p1, w0 ) ;
            f = std::exp( g );
            dfdx = f * dgdw / x ;
            d2fdx2 = f * ( dgdx * dgdx + d2dgw2 ) ;
            p0.set_size( 4 );
            p0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            q0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            q0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            q0( 3 ) = 0.0 ;

            mBSKohlerSwitch( 0 ) = x;


            // compute transition point
            real w1 = -0.5*p1(1)/p1(0);

            mBSKohlerSwitch( 1 ) = std::exp( w1 );

            // longitudinal extension (constant)
            p2.set_size( 1, std::exp(polyval( p1, w1 ) ));
            q2.set_size( 2 );

            // transversal extension (linear)

            q2(0) = dpolyval( q1, w1 );
            q2(1) = polyval( q1, w1 ) - q2(0) * w1 ;

            this->set_kohler_dependencies();
        }

        real Silver::kohler( const real B, const real S, const real beta ) const
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
                Along  = std::exp( polyval( mKohlerLongPolys(1), w ) );
                Atrans =  std::exp( polyval( mKohlerTransPolys(1), w ) );
            }
            else
            {
                Along  =  mKohlerLongPolys(2)(0);
                Atrans = std::exp( polyval( mKohlerTransPolys(2), std::log( BS ) ) );
            }

            real c = std::cos( beta );
            real c2 = c*c ;
            real s2 = 1.0 - c2 ;

            // Pippard's angular interpolation formula
            return Along * c2 + Atrans * s2 ;
        }

        real
        Silver::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        real Silver::debye_custom( const real T ) const
        {
            if ( T <= 0 ) return this->constant_property( MaterialProperty::debye0K );

            if ( T < mTDebyeSwitch( 0 ) ) return polyval( mDebyePolys(0), T );
            if ( T < mTDebyeSwitch( 1 ) ) return polyval( mDebyePolys(1), T );
            if ( T < mTDebyeSwitch( 2 ) ) return polyval( mDebyePolys(2), T );

            return polyval( mDebyePolys( 3 ), T );
        }

    }


}
