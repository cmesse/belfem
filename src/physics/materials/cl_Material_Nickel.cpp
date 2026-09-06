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

#include "cl_Material_Nickel.hpp"
#include "cl_Material_Abundance.hpp"

namespace belfem
{
    namespace material
    {
        Nickel::Nickel( const real RRR, const bool aBuildTables ) :
            Ferromagnetic( "Nickel", aBuildTables )
        {
            this->create_magnetization_curve();
            this->set_constants();

            this->create_cp();
            this->create_alpha();
            this->create_debye_and_rho();
            this->create_kohler();

            // fitted against Touloukian dataset, Curves 5, 6 and reference
            this->set_lambda_coefficients(
{  0.0, 9.17539017e-7, 2.04951604,  1.08209240e3, -0.438026164,
             78.1833995,  1.77403741,  0.0  } );

            // Nickel needs two Wachtman curves and a transition
            this->create_young();

            // from Wolfram Cloud, skipping the Wachtman data
            this->create_mech( 0., 0., 0., 293.15, 0.31 );

            // must run after create_debye_and_rho(), which sets rho_i_ref —
            // with rho_i_ref unset, the regula falsi in set_RRR() exits
            // silently on NaN and leaves rho_0 as NaN
            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }
        }

        Nickel::~Nickel()
        {
            delete mReducedMagnetization;
            delete mThermalExpansion;
            if ( mKohlerLongBezier != nullptr ) delete mKohlerLongBezier ;
            if ( mYoungBezier != nullptr ) delete mYoungBezier ;
        }

        void
        Nickel::create_magnetization_curve()
        {
            // fitted against
            // J. Crangle; G. M. Goodman, 1971
            // The magnetization of pure iron and nickel
            // 10.1098/rspa.1971.0044
            mReducedMagnetization = new Bezier(
                { 0.0, 0.937409184, 1.0, 1.0},
                { 1.0, 1.0,0.485425638 , 0.0 } );
        }

        void Nickel::set_constants()
        {

            mSplineJn = this->create_J_spline( 4.5 );
            this->set_bloch_gruen_parameter( 4.5 );

            // Spin-disorder amplitude ( high-T limit ): 30 nΩ·m for nickel
            this->set_constant( MaterialProperty::A_spin_disorder, 3e-8 );

            // ferromagnetic Curie point of nickel. compute_mred() scales the
            // reduced magnetization curve by this value; the Bezier in
            // create_magnetization_curve() is calibrated against Crangle &
            // Goodman with this true Tc, verified at 600 K: m = 0.47 vs the
            // measured ~0.45-0.50
            this->set_constant( MaterialProperty::Tcurie, 631. );
            this->set_constant( MaterialProperty::T_max, 631 );

            this->set_angular_momentum( 0.5 );

            // White & Woods, 1959, Phil. Trans. R. Soc. A 251, 273, p. 288;
            // DOI: 10.1098/rsta.1959.0004
            this->set_constant( MaterialProperty::A_electron_magnon, 1.6e-13 );

            Abundance tAbundance ;

            real M = tAbundance.compute_molar_mass( "Ni" );

            this->set_constant( MaterialProperty::M, M );

            this->set_constant( MaterialProperty::ref_density, 8908. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );

            // residual resistivity of the ideal crystal; rho_custom and
            // lambda_custom read this, a real sample gets its value via set_RRR()
            this->set_constant( MaterialProperty::rho_0, 0. );


            // fitted against TPPM dataset
            this->set_constant( MaterialProperty::beta, 0.0003491874 );
            this->set_constant( MaterialProperty::gamma, 0.1271067 );

            // TPPM vol. 1, p. 237
            this->set_constant( MaterialProperty::rho0_pure, 0.0384e-8 );

        }

        void
        Nickel::create_alpha()
        {
            // mapped against Touloukian dataset
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 110.2485, 293.15, 600. };
            mThermalExpansion->basis_y() = { -0.24179003e-2,  -0.24179003e-2,-0.04203828e-2, 0.44788604e-2 };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }


        void
        Nickel::create_cp()
        {
            // fitted against Touloukian dataset
            Metal::create_cp( { 2.80562198609, 3.09174044467, 3.58502987240, 3.76983270768 },
                              { 1.28156678044, 1.80227683292, 3.36266505280, 3.88855007170 },
                              { 3.76983270768, 4.65779557577, 5.73266272111, 6.40258032700 },
                              { 3.88855007170, 6.41538559867, 5.85730691866, 6.39298025828 } );
        }

        void
        Nickel::create_debye_and_rho()
        {
            // calibrated against White and Woods, 1959, 0.1098/rsta.1959.0004, and Farrell and Greig, 1968 10.1088/0022-3719/1/5/326
            mDebyePoly = {- 1.57719731763856E-15 , 0, 1.21020850554333E-09, 0.,- 8.52190158963375E-04, 0., this->constant_property( MaterialProperty::debye0K ) };
            this->set_custom( MaterialProperty::debye );

            // value from0.1098/rsta.1959.0004
            real T = 295. ;
            real theta = polyval( mDebyePoly, T ); // ca. 390.
            real rho = 7.0e-8 - this->rho_mag( T );
            real n = this->constant_property( MaterialProperty::n_bloch_gruen ) ;
            real Z = theta / T ;
            real A = rho * std::pow( Z, n ) / this->Jn( Z );

            T = 273.15 ;
            theta = polyval( mDebyePoly, T );
            Z = theta / T ;
            rho = A * std::pow( Z, -n ) * this->Jn( Z );
            this->set_rho_i_ref( T, rho, theta );

            this->create_spline( MaterialProperty::rho_i, 0.0 );
            this->set_custom( MaterialProperty::rho );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho_i, MaterialDependency::T );
        }

        real
        Nickel::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        void
        Nickel::create_kohler()
        {

            // based on Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //                      Dissertation, ETH Zürich, 1960
            mKohlerTransPolys.set_size( 3 );
            Vector< real > & q1 = mKohlerTransPolys( 1 );
            q1 = { 8.350940372517E-02, 6.472680585001E-01, - 6.232441769620E+00 };

            mKohlerLongBezier = new Bezier(  {  3.000000000000, 5.901022203328, 6.079578875969, 6.237686705000 },
                      { -3.935015144441, -1.091459238608, -1.091459238608, -1.091459238608 } );

            // polynomial for low field area
            mBSKohlerSwitch.set_size( 2 );

            real w = mKohlerLongBezier->basis_x()(0);
            real x = std::exp( w );

            // transition points
            real w1 = mKohlerLongBezier->basis_x()(3);
            real x1 = std::exp( w1 );
            mBSKohlerSwitch( 0 ) = x ;
            mBSKohlerSwitch( 1 ) = x1; // field saturates here


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
            p1.set_size( 1, std::exp(mKohlerLongBezier->y( w1 )) );


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
            q2(0) = dpolyval( q1, w1 ) ;
            q2(1) = polyval( q1, w1 ) - q2(0)*w1 ;

            this->set_dependency( MaterialProperty::rho, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::normB);
            this->set_dependency( MaterialProperty::rho, MaterialDependency::angleBxJ );

            this->set_dependency( MaterialProperty::lambda, MaterialDependency::T );
            this->set_dependency( MaterialProperty::lambda, MaterialDependency::normB);
            this->set_dependency( MaterialProperty::lambda, MaterialDependency::angleBxJ );
        }

        real
        Nickel::kohler( const real B, const real S, const real beta ) const
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

        void Nickel::create_young()
        {
            mYoungData.set_size( 2, {} );

            Vector< real > & p0 = mYoungData( 0 );
            Vector< real > & p1 = mYoungData( 1 );

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Wolfram Cloud

            p0 = { 213.721e9, 0.51216e9, 474.73 };
            p1 = { 204.627e9, 0.17663e9, 1664.71 };



            Vector< real > x = { 465.4, 527.3, 631.8, 637.};


            real T = x(0);
            real E0 = p0( 0 );
            real  b = p0( 1 );
            real T0 = p0( 2 );

            real z = b * std::exp( - T0 / T ) ;
            real y0  = E0 - T * z ;
            real dy0 = -z*( T + T0 ) / T ;
            real y1 = y0 + dy0 * ( x( 1 ) - x( 0 ) );

            T = x(3);
            E0 = p1( 0 );
            b = p1( 1 );
            T0 = p1( 2 );

            z = b * std::exp( - T0 / T ) ;
            real y3  = E0 - T * z ;
            real dy3 = -z*( T + T0 ) / T ;

            real y2 = y3 - dy3 * ( x( 3 ) - x( 2 ) );

            Vector< real > y( { y0, y1, y2, y3 } );

            mYoungBezier = new Bezier( x, y );

            // the E spline itself is built by Metal::create_mech(), which also
            // supplies the 0 K slope
            this->set_have( MaterialProperty::E );
        }

    }
}
