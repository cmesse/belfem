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

#include "cl_Material_Iron.hpp"
#include "cl_Material_Abundance.hpp"
#include "../../linalg/lapack/fn_gesv.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
        Iron::Iron( const real RRR, const bool aBuildTables ) :
            Ferromagnetic( "Iron", aBuildTables )
        {

            this->set_constants();
            this->create_magnetization_curve();

            this->create_debye_and_rho();
            this->create_cp();

            this->create_alpha();

            // from NBSIR 84-3007
            this->set_lambda_coefficients( {
                0., 166.9e-8, 1.868, 1.503e5, -1.22, 238.6, 1.392, 0. } );

            this->create_kohler();

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Wolfram Cloud
           this->create_mech( 214.886, 0.18823, 998.94, 293.15, 0.29 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }
        }

        Iron::~Iron()
        {
            delete mReducedMagnetization;
            delete mThermalExpansion;
            delete mDebyeBezier;
            if ( mKohlerTransBezier != nullptr ) delete mKohlerTransBezier ;
        }

        void
        Iron::create_magnetization_curve()
        {
            // fitted against
            // J. Crangle; G. M. Goodman, 1971
            // The magnetization of pure iron and nickel
            // 10.1098/rspa.1971.0044

            mReducedMagnetization = new Bezier(
                { 0.0, 0.998609292, 1.0, 1.0},
                { 1.0, 1.0,0.512733285, 0.0 } );
        }

        void
        Iron::create_cp()
        {
            // fitted against Touloukian dataset, smoothening Curie Anomaly
            Metal::create_cp( { 2.73540715545, 3.47358887137, 3.57594573420, 3.89691502694 },
                { 0.93884823153, 2.36862704577, 3.04957208559, 3.94857248473 },
                { 3.89691502694, 4.86719614149, 5.82580172107, 6.79458658100 },
                { 3.94857248473, 6.66622532259, 5.74833734860, 6.62345775080 },
                { 6.794586581, 7.021437538311806, 7.01602489, 7.12 },
                { 6.6234577508, 6.828376210215755, 6.60047746647, 6.6142704735 } );
        }

        void Iron::create_debye_and_rho()
        {
            mSplineJn = this->create_J_spline( 4.5 );

            real y0 = 470.0 ;

            // well known literature value; beta and gamma are set in
            // set_constants() — do not re-set them here, create_cp() reads them
            this->set_constant( MaterialProperty::debye0K, y0 );

            // values < 295 K fitted against data from
            // J. Crangle; G. M. Goodman, 1971
            // "The Magnetization of Pure Iron and Nickel"
            // 10.1098/rspa.1971.0044
            // cross checked against NBSIR 84-3007
            // fitted so that theta = 400 K @ T = 295 K

            // values > 295 K fitted against
            // Arays and Colvin 1964
            // "Electrical Resistivity of High Purity Iron from 300 to 1300 °K"
            // 10.1002/pssb.19640060317


            real T1 = 39.99983608937119 ;
            real y1 = y0 ;

            real T2 = 922.8395244075266 ;
            real y2 = 156.1427835433610 ;

            mDebyeBezier = new Bezier(
                { std::log(T1), 5.6262856757, 5.6262856757, std::log( T2 )},
                { y1, y1, 439.044319081386, y2 } );


            real dy2 = mDebyeBezier->dydx( std::log( T2 )  ) / T2 ;

            real f = mDebyeBezier->dydx( std::log( T2 )  ) ;
            real g = 1./T2 ;

            real df = mDebyeBezier->d2ydx2( std::log( T2 )  ) / T2 ;
            real dg = -g*g ;
            real ddy2 = f*dg + g*df ;

            Matrix< real > V( 4, 4 );

            real x0 = T2 * 0.001 ;
            real x1 = this->constant_property( MaterialProperty::Tcurie ) * 0.001 ;


            V( 0, 0 ) = x0*x0*x0 ;
            V( 0, 1 ) = x0*x0 ;
            V( 0, 2 ) = x0 ;
            V( 0, 3 ) = 1.0 ;

            V( 1, 0 ) = 3. * x0*x0 ;
            V( 1, 1 ) = 2. * x0 ;
            V( 1, 2 ) = 1.0 ;
            V( 1, 3 ) = 0.0 ;

            V( 2, 0 ) = 6. * x0 ;
            V( 2, 1 ) = 2. ;
            V( 2, 2 ) = 0.0 ;
            V( 2, 3 ) = 0.0 ;

            V( 3, 0 ) = 3. * x1*x1 ;
            V( 3, 1 ) = 2. * x1 ;
            V( 3, 2 ) = 1.0 ;
            V( 3, 3 ) = 0.0 ;

            mDebyePolys.set_size( 2, {} );
            Vector< real > & p = mDebyePolys( 0 ) ;
            p = { y2, dy2 * 1e3, ddy2 * 1e6, 0.0 };

            Vector< int_t > q( 4 );
            gesv( V, p, q );

            p( 0 ) *= 1.e-9 ;
            p( 1 ) *= 1.e-6 ;
            p( 2 ) *= 1.e-3 ;

            mTDebyeSwitch = { T1, T2, this->constant_property( MaterialProperty::Tcurie ) };
            mDebyePolys( 1 ) = { polyval( p, this->constant_property( MaterialProperty::Tcurie ) ) };

            this->set_custom( MaterialProperty::debye );


            // value from0.1098/rsta.1959.0004
            real T = 295. ;
            real theta = this->debye( T );
            real rho = 9.8e-8 - this->rho_mag( T );
            real n = this->constant_property( MaterialProperty::n_bloch_gruen ) ;
            real Z = theta / T ;
            real A = rho * std::pow( Z, n ) / this->Jn( Z );

            T = 273.15 ;
            theta = this->debye( T );
            Z = theta / T ;
            rho = A * std::pow( Z, -n ) * this->Jn( Z );
            this->set_rho_i_ref( T, rho, theta );
            this->create_spline( MaterialProperty::debye, 0.0 );

            this->create_spline( MaterialProperty::rho_i, 0.0 );
            this->set_custom( MaterialProperty::rho );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho_i, MaterialDependency::T );
        }

        void
        Iron::set_constants()
        {
            // deliberately using n = 4.5 ( Jn spline ), not the textbook 3 or 5
            this->set_bloch_gruen_parameter( 4.5 );

            // Spin-disorder amplitude ( high-T limit )
            // Typical value for iron: ~ 80–120 nΩ·m
            this->set_constant( MaterialProperty::A_spin_disorder, 1e-7 );
            this->set_constant( MaterialProperty::Tcurie, 1043. );

            // we cap the maximum temperature at 860 K because the resistivity is wrong above here
            this->set_constant( MaterialProperty::T_max, 860. );

            this->set_angular_momentum( 1. );

            // White & Woods, 1959, Phil. Trans. R. Soc. A 251, 273, p. 288;
            // derived from polycrystalline Fe specimens with RRR ≈ 40–104
            // DOI: 10.1098/rsta.1959.0004
            this->set_constant( MaterialProperty::A_electron_magnon, 1.3e-13 );

            Abundance tAbundance ;

            real M = tAbundance.compute_molar_mass( "Fe" );

            this->set_constant( MaterialProperty::M, M );

            this->set_constant( MaterialProperty::ref_density, 7870. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );

            // residual resistivity of the ideal crystal; rho_custom and
            // lambda_custom read this, a real sample gets its value via set_RRR()
            this->set_constant( MaterialProperty::rho_0, 0. );

            // fitted against Touloukian dataset
            this->set_constant( MaterialProperty::beta, 0.000353438392636 );
            this->set_constant( MaterialProperty::gamma, 0.089076406330672 );

            // TPPM vol. 1
            this->set_constant( MaterialProperty::rho0_pure, 0.0327e-8 );
        }

        real
        Iron::debye_custom( const real T ) const
        {
            if ( T < mTDebyeSwitch( 0 ) ) return this->constant_property( MaterialProperty::debye0K );
            if ( T < mTDebyeSwitch( 1 ) ) return mDebyeBezier->y( std::log( T ) ) ;
            if ( T < mTDebyeSwitch( 2 ) ) return polyval( mDebyePolys( 0 ), T );
            return mDebyePolys( 1 )( 0 );
        }

        real
        Iron::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        void
        Iron::create_alpha()
        {
            // mapped against Touloukian dataset
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 163.007071576, 293.15, 1185. };
            mThermalExpansion->basis_y() = { -0.202550622361, -0.202550622361, -0.092309804530, 1.345115642917 };
            mThermalExpansion->basis_y() *= 0.01 ;

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }

        void
        Iron::create_kohler()
        {
            // fitted against Klaffky and Coleman, 1974, 10.1103/PhysRevB.10.2915
            // and Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //            Dissertation, ETH Zürich, 1960

            mKohlerLongPolys.set_size( 3, {} );
            mKohlerTransPolys.set_size( 2, {} );

            // longitudinal branch: quadratic in w = ln( B*S ), saturating at
            // its vertex w1
            Vector< real > & p1 = mKohlerLongPolys( 1 );
            p1 = {-7.049243881986E-02, + 2.002713598121E+00, - 1.068149544139E+01 };

            real w = 3.0 ;
            real x = std::exp(w);

            // saturation value, continuous with the mid field branch, which
            // evaluates p1 in w — not in B*S
            real w1 = -0.5*p1(1)/p1(0);
            real x1 = std::exp(w1);
            mKohlerLongPolys( 2 ) = { std::exp( polyval( p1, w1 ) ) };

            // transversal branch: Bezier over ( w, ln( delta rho / rho ) ) with
            // a flat end ( y2 == y3 ) at the peak the dataset itself implies,
            // the vertex of its least-squares quadratic. Bezier::xi_by_x clamps
            // past the last control point, so the curve holds its peak from
            // w = 11.17 up to the longitudinal switch at w1 = 14.21.
            mKohlerTransBezier = new Bezier(
                { 3.0,       7.132373, 10.138093, 11.172802 },
                { -2.324512, 4.549139,  4.557684,  4.557684 } );

            // transversal saturation
            mKohlerTransPolys( 1 ) = { std::exp( mKohlerTransBezier->basis_y()( 3 ) ) };

            // compute longitudinal polynomial for low field
            // let f = exp(g) and w = ln(x)
            real g      = polyval( p1, w );
            real dgdw   = dpolyval(p1, w );
            real d2dgw2 = ddpolyval( p1, w );
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

            // compute transversal polynomial for low field, from the Bezier
            g      = mKohlerTransBezier->y( w );
            dgdw   = mKohlerTransBezier->dydx( w );
            d2dgw2 = mKohlerTransBezier->d2ydx2( w );
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

            mBSKohlerSwitch = { x, x1 };

            this->set_kohler_dependencies();
        }

        real Iron::kohler( const real B, const real S, const real beta ) const
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
                Along  = std::exp( polyval( mKohlerLongPolys( 1 ), w ) );

                // clamps at the peak for w > 11.17, see create_kohler()
                Atrans = std::exp( mKohlerTransBezier->y( w ) );
            }
            else
            {
                Along  = mKohlerLongPolys(2)(0);
                Atrans = mKohlerTransPolys(1)(0);
            }

            real c = std::cos( beta );
            real c2 = c*c ;
            real s2 = 1.0 - c2 ;

            // Pippard's angular interpolation formula
            return Along * c2 + Atrans * s2 ;
        }

    }
}
