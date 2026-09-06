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

#include <cmath>

#include "assert.hpp"
#include "constants.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_Material_Abundance.hpp"
#include "cl_Material_Chromium.hpp"


namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        Chromium::Chromium( const real RRR, const bool aBuildTables ) :
            Metal( "Chromium", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();
            this->create_cp();
            this->create_alpha();
            this->create_debye();

            // White and Woods, 1959, 10.1098/rsta.1959.0004
            this->set_rho_i_ref( 295, 12.9e-8 );
            this->set_rho_i_ref( 273.15, this->rho_i_custom( 273.15 ) );

            // fitted against Touloukian dataset, see also 10.1080/14786435708242698
            this->set_lambda_coefficients( { 0.0,
                4.054361011e-06,
                1.659582031,
                178.028422,
                -0.1474303106,
                89.58673092,
                3.820930427,
                0.0});

            this->create_kohler();

            // dynamic modulus taken from paper, scaled to the literature
            // 279 GPa at 293.15 K, ignoring the hump at room temperature

            // P. E. Armstrong and H. L. Brown
            // "Dynamic Young's Modulus Measurements Above 1000°C
            //  on Some Pure Polycrystalline Metals and Commercial Graphites"
            // Trans. AIME, 1964

            // poisson value at room temperature from Wolfram Cloud
            this->create_mech(279.837, 0.08741, 900.83, 293.15, 0.21 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR );
            }
        }

//------------------------------------------------------------------------------

        Chromium::~Chromium()
        {
            if ( mThermalExpansion != nullptr ) delete mThermalExpansion ;

            for ( Bezier * tBezier : mDebyeBeziers )
            {
                delete tBezier ;
            }
        }

//------------------------------------------------------------------------------

        void
        Chromium::set_constants()
        {
            // melting temperature
            this->set_constant( MaterialProperty::T_max, 2180. );

            // reference density at room temperature
            this->set_constant( MaterialProperty::ref_density, 7150. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );

            Abundance tAbundance ;
            real M = tAbundance.compute_molar_mass( "Cr" );
            this->set_constant( MaterialProperty::M, M );

            // Sommerfeld and Debye coefficients of cp ~ cv = gamma * T + beta * T^3.
            // Must be set after M, because writing beta derives debye0K from the
            // molar mass - it comes out at 606.0 K, at the lower edge of the
            // literature spread of 606 - 630 K. The Sommerfeld coefficient
            // corresponds to 1.74 mJ/(mol K^2), inside the spread of 1.40 - 1.78.
            this->set_constant( MaterialProperty::gamma,  0.0335375739423 );
            this->set_constant( MaterialProperty::beta,  0.000167979496326 );

            this->set_constant( MaterialProperty::n_bloch_gruen, 5.0 );

            // TPPM vol. 1, p. 60
            this->set_constant( MaterialProperty::rho0_pure, 0.0609e-8 );
        }
//------------------------------------------------------------------------------

        void
        Chromium::create_alpha()
        {
            // dL/L referred to 293.15 K, fitted primarily against
            // 10.1051/e3sconf/202459202009, with the Touloukian dataset as the
            // secondary source. dL/L( 293.15 K ) = 0 holds exactly, and
            // basis_y(0) == basis_y(1) gives a zero derivative at 0 K.
            //
            // Below the split temperature this curve is NOT used - alpha comes
            // from the cryogenic branch instead. That matters here more than for
            // any other material: a cubic Bezier with basis_y(0) == basis_y(1) is
            // linear in T near the origin, whereas alpha must follow cp and fall
            // off far more steeply. Measured against the fitted cp, this curve
            // overstates alpha by 2.1x at 100 K, 4.1x at 77 K and 12.4x at 50 K.
            // Chromium has the highest Debye temperature in the roster, so the
            // defect reaches further up in temperature than anywhere else.
            //
            // The anomalies at the Neel point ( 311 K ) and the spin flip
            // ( 123 K ) are not represented - see the class documentation.
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 100., 1024.15381264, 1700. };
            mThermalExpansion->basis_y() = { -0.13338499e-2, -0.13338499e-2,
                                              0.27598441e-2,  1.72437789e-2 };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }

//------------------------------------------------------------------------------

        void
        Chromium::create_cp()
        {
            // The Sommerfeld-Debye form behind Metal::create_cp() has no magnetic
            // channel, so it does not carry the lambda anomaly at the Neel point.
            // cp( 293 K ) comes out at 451.0 J/(kg K) against a literature 449.
            // The third segment carries the curve to 2103 K; the linear tail
            // covers the last stretch to the melting point at 2180 K. The
            // expansion Bezier ends at 1700 K and saturates above it
            // ( Bezier::xi_by_x clamps ), so alpha is held constant from there up.
            Metal::create_cp( { 2.64617543054, 3.18408331589, 3.84458683165, 4.43839211142 },
                              { -0.04417315023, 0.90468784886, 3.80265419805, 4.98362873296 },
                              { 4.43839211142, 4.97747978397, 5.68163244654, 6.39404801601 },
                              { 4.98362873296, 6.05577954580, 6.15040781077, 6.24614644164 },
                              { 6.39404801601, 6.48847862122, 7.03352159527, 7.65112017600 },
                              { 6.24614644164, 6.25883658661, 6.33208269419, 7.02635078400 });
        }

//------------------------------------------------------------------------------

        real
        Chromium::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

//------------------------------------------------------------------------------


        void Chromium::create_debye()
        {
            // for T > 140 K, fitted against  White and Woods, 1959, 10.1098/rsta.1959.0004,
            //                and Anderson, Stewart, Ramsay, 1970, 10.1002/pssb.19700370137

            mDebyeBeziers.set_size( 5, nullptr );

            real theta0 = this->constant_property( MaterialProperty::debye0K );


            mDebyeBeziers( 0 ) = new Bezier( {0, 9.5179, 10.2759, 16.8048}, {theta0, theta0, 332.8282, 332.8282} );
            mDebyeBeziers( 1 ) = new Bezier( {16.8048, 20.5192, 25.2865, 89.9985}, {332.8282, 332.8282, 498.1866, 498.1866} );
            mDebyeBeziers( 2 ) = new Bezier( {89.9985, 112.7902, 172.3140, 250.5671}, {498.1866, 498.1866, 463.4315, 463.4315} );
            mDebyeBeziers( 3 ) = new Bezier( {250.5671, 369.5911, 198.0318, 623.8304}, {463.4315, 463.4315, 612.9950, 612.9950} );
            mDebyeBeziers( 4 ) = new Bezier( {623.8304,  971.1774,  2102.1915,  2180.0000}, {612.9950,  612.9950, 417.4560, 377.4521} );

            this->set_have( MaterialProperty::debye );
            this->create_spline( MaterialProperty::debye );
        }

        real
        Chromium::debye_custom( const real T ) const
        {
            if ( T < BELFEM_EPSILON ) return this->constant_property( MaterialProperty::debye0K ) ;

            for ( Bezier * tBezier : mDebyeBeziers )
            {
                if ( T < tBezier->basis_x()(3) )
                {
                    return tBezier->y( T );
                }
            }
            return mDebyeBeziers.last()->basis_y()(3);
        }

        void
        Chromium::create_kohler()
        {
            mKohlerLongPolys.set_size( 3, {} );

            Vector< real > & p1 = mKohlerLongPolys( 1 );

            // Kozlova and Kondorskii: Dependence of the Electric and Magnetic Properties
            // of Chromium on the Magnetic Field Strength and Temperature,
            // Soviet Physics JETP, 1963
            p1 = { -9.449966e-2, + 1.211754, - 3.719104, - 2.299049 };

            real a = 3.*p1(0);
            real b = 2.*p1(1);
            real c = p1(2);
            real w1 = (-b - std::sqrt(b*b-4.*a*c))/(2.*a);
            real x1 = std::exp(w1);

            real w = -b/(2.*a);     // ca. 4.2
            real x = std::exp(w);

            Vector< real > & p2 = mKohlerLongPolys( 2 );
            p2 = { std::exp( polyval( p1, w1 ) ) };

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

            mKohlerTransPolys.set_size( 3, {} );

            Vector< real > & q1 = mKohlerTransPolys( 1 );

            // Arajs and Dunmyre, 1965, 10.1063/1.1703039
            q1 = { -2.198538e-2,  + 9.910930e-2, + 1.630813, - 7.555265 };

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

            mBSKohlerSwitch = { x, x1 };

            this->set_kohler_dependencies();
        }

        real
        Chromium::kohler( const real B, const real S, const real beta ) const
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
    }
}
