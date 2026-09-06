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
#include "cl_Material_Aluminum.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        Aluminum::Aluminum( const real RRR, const bool aBuildTables ) :
            Metal( "Aluminum", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();

            this->create_cp();
            this->create_alpha();

            this->create_debye() ;
            this->set_rho_i_ref( 298.15, 2.71e-8 );
            this->set_rho_i_ref( 273.15, this->rho_i_custom( 273.15 ) );

            // Thermal conductivity coefficients from Hust 1984, 10.6028/nbs.ir.84-3007
            this->set_lambda_coefficients( {0., 4.716e-8, 2.446, 623.6, -0.16, 130.9, 2.5, 0.8168} );

            this->create_kohler();

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Wolfram Cloud
            this->create_mech( 75.706, 0.07071, 649.62, 293.15, 0.35 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR );
            }
        }

//------------------------------------------------------------------------------

        Aluminum::~Aluminum()
        {
            if ( mThermalExpansion != nullptr ) delete mThermalExpansion ;

            for ( Bezier * tBezier : mDebyeBeziers )
            {
                delete tBezier ;
            }
        }

//------------------------------------------------------------------------------

        void
        Aluminum::set_constants()
        {
            // melting temperature
            this->set_constant( MaterialProperty::T_max, 933.47 );

            // reference density at room temperature
            this->set_constant( MaterialProperty::ref_density, 2700. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );

            Abundance tAbundance ;
            real M = tAbundance.compute_molar_mass( "Al" );
            this->set_constant( MaterialProperty::M, M );

            // Sommerfeld and Debye coefficients of cp ~ cv = gamma * T + beta * T^3,
            // fitted against the Touloukian dataset. Must be set after M, because
            // writing beta derives debye0K from the molar mass - it comes out at
            // 417.9 K, against a literature Debye temperature of 428 K.
            this->set_constant( MaterialProperty::gamma, 5.39935764122e-2 );
            this->set_constant( MaterialProperty::beta,  9.87430801769e-4 );

            // TPPM vol. 1, p. 9
            this->set_constant( MaterialProperty::rho0_pure, 0.000593e-8 );
        }

//------------------------------------------------------------------------------

        void
        Aluminum::create_alpha()
        {
            // mapped against Touloukian dataset. dL/L is referred to 293.15 K and
            // is zero there exactly; basis_y(0) == basis_y(1) gives a zero
            // derivative at 0 K. The upper control point is the melting temperature.
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 184.2839196154, 686.6437602561, 933.45 };
            mThermalExpansion->basis_y() = { -0.44257783e-2, -0.44257783e-2,
                                              0.92163458e-2,  1.90023027e-2 };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }

//------------------------------------------------------------------------------

        void
        Aluminum::create_cp()
        {
            // fitted against Touloukian dataset. The upper control point sits at
            // 932.2 K, essentially the melting temperature. Residuals are below
            // 2 % from 20 K to 300 K and run about 5 % low above 400 K.
            Metal::create_cp( { 2.31871750933, 3.24713918683, 3.54296959574, 4.20022904346 },
                              { 0.47172176910, 2.49035257537, 4.35612664303, 5.57401591925 },
                              { 4.20022904346, 5.08051627224, 6.02946712178, 6.83750447000 },
                              { 5.57401591925, 7.20517143878, 6.60003178609, 7.10422754837 } );
        }

//------------------------------------------------------------------------------

        real
        Aluminum::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

//------------------------------------------------------------------------------

        void
        Aluminum::create_debye()
        {
            // fitted against Desai et al, 1984 doi.org/10.1063/1.555725
            //                and Cook et al, 1975, ORNL-5079

            mDebyeBeziers.set_size( 4, nullptr );

            real theta0 = this->constant_property( MaterialProperty::debye0K );

            mDebyeBeziers( 0 ) = new Bezier(  {0.000000, 7.000000, 11.000000, 20.000000}, {theta0, theta0, 383.438965, 383.438965 } );
            mDebyeBeziers( 1 ) = new Bezier( {20.000000, 31.913256, 6.281909, 80.600100},{383.438965, 383.438965, 403.615200, 403.615200} );
            mDebyeBeziers( 2 ) = new Bezier(  {80.600100, 131.561868, 361.672012, 400.000000}, {403.615200, 403.615200, 385.042453, 382.091610} );
            mDebyeBeziers( 3 ) = new Bezier({400.000000, 516.080221, 640.619091, 933.520000}, {382.091610, 373.154680, 383.730010, 331.351000} );

            this->set_have( MaterialProperty::debye );
            this->create_spline( MaterialProperty::debye );
        }

        real
        Aluminum::debye_custom( const real T ) const
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
        Aluminum::create_kohler()
        {
            mKohlerLongPolys.set_size( 3, {} );

            Vector< real > & p1 = mKohlerLongPolys( 1 );

            // based on Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //                      Dissertation, ETH Zürich, 1960
            p1 = { -9.017754E-02, 1.483167E+00, - 6.663454E+00 };

            real w = 3.0 ;
            real x = std::exp(w);

            // saturation value, continuous with the mid field branch, which
            // evaluates p1 in w — not in B*S
            real w1 = -0.5*p1(1)/p1(0);
            real x1 = std::exp(w1);
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

            // based on Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //                      Dissertation, ETH Zürich, 1960
            // and Fickett: Magnetoresistivity of Copper and Aluminum at Cryogenic Temperatures, 1972, FNAL
            q1 = { 1.100949E-02, - 3.229136E-01, + 3.443542E+00, - 1.138694E+01 };

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

        real Aluminum::kohler( const real B, const real S, const real beta ) const
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
