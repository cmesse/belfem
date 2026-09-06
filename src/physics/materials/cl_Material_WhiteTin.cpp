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


#include "cl_Material_WhiteTin.hpp"
#include "cl_Material_Abundance.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_create_fifth_order_beam_poly.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "fn_hust.hpp"
#include "fn_cardano.hpp"

namespace belfem
{
    namespace material
    {
        WhiteTin::WhiteTin( const real RRR, const bool aBuildTables ) :
            Metal( "WhiteTin", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();

            this->create_cp();
            this->create_alpha();

            this->create_debye();

            this->create_rho();

            // fitted against recommended values from TPRC Vol 1, p. 408
            this->set_lambda_coefficients( {
                0., 9.7309e-07 , 2.7040, 170.25,  -0.1641, 23.742, 5.5249, 0.
            } );

            this->create_kohler();

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Wolfram Cloud
            this->create_mech( 63.532, 0.21789,600.08, 293.15, 0.36 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }
        }

        WhiteTin::~WhiteTin()
        {
            if ( mThermalExpansion != nullptr ) delete mThermalExpansion ;
            if ( mDebyeTemperature != nullptr ) delete mDebyeTemperature ;
        }

        void WhiteTin::set_constants()
        {
            // set the melting temperature
            this->set_constant( MaterialProperty::T_max, 505.08 );

            this->set_constant( MaterialProperty::ref_density, 7289. );
            this->set_constant( MaterialProperty::T_ref_density, 293.15 );

            Abundance tAbundance ;
            this->set_constant( MaterialProperty::M, tAbundance.compute_molar_mass( "Sn" ) );


            real M = this->constant_property( MaterialProperty::M ) ;
            real gamma = 1.78e-3 / M ; // J/(kg K²) // O'Neil et al, 1964 10.1103/PhysRev.137.A748, Fig. 10
            real beta = 0.246e-3 / M ; // J/(kg K⁴) // O'Neil et al, 1964 10.1103/PhysRev.137.A748, Fig. 10

            this->set_constant( MaterialProperty::gamma, gamma );
            this->set_constant( MaterialProperty::beta, beta );
            real R = constant::Rm / M ;

            real theta0 = std::pow( ( 12. * std::pow( constant::pi, 4 ) * R  / ( 5. * beta ) ) , 1./3. ) ;
            this->set_constant( MaterialProperty::debye0K, theta0);

            // referene value from White: "Experimental Techniques in Low Temperature Physics", 1968
            //                     via Hariharan,1979 10.1007/BF02872130
            this->set_rho_i_ref( 273.15, 10.1e-8, 160. );

            // TPPM vol. 1
            this->set_constant( MaterialProperty::rho0_pure, 0.000133e-8  );

        }

        void
        WhiteTin::create_alpha()
        {
            // mapped against Touloukian dataset
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 50, 293.15, 600. };
            mThermalExpansion->basis_y() = { -0.47890123e-2,  -0.47890123e-2, -0.12758851e-2, 0.80820449e-2 };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }

        void WhiteTin::create_cp()
        {
            // fitted against Touloukian dataset
            Metal::create_cp( { 2.94566860085, 3.38256637736, 3.39712979715, 3.52819922546 },
                             {3.57415211848, 4.33697081305, 4.37073814112, 4.50872681872  },
                             { 3.52819922546, 4.53106921723, 5.25024866485, 6.22455842900  },
                             { 4.50872681872,  5.56453903787, 5.24241017537,  5.54850299467 } );

        }

        void WhiteTin::create_debye()
        {
            mDebyePolys.set_size( 5, {} );
            mTDebyeSwitch.set_size( 4, 0.0 );

            real & T0 = mTDebyeSwitch( 0 );
            real & T1 = mTDebyeSwitch( 1 );
            real & T2 = mTDebyeSwitch( 2 );
            real & T3 = mTDebyeSwitch( 3 );

            Vector< real > & p0 = mDebyePolys(0);
            Vector< real > & p1 = mDebyePolys(1);
            //Vector< real > & p2 = mDebyePolys(2);
            Vector< real > & p3 = mDebyePolys(3);
            Vector< real > & p4 = mDebyePolys(4);

            // fitted against cp data from TPRC
            p1 = { 1.1196e-2, -3.8340e-1, 5.0489, - 3.0290e1, 1.9914e2 };

            // continuity for low temperatures
            T0 = 4.0 ;
            real y0 = polyval( p1, T0 );
            real dy0 = 0.5*dpolyval( p1, T0 ) / std::sqrt( T0 ) ;

            create_beam_poly( 0., this->constant_property( MaterialProperty::debye0K  ), 0.0,
                              T0 * T0, y0, dy0, p0 );

            // find local minumum
            real a =  4.*p1(0);
            real b =  3.*p1(1);
            real c =  2*p1(2);
            real d =    p1(3);

            Vector< real > tT ;
            cardano( { a, b, c, d }, tT ) ;

            T1 = tT( 0 );
            real y1 = polyval( p1, T1 );


            // fitted against data from Meaden, "Electrical Resistance of Metals", 1965
            //                          via Hall, NBS TN 365, 1968
            p3 = { 2.1865e-6, -1.2669e-3, 1.49358e-1, 0.0 };

            a = 3. * p3(0);
            b = 2. * p3(1);
            c = p3(2);
            d = std::sqrt( b*b-4*a*c );
            real x1 = ( -b + d ) / ( 2*a );
            real x2 = ( -b - d ) / ( 2*a );

            T2 = std::min( x1, x2 );

            T3 = this->constant_property( MaterialProperty::T_ref_rho_i );
            p3( 3 ) = this->constant_property( MaterialProperty::debye ) - polyval( p3, T3 );

            real y2 = polyval( p3, T2 );

            mDebyeTemperature = new Bezier( );

            mDebyeTemperature->basis_x() = { T1 , 15.5 , 15.5, T2 };
            mDebyeTemperature->basis_y() = { y1, y1, y2, y2 };

            a = dpolyval( p3, T3 );
            b = this->constant_property( MaterialProperty::debye )  - a * T3 ;

            p4 = { a, b };

            this->set_custom( MaterialProperty::debye );
        }

        real
        WhiteTin::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        real
        WhiteTin::debye_custom( const real T ) const
        {
            if ( T < mTDebyeSwitch( 0 ) ) return polyval( mDebyePolys( 0 ), T * T );
            if ( T < mTDebyeSwitch( 1 ) ) return polyval( mDebyePolys( 1 ), T );
            if ( T < mTDebyeSwitch( 2 ) ) return mDebyeTemperature->y( T );
            if ( T < mTDebyeSwitch( 3 ) ) return polyval( mDebyePolys( 3 ), T );
            return polyval( mDebyePolys( 4 ), T );
        }

        void
        WhiteTin::create_kohler()
        {
            mKohlerLongPolys.set_size( 3, {} );
            mBSKohlerLongSwitch.set_size( 3 );
            mKohlerTransPolys.set_size( 2, {} );
            mBSKohlerTransSwitch.set_size( 1 );

            Vector< real > & p0 = mKohlerLongPolys( 0 );
            Vector< real > & p1 = mKohlerLongPolys(1);
            Vector< real > & p2 = mKohlerLongPolys( 2 );

            Vector< real > & q0 = mKohlerTransPolys( 0 );
            Vector< real > & q1 = mKohlerTransPolys(1);

            // based on data Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //                      Dissertation, ETH Zürich, 1960

            p1 = { 1.4158e-2, - 4.4977e-1, + 4.5121, - 1.3298e1 };
            q1 = { 1.8073, - 1.2545e1 };

            // compute low field polynomial
            real w1 = 2.5; // w = ln(x) w' = 1/x w'' = -1/x^2
            real x = std::exp( w1 );

            real g      = polyval( p1, w1 );
            real dgdw   = dpolyval( p1, w1 );
            real d2dgw2 = ddpolyval( p1, w1 );
            real dgdx = dgdw / x ;
            real d2gdx2 = ( d2dgw2 - dgdw)/(x*x);

            real f = std::exp( g );
            real dfdx = f * dgdw / x ;
            real d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;


            p0.set_size( 4 );
            p0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            p0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            p0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            p0( 3 ) = 0.0 ;

            // compute high field saturation

            real a = 3. * p1(0);
            real b = 2. * p1(1);
            real c = p1(2);
            real d = std::sqrt( b*b-4*a*c );

            real x1 = ( -b + d ) / ( 2*a );
            real x2 = ( -b - d ) / ( 2*a );

            real w2 = std::min( x1, x2 );
            p2 = { std::exp( polyval( p1, w2 ) ) };

            mBSKohlerLongSwitch( 0 ) = std::exp(w1);
            mBSKohlerLongSwitch( 1 ) = std::exp(w2);

            real w = 8.4 ;
                 x = std::exp( w );

            // compute transversal polynomial for low field
            g = polyval( q1, w );
            dgdw = dpolyval( q1, w );
            d2dgw2 = ddpolyval( q1, w );
            dgdx = dgdw / x ;
            d2gdx2 = ( d2dgw2 - dgdw)/(x*x);

            f = std::exp( g );
            dfdx = f * dgdw / x ;
            d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;
            q0.set_size( 4 );
            q0( 0 ) = ( f  + x * ( 0.5 * d2fdx2*x- dfdx) ) / ( x*x*x );
            q0( 1 ) = ( x*(3.*dfdx-d2fdx2*x) - 3.*f ) / ( x*x );
            q0( 2 ) = ( 3.*f + x * ( 0.5 * d2fdx2 * x - 2. * dfdx ) ) / x ;
            q0( 3 ) = 0.0 ;

            this->set_kohler_dependencies();
        }
    }
}
