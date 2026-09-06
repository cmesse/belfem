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

#include "cl_Material_Lead.hpp"
#include "cl_Material_Abundance.hpp"
#include "fn_dpolyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_cardano.hpp"
#include "fn_polyfit.hpp"

namespace belfem
{
    namespace material
    {
        Lead::Lead( const real RRR, const bool aBuildTables ) :
            Metal( "Lead", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();
            this->create_cp();

            this->create_alpha();

            this->create_debye();


            // we need to cap the spline, the values become unrealistic otherwise
            real Tmax = this->constant_property( MaterialProperty::T_max );
            this->set_constant( MaterialProperty::T_max, 273.15 );
            this->create_rho();
            this->set_constant( MaterialProperty::T_max, Tmax );

            // fitted against recommended values from TPRC Vol 1, p. 191
            this->set_lambda_coefficients( {
                0., 8.0921e-6 , 2.7644, 51.31,  -6.8696e-2, 15.663, 4.6127, 0.
            } );

            this->create_kohler();

            // Wachtman data fit against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // poisson ratio at room temperature from Wolfram Cloud
            this->create_mech( 19.011, 0.03196, 383.6, 293.15, 0.44 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }
        }

        Lead::~Lead()
        {
            if ( mThermalExpansion != nullptr ) delete mThermalExpansion ;
            if ( mDebyeTemperature != nullptr ) delete mDebyeTemperature ;
        }

        void
        Lead::set_constants()
        {
            // set the melting temperature
            this->set_constant( MaterialProperty::T_max, 600.61 );

            this->set_constant( MaterialProperty::ref_density, 11340. );
            this->set_constant( MaterialProperty::T_ref_density, 293.15 );

            Abundance tAbundance ;
            this->set_constant( MaterialProperty::M, tAbundance.compute_molar_mass( "Pb" ) );

            // extracted from touloukian dataset TPRC
            real gamma = 1.4718e-2 ;
            real beta  = 8.3156e-3 ;

            this->set_constant( MaterialProperty::gamma, gamma );
            this->set_constant( MaterialProperty::beta, beta );


            real M = this->constant_property( MaterialProperty::M ) ;

            real R = constant::Rm / M ;

            real theta0 = std::pow( ( 12. * std::pow( constant::pi, 4 ) * R  / ( 5. * beta ) ) , 1./3. ) ;
            this->set_constant( MaterialProperty::debye0K, theta0 );
            this->set_constant( MaterialProperty::rho_0, 0 );

            // referene value from White: "Experimental Techniques in Low Temperature Physics", 1968
            //                     via Hariharan,1979 10.1007/BF02872130
            this->set_rho_i_ref( 273.15, 19.3e-8, 90);

            // TPPM vol. 1, p. 237
            this->set_constant( MaterialProperty::rho0_pure, 0.00088e-8 );

        }

        void
        Lead::create_alpha()
        {
            // mapped against Touloukian dataset
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 50, 293.15, 600. };
            mThermalExpansion->basis_y() = { -0.7134952e-2,  -0.71349527e-2, -0.06237054e-2, 0.98045527e-2 };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }


        void
        Lead::create_cp()
        {
            // fitted against Touloukian dataset
            Metal::create_cp( { 2.02641053843, 2.92519047395, 3.74214498125, 4.85578538678 },
                             { 1.92844531047, 4.61831528168, 4.70632666703, 4.79183717770 },
                             { 4.85578538678, 5.94375160754, 5.97906703086, 6.38519439900 },
                             { 4.79183717770, 4.87537630456, 4.85395845489, 4.96559014646 } );
        }

        void
        Lead::create_debye()
        {
            mDebyePolys.set_size( 4, {} );
            mTDebyeSwitch.set_size( 3, 0.0 );

            real & T0 = mTDebyeSwitch( 0 );
            real & T1 = mTDebyeSwitch( 1 );
            real & T2 = mTDebyeSwitch( 2 );

            Vector< real > & p0 = mDebyePolys(0);
            //Vector< real > & p1 = mDebyePolys(1);
            Vector< real > & p2 = mDebyePolys(2);
            Vector< real > & p3 = mDebyePolys(3);

            // we need to be careful that the spline does not become negative
            // this value has been found by trial and error
            T0 = 4.0 ;
            real y0 = 98.5 ;

            // we grab the data for 12 K from the cp function
            T1 = 12.0 ;
            real y1 = this->compute_debye_from_cv( T1, this->cp_custom( T1 ) ) ;

            real x0 = T0*T0 ;
            real x1 = T1*T1 ;

            Vector< real > x = { 0.0, x0, x1 };
            Vector< real > y = {this->constant_property( MaterialProperty::debye0K ), y0, y1 };

            polyfit( x, y, 2, p0 );

            x0 = -0.5 * p0(1)/p0(0) ;
            T0 = std::sqrt(  x0 );
            y0 = polyval( p0, x0 );

            // fitted against data from Meaden, "Electrical Resistance of Metals", 1965
            //                          via Hall, NBS TN 365, 1968
            p2 = { -6.3504e-9, 5.1659e-6, -1.3904e-3, 1.1028e-1, 0. } ;

            real a = 4. * p2(0);
            real b = 3. * p2(1);
            real c = 2. * p2(2);
            real d = p2(3);

            Vector< real > tT ;
            cardano( { a, b, c, d }, tT ) ;
            T1 = tT( 0 );
            T2 = this->constant_property( MaterialProperty::T_ref_rho_i );
            p2( 4 ) = this->constant_property( MaterialProperty::debye ) - polyval( p2, T2 );

            y1 = polyval( p2, T1 );

            mDebyeTemperature = new Bezier( );

            mDebyeTemperature->basis_x() = { T0 , 22.5 , 22.5, T1 };

            mDebyeTemperature->basis_y() = { y0, y0, y1, y1 };

            a = dpolyval( p2, T2 );
            b = this->constant_property( MaterialProperty::debye ) - a * T2 ;
            p3 = { a, b };

            this->set_custom( MaterialProperty::debye );
        }

        void
        Lead::create_kohler()
        {
            mKohlerLongPolys.set_size( 3, {} );
            mBSKohlerLongSwitch.set_size( 2 );
            mKohlerTransPolys.set_size( 2, {} );
            mBSKohlerTransSwitch.set_size( 1 );

            Vector< real > & p0 = mKohlerLongPolys( 0 );
            Vector< real > & p1 = mKohlerLongPolys( 1 );
            Vector< real > & p2 = mKohlerLongPolys( 2 );

            Vector< real > & q0 = mKohlerTransPolys( 0 );
            Vector< real > & q1 = mKohlerTransPolys(1);

            // based on data Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //                      Dissertation, ETH Zürich, 1960
            p1 = { -5.0939e-2, 1.273, - 7.0351 };
            q1 = { 1.458, - 9.7936 };

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
            real w2 = -0.5*p1(1)/p1(0) ;

            p2 = { std::exp( polyval( p1, w2 ) ) };

            mBSKohlerLongSwitch( 0 ) = std::exp(w1);
            mBSKohlerLongSwitch( 1 ) = std::exp(w2);

            real w = 3.0 ;
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

        real
        Lead::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        real
        Lead::debye_custom( const real T ) const
        {
            if ( T < mTDebyeSwitch( 0 ) ) return polyval( mDebyePolys( 0 ), T*T );
            if ( T < mTDebyeSwitch( 1 ) ) return mDebyeTemperature->y( T );
            if ( T < mTDebyeSwitch( 2 ) ) return polyval( mDebyePolys( 2 ), T );
            return polyval( mDebyePolys( 3 ), T );
        }

    }
}