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
#include "cl_Material_Magnesia.hpp"

#include "fn_cardano.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace material
    {
        Magnesia::Magnesia() :
            SplineLookupTable( MaterialType::NonMetal, true )
        {
            // Set the material label
            this->set_label( "Magnesia" );

            this->set_constants();

            this->create_cp();
            this->create_mech();           // K( T ) is needed by the cryogenic alpha branch
            this->create_alpha();
            this->create_lambda();

        }

        Magnesia::~Magnesia()
        {
            if ( mThermalExpansion != nullptr ) delete mThermalExpansion ;
        }

        void
        Magnesia::set_constants()
        {
            // The maximum service temperature of magnesia is actually around 2000 K,
            // but we don't expect to need the full band
            this->set_constant( MaterialProperty::T_max, 600.0 );

            this->set_constant( MaterialProperty::ref_density, 3580. );
            this->set_constant( MaterialProperty::T_ref_density, 293.15 );

            // MgO: molar mass and two atoms per formula unit, so that beta below
            // derives the 0 K Debye temperature ( 939 K, literature 940 - 946 K )
            this->set_constant( MaterialProperty::M, 0.0403044 );
            this->set_constant( MaterialProperty::q, 2.0 );

            // fitted against curve 4 of Touloukian TPRC dataset
            this->set_constant( MaterialProperty::beta, 1.1674e-04 );

            // a very small value for gamma that is > 0 prevents cp from becoming negative at very low temperatures
            // the value is chosen so that we get a minimal error for 0.1 K
            this->set_constant( MaterialProperty::gamma, 0.5426e-6 );
        }

        void
        Magnesia::create_alpha()
        {
            // mapped against Simon, 1994 : NISTIR-5030 and Durand, 1936: 10.1063/1.1745396
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 50, 293.15, 600. };
            mThermalExpansion->basis_y() = { -0.15166635e-2, -0.15166635e-2, 0., 0.4507889e-2,  };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            // Below the split temperature alpha follows cp ( Grueneisen ). The
            // Bezier is anchored to data at room temperature but too flat below
            // ~200 K ( its alpha rises slower than cp everywhere, dln(C)/dT < 0 ),
            // so only its VALUE at the split is used; the derivatives come from
            // the Grueneisen relation itself, see create_cryo_expansion_anchored()
            real T = this->set_alpha_switch_temperature();
            this->create_cryo_expansion_anchored(
                SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T ),
                mThermalExpansionCryo );

        }

        void
        Magnesia::create_mech()
        {
            // polynomial for poisson, fitted against Simon, 1994 : NISTIR-5030
            // ν = a*T^4 + b * T^2 + c

            real a = -4.5534e-13 ;
            real b = 1.6122e-7 ;
            real c = 2.2444e-1 ;

            mPoissonPolys.set_size( 2, {} );
            mPoissonPolys( 0 ) = { a, b, c };

            // find point where curvature is zero (ca. 242.92 K)
            real T  = std::sqrt( -b / ( 6. * a ) );

            real f  =  T*T*( a *T*T +  b ) + c ;
            real df = T*( 4. * a *T*T +  2. * b );

            // continuity
            real d = f - df * T ;
            mPoissonPolys( 1 ) = { df, d };
            mTPoissonSwitch = T;
            mYoungPolys.set_size( 2, {} );

            // polynomial for poisson, fitted against Blanke, Thermophysikalische Stoffgrößen, Springer 1989
            // note: Simon's data seem too low for HTS tapes
            // E = a* + b ;
            a = -4.2695E+07 ;
            b = 3.1195E+11 ;

            // we just use the same temperature. Values for very high temperatures, however, might be off by a lot
            mTYoungSwitch = T ;
            mYoungPolys( 1 ) = { a, b };
            f = a * T + b ;
            df = a ;

            // now we create a continuous polynomial that has dE/dT|T=0 = 0
            a = df / ( 2. * T );
            b = f - 0.5 * T * df ;
            mYoungPolys( 0 ) = { a, b };

            this->set_custom( MaterialProperty::E );
            this->set_custom( MaterialProperty::nu );
        }

        void
        Magnesia::create_cp()
        {

            Vector< real > theta = { 2.0, 3.0, 4.0, 5.0, 6.0, 6.5 };
            uint n = theta.length();
            mTCpSwitch.set_size( n );
            for ( uint t=0; t<n; ++t )
            {
                mTCpSwitch( t ) = std::exp( theta( t ) );
            }

            mCpPolys.set_size( 6, {} );
            // fitted against Touloukian TPRC and  imon, 1994 : NISTIR-5030

            // 2 < theta < 4.5
            mCpPolys( 1 ) = { 1.5577e-1, + 2.2610, - 8.1909};

            // 4.5 < theta < 6.5
            mCpPolys( 3 ) = { 2.6460e-1, - 5.0056, + 3.1799e1, - 6.0778e1 };

            // 6 < theta < 7.2
            mCpPolys( 5 ) = { 7.3150e-2, - 1.5197, + 1.0641e1, - 1.7963e1 };

            // poly 0
            real gamma = this->constant_property( MaterialProperty::gamma );
            real beta = this->constant_property( MaterialProperty::beta );

            real T0 = mTCpSwitch( 0 );
            real x0 = theta( 0 );
            real x1 = theta( 1 );

            real cp0 = ( gamma + beta * T0 * T0 ) * T0 ;
            real y0  = std::log( cp0 );
            real dydx0 = ( gamma  + 3. * beta * T0 * T0 ) * T0 / cp0 ;

            real y1 = polyval( mCpPolys(1), x1 );
            real dydx1 = dpolyval( mCpPolys(1), x1 );

            create_beam_poly( x0, y0, dydx0, x1, y1, dydx1, mCpPolys(0) );

            // poly 2
            x0 = theta( 2 );
            y0 = polyval( mCpPolys(1), x0 );
            dydx0 = dpolyval( mCpPolys(1), x0 );

            x1 = theta( 3 );
            y1 = polyval( mCpPolys(3), x1 );
            dydx1 = dpolyval( mCpPolys(3), x1 );

            create_beam_poly( x0, y0, dydx0, x1, y1, dydx1, mCpPolys(2) );

            // poly 4
            x0 = theta( 4 );
            y0 = polyval( mCpPolys(3), x0 );
            dydx0 = dpolyval( mCpPolys(3), x0 );

            x1 = theta( 5 );
            y1 = polyval( mCpPolys(5), x1 );
            dydx1 = dpolyval( mCpPolys(5), x1 );

            create_beam_poly( x0, y0, dydx0, x1, y1, dydx1, mCpPolys(4) );

            this->set_custom( MaterialProperty::cp );
            this->create_spline( MaterialProperty::cp, gamma );
        }

        void
        Magnesia::create_lambda()
        {
            mLambdaPolys.set_size( 3, {} );

            Vector< real > & p = mLambdaPolys( 1 );

            // fitted against recommended values from touloukian
            p = {  -5.3749e-2, + 1.0666, - 7.8520, + 2.5754e1, -3.5562e1, + 2.1733e1 };
            real T0 = 0.408 ; // magic number that causes mimimal error for low temperatures

            // compute second derivative
            Vector< real > q = { 20. * p( 0 ), 12. * p( 1 ), 6. * p( 2 ), 2. * p( 3 )};
            Vector< real > x ;
            cardano( q, x  );

            // linear extrapolation to very low temperatures
            real y = polyval( p, x( 0 ) );
            real a = dpolyval( p, x( 0 ) );
            real b = y - a * x( 0 );
            mLambdaPolys( 0 ) = { a, b };

            // linear extrapolation to higher temperatures
            y = polyval( p, x( 2 ) );
            a = dpolyval( p, x( 2 ) );
            b = y - a * x( 2 );
            mLambdaPolys( 2 ) = { a, b };

            mTLambdaSwitch = { std::exp( x( 0 ) ), std::exp( x( 2 ) ) };

            this->set_custom( MaterialProperty::lambda );


            real dkdT = std::exp( polyval( mLambdaPolys( 0 ), std::log( T0 ) ) ) / T0
                         * dpolyval( mLambdaPolys( 0 ), std::log( T0 ) );

            this->create_spline( MaterialProperty::lambda, dkdT );
        }

        real
        Magnesia::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        real
        Magnesia::cp_custom( const real T ) const
        {
            if ( T < mTCpSwitch( 0 ) )
            {
                return ( this->constant_property( MaterialProperty::gamma )
                    +    this->constant_property( MaterialProperty::beta ) * T * T ) * T ;
            }

            uint n = mTCpSwitch.length();
            for ( uint t=1; t<n; ++t )
            {
                if ( T < mTCpSwitch( t ) )
                {
                    return std::exp( polyval( mCpPolys( t-1 ), std::log( T ) ) );
                }
            }

            return std::exp( polyval( mCpPolys( n-1 ), std::log( T ) ) );
        }

        real
        Magnesia::lambda_custom( const real T ) const
        {
            if ( T < BELFEM_EPSILON ) return 0. ;
            if ( T < mTLambdaSwitch( 0 ) ) return std::exp( polyval( mLambdaPolys( 0 ), std::log( T ) ) );
            if ( T < mTLambdaSwitch( 1 ) ) return std::exp( polyval( mLambdaPolys( 1 ), std::log( T ) ) );
            return std::exp( polyval( mLambdaPolys( 2 ), std::log( T ) ) );
        }

    }
}
