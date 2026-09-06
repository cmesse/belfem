//
// Created by Christian Messe on 5/28/26.
//

#include "cl_Material_Indium.hpp"
#include "cl_Material_Abundance.hpp"
#include "fn_create_beam_poly.hpp"

namespace belfem
{
    namespace material
    {
        Indium::Indium( const real RRR,
                        const bool aBuildTables ) :
            Metal( "Indium", MaterialType::PureMetal, aBuildTables )
        {
            this->set_constants();
            this->create_cp();
            this->create_alpha();


            this->create_debye();

            this->create_rho();

            // fitted against recommended values from TPRC Vol 1, p. 151
            this->set_lambda_coefficients({ 0., 5.7641e-06, 2.5263, 138.5, -0.09006, 9.5919, 6.8058, 0. } );

            this->create_kohler();

            // fitted against Kim and Ledbetter, 1998, 10.1016/S0921-5093(98)00490-0
            this->create_mech( 19.603, 0.02701, 44.60, 295, 0.4498 );

            if ( ! std::isnan( RRR ) )
            {
                this->set_RRR( RRR ) ;
            }
        }

        Indium::~Indium()
        {
            if ( mThermalExpansion != nullptr ) delete mThermalExpansion ;
            if ( mDebyeBezier != nullptr ) delete mDebyeBezier ;

            if ( mKohlerLongBezier != nullptr ) delete mKohlerLongBezier ;
        }

        void
        Indium::set_constants()
        {
            // set the melting temperature
            this->set_constant( MaterialProperty::T_max, 	429.7485 );

            // reference density at room temperature
            this->set_constant( MaterialProperty::ref_density, 7290. );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );


            Abundance tAbundance ;

            real M = tAbundance.compute_molar_mass( "In" );

            this->set_constant( MaterialProperty::M, M );

            this->set_constant( MaterialProperty::rho_0, 0 );

            // extracted from 10.1103/PhysRevB.23.3845
            this->set_rho_i_ref( 273.15, 7.9873e-8, 112.);

            // fitted against TPRC data, curve 3 up to 2 K
            real gamma = 0.014214525 ;
            real beta = 0.0141009802706 ;

            this->set_constant( MaterialProperty::beta, beta );
            this->set_constant( MaterialProperty::gamma, gamma );

            // TPPM vol. 1
            this->set_constant( MaterialProperty::rho0_pure, 0.000828e-8 );

        }

        void
        Indium::create_alpha()
        {
            // Fitted against the TPRC dL/L dataset ( 23 points, 6 K .. 374 K ),
            // rms residual 0.0053 %. Constraints: dL/L( 293.15 K ) = 0 exactly,
            // zero derivative at 0 K ( y0 == y1 ), and a monotone control
            // polygon in y so that alpha stays positive everywhere.
            mThermalExpansion = new Bezier( );
            mThermalExpansion->basis_x() = { 0., 30.824929466, 69.702017775, 530. };
            mThermalExpansion->basis_y() = { -0.7062524808e-2, -0.7062524808e-2, -0.7062522006e-2, 0.7189422310e-2 };

            this->set_custom( MaterialProperty::alpha );
            this->set_custom( MaterialProperty::density );

            this->create_cryo_expansion( mThermalExpansion, mThermalExpansionCryo );
        }


        real
        Indium::alpha_custom( const real T ) const
        {
            return SplineLookupTable::alpha_composite( mThermalExpansion, mThermalExpansionCryo, T );
        }

        void
        Indium::create_cp()
        {
            // fitted against Touloukian dataset
            Metal::create_cp( { 2.58503846426, 2.59202348982, 2.88876390637, 2.96201718933 },
                              { 3.47186315293, 3.48305438663, 3.95800602737, 4.07415081171 },
                              { 2.96201718933, 3.57338489797, 3.93823568857, 5.60844562800 },
                              { 4.07415081171, 5.04348853709, 5.31111222017, 5.44278926752 } );
        }

        void
        Indium::create_debye()
        {
            Bezier * b1 = new Bezier( );

            b1->basis_x() = { 2.091864062, 3.07861619482, 3.44571421673, 4.441723768 };
            b1->basis_y() = { 101.658, 101.658, 131.077,  131.077 };
            mDebyeBezier = b1 ;
            mTDebyeSwitch.set_size( 2 );
            mDebyePolys.set_size( 2, {} );
            Vector< real > & p0 = mDebyePolys( 0 );
            Vector< real > & p2 = mDebyePolys( 1 );
            real x1 = b1->basis_x()( 0 );
            real x2 = b1->basis_x()( 3 );
            real x3 = std::log( this->constant_property( MaterialProperty::T_ref_rho_i ) );

            mTDebyeSwitch( 0 ) = std::exp( x1 );
            mTDebyeSwitch( 1 ) = std::exp( x2 );

            real f0 = this->constant_property( MaterialProperty::debye0K );
            real f1 = b1->basis_y()( 0 );
            real f2 = b1->basis_y()( 3 );
            real f3 = this->constant_property( MaterialProperty::debye );

            create_beam_poly( 0.0, f0, 0.0, mTDebyeSwitch( 0 ), f1, 0.0, p0 );

            // this polynomial is in the log-space
            // in continues with zero derivative at x2 and hits the reference value of theta at x3
            p2.set_size( 3 );
            p2( 0 ) = f3 - f2;
            p2( 1 ) = 2.*(x2*f2 - x2*f3);
            p2( 2 ) = f3*x2*x2 - 2*f2*x2*x3 + f2*x3*x3;
            p2 /= (x3-x2)*(x3-x2);

            this->set_custom( MaterialProperty::debye );
            this->create_spline( MaterialProperty::debye, 0 );
        }

        real
        Indium::debye_custom( const real T ) const
        {
            if ( T < mTDebyeSwitch( 0 ) )
            {
                return polyval( mDebyePolys( 0 ), T );
            }
            if ( T < mTDebyeSwitch( 1 ) )
            {
                return mDebyeBezier->y( std::log( T ) );
            }
            return polyval( mDebyePolys( 1 ), std::log( T ) );
        }

        void
        Indium::create_kohler()
        {
            mKohlerLongBezier = new Bezier();

            // based on data Lüthi: Widerstandsänderung von Metallen in hohen Magnetfeldern,
            //                      Dissertation, ETH Zürich, 1960
            // and data from 10.1103/PhysRev.120.1167

            mKohlerLongBezier->basis_x() = { 4.41376149, 5.97953263172, 7.31219750714, 9.631883912 };
            mKohlerLongBezier->basis_y() = { -4.735962476, -1.09998491267, -0.138209009,  -0.138209009 };

            real w0 =  mKohlerLongBezier->basis_x()( 0 );
            real w1 =  mKohlerLongBezier->basis_x()( 3 );

            real x0 = std::exp( w0 );
            real x1 = std::exp( w1 );

            mBSKohlerSwitch = { x0, x1 };

            // compute longitudinal polynomial for low field
            // let f = exp(g) and w = ln(x)
            real g      = mKohlerLongBezier->y( w0 );
            real dgdw   = mKohlerLongBezier->dydx( w0 ) ;
            real d2dgw2 = mKohlerLongBezier->d2ydx2( w0 ) ;
            real dgdx = dgdw / x0 ;
            real d2gdx2 = ( d2dgw2 - dgdw)/(x0*x0);

            real f = std::exp( g );
            real dfdx = f * dgdw / x0 ;
            real d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;

            mKohlerLongPolys.set_size( 2 );

            Vector< real > & p0 = mKohlerLongPolys( 0 );
            p0.set_size( 4 );
            p0( 0 ) = ( f  + x0 * ( 0.5 * d2fdx2*x0- dfdx) ) / ( x0*x0*x0 );
            p0( 1 ) = ( x0*(3.*dfdx-d2fdx2*x0) - 3.*f ) / ( x0*x0 );
            p0( 2 ) = ( 3.*f + x0 * ( 0.5 * d2fdx2 * x0 - 2. * dfdx ) ) / x0 ;
            p0( 3 ) = 0.0 ;


            Vector< real > & p1 = mKohlerLongPolys( 1 );
            p1.set_size( 1, std::exp(mKohlerLongBezier->y( w1 )) );

            mKohlerTransCoeffs.set_size( 2 );
            Vector< real > & q0 = mKohlerTransCoeffs( 0 );
            Vector< real > & q1 = mKohlerTransCoeffs( 1 );

            // also based on Lüthi
            q1 = { 0.7392, 0.4963 };

            real a = q1( 0 );
            real b = q1( 1 );

            g = a * std::pow( w0, b );
            dgdw = b * g / w0 ;
            d2dgw2 = dgdw * ( b - 1. ) / w0 ;

            dgdx = dgdw / x0 ;
            d2gdx2 = ( d2dgw2 - dgdw)/(x0*x0);

            f = std::exp( g );
            dfdx = f * dgdw / x0 ;
            d2fdx2 = f * ( dgdx * dgdx + d2gdx2 ) ;

            q0.set_size( 4 );
            q0( 0 ) = ( f  + x0 * ( 0.5 * d2fdx2*x0- dfdx) ) / ( x0*x0*x0 );
            q0( 1 ) = ( x0*(3.*dfdx-d2fdx2*x0) - 3.*f ) / ( x0*x0 );
            q0( 2 ) = ( 3.*f + x0 * ( 0.5 * d2fdx2 * x0 - 2. * dfdx ) ) / x0 ;
            q0( 3 ) = 0.0 ;

            this->set_kohler_dependencies();
        }

        inline real Indium::kohler( const real B, const real S, const real beta ) const
        {
            if ( B < BELFEM_EPSILON ) return 0.0 ;

            real BS = B * S ;

            real Along ;
            real Atrans ;
            if ( BS < mBSKohlerSwitch( 0 ) )
            {
                Along  = polyval( mKohlerLongPolys(0), BS );
                Atrans = polyval( mKohlerTransCoeffs(0), BS );
            }
            else if ( BS < mBSKohlerSwitch( 1 ) )
            {
                real w = std::log( BS ) ;
                Along  = std::exp( mKohlerLongBezier->y( w ) ) ;
                Atrans = std::exp( mKohlerTransCoeffs( 1 )( 0 ) * std::pow( w, mKohlerTransCoeffs( 1 )( 1 ) ) );
            }
            else
            {
                real w  = std::log( BS ) ;
                Along   =  mKohlerLongPolys(1)(0);
                Atrans  = std::exp( mKohlerTransCoeffs( 1 )( 0 ) * std::pow( w, mKohlerTransCoeffs( 1 )( 1 ) ) );

            }

            real c = std::cos( beta );
            real c2 = c*c ;
            real s2 = 1.0 - c2 ;

            // Pippard's angular interpolation formula
            return Along * c2 + Atrans * s2 ;
        }

    }
}