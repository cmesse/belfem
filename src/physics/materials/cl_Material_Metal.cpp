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

#include "debye.hpp"
#include <algorithm>
#include "cl_Material_Metal.hpp"

#include "commtools.hpp"

#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "cl_Mesh_Distributor.hpp"
#include "cl_Timer.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_create_database_mesh.hpp"
#include "fn_create_fifth_order_beam_poly.hpp"
#include "fn_hust.hpp"
#include "fn_polyval.hpp"
#include "fn_rho_database_is_current.hpp"
#include "fn_embed_python_guide.hpp"

namespace belfem
{
    namespace material
    {
        Metal::Metal( const string & aLabel,
                      const MaterialType aType,
                      const bool aBuildTables ) :
            SplineLookupTable( aType )
        {
            mComputeTables = aBuildTables ;
            this->set_label( aLabel );

            mSplineJ3 = this->create_J_spline( 3 );
            mSplineJ4 = this->create_J_spline( 4 );
            mSplineJ5 = this->create_J_spline( 5 );

            // default value is 5
            this->set_bloch_gruen_parameter( 5 );

            mFunctionRhoKohler    = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::rho_kohler );
            mFunctiondRhoKohlerdT = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::drhodT_kohler );
            mFunctiondRhoKohlerdB = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::drhodB_kohler );
            mFunctiondRhoKohlerdbeta = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::drhodbeta_kohler );

            this->set_constant( MaterialProperty::mu, constant::mu0 );
        }

        Metal::~Metal()
        {
            if ( mSplineJ3 != nullptr ) delete mSplineJ3 ;
            if ( mSplineJ4 != nullptr ) delete mSplineJ4 ;
            if ( mSplineJ5 != nullptr ) delete mSplineJ5 ;
            if ( mSplineJn != nullptr ) delete mSplineJn ;
            if ( mRhoData != nullptr ) delete mRhoData ;
            if ( mBhCurve != nullptr ) delete mBhCurve ;
            if ( mCpBezierLow != nullptr ) delete mCpBezierLow ;
            if ( mCpBezierMedium != nullptr ) delete mCpBezierMedium ;
            if ( mCpBezierHigh != nullptr ) delete mCpBezierHigh ;
        }

        void
        Metal::set_bloch_gruen_parameter( const real n )
        {
            this->set_constant( MaterialProperty::n_bloch_gruen, n );

            if ( n == 3 ) mFunJ = & Metal::J3 ;
            else if ( n == 4 ) mFunJ = & Metal::J4 ;
            else if ( n == 5 ) mFunJ = & Metal::J5 ;
            else mFunJ = & Metal::Jn ;

        }

        void Metal::create_cp(
            const Vector< real > & Px,
            const Vector< real > & Py,
            const Vector< real > & Qx,
            const Vector< real > & Qy,
            const Vector< real > & Rx,
            const Vector< real > & Ry )
        {
            BELFEM_ERROR( this->have( MaterialProperty::beta ) ,  "parameter beta not set for %s", this->label().c_str() );
            BELFEM_ASSERT( this->have( MaterialProperty::gamma ) , "parameter gamma not set for %s", this->label().c_str()  );
            BELFEM_ERROR( this->have( MaterialProperty::debye0K ) , "parameter debye0K not set for %s", this->label().c_str()  );

            BELFEM_ERROR( mCpBezierLow == nullptr, "Bezier curve CpBezierLow already set for %s", this->label().c_str() );
            BELFEM_ERROR( mCpBezierMedium == nullptr, "Bezier curve CpBezierMedium already set for %s", this->label().c_str() );
            BELFEM_ERROR( mCpBezierHigh == nullptr, "Bezier curve CpBezierHigh already set for %s", this->label().c_str() );

            // consistency check

            BELFEM_ERROR( Px.length() == 4 && Py.length() == 4 && Qx.length() == 4 && Qy.length() == 4,
                           "Bezier bases for cp of %s must have four control points each", this->label().c_str() );

            // the bisection in Bezier::xi_by_x needs a monotonic x
            BELFEM_ERROR( Px( 0 ) < Px( 1 ) && Px( 1 ) < Px( 2 ) && Px( 2 ) < Px( 3 ),
                           "Bezier basis Px for cp of %s is not monotonic", this->label().c_str() );
            BELFEM_ERROR( Qx( 0 ) < Qx( 1 ) && Qx( 1 ) < Qx( 2 ) && Qx( 2 ) < Qx( 3 ),
                           "Bezier basis Qx for cp of %s is not monotonic", this->label().c_str() );

            // the two curves must share their connection point
            BELFEM_ERROR( Px( 3 ) == Qx( 0 ) && Py( 3 ) == Qy( 0 ),
                           "Inconsistent connection point for P and Q of %s", this->label().c_str() );

            //  cp ≈ cv = γ * T + β * T³
            real beta  = this->constant_property( MaterialProperty::beta );
            real gamma = this->constant_property( MaterialProperty::gamma );
            real theta = this->constant_property( MaterialProperty::debye0K );

            uint n = Rx.length() == 4 && Ry.length() == 4 ? 5 : 4 ;

            mTCpSwitch.set_size( n );

            real x1 = Px( 0 );
            real x2 = Qx( 0 );
            real x3 = Qx( 3 );
            real & T0 = mTCpSwitch( 0 );
            real & T1 = mTCpSwitch( 1 );
            real & T2 = mTCpSwitch( 2 );
            real & T3 = mTCpSwitch( 3 );

            // the low polynomial is valid for T < theta/50
            T0 = theta * 0.02 ;
            T1 = std::exp( x1 ) ;
            T2 = std::exp( x2 );
            T3 = std::exp( x3 );

            // now the splines
            mCpBezierLow = new Bezier();
            mCpBezierLow->basis_x() = Px ;
            mCpBezierLow->basis_y() = Py ;

            mCpBezierMedium = new Bezier();
            mCpBezierMedium->basis_x() = Qx ;
            mCpBezierMedium->basis_y() = Qy ;

            if ( n == 5 )
            {
                // the two curves must share their connection point
                BELFEM_ERROR( Qx( 3 ) == Rx( 0 ) && Qy( 3 ) == Ry( 0 ),
                    "Inconsistent connection point for Q and R of %s", this->label().c_str() );

                // the bisection in Bezier::xi_by_x needs a monotonic x( t ), but
                // unlike Px and Qx the control polygon of Rx may dip ( iron does ),
                // so check the sharp condition instead: x'( t ) is the Bernstein
                // quadratic with coefficients d0, d1, d2, which stays positive on
                // [0,1] iff d0 > 0, d2 > 0 and d1 > -sqrt( d0 * d2 )
                real d0 = Rx( 1 ) - Rx( 0 );
                real d1 = Rx( 2 ) - Rx( 1 );
                real d2 = Rx( 3 ) - Rx( 2 );
                BELFEM_ERROR( d0 > 0.0 && d2 > 0.0 && d1 > -std::sqrt( d0 * d2 ),
                    "Bezier basis Rx for cp of %s folds back ( x(t) not monotonic )",
                    this->label().c_str() );

                mCpBezierHigh = new Bezier ;
                mCpBezierHigh->basis_x() = Rx ;
                mCpBezierHigh->basis_y() = Ry ;
            }

            // the Sommerfeld-Debye cubic must not reach past the first control point
            BELFEM_ERROR( T0 < T1,
                          "Debye temperature of %s places the low cp polynomial at %f K, "
                          "which is above the first Bezier control point at %f K",
                          this->label().c_str(), ( double ) T0, ( double ) T1 );

            // computing the supporting point and the temperature
            real cp0      = T0 * ( gamma + beta * T0 * T0 );
            real dcpdT0   = gamma + 3.0 * beta * T0 * T0 ;
            real d2cpdT20 = 6.0 * beta * T0 ;

            // expressing the derivatives in log-log space
            real x0 = std::log( T0 );
            real y0 = std::log( cp0 );
            real dydx0 = T0*dcpdT0/cp0 ;
            real d2ydx20 = T0*T0* d2cpdT20 / cp0 - dydx0*dydx0 + dydx0 ;

            mCpPolys.set_size( 3, {} );
            mCpPolys( 0 ) = { beta, 0., gamma, 0. };
            mCpPolys( 1 ).set_size( 6 );
            create_fifth_order_beam_poly(
                x0,
                y0,
                dydx0,
                d2ydx20,
                x1,
                mCpBezierLow->y( x1 ),
                mCpBezierLow->dydx( x1 ),
                mCpBezierLow->d2ydx2( x1 ),
                mCpPolys( 1 ) );

            // The curve is fitted in log-log space, but the extrapolation above T3
            // is a plain linear polynomial in T, so both the value and the slope
            // have to be carried back: cp = exp( y ), dcp/dT = cp * dy/dx / T .
            if ( mCpBezierHigh == nullptr )
            {
                real cp3    = std::exp( Qy( 3 ) );
                real dcpdT3 = cp3 * mCpBezierMedium->dydx( x3 ) / T3 ;
                mCpPolys( 2 ) = { dcpdT3, cp3 - dcpdT3 * T3 };
            }
            else
            {
                real x4 = Rx( 3 );
                real T4 = std::exp( x4 );
                real cp4    = std::exp( Ry( 3 ) );
                real dcpdT4 = cp4 * mCpBezierHigh->dydx( x4 ) / T4 ;
                mCpPolys( 2 ) = { dcpdT4, cp4 - dcpdT4 * T4 };
                mTCpSwitch( 4 ) = T4 ;
            }
            this->set_have( MaterialProperty::cp );
            this->create_spline( MaterialProperty::cp, gamma );
        }

        void
        Metal::set_RRR( const real RRR )
        {

            this->set_constant( MaterialProperty::RRR, RRR );

            // inner resistivity
            real rho_i = this->constant_property( MaterialProperty::rho_i_ref );

            // initial guess
            real rho_0 = rho_i / ( RRR - 1.0 );
            real T = this->constant_property( MaterialProperty::T_ref_rho_i );

            real x0 = 0.9 * rho_0 ;
            this->set_constant( MaterialProperty::rho_0, x0 );
            real f0 = ( this->rho_i_custom( T ) + x0 )/ x0 - RRR ;

            real x1 = 1.1 * rho_0 ;
            this->set_constant( MaterialProperty::rho_0, x1 );
            real f1 =  ( this->rho_i_custom( T ) + x1 )/ x1 - RRR ;


            real f = f1 ;
            bool s = false ;
            real x = 0.5 * ( x0 + x1 ) ;

            uint it = 0 ;
            while (  std::abs( f ) > 1e-12 )
            {
                x = ( x0 * f1 - x1 * f0 ) / ( f1 - f0 ) ;
                this->set_constant( MaterialProperty::rho_0, x );
                f = ( this->rho_i_custom( T ) + x ) / x - RRR ;
                if ( f0 * f > 0 )
                {
                    if ( ! s ) f1 *= 0.5 ;

                    x0 = x ;
                    f0 = f ;
                    s = false ;
                }
                else
                {
                    if ( s ) f0 *= 0.5 ;
                    x1 = x ;
                    f1 = f ;
                    s = true ;
                }

                BELFEM_ERROR( it++ < 100, "Created infinite loop while trying to iterate rho_0" );
            }


            this->set_constant( MaterialProperty::rho_0, x );

            if ( this->have( MaterialProperty::lambda ) )
            {
                // start tangent: lambda -> L0 * T / rho_0 as T -> 0, so the
                // slope at the origin is L0 / rho_0, not its reciprocal
                this->create_spline( MaterialProperty::lambda,
                           constant::L0 / x );
            }


            // lookup tables for temperature, field and angle dependent
            // electric resistivity and thermal conductivity
            if ( this->depends( MaterialProperty::rho, MaterialDependency::angleBxJ ) && mComputeTables )
            {
                this->populate_rho_database();
            }


            this->set_have( MaterialProperty::rho );
            this->set_have( MaterialProperty::lambda );

            // set_spline() resets the lambda dependencies to T-only, which
            // drops the field flags create_kohler() registered before this
            // call; without them the tool and the calculator never select
            // lambda( T, B, beta ). rho keeps its flags, so use them as the
            // marker that this metal carries a magnetoresistance model
            if ( this->depends( MaterialProperty::rho, MaterialDependency::normB ) )
            {
                this->set_kohler_dependencies();
            }
        }

        Spline *
        Metal::create_J_spline( const real aExponent )
        {
            Vector< real > tZ ;
            Vector< real > tY ;

            double n = aExponent ;

            // number of integration points per interval
            int m = 10 ;

            // number of points to peak
            int p = 20 ;

            // number of points
            int q  = static_cast< int >( std::ceil( ( 1271.736 / n + 329.095 ) / n + 78.162 ) ) ;


            tZ.set_size( q );
            tY.set_size( q );
            debye_table( &n, &m, &p, &q, tZ.data(), tY.data() );

            if ( n == 3. )
            {
                mZ3Max = tZ( q-1 ) ;
                mJ3Max = tY( q-1 ) ;
            }
            else if ( n == 4. )
            {
                mZ4Max = tZ( q-1 ) ;
                mJ4Max = tY( q-1 ) ;
            }
            else if ( n == 5. )
            {
                mZ5Max = tZ( q-1 ) ;
                mJ5Max = tY( q-1 ) ;
            }
            else
            {
                mZnMax = tZ( q-1 ) ;
                mJnMax = tY( q-1 ) ;
            }

            SpMatrix * tA = new SpMatrix ;
            spline::create_helpmatrix(
                static_cast< index_t >( q ),
                tZ(1),
                *tA,
                spline::SplineBC::Tangent,
                spline::SplineBC::Tangent );

            Spline * aSpline = new Spline(
                tZ,
                tY,
                *tA,
                spline::SplineBC::Tangent,
                spline::SplineBC::Tangent );
            delete tA ;

            return aSpline ;
        }

        real
        Metal::l( const real T ) const
        {
            if ( std::abs( T - this->constant_property( MaterialProperty::T_ref_density ) ) < BELFEM_EPSILON )
            {
                return 1.0 ;
            }

            return std::exp(
                mSplines( static_cast< uint >( MaterialProperty::alpha ) )->integrate( T ) );
        }

        real
        Metal::spline_property( const MaterialProperty aProperty, const real aX ) const
        {
            BELFEM_ASSERT( mSplines( static_cast< size_t >( aProperty ) ) != nullptr,
                           "Spline for property is not set" );
            return mSplines( static_cast< size_t >( aProperty ) )->eval( aX ) ;
        }

        void
        Metal::set_rho_i_ref( const real T_ref, const real rho_i_ref, const real theta )
        {
            real Z = std::isnan( theta ) ? this->debye( T_ref ) / T_ref : theta / T_ref ;

            if ( ! std::isnan( theta ) )
            {
                this->set_constant( MaterialProperty::debye, theta );
            }

            real n = this->constant_property( MaterialProperty::n_bloch_gruen );
            real A = rho_i_ref * std::pow( Z, n ) / this->J( Z );

            this->set_constant( MaterialProperty::rho_i_ref, rho_i_ref );
            this->set_constant( MaterialProperty::T_ref_rho_i, T_ref );
            this->set_constant( MaterialProperty::A_bloch_gruen, A );

            this->set_have( MaterialProperty::rho_i_ref );
            this->set_have( MaterialProperty::T_ref_rho_i );
            this->set_have( MaterialProperty::A_bloch_gruen );
            this->set_custom( MaterialProperty::rho );
            this->set_custom( MaterialProperty::rho_i ) ;
        }


        void
        Metal::create_rho()
        {
            this->set_dependency( MaterialProperty::rho_i, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho,   MaterialDependency::T );
            this->set_have( MaterialProperty::rho_i );
            this->set_have( MaterialProperty::rho );
            this->create_spline( MaterialProperty::rho_i, 0.0 );
            this->set_custom( MaterialProperty::rho );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho_i, MaterialDependency::T );
        }


        real Metal::rho_i_custom( const real T ) const
        {
            real Z = this->debye( T ) / T ;

            return this->constant_property( MaterialProperty::A_bloch_gruen ) * this->J( Z ) /
                std::pow( Z, this->constant_property( MaterialProperty::n_bloch_gruen ) );
        }

        real
        Metal::invert_debye( const real T, const real y, const real theta_guess,
                             real ( Metal::*aFunction )( const real, const real ) const )
        {
            real x  =  std::isnan( theta_guess ) ? this->constant_property( MaterialProperty::debye0K ) : theta_guess ;

            real x0 = x - 1.0 ;
            real x1 = x + 1.0 ;
            real f0 = ( this->*aFunction )( T, x0 ) - y;
            real f1 = ( this->*aFunction )( T, x1 ) - y;
            real f  = ( this->*aFunction )( T, x )  - y;

            // check gradient
            real dT = 1.0 ;
            real xi = 1.1 ;
            if ( f0 * f1 > 0 )
            {
                if ( std::abs(f0) < std::abs(f) )
                {
                    x1 = x ;
                    f1 = f ;
                    // loop backwards
                    while ( f0*f1 > 0 && x0 > 10.0 )
                    {
                        x1 = x0 ;
                        f1 = f0 ;
                        x0 -= dT ;
                        dT *= xi ;
                        f0 = ( this->*aFunction )( T, x0 ) - y;
                    }
                }
                else if ( std::abs(f1) < std::abs(f) )
                {
                    x0 = x ;
                    f0 = f ;
                    while ( f1*f0 > 0 && x1 < 1000.0 )
                    {
                        x0 = x1 ;
                        f0 = f1 ;
                        x1 += dT ;
                        dT *= xi ;
                        f1 = ( this->*aFunction )( T, x1 ) - y;
                    }
                }
                else
                {
                    BELFEM_ERROR( false, "error in computing debye temperature from cp" );
                }
            }
            else
            {
                x0 = 10 ;
                x1 = this->constant_property( MaterialProperty::debye0K ) ;
                f0 = ( this->*aFunction )( T, x0 ) - y;
                f1 = ( this->*aFunction )( T, x1 ) - y;

                BELFEM_ERROR( f0*f1 < 0, "Couldn't detect debye temperature interval" );


            }

            uint k = 0 ;
            x = x0 ;
            f = f0 ;

            while ( k++ < 100 && std::abs( f ) > 1e-12 )
            {
                x -= 0.95*f0 *( x1-x0)/( f1-f0 );
                if ( x < x0 || x > x1 )
                {
                    x = 0.5 * ( x0 + x1 );
                }
                f = ( this->*aFunction )( T, x) - y;
                if ( f0 * f > 0.0 )
                {
                    x0 = x ;
                    f0 = f ;
                }
                else
                {
                    x1 = x ;
                    f1 = f ;
                }

            }
            BELFEM_ERROR( k < 100, "Created infinite loop while trying to iterate debye temperature" );

            return x ;
        }


        real
        Metal::compute_debye_from_cp( const real T, const real cp, const real theta_guess )
        {
            // cp_from_debye() carries the dilation term 9 alpha^2 K T / rho, so this
            // inversion needs the elastic data. Below ~ 30 K use compute_debye_from_cv().
            BELFEM_ERROR( this->have( MaterialProperty::E ) && this->have( MaterialProperty::nu ),
                "%s: compute_debye_from_cp() needs E and nu - call create_mech() first, "
                "or use compute_debye_from_cv() at cryogenic temperatures",
                this->label().c_str() );

            real y = std::isnan( cp ) ? this->cp( T ) : cp ;

            return this->invert_debye( T, y, theta_guess, & Metal::cp_from_debye );
        }

        real
        Metal::compute_debye_from_cv( const real T, const real cv, const real theta_guess )
        {
            // The neglected dilation term is ( cp - cv ) / cp = alpha_V gamma_G T
            // ( Grueneisen identity ), of order 1e-5 at the 12 - 14 K where the
            // Debye curves are anchored. With gamma_G ~ 2 the bound below rejects
            // any call above roughly 100 K, and needs only alpha - which every
            // metal has before create_debye(), unlike the elastic data.
            real tNeglected = 6.0 * this->alpha( T ) * T ;

            BELFEM_ERROR( tNeglected < 1e-3,
                "%s: compute_debye_from_cv() called at %g K, where cp - cv is %g of cp - "
                "use compute_debye_from_cp() there",
                this->label().c_str(), ( double ) T, ( double ) tNeglected );

            real y = std::isnan( cv ) ? this->cp( T ) : cv ;

            return this->invert_debye( T, y, theta_guess, & Metal::cv_from_debye );
        }

        real
        Metal::compute_debye_from_rho( const real T, const real rho, const real theta_guess )
        {
            real x = std::isnan( theta_guess ) ? this->constant_property( MaterialProperty::debye0K ) : theta_guess ;
            real x0 = x - 10 ;
            real x1 = x + 10 ;

            real f0 = this->fun_debye_from_rho( T, rho, x0 ) ;
            real f  = this->fun_debye_from_rho( T, rho, x1 ) ;

            while ( f0 * f > 0 )
            {
                x0 -= 1.0 ;
                x1 += 1.0 ;
                f0 = this->fun_debye_from_rho( T, rho, x0 ) ;
                f = this->fun_debye_from_rho( T, rho, x1 ) ;
            }

            uint tCount = 0 ;

            while ( std::abs( f ) > 1e-6 )
            {
                x = 0.5 * ( x0 + x1 );
                f = this->fun_debye_from_rho( T, rho, x ) ;
                if ( f0 * f > 0.0 )
                {
                    x0 = x ;
                    f0 = f ;
                }
                else
                {
                    x1 = x ;
                }
                BELFEM_ERROR( tCount++ < 100, "Created infinite loop while trying to iterate debye temperature" );
            }
            return x ;
        }

        real
        Metal::cv_from_debye( const real T, const real theta ) const
        {
            // electron contribution
            real aValue = this->constant_property( MaterialProperty::gamma ) * T ;

            // phonon contribution
            real Z = theta / T ;
            aValue += 9.0 * this->constant_property( MaterialProperty::R ) * this->J4( Z ) / ( Z * Z * Z ) ;

            return aValue ;
        }

        real
        Metal::cp_from_debye( const real T, const real theta ) const
        {
            BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", this->label().c_str() ) ;

            // dilation term cp - cv = 9 alpha^2 K T / rho, assuming isotropy
            real alpha = this->alpha( T ) ;

            return this->cv_from_debye( T, theta )
                 + 9.0 * alpha * alpha *  this->K( T )/this->density( T ) * T ;
        }

        real
        Metal::fun_debye_from_rho( const real T, const real rho, const real theta )
        {
            BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", this->label().c_str() ) ;

            real A = this->constant_property( MaterialProperty::A_bloch_gruen ) ;
            real n = this->constant_property( MaterialProperty::n_bloch_gruen ) ;
            real Z = theta / T ;
            real rho_i = A * this->J( Z ) / std::pow( Z, n ) ;

            real rho_0 = this->constant_property( MaterialProperty::rho_0 ) ;
            return ( rho_i + rho_0 - rho ) / rho ;
        }

        void
        Metal::set_lambda_coefficients( const Vector< real > & aCoeffs )
        {
            BELFEM_ERROR( this->type() == MaterialType::PureMetal, "lambda coefficients can only be used in conjuction with pure metals" );
            mLambdaCoefficients = aCoeffs ;
            this->set_have( MaterialProperty::lambda );

        }

        real
        Metal::lambda_custom( const real T ) const
        {
            BELFEM_ASSERT( this->type() == MaterialType::PureMetal, "lambda coefficients can only be used in conjunction with pure metals" );

            if ( T <= 0.0 ) return 0.0 ;

            const Vector< real > & p = mLambdaCoefficients ;

            real rho_0  = this->constant_property( MaterialProperty::rho_0 );
            real beta   = rho_0 / constant::L0 ;
            real w_0 = beta / T ;

            real w_i = hust( p, T );


            real C = p(7) / std::pow( beta / 0.0003, p(0) );
            real w_i0 = C * w_i * w_0 / ( w_i + w_0 );

            return 1.0 / ( w_0 + w_i + w_i0 );
        }

        real
        Metal::group_velocity( const real T ) const
        {
            real M  = this->constant_property( MaterialProperty::M ) ;
            real kB = constant::kB ;
            real q = this->constant_property( MaterialProperty::q ) ;
            real theta = this->debye( T ) ;
            real NA = constant::NA ;
            real h = constant::h ;
            real pi = constant::pi ;
            real rho = this->density( T );

            return std::pow( 4. * pi * M / ( 3. * q * NA * rho ), 1.0 / 3.0 )
                   * kB * theta / h ;
        }

        real Metal::grueneisen( const real T ) const
        {
            // sepecific heat at constant pressure
            real cp = this->cp( T ) ;

            // density
            real rho = this->density( T ) ;

            // volumentric expansion coefficient
            real alpha = 3. * this->alpha( T ) ;

            // bulk modulus
            real Kt = this->K( T );

            // specific heat capacity at constant volume
            real cv = cp - alpha * alpha * T * Kt / rho ;

            // Grüneisen parameter
            return alpha * Kt / ( cv * rho );
        }

        void
        Metal::grueneisen( const real T,real & gamma, real & c ) const
        {
            // sepecific heat at constant pressure
            real cp = this->cp( T ) ;

            // density
            real rho = this->density( T ) ;

            // volumentric expansion coefficient
            real alpha = 3. * this->alpha( T ) ;

            // bulk modulus
            real Kt = this->K( T );

            // specific heat capacity at constant volume
            real cv = cp - alpha * alpha * T * Kt / rho ;

            // Grüneisen parameter
            gamma = alpha * Kt / ( cv * rho );

            real Ks = Kt * cp/cv ;

            // speed of sound
            c = std::sqrt( Ks / rho );
        }

        real
        Metal::kohler( const real B, const real S, const real beta ) const
        {
            BELFEM_ERROR( false, "the material %s doesn't have a kohler curve",
                this->label().c_str() );
            return BELFEM_QUIET_NAN ;
        }

        real
        Metal::rho_kohler( const real T, const real B, const real beta  ) const
        {

            real rho_i  = this->constant_property( MaterialProperty::rho_i_ref );
            real rho_0K = this->constant_property( MaterialProperty::rho_0 ) ;
            real rho_ref = rho_i + rho_0K ;
            real rho_0T = T < BELFEM_EPSILON ? rho_0K : this->rho_i_custom( T ) + rho_0K ;

            if ( B < 0.001 ) return rho_0T ;

            real S = rho_ref / rho_0T ;

            return ( this->kohler( B, S, beta ) + 1.) * rho_0T ;
        }

        real
        Metal::drhodT_kohler( const real T, const real B, const real beta ) const
        {
            if ( T < gFinDiffDeltaT + BELFEM_EPSILON ) return 0 ;

            real rho_ia = this->rho_i_custom( T - gFinDiffDeltaT );
            real rho_ib = this->rho_i_custom( T + gFinDiffDeltaT );

            if ( B < 0.001 )
            {
                return 0.5 * ( rho_ib - rho_ia ) / gFinDiffDeltaT ;
            }

            real rho_i  = this->constant_property( MaterialProperty::rho_i_ref );
            real rho_0K = this->constant_property( MaterialProperty::rho_0 ) ;

            real rho_ref = rho_i + rho_0K ;

            real g = this->rho_i_custom( T ) + rho_0K ;
            real dg = 0.5 * ( rho_ib - rho_ia ) / gFinDiffDeltaT ;

            real Sa = rho_ref / ( rho_ia + rho_0K );
            real Sb = rho_ref / ( rho_ib + rho_0K );
            real S  = rho_ref / g ;

            real f = this->kohler( B, S, beta ) + 1. ;
            real df = 0.5 * ( this->kohler( B, Sb, beta ) - this->kohler( B, Sa, beta ) ) / gFinDiffDeltaT ;

            return g * df + dg * f ;
        }

        real
        Metal::drhodB_kohler( const real T, const real B, const real beta ) const
        {
            if ( T < gFinDiffDeltaT + BELFEM_EPSILON ) return 0 ;

            if ( B < gFinDiffDeltaB + BELFEM_EPSILON ) return 0 ;

            real rho_i  = this->constant_property( MaterialProperty::rho_i_ref );
            real rho_0K = this->constant_property( MaterialProperty::rho_0 ) ;

            real rho_ref = rho_i + rho_0K ;

            real g = this->rho_i_custom( T ) + rho_0K ;
            real S  = rho_ref / g ;

            real fa = this->kohler( B - gFinDiffDeltaB, S, beta ) ;
            real fb = this->kohler( B + gFinDiffDeltaB, S, beta ) ;

            real df = 0.5 * ( fb-fa ) / gFinDiffDeltaB ;

            return g * df ;
        }

        real
        Metal::drhodbeta_kohler( const real T, const real B, const real beta ) const
        {
            if ( T < gFinDiffDeltaT + BELFEM_EPSILON ) return 0 ;
            if ( B < gFinDiffDeltaB + BELFEM_EPSILON ) return 0 ;

            real rho_i  = this->constant_property( MaterialProperty::rho_i_ref );
            real rho_0K = this->constant_property( MaterialProperty::rho_0 ) ;

            real rho_ref = rho_i + rho_0K ;

            real g = this->rho_i_custom( T ) + rho_0K ;
            real S  = rho_ref / g ;

            real fa = this->kohler( B, S, beta - gFinDiffDeltaAngle ) ;
            real fb = this->kohler( B, S, beta + gFinDiffDeltaAngle ) ;

            real df = 0.5 * ( fb-fa ) / gFinDiffDeltaAngle ;

            return g * df ;
        }

        void
        Metal::set_kohler_dependencies()
        {
            this->set_dependency( MaterialProperty::rho, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::normB);
            this->set_dependency( MaterialProperty::rho, MaterialDependency::angleBxJ );

            this->set_dependency( MaterialProperty::lambda, MaterialDependency::T );
            this->set_dependency( MaterialProperty::lambda, MaterialDependency::normB);
            this->set_dependency( MaterialProperty::lambda, MaterialDependency::angleBxJ );
        }

        real
        Metal::kohler_find_bs_crit( const real delta_rho_res )
        {
            real x0 = -10 ;
            real x1 = 10 ;
            real df = std::log( delta_rho_res );

            real beta = 0.5 * constant::pi ;
            real f0 = std::log(this->kohler( 1.0, std::exp(x0), beta ) ) - df ;
            real f1 = std::log( this->kohler( 1.0, std::exp(x1), beta ) ) - df ;

            real x = x0 ;
            real f = f0 ;

            index_t tCount = 0 ;

            while ( std::abs( f ) > 1e-12 && tCount++ < 100 )
            {
                x -= 0.95 * f0 *( x1-x0)/( f1-f0 );
                if ( x < x0 || x > x1 )
                {
                    x = 0.5 * ( x0 + x1 );
                }


                f = std::log( this->kohler( 1.0, std::exp(x), beta ) ) - df ;
                if ( f0 * f > 0.0 )
                {
                    x0 = x ;
                    f0 = f ;
                }
                else
                {
                    x1 = x ;
                    f1 = f ;
                }
            }
            return std::exp(x);
        }

        void
        Metal::populate_rho_database()
        {
            if ( mRhoData != nullptr )
            {
                delete mRhoData ;
            }

            string tFile = sprint( "%s_RRR%u.hdf5" ,
                this->label().c_str(),
                ( uint ) this->constant_property( MaterialProperty::RRR ) );

            // A cached database carries no format version, so an old file would be
            // loaded and then fail deep inside load_rho_database with an opaque
            // "Dataset RRR ... does not exist". Probe for the marker instead and
            // rebuild. The decision must be identical on every rank ( all ranks run
            // this branch ), so the master probes and broadcasts.
            bool tUsable = file_exists( tFile );

            if ( tUsable && comm_rank() == 0 )
            {
                tUsable = rho_database_is_current( tFile );

                if ( ! tUsable )
                {
                    message( InfoLevel::Default,
                        "    Warning: %s predates the current rho-database format "
                        "( no RRR entry ) and is being rebuilt",
                        tFile.c_str() );
                }
            }
            if ( comm_size() > 1 )
            {
                broadcast( tUsable );
            }

            if ( ! tUsable )
            {
                // every rank takes this branch in lockstep: the master holds
                // the real work grid and evaluates every node ( the grid is
                // never partitioned, so all nodes carry owner zero ), while
                // the other ranks hold an empty grid copy and join the
                // projection's collective pair inside the Database
                // constructor, receiving the finished values there. The
                // former distributor-based build is retired -- it died on its
                // own empty worker partitions before ever completing.
                this->populate_rho_database_serial();

                this->save_rho_database( tFile );
            }
            else
            {
                this->load_rho_database( tFile );
            }

            mFunctionRhoKohler       = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::rho_table );
            mFunctiondRhoKohlerdT    = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::drhodT_table );
            mFunctiondRhoKohlerdB    = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::drhodB_table );
            mFunctiondRhoKohlerdbeta = static_cast< real ( Material::* ) ( const real, const real, const real ) const >( & Metal::drhodbeta_table );
            comm_barrier() ;
        }

        void
        Metal::populate_rho_database_serial()
        {
            Mesh * tMesh ;

            // create a tensor mesh
            tMesh = create_database_mesh();

            mDatabaseTmin = tMesh->tensorconf()->min( 0 ) ;
            mDatabaseTmax = tMesh->tensorconf()->max( 0 ) ;
            mDatabaseBmin = std::pow( 10, tMesh->tensorconf()->min( 1 )  );
            mDatabaseBmax = std::pow( 10, tMesh->tensorconf()->max( 1 )  );


            Vector< real > & tRho = tMesh->create_field( "rho" );
            this->populate_rho_database( tMesh, tRho );

            string tDisplayLabel = sprint( "%s RRR %u",
                this->label().c_str(),
                ( uint ) this->constant_property( MaterialProperty::RRR ) );
            mRhoData = new Database( tMesh, "rho", true, tDisplayLabel );

            // the Database keeps its own copy of the values and config;
            // the work grid is ours to free
            delete tMesh ;
        }

        void
        Metal::populate_rho_database( Mesh * aMesh, Vector< real > & aRho )
        {
            proc_t tRank = comm_rank() ;

            real T0 = mDatabaseTmin ;

            Cell< mesh::Node * > & tNodes = aMesh->nodes();

            // reference value at 273.15 K and 0 T
            real rho_i  = this->constant_property( MaterialProperty::rho_i_ref );
            real rho_0K = this->constant_property( MaterialProperty::rho_0 ) ;
            real rho_ref = rho_i + rho_0K ;

            // populate dataset
            index_t tCount = 0 ;

            real rho_0T  = T0 < BELFEM_EPSILON ? rho_0K : this->rho_i_custom( T0 ) + rho_0K ;

            real rho ;

            for ( mesh::Node * tNode : tNodes )
            {
                if ( tNode->owner() != tRank ) continue ;

                real T = tNode->x();

                if ( T != T0 )
                {
                    // reference at 0 T
                    rho_0T = this->rho_i_custom( T ) + rho_0K ;
                    T0 = T ;
                }

                real B = std::pow( 10., tNode->y() );
                real beta = tNode->z();


                // similarity parameter
                real S = rho_ref / rho_0T ;

                // projected value for temperature, field and angle
                if ( T < BELFEM_EPSILON )
                {
                    rho = rho_0K ;
                }
                else if ( B < 0.001 )
                {
                    rho = rho_0T ;
                }
                else
                {
                    rho = ( this->kohler( B, S, beta )  + 1. ) * rho_0T ;
                }

                // store data in field
                aRho( tCount++ ) = std::log(rho) ;
            }
        }

        void
        Metal::save_rho_database( const std::string & aPath )
        {
            if ( comm_rank() != 0 ) return ;

            HDF5 tFile( aPath, FileMode::NEW );
            tFile.save_data( "label" , this->label() );
            tFile.save_data( "RRR",  this->constant_property( MaterialProperty::RRR ) );
            tFile.save_data( "Tmin", mDatabaseTmin );
            tFile.save_data( "Tmax", mDatabaseTmax );
            tFile.save_data( "Bmin", mDatabaseBmin );
            tFile.save_data( "Bmax", mDatabaseBmax );

            tFile.create_group( "rho" );
            mRhoData->save( tFile.active_group() );
            tFile.close_active_group();

            // ship the reference reader with the table, if it is available
            material::embed_python_guide( tFile );

            tFile.close();
        }

        void
        Metal::load_rho_database( const std::string & aPath )
        {
            if ( comm_rank() == 0 )
            {
                HDF5 tFile( aPath, FileMode::OPEN_RDONLY );

                string tLabel ;
                tFile.load_data( "label", tLabel );

                BELFEM_ERROR( tLabel == this->label(), "Found material %s but expect %s",
                    tLabel.c_str(), this->label().c_str() );

                real tRRR ;
                tFile.load_data( "RRR", tRRR );
                BELFEM_ERROR( tRRR == this->constant_property( MaterialProperty::RRR ), "Found RRR %g for material %s but expect %g",
                   ( double ) tRRR, this->label().c_str(), ( double ) this->constant_property( MaterialProperty::RRR ) );

                tFile.load_data( "Tmin", mDatabaseTmin );
                tFile.load_data( "Tmax", mDatabaseTmax );
                tFile.load_data( "Bmin", mDatabaseBmin );
                tFile.load_data( "Bmax", mDatabaseBmax );

                tFile.select_group( "rho" );
                if ( mRhoData != nullptr ) delete mRhoData ;
                mRhoData = new Database( tFile.active_group() , "rho");
                tFile.close_active_group();
                tFile.close();

                Vector< real > tLimits( 4 );
                tLimits(0) = mDatabaseTmin ;
                tLimits(1) = mDatabaseTmax ;
                tLimits(2) = mDatabaseBmin ;
                tLimits(3) = mDatabaseBmax ;
                share( tLimits );
            }
            else
            {
                mRhoData    = new Database( ( hid_t ) 0 , "rho");

                Vector< real > tLimits( 4 );
                receive( tLimits );
                mDatabaseTmin = tLimits(0) ;
                mDatabaseTmax = tLimits(1) ;
                mDatabaseBmin = tLimits(2) ;
                mDatabaseBmax = tLimits(3) ;
            }


        }

        void
        Metal::set_bh_curve( const BhCurve * aCurve )
        {
            BELFEM_ERROR( mBhCurve == nullptr, "Material %s already has a BH curve", this->label().c_str() );

            mBhCurve = aCurve ;
        }

        void
        Metal::set_table_flags( const bool aFlag )
        {
            mComputeTables = aFlag ;
        }

        void
        Metal::create_mech( const real E2, const real nu2, const real T2, const real deltaK, const real deltaG )
        {
            // once-per-material setup: the checks stay active in release builds
            BELFEM_ERROR( this->have( MaterialProperty::alpha ),
                "Need alpha table for %s", this->label().c_str() );
            BELFEM_ERROR( this->have( MaterialProperty::cp ),
                "Need cp table for %s", this->label().c_str() );
            BELFEM_ERROR( this->have( MaterialProperty::ref_density ),
                "Need reference density for %s", this->label().c_str() );
            BELFEM_ERROR( this->have( MaterialProperty::T_max ),
                "Need maximum temperature for %s", this->label().c_str() );

            // the strain integral must start on the cryogenic alpha branch ( alpha = C cp ),
            // otherwise eps( T ) at low T inherits a fitted Bezier slope and dE/dT( 0 ),
            // dnu/dT( 0 ) are no longer zero
            BELFEM_ERROR( this->alpha_switch_temperature() > 0.0,
                "%s: create_mech() needs the cryogenic alpha branch, but the split temperature is %g K",
                this->label().c_str(), ( double ) this->alpha_switch_temperature() );

            BELFEM_ERROR( std::isfinite( E2 ) && E2 > 0.0,
                "%s: Young's modulus at the anchor must be positive, reads %g GPa",
                this->label().c_str(), ( double ) E2 );
            BELFEM_ERROR( nu2 > -1.0 && nu2 < 0.5,
                "%s: Poisson ratio at the anchor must lie in ( -1, 1/2 ), reads %g",
                this->label().c_str(), ( double ) nu2 );
            BELFEM_ERROR( std::isfinite( deltaK ) && deltaK > 0.0 && std::isfinite( deltaG ) && deltaG > 0.0,
                "%s: softening constants must be positive, read deltaK = %g, deltaG = %g",
                this->label().c_str(), ( double ) deltaK, ( double ) deltaG );

            const real Tmax = this->constant_property( MaterialProperty::T_max );
            BELFEM_ERROR( std::isfinite( T2 ) && T2 > 0.0 && T2 <= Tmax,
                "%s: anchor temperature %g K is outside ( 0, %g ]",
                this->label().c_str(), ( double ) T2, ( double ) Tmax );

            // reference length at 0 K, cached: eps( T ) = 3 ln( l( T ) / l( 0 ) )
            mL0 = this->l( 0.0 );
            BELFEM_ERROR( std::isfinite( mL0 ) && mL0 > 0.0,
                "%s: relative length at 0 K reads %g - the alpha integral is not usable",
                this->label().c_str(), ( double ) mL0 );

            mDeltaK = deltaK ;
            mDeltaG = deltaG ;

            // anchor: isotropic identities, stored in Pa, then back to 0 K along the exponential
            const real E  = E2 * 1e9 ;
            const real K2 = E / ( 3.0 - 6.0 * nu2 );
            const real G2 = E / ( 2.0 + 2.0 * nu2 );
            const real e2 = this->eps( T2 );

            mK0 = K2 * std::exp( deltaK * e2 );
            mG0 = G2 * std::exp( deltaG * e2 );

            BELFEM_ERROR( std::isfinite( mK0 ) && mK0 > 0.0 && std::isfinite( mG0 ) && mG0 > 0.0,
                "%s: elastic constants at 0 K read K0 = %g Pa, G0 = %g Pa",
                this->label().c_str(), ( double ) mK0, ( double ) mG0 );

            // the report constant: the Grueneisen parameter as a DIAGNOSTIC, gamma = alpha_V K_S / ( rho cp ),
            // evaluated at the reference temperature of the printed density ( or at the anchor if the
            // material does not reach it ), with the adiabatic bulk modulus recovered exactly from the
            // isothermal one, K_S = K_T / ( 1 - alpha_V^2 T K_T / ( rho cp ) ). It is not prescribed
            // anywhere: it is what the served K, alpha, cp and rho imply, so an inconsistent set shows
            // up here as a value outside the O( 2 ) the metals share
            {
                const real Tg     = BELFEM_TREF <= Tmax ? BELFEM_TREF : T2 ;
                const real Kg     = mK0 * std::exp( -mDeltaK * this->eps( Tg ) );
                const real alphaV = 3.0 * this->alpha( Tg );
                const real rho    = this->density( Tg );
                const real cp     = this->cp( Tg );
                const real q      = alphaV * alphaV * Tg * Kg / ( rho * cp );

                BELFEM_ERROR( q < 1.0,
                    "%s: thermoelastic correction %g at %g K is not physical - check alpha, cp and rho there",
                    this->label().c_str(), ( double ) q, ( double ) Tg );

                const real Ks = Kg / ( 1.0 - q );
                this->set_constant( MaterialProperty::grueneisen, alphaV * Ks / ( rho * cp ) );
            }

            // the served splines sample E_custom and nu_custom on the 4 K grid; both slopes
            // at 0 K are exactly zero because alpha( 0 ) = 0
            this->set_have( MaterialProperty::E );
            this->set_have( MaterialProperty::nu );
            this->create_spline( MaterialProperty::E, 0.0 );
            this->create_spline( MaterialProperty::nu, 0.0 );

            // sampled guards over the whole served range: alpha >= 0, nu inside ( -1, 1/2 ), and nu
            // non-decreasing within data accuracy. The analytic curve is checked at the spline knots,
            // the served spline half-way between them ( a cubic reproduces its knots, so only the
            // midpoints can reveal an overshoot ). With pure exponentials a drop of 1e-3 is a sign
            // error in the constants; the 1e-4 drift of a metal whose isothermal Poisson ratio is
            // flat ( iron ) passes
            const real   dT  = 4.0 ;
            const real   tol = 1e-3 ;
            const size_t n   = static_cast< size_t >( std::ceil( Tmax / dT ) ) + 1 ;
            real numax = -1.0 ;
            real T = 0.0 ;
            for ( size_t k = 0; k < n; ++k )
            {
                BELFEM_ERROR( this->alpha( T ) >= 0.0,
                    "%s: negative thermal expansion %g at %g K - the quasi-harmonic closure assumes alpha >= 0",
                    this->label().c_str(), ( double ) this->alpha( T ), ( double ) T );

                const real nu = this->nu_custom( T );

                BELFEM_ERROR( nu > -1.0 && nu < 0.5,
                    "Poisson ratio of %s reads %g at %g K - check the anchor and the softening constants",
                    this->label().c_str(), ( double ) nu, ( double ) T );

                numax = std::max( numax, nu );

                BELFEM_ERROR( numax - nu <= tol,
                    "Poisson ratio of %s falls with temperature ( %g at %g K, %g before ) - "
                    "the shear modulus must not soften more slowly than the bulk modulus",
                    this->label().c_str(), ( double ) nu, ( double ) T, ( double ) numax );

                const real Tmid = T + 0.5 * dT ;
                if ( Tmid <= Tmax )
                {
                    const real numid = this->nu( Tmid );   // the served spline

                    BELFEM_ERROR( numid > -1.0 && numid < 0.5 && numax - numid <= tol,
                        "served Poisson ratio of %s reads %g at %g K between the spline knots ( %g before )",
                        this->label().c_str(), ( double ) numid, ( double ) Tmid, ( double ) numax );
                }

                T += dT ;
            }
        }

    }
}
