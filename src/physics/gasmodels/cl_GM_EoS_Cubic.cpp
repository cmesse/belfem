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

#include "constants.hpp"
#include "cl_Gas.hpp"
#include "cl_GM_EoS_Cubic.hpp"
#include "cl_GM_EoS_AlphaFunction.hpp"
#include "cl_GM_EoS_AlphaFunctionFactory.hpp"

#include "fn_dot.hpp"
#include "fn_cardano.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"
#include "fn_linspace.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        EoS_Cubic::EoS_Cubic( Gas & aParent, const GasModel aGasModel ) :
                EoS( aParent ),
                mDepartureSpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax )
        {
            switch( aGasModel )
            {
                case( GasModel::SRK ) :
                {
                    this->init_srk();
                    break;
                }
                case( GasModel::PR ) :
                {
                    this->init_pr();
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "unknown gas model" );
                }
            }

            this->init_common();
            this->init_departure_splines();
        }

//----------------------------------------------------------------------------

        EoS_Cubic::~EoS_Cubic()
        {
            for( auto tAlpha : mAlpha )
            {
                delete tAlpha;
            }
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::remix()
        {
            // update the B-Constant for Cubic gas
            if( mParent.number_of_components() == 1 )
            {
                mB( 0 ) = mBc( 0 ) / mM ;
            }
            else
            {
                mB( 0 ) = dot( mParent.molar_fractions(), mBc ) / mM;
            }

            mB( 1 ) = mB( 0 ) * mR1;
            mB( 2 ) = mB( 0 ) * mR2;
            mB( 3 ) = mB( 2 ) - mB( 1 ) ; // for hdep

            // reset values for departure function
            mDepartureV = BELFEM_QUIET_NAN;

            // the a(T) cache is composition dependent: invalidate it,
            // otherwise a query at a previously seen temperature returns
            // the attraction term of the old mixture
            mCubicTemperatures.fill( BELFEM_QUIET_NAN );

            // get deperture coefficients
            Matrix< real > & tDeparture = mDepartureSpline.matrix_data();

            const Vector< real > & tY = mParent.mass_fractions();

            if( tDeparture.n_cols() == 0 )
            {
                // bugfix
                tDeparture.set_size( mDepartureCoefficients( 0 ).n_rows(),
                                     mDepartureCoefficients( 0 ).n_cols(),
                                     0.0 );
            }
            else
            {
                // reset
                tDeparture.fill( 0.0 );
            }



            for( uint k=0; k<mParent.number_of_components(); ++k )
            {
                tDeparture.matrix_data() += tY( k ) * mDepartureCoefficients( k ).matrix_data();
            }
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::p( const real T, const real v ) const
        {
            return mR * T / ( v - mB( 0 ) ) - this->a( T ) /
                    ( ( v - mB( 1 ) ) * ( v - mB( 2 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::v( const real T, const real p ) const
        {
            // coefficient vectors
            // f = a*Z^3 + b*Z^2 + c*Z + d

            // a
            mWorkA( 0 ) = std::pow(-(mR*T), 3 );

            // b
            mWorkA( 1 ) = std::pow(mR*T, 2 )
                          *((mB(0) + mB(1) + mB(2))*p + mR*T);

            // c
            mWorkA( 2 ) =-(p*mR*T*(this->a(T) + mB(1)*mB(2)*p
                       + mB(0)*(mB(1) + mB(2))*p
                       + (mB(1) + mB(2))*mR*T));

            // d
            mWorkA( 3 ) =  std::pow(p, 2 )*(this->a(T)*mB(0)
                       + mB(1)*mB(2)*(mB(0)*p + mR*T));




            /* the cubic is always solved by Cardano; the Newton
             * alternative that used to sit behind mUseAlwaysCardano was
             * unreachable and had not been kept in step - it carried
             * neither the liquid root selection nor the Z > B filter */
            cardano( mWorkA, mWorkZ );

            // cardano may return one, two or three real roots
            const uint tNumRoots = mWorkZ.length();

            if ( mParent.is_liquid() )
            {
                sort( mWorkZ );

                /* the physical liquid root is the smallest one above the
                 * covolume, Z > B ( i.e. v > b ); a root below it sits on
                 * the unphysical branch and would poison chi() with a
                 * negative log argument */
                const real tB = p * mB( 0 ) / ( mR * T );

                for( uint k=0; k<tNumRoots; ++k )
                {
                    if( mWorkZ( k ) > tB )
                    {
                        return mWorkZ( k ) * mR * T / p;
                    }
                }

                BELFEM_ERROR( false,
                        "No physical liquid root for T=%f and p=%f",
                        ( float ) T, ( float ) p );

                return BELFEM_QUIET_NAN;
            }
            else
            {
                // largest root is the gas one; it must be positive
                real tZ = max( mWorkZ );

                BELFEM_ERROR( tZ > 0.0,
                        "No physical gas root for T=%f and p=%f",
                        ( float ) T, ( float ) p );

                return tZ * mR * T / p;
            }
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::T( const real p, const real v ) const
        {
            // a silent NaN here would be cached as a valid stateval by Gas::T
            BELFEM_ERROR( false,
                    "T( p, v ) is not implemented for the cubic equation of state" );

            return BELFEM_QUIET_NAN;
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::dpdT( const real T, const real v ) const
        {
            return  mR / ( v - mB( 0 ) ) - this->dadT( T )  /
                    ( ( v - mB( 1 ) ) * ( v - mB( 2 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::d2pdT2( const real T, const real v ) const
        {
            return -this->d2adT2( T )  /
                   ( ( v - mB( 1 ) ) * ( v - mB( 2 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::dpdv( const real T, const real v ) const
        {
            return -mR*T/std::pow( ( v - mB( 0 ) ), 2 )
                   + this->a( T ) * ( 2.0 * v  - ( mB( 1 ) + mB( 2 ) ) ) /
                   std::pow( ( v - mB( 1 ) ) * ( v - mB( 2 )), 2 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::d2pdv2( const real T, const real v ) const
        {
            return 2.0 * ( mR*T/std::pow( ( v - mB( 0 ) ), 3 )
                           - this->a( T ) / mB( 3 )
                             * (   std::pow( mB( 1 ) -v , -3 )
                                 + std::pow( v -  mB( 2 ), -3 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::alpha( const real T, const real p ) const
        {
            mStatevals.update_Tp( T, p );

            if( ! mStatevals.test( BELFEM_STATEVAL_ALPHA ) )
            {
                mStatevals.set( BELFEM_STATEVAL_ALPHA,
                p * this->beta( T, p ) * this->kappa( T, p ) );
            }

            return mStatevals.get( BELFEM_STATEVAL_ALPHA );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::beta( const real T, const real p ) const
        {
            mStatevals.update_Tp( T, p );

            if( ! mStatevals.test( BELFEM_STATEVAL_BETA ) )
            {
                // 10.18419/opus-9381 ( 2.6 )
                mStatevals.set( BELFEM_STATEVAL_BETA,
                this->dpdT( T, mParent.v( T, p ) ) / p );
            }

            return mStatevals.get( BELFEM_STATEVAL_BETA );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::kappa( const real T, const real p ) const
        {
            mStatevals.update_Tp( T, p );

            if( ! mStatevals.test( BELFEM_STATEVAL_KAPPA ) )
            {
                // 10.18419/opus-9381 ( 2.7 )
                real tV = mParent.v( T, p );

                mStatevals.set( BELFEM_STATEVAL_KAPPA,
                                -1.0 / ( tV * this->dpdv( T, tV ) ) );
            }

            return mStatevals.get( BELFEM_STATEVAL_KAPPA );
        }

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

        real
        EoS_Cubic::hdep( const real T, const real p ) const
        {
            real tV = mParent.v( T, p );

            return p * tV - mR * T
                + ( T * this->dadT( T ) - this->a( T )) *  this->chi( tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::cpdep( const real T, const real p ) const
        {
            // get specific volume from parent
            real tV =  mParent.v( T, p );
            real tdVdT = this->alpha( T, p ) * tV;

            return p*tdVdT - mR + ( T * this->dadT( T )
                    - this->a( T )) * this->dchidT( T, p, tV )
                    + T * this->d2adT2( T )*this->chi( tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::sdep( const real T, const real p ) const
        {
            // get specific volume from parent
            real tV = mParent.v( T, p );

            return mR*std::log( p  * ( tV - mB( 0 ) )/ ( mR * T ) )
                        + this->dadT( T )
                        * this->chi( tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dhdepdp( const real T, const real p ) const
        {
            // get specific volume from parent
            real tV = mParent.v( T, p );

            real tdVdP =  - mParent.kappa( T, p ) * tV ;

            // p * tV - mR * T + ( T * this->dadT( T ) - this->a( T )) *
            //                    this->chi( tV );

            return  tV + p * tdVdP + ( T * this->dadT( T ) - this->a( T )) *
                this->dchidp( T, p , tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dsdepdT( const real T, const real p ) const
        {
            real tV = mParent.v( T, p );
            real tdVdT = this->alpha( T, p ) * tV;

            return mR * ( - tdVdT/( mB( 0 ) - tV ) - 1.0 / T )
                + this->d2adT2( T ) * this->chi( tV )
                + this->dadT( T ) * this->dchidT(  T, p, tV ) ;
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dsdepdp( const real T, const real p ) const
        {
            real tV = mParent.v( T, p );

            //   R (1/p + Dt[v, p]/v)
            //  + (b^3 R r1 r2 + b (dadT - b R (r1 + r2)) v + (-dadT + b R) v^2)
            //  * Dt[v,p])/(v (-b + v) (-b r1 + v) (-b r2 + v))

            // ( 2.7 )
            real tdVdP = -  mParent.kappa( T, p ) * tV ;

            // v^2 \, \left( b \, R - \frac{\partial a}{\partial T}\right)
            real tResult = tV * ( mB( 0 ) * mR - this->dadT( T )  ) ;

            // v \, b \, \left[ \frac{\partial a}{\partial T} - b \, R \, \left( r_1 + r_2 \right) \right]
            tResult += mB( 0 ) * ( this->dadT( T ) - mR * ( mB( 1 ) + mB( 2 ) ) ) ;

            tResult *= tV;

            // R \, b^3 \, r_1 \, r_2
            tResult +=  mR * mB( 0 ) * mB( 1 ) * mB( 2 );

            tResult *= tdVdP;

            tResult /= tV * ( tV - mB( 0 ) ) * ( tV - mB( 1 ) ) * ( tV - mB( 2 ) ) ;

            tResult += mR * ( 1.0 / p + tdVdP / tV );

            return  tResult;
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::hdep0( const real T ) const
        {
            return mDepartureSpline.eval( T );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::cpdep0( const real T ) const
        {
            return mDepartureSpline.deval( T );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::sdep0( const real T ) const
        {
            return mDepartureSpline.entropy( T );
        }

//------------------------------------------------------------------------------

        // temperature derivative of entropy departure
        real
        EoS_Cubic::dsdepdT0( const real T ) const
        {
            return mDepartureSpline.dentropy( T );
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::init_srk()
        {
            // parameters for Soave-Redlich-Kwong
            real tX   = std::pow( 2.0, 1.0/3.0 ) - 1.0;

            mR1     = 0.0;
            mR2     = -1.0;
            mOmegaA = 1.0/( 9.0 * tX );
            mOmegaB = tX/3.0;

            mAlpha.set_size( mParent.number_of_components(), nullptr );

            // create factory
            AlphaFunctionFactory tFactory;

            for( uint k=0; k<mParent.number_of_components(); ++k )
            {
                const gastables::GasData * tData = mParent.data( k );

                if( tData->has_cubic() )
                {
                    mAlpha( k ) = tFactory.create_pm_srk( tData );
                }
                else if( tData->has_crit() )
                {
                    mAlpha( k ) = tFactory.create_ccr_mc_srk( tData );
                }
                else
                {
                    mAlpha( k ) = tFactory.create_empty();
                }
            }
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::init_pr()
        {
            real tX      =     std::pow( 6.0*std::sqrt( 2.0 ) + 8.0, 1.0/3.0 )
                             - std::pow( 6.0*std::sqrt( 2.0 ) - 8.0, 1.0/3.0 )
                             - 1.0;

            // parameters for Peng-Robinson
            mR1  = -1.0 - std::sqrt( 2.0 );
            mR2  = -1.0 + std::sqrt( 2.0 );

            mOmegaA = ( 40.0*tX + 24.0 )/( 147.0 - 37.0 * tX );

            mOmegaB = tX / ( tX + 9.0 );

            mAlpha.set_size( mParent.number_of_components(), nullptr );

            // create factory
            AlphaFunctionFactory tFactory;

            for( uint k=0; k< mParent.number_of_components(); ++k )
            {

                const gastables::GasData * tData = mParent.data( k );

                if( tData->has_cubic() )
                {
                    mAlpha( k ) = tFactory.create_pm_pr( tData );
                }
                else if( tData->has_crit() )
                {
                    mAlpha( k ) = tFactory.create_ccr_pr( tData );
                }
                else
                {
                    mAlpha( k ) = tFactory.create_empty();
                }
            }
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::init_common()
        {
            // save a factors
            mAc.set_size( mParent.number_of_components() );
            mBc.set_size( mParent.number_of_components() );

            for( uint k=0; k< mParent.number_of_components(); ++k )
            {
                gastables::GasData * tData = mParent.data( k );

                if( tData->has_crit() )
                {
                    mAc( k ) = mOmegaA * std::pow( constant::Rm
                                                   * tData->T_crit(), 2 ) /
                               tData->p_crit();

                    mBc( k ) = ( mOmegaB * constant::Rm
                                 * tData->T_crit() ) / tData->p_crit();
                }
                else
                {
                    mAc( k ) = 0.0;
                    mBc( k ) = 0.0;
                }

            }

            mA.set_size( mParent.number_of_components(), 0.0 );
            mB.set_size( 4, 0.0 );

            mdAdT.set_size( mParent.number_of_components(), 0.0 );
            md2AdT2.set_size( mParent.number_of_components(), 0.0 );

            mCubicStatevals.set_size( 3, 0 );
            mCubicTemperatures.set_size( 3, 0 );

            mWorkA.set_size( 4 );

            // mVM.set_size( mParent.number_of_components(), 0.0 );

            mComponentHDEP.set_size( mParent.number_of_components(), 0.0 );
            mComponentCPDEP.set_size( mParent.number_of_components(), 0.0 );
            mComponentV.set_size( mParent.number_of_components(), 0.0 );

        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::init_departure_splines()
        {
            // create a help matrix for the splines
            SpMatrix tHelpMatrix;
            spline::create_helpmatrix(
                    gastables::gNumberOfSplinePoints,
                    gastables::gDeltaT,
                    tHelpMatrix );

            Matrix< real > tEmpty;

            // allocate memory for cell
            mDepartureCoefficients.set_size( mParent.number_of_components(), tEmpty );

            // temperature steps
            Vector< real > tT;

            linspace(
                    0.0,
                    gastables::gTmax,
                    gastables::gNumberOfSplinePoints,
                    tT );

            Vector< real > tValues( gastables::gNumberOfSplinePoints );

            // loop over all gases
            for( uint g=0; g<mParent.number_of_components(); ++g )
            {
                // get data pointer of component
                real tV;
                real tCPDEP;

                // mDepartureCoefficients( g ) = new Matrix<real>
                if ( mParent.data( g )->has_crit() )
                {
                    /**
                     * the following lines create the departure function for each gas
                     * at reference pressure.
                     */

                    // loop over all temperature steps
                    for ( uint k = 1; k < gastables::gNumberOfSplinePoints; ++k )
                    {
                        this->component_wise_parameters(
                                g,
                                tT( k ),
                                gastables::gPref,
                                tV,
                                tValues( k ),
                                tCPDEP );
                    }

                    // extrapolate first step
                    tValues( 0 ) = tValues( 1 ) - ( tValues( 2 ) - tValues( 1 ));


                    /* The entropy row of this spline is sdep0, an entropy
                     * DEPARTURE, so it is anchored at zero. Anchoring it at the
                     * standard state entropy of the component would add that
                     * constant to every sdep0 call and, through
                     * Gas::realgas_s, subtract the standard state entropy from
                     * the absolute entropy of every cubic gas. */
                    mDepartureSpline.update_data( tHelpMatrix, tValues,
                            spline::SplineBC::NoCurvature, spline::SplineBC::NoCurvature,
                            0.0, 0.0,
                            gastables::gTref,
                            0.0 );

                    // copy date into matrix
                    mDepartureCoefficients( g ) = mDepartureSpline.matrix_data();
                }
                else
                {
                    mDepartureCoefficients( g ).set_size( 5, gastables::gNumberOfSplinePoints, 0.0 );
                }
            }
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::eval_a( const real T, const int aDeriv ) const
        {

            uint tNumberOfComponents = mParent.number_of_components();
            const Vector< real > & tMolarFractions = mParent.molar_fractions();

            if( mCubicTemperatures( 0 ) != T )
            {
                // remember this temperature
                mCubicTemperatures( 0 ) = T;

                // evaluate A
                for ( uint k = 0; k < tNumberOfComponents; ++k )
                {
                    mA( k ) = mAc( k ) * mAlpha( k )->alpha( T );
                }

                // store value
                mCubicStatevals( 0 ) = 0.0;

                for ( uint i = 0; i < tNumberOfComponents; ++i )
                {
                    for ( uint j = 0; j < tNumberOfComponents; ++j )
                    {
                        mCubicStatevals( 0 ) +=
                                tMolarFractions( i ) * tMolarFractions( j )
                                * std::sqrt( mA( i ) * mA( j ) );
                    }
                }

                mCubicStatevals( 0 ) /= std::pow( mM, 2 );
            }

            if( aDeriv >= 1 )
            {
                if( mCubicTemperatures( 1 ) != T )
                {
                    // remember this temperature
                    mCubicTemperatures( 1 ) = T;

                    // evaluate dAdT
                    for ( uint k = 0; k < tNumberOfComponents; ++k )
                    {
                        mdAdT( k ) = mAc( k ) * mAlpha( k )->dalphadT( T );
                    }

                    mCubicStatevals( 1 ) = 0.0;

                    for ( uint i = 0; i < tNumberOfComponents; ++i )
                    {
                        if( std::abs( mA( i ) ) > BELFEM_EPSILON )
                        {
                            for ( uint j = 0; j < tNumberOfComponents; ++j )
                            {
                                if( std::abs( mA( j ) ) > BELFEM_EPSILON )
                                {
                                    mCubicStatevals( 1 ) +=
                                            tMolarFractions( i ) * tMolarFractions( j )
                                            * ( mdAdT( i ) * mA( j ) + mA( i ) * mdAdT( j ) ) /
                                            std::sqrt( mA( i ) * mA( j ) );
                                }
                            }
                        }

                    }

                    mCubicStatevals( 1 ) /= 2.0 * std::pow( mM, 2 );
                }

                if( aDeriv == 2 )
                {
                    if( mCubicTemperatures( 2 ) != T )
                    {
                        // remember this temperature
                        mCubicTemperatures( 2 ) = T;

                        // evaluate dAdT
                        for ( uint k = 0; k < tNumberOfComponents; ++k )
                        {
                            md2AdT2( k ) = mAc( k ) * mAlpha( k )->d2alphadT2( T );
                        }

                        mCubicStatevals( 2 ) = 0.0;

                        for ( uint i = 0; i < tNumberOfComponents; ++i )
                        {
                            if( std::abs( mA( i ) ) > BELFEM_EPSILON )
                            {
                                for ( uint j = 0; j < tNumberOfComponents; ++j )
                                {
                                    if( std::abs( mA( j ) ) > BELFEM_EPSILON )
                                    {
                                        mCubicStatevals( 2 ) +=
                                                (-std::pow(mA(j)*mdAdT(i) + mA(i)*mdAdT(j),2) +
                                                 2.0*mA(i)*mA(j)*(
                                                         mA(j)*md2AdT2(i) +2.0*mdAdT(i)*mdAdT(j)
                                                         + mA(i)*md2AdT2(j)) ) *
                                                (tMolarFractions(i)*tMolarFractions(j))/std::pow(mA(i)*mA(j),1.5);
                                    }
                                }
                            }

                        }

                        mCubicStatevals( 2 ) /= 4.0 * std::pow( mM, 2 );
                    }
                }
            }
        }

        real
        EoS_Cubic::a( const real T ) const
        {
            this->eval_a( T, 0 );

            // return value
            return mCubicStatevals( 0 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::dadT( const real T ) const
        {

            this->eval_a( T, 1 );

            // return value
            return mCubicStatevals( 1 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::d2adT2( const real T ) const
        {
            this->eval_a( T, 2 );

            // return value
            return mCubicStatevals( 2 );
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::eval_critical_point( real & T, real & p, real & v ) const
        {
            // guess T_crit using Kay's method
            T = 0.0;

            const Vector< real > & tY = mParent.molar_fractions();

            for( uint k=0; k<tY.length(); ++k )
            {
                T += tY( k ) * mParent.data( k )->T_crit();
            }

            real tT = 0;
            real tdPdV = 1e12;

            // Run Newton Loop
            uint tCount = 0;
            while( std::abs( tT - T ) > BELFEM_EPSILON_T && std::abs( tdPdV ) > 1 )
            {
                // see ( 2. 4 )
                p = mOmegaB * mR * T / mB( 0 );
                real tF  = this->a( T ) - mOmegaA * std::pow( mR * T, 2 ) / p;
                real tdF = this->dadT( T ) - mB( 0 )*mR * mOmegaA / mOmegaB;

                tT = T;
                T -= 0.9 * tF/tdF;

                // evaluate v
                v = this->v( T, p );

                // perform funcition test
                tdPdV = this->dpdv( T, v );

                ++tCount;
                BELFEM_ERROR( tCount < 1000, "Too many Iterations while trying to calculate T_crit" );
            }

            // the loop can exit on the dpdv criterion one step after the last
            // temperature update: make p and v consistent with the final T
            p = mOmegaB * mR * T / mB( 0 );
            v = this->v( T, p );
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::component_wise_parameters(
                const uint aIndex,
                const real T,
                const real p,
                real & v,
                real & aHDEP,
                real & aCPDEP ) const
        {
            const real & tM = mParent.data( aIndex )->M();
            const real & tR = mParent.data( aIndex )->R();

            real tB = mBc( aIndex ) / tM;
            real tB1 = tB * mR1;
            real tB2 = tB * mR2;

            // alpha function and derivatives
            real tA = mAc( aIndex ) * mAlpha( aIndex )->alpha( T ) / ( tM * tM );

            real tdAdT = mAc( aIndex ) * mAlpha( aIndex )->dalphadT( T ) / ( tM * tM );

            real td2AdT2 =  mAc( aIndex ) * mAlpha( aIndex )->d2alphadT2( T ) / ( tM * tM );

            // step 1: calculate volume

            // a
            mWorkA( 0 ) = std::pow( -( tR * T ), 3 );

            // b
            mWorkA( 1 ) = std::pow( tR * T, 2 )  * ( ( tB + tB1 + tB2 ) * p + tR * T );

            // c
            mWorkA( 2 ) = -( p * tR * T * ( tA + tB1 * tB2 * p
                            + tB * ( tB1 + tB2 ) * p
                            + ( tB1 + tB2 ) * tR * T ));

            // d
            mWorkA( 3 ) =   std::pow( p, 2 ) * ( tA * tB
                          + tB1 * tB2 * ( tB * p + tR * T ));

            cardano( mWorkA, mWorkZ );

            // density
            v = max( mWorkZ ) * tR * T / p;

            // derivatives
            real tdPdT = tR / ( v - tB ) - tdAdT  /
                       ( ( v - tB1 ) * ( v - tB2 ) );

            real tdPdV = -tR*T/std::pow( ( v - tB ), 2 )
                        + tA * ( 2.0 * v  - ( tB1 + tB2 ) ) /
                         std::pow( ( v - tB1 ) * ( v - tB2), 2 );

            real tdVdT = -tdPdT / tdPdV;

            // chi function
            real tChi = std::log( ( v - tB1 ) / ( v - tB2 ) ) / ( tB2 - tB1 );

            // derivative of chi function
            real tdChidT = -tdVdT / ( ( v - tB1 ) * ( v - tB2 ) );


            // calculate enthalpy departure
            aHDEP = p * v - tR * T + ( T  * tdAdT - tA ) * tChi;

            // departure of specific heat
            aCPDEP = p * tdVdT - tR + ( T * tdAdT - tA ) * tdChidT + T * tChi * td2AdT2 ;

        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::update_component_parameters(
                const real T,
                const real p ) const
        {
            if( T != mComponentT || p != mComponentP )
            {
                for ( uint k = 0; k < mParent.number_of_components(); ++k )
                {
                    if ( mParent.data( k )->has_crit() )
                    {
                        this->component_wise_parameters(
                                k,
                                T,
                                p,
                                mComponentV( k ),
                                mComponentHDEP( k ),
                                mComponentCPDEP( k ));
                    }
                }

                mComponentT = T;
                mComponentP = p;

                mComponentCol = mDepartureSpline.find_col( T );
            }
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::v( const uint aIndex, const real T, const real p ) const
        {
            this->update_component_parameters( T, p );
            return mComponentV( aIndex );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::hdep( const uint aIndex, const real T, const real p ) const
        {
            this->update_component_parameters( T, p );

            const Matrix< real > & tData = mDepartureCoefficients( aIndex );

            return mComponentHDEP( aIndex ) -
                ( (   ( tData( 0, mComponentCol )   * T
                      + tData( 1, mComponentCol ) ) * T
                      + tData( 2, mComponentCol ) ) * T
                      + tData( 3, mComponentCol ) );

        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::cpdep( const uint aIndex, const real T, const real p ) const
        {
            this->update_component_parameters( T, p );

            const Matrix< real > & tData = mDepartureCoefficients( aIndex );

            return mComponentCPDEP( aIndex ) -
                      ( ( ( 3.0 * tData( 0, mComponentCol )   * T
                          + 2.0 * tData( 1, mComponentCol ) ) * T
                          +       tData( 2, mComponentCol ) ) );

        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::chi( const real v ) const
        {
            if ( mDepartureV != v )
            {
                mDepartureValue = std::log( ( v - mB( 1 ) )/ ( v -  mB( 2 ) )) / mB( 3 );
                mDepartureV = v;
            }
            return mDepartureValue;
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dchidT( const real T, const real p, const real v  ) const
        {
            return -this->alpha( T, p ) * v / ( ( v - mB( 1 ) ) * ( v - mB( 2 ) ) );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dchidp( const real T, const real p, const real v  ) const
        {

            real tdVdP =  - mParent.kappa( T, p ) * v ;

            return -tdVdP / ( ( v - mB( 1 ) ) * ( v - mB( 2 ) ) );
        }

//------------------------------------------------------------------------------
    }
}