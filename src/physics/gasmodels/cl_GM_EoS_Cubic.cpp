//
// Created by Christian Messe on 13.09.19.
//

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
        EoS_Cubic::p( const real & aT, const real & aV )
        {
            return mR * aT / ( aV - mB( 0 ) ) - this->a( aT ) /
                    ( ( aV - mB( 1 ) ) * ( aV - mB( 2 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::v( const real & aT, const real & aP )
        {
            // coefficient vectors
            // f = a*Z^3 + b*Z^2 + c*Z + d

            // a
            mWorkA( 0 ) = std::pow(-(mR*aT), 3 );

            // b
            mWorkA( 1 ) = std::pow(mR*aT, 2 )
                          *((mB(0) + mB(1) + mB(2))*aP + mR*aT);

            // c
            mWorkA( 2 ) =-(aP*mR*aT*(this->a(aT) + mB(1)*mB(2)*aP
                       + mB(0)*(mB(1) + mB(2))*aP
                       + (mB(1) + mB(2))*mR*aT));

            // d
            mWorkA( 3 ) =  std::pow(aP, 2 )*(this->a(aT)*mB(0)
                       + mB(1)*mB(2)*(mB(0)*aP + mR*aT));




            if( mUseAlwaysCardano || mParent.is_liquid() )
            {
                cardano( mWorkA, mWorkZ );

                if ( mParent.is_liquid() )
                {
                    sort( mWorkZ );

                    if( mWorkZ( 0 ) < 1e-6 )
                    {
                        if( mWorkZ( 1 ) < 1e-6 )
                        {
                            return mWorkZ( 2 ) * mR * aT / aP;
                        }
                        else
                        {
                            return mWorkZ( 1 ) * mR * aT / aP;
                        }
                    }
                    else
                    {
                        return mWorkZ( 0 ) * mR * aT / aP;
                    }
                }
                else
                {
                    return max( mWorkZ ) * mR * aT / aP;
                }
            }
            else
            {
                real tZ0 = 0.0;

                // initial guess
                real tZ = 1.0;

                uint tCount = 0;

                while( std::abs( tZ - tZ0 ) > 1e-6 )
                {
                    tZ0 = tZ;
                    real tF = ( (  mWorkA( 0 ) * tZ
                                 + mWorkA( 1 ) ) * tZ
                                 + mWorkA( 2 ) ) * tZ
                                 + mWorkA( 3 );

                    real tdFdZ =
                              ( 3.0 * mWorkA( 0 ) * tZ
                              + 2.0 * mWorkA( 1 )) * tZ
                                    + mWorkA( 2 );

                    tZ -= tF / tdFdZ;

                    tCount++;

                    BELFEM_ERROR( tCount < 100, "Created infinite loop for T=%f and p=%f",
                                 ( float ) aT, ( float ) aP );
                }

                return tZ * mR * aT / aP;
            }
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::T( const real & aP, const real & aV )
        {
            return BELFEM_QUIET_NAN;
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::dpdT( const real & aT, const real & aV )
        {
            return  mR / ( aV - mB( 0 ) ) - this->dadT( aT )  /
                    ( ( aV - mB( 1 ) ) * ( aV - mB( 2 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::d2pdT2( const real & aT, const real & aV )
        {
            return -this->d2adT2( aT )  /
                   ( ( aV - mB( 1 ) ) * ( aV - mB( 2 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::dpdv( const real & aT, const real & aV )
        {
            return -mR*aT/std::pow( ( aV - mB( 0 ) ), 2 )
                   + this->a( aT ) * ( 2.0 * aV  - ( mB( 1 ) + mB( 2 ) ) ) /
                   std::pow( ( aV - mB( 1 ) ) * ( aV - mB( 2 )), 2 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::d2pdv2( const real & aT, const real & aV )
        {
            return 2.0 * ( mR*aT/std::pow( ( aV - mB( 0 ) ), 3 )
                           - this->a( aT ) / mB( 3 )
                             * (   std::pow( mB( 1 ) -aV , -3 )
                                 + std::pow( aV -  mB( 2 ), -3 ) ) );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::alpha( const real & aT, const real & aP )
        {
            mStatevals.update_Tp( aT, aP );

            if( ! mStatevals.test( BELFEM_STATEVAL_KAPPA ) )
            {
                mStatevals.set( BELFEM_STATEVAL_ALPHA,
                aP * this->beta( aT, aP ) * this->kappa( aT, aP ) );
            }

            return mStatevals.get( BELFEM_STATEVAL_ALPHA );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::beta( const real & aT, const real & aP )
        {
            mStatevals.update_Tp( aT, aP );

            if( ! mStatevals.test( BELFEM_STATEVAL_BETA ) )
            {
                // 10.18419/opus-9381 ( 2.6 )
                mStatevals.set( BELFEM_STATEVAL_BETA,
                this->dpdT( aT, mParent.v( aT, aP ) ) / aP );
            }

            return mStatevals.get( BELFEM_STATEVAL_BETA );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::kappa( const real & aT, const real & aP )
        {
            mStatevals.update_Tp( aT, aP );

            if( ! mStatevals.test( BELFEM_STATEVAL_KAPPA ) )
            {
                // 10.18419/opus-9381 ( 2.7 )
                real tV = mParent.v( aT, aP );

                mStatevals.set( BELFEM_STATEVAL_KAPPA,
                                -1.0 / ( tV * this->dpdv( aT, tV ) ) );
            }

            return mStatevals.get( BELFEM_STATEVAL_KAPPA );
        }

//------------------------------------------------------------------------------
// Departure Functions
//------------------------------------------------------------------------------

        real
        EoS_Cubic::hdep( const real & aT, const real & aP )
        {
            real tV = mParent.v( aT, aP );

            return aP * tV - mR * aT
                + ( aT * this->dadT( aT ) - this->a( aT )) *  this->chi( tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::cpdep( const real & aT, const real & aP )
        {
            // get specific volume from parent
            real tV =  mParent.v( aT, aP );
            real tdVdT = this->alpha( aT, aP ) * tV;

            return aP*tdVdT - mR + ( aT * this->dadT( aT )
                    - this->a( aT )) * this->dchidT( aT, aP, tV )
                    + aT * this->d2adT2( aT )*this->chi( tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::sdep( const real & aT, const real & aP )
        {
            // get specific volume from parent
            real tV = mParent.v( aT, aP );

            return mR*std::log( aP  * ( tV - mB( 0 ) )/ ( mR * aT ) )
                        + this->dadT( aT )
                        * this->chi( tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dhdepdp( const real & aT, const real & aP )
        {
            // get specific volume from parent
            real tV = mParent.v( aT, aP );

            real tdVdP =  - mParent.kappa( aT, aP ) * tV ;

            // aP * tV - mR * aT + ( aT * this->dadT( aT ) - this->a( aT )) *
            //                    this->chi( tV );

            return  tV + aP * tdVdP + ( aT * this->dadT( aT ) - this->a( aT )) *
                this->dchidp( aT, aP , tV );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dsdepdT( const real & aT, const real & aP )
        {
            real tV = mParent.v( aT, aP );
            real tdVdT = this->alpha( aT, aP ) * tV;

            return mR * ( - tdVdT/( mB( 0 ) - tV ) - 1.0 / aT )
                + this->d2adT2( aT ) * this->chi( tV )
                + this->dadT( aT ) * this->dchidT(  aT, aP, tV ) ;
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dsdepdp( const real & aT, const real & aP )
        {
            real tV = mParent.v( aT, aP );

            //   R (1/p + Dt[v, p]/v)
            //  + (b^3 R r1 r2 + b (dadT - b R (r1 + r2)) v + (-dadT + b R) v^2)
            //  * Dt[v,p])/(v (-b + v) (-b r1 + v) (-b r2 + v))

            // ( 2.7 )
            real tdVdP = -  mParent.kappa( aT, aP ) * tV ;

            // v^2 \, \left( b \, R - \frac{\partial a}{\partial T}\right)
            real aResult = tV * ( mB( 0 ) * mR - this->dadT( aT )  ) ;

            // v \, b \, \left[ \frac{\partial a}{\partial T} - b \, R \, \left( r_1 + r_2 \right) \right]
            aResult += mB( 0 ) * ( this->dadT( aT ) - mR * ( mB( 1 ) + mB( 2 ) ) ) ;

            aResult *= tV;

            // R \, b^3 \, r_1 \, r_2
            aResult +=  mR * mB( 0 ) * mB( 1 ) * mB( 2 );

            aResult *= tdVdP;

            aResult /= tV * ( tV - mB( 0 ) ) * ( tV - mB( 1 ) ) * ( tV - mB( 2 ) ) ;

            aResult += mR * ( 1.0 / aP + tdVdP / tV );

            return  aResult;
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::hdep0( const real & aT )
        {
            return mDepartureSpline.eval( aT );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::cpdep0( const real & aT )
        {
            return mDepartureSpline.deval( aT );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::sdep0( const real & aT )
        {
            return mDepartureSpline.entropy( aT );
        }

//------------------------------------------------------------------------------

        // temperature derivative of entropy departure
        real
        EoS_Cubic::dsdepdT0( const real & aT )
        {
            return mDepartureSpline.dentropy( aT );
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


                    mDepartureSpline.update_data( tHelpMatrix, tValues, gastables::gTref,
                            mParent.data( g )->Sref() / mParent.data( g )->M() );

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
        EoS_Cubic::eval_a( const real & aT, const int aDeriv )
        {

            uint tNumberOfComponents = mParent.number_of_components();
            const Vector< real > & tMolarFractions = mParent.molar_fractions();

            if( mCubicTemperatures( 0 ) != aT )
            {
                // remember this temperature
                mCubicTemperatures( 0 ) = aT;

                // evaluate A
                for ( uint k = 0; k < tNumberOfComponents; ++k )
                {
                    mA( k ) = mAc( k ) * mAlpha( k )->alpha( aT );
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
                if( mCubicTemperatures( 1 ) != aT )
                {
                    // remember this temperature
                    mCubicTemperatures( 1 ) = aT;

                    // evaluate dAdT
                    for ( uint k = 0; k < tNumberOfComponents; ++k )
                    {
                        mdAdT( k ) = mAc( k ) * mAlpha( k )->dalphadT( aT );
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
                    if( mCubicTemperatures( 2 ) != aT )
                    {
                        // remember this temperature
                        mCubicTemperatures( 2 ) = aT;

                        // evaluate dAdT
                        for ( uint k = 0; k < tNumberOfComponents; ++k )
                        {
                            md2AdT2( k ) = mAc( k ) * mAlpha( k )->d2alphadT2( aT );
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
        EoS_Cubic::a( const real & aT )
        {
            this->eval_a( aT, 0 );

            // return value
            return mCubicStatevals( 0 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::dadT( const real & aT )
        {

            this->eval_a( aT, 1 );

            // return value
            return mCubicStatevals( 1 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::d2adT2( const real & aT )
        {
            this->eval_a( aT, 2 );

            // return value
            return mCubicStatevals( 2 );
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::eval_critical_point( real & aT, real & aP, real & aV )
        {
            // guess T_crit using Kay's method
            aT = 0.0;

            const Vector< real > & tY = mParent.molar_fractions();

            for( uint k=0; k<tY.length(); ++k )
            {
                aT += tY( k ) * mParent.data( k )->T_crit();
            }

            real tT = 0;
            real tdPdV = 1e12;

            // Run Newton Loop
            uint tCount = 0;
            while( std::abs( tT - aT ) > BELFEM_EPSILON_T && std::abs( tdPdV ) > 1 )
            {
                // see ( 2. 4 )
                aP = mOmegaB * mR * aT / mB( 0 );
                real tF  = this->a( aT ) - mOmegaA * std::pow( mR * aT, 2 ) / aP;
                real tdF = this->dadT( aT ) - mB( 0 )*mR * mOmegaA / mOmegaB;

                tT = aT;
                aT -= 0.9 * tF/tdF;

                // evaluate v
                aV = this->v( aT, aP );

                // perform funcition test
                tdPdV = this->dpdv( aT, aV );

                ++tCount;
                BELFEM_ERROR( tCount < 1000, "Too many Iterations while trying to calculate T_crit" );
            }
        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::component_wise_parameters(
                const uint & aIndex,
                const real & aT,
                const real & aP,
                real & aV,
                real & aHDEP,
                real & aCPDEP )
        {
            const real & tM = mParent.data( aIndex )->M();
            const real & tR = mParent.data( aIndex )->R();

            real tB = mBc( aIndex ) / tM;
            real tB1 = tB * mR1;
            real tB2 = tB * mR2;

            // alpha function and derivatives
            real tA = mAc( aIndex ) * mAlpha( aIndex )->alpha( aT ) / ( tM * tM );

            real tdAdT = mAc( aIndex ) * mAlpha( aIndex )->dalphadT( aT ) / ( tM * tM );

            real td2AdT2 =  mAc( aIndex ) * mAlpha( aIndex )->d2alphadT2( aT ) / ( tM * tM );

            // step 1: calculate volume

            // a
            mWorkA( 0 ) = std::pow( -( tR * aT ), 3 );

            // b
            mWorkA( 1 ) = std::pow( tR * aT, 2 )  * ( ( tB + tB1 + tB2 ) * aP + tR * aT );

            // c
            mWorkA( 2 ) = -( aP * tR * aT * ( tA + tB1 * tB2 * aP
                            + tB * ( tB1 + tB2 ) * aP
                            + ( tB1 + tB2 ) * tR * aT ));

            // d
            mWorkA( 3 ) =   std::pow( aP, 2 ) * ( tA * tB
                          + tB1 * tB2 * ( tB * aP + tR * aT ));

            cardano( mWorkA, mWorkZ );

            // density
            aV = max( mWorkZ ) * tR * aT / aP;

            // derivatives
            real tdPdT = tR / ( aV - tB ) - tdAdT  /
                       ( ( aV - tB1 ) * ( aV - tB2 ) );

            real tdPdV = -tR*aT/std::pow( ( aV - tB ), 2 )
                        + tA * ( 2.0 * aV  - ( tB1 + tB2 ) ) /
                         std::pow( ( aV - tB1 ) * ( aV - tB2), 2 );

            real tdVdT = -tdPdT / tdPdV;

            // chi function
            real tChi = std::log( ( aV - tB1 ) / ( aV - tB2 ) ) / ( tB2 - tB1 );

            // derivative of chi function
            real tdChidT = -tdVdT / ( ( aV - tB1 ) * ( aV - tB2 ) );


            // calculate enthalpy departure
            aHDEP = aP * aV - tR * aT + ( aT  * tdAdT - tA ) * tChi;

            // departure of specific heat
            aCPDEP = aP * tdVdT - tR + ( aT * tdAdT - tA ) * tdChidT + aT * tChi * td2AdT2 ;

        }

//----------------------------------------------------------------------------

        void
        EoS_Cubic::update_component_parameters(
                const real & aT,
                const real & aP )
        {
            if( aT != mComponentT || aP != mComponentP )
            {
                for ( uint k = 0; k < mParent.number_of_components(); ++k )
                {
                    if ( mParent.data( k )->has_crit() )
                    {
                        this->component_wise_parameters(
                                k,
                                aT,
                                aP,
                                mComponentV( k ),
                                mComponentHDEP( k ),
                                mComponentCPDEP( k ));
                    }
                }

                mComponentT = aT;
                mComponentP = aP;

                mComponentCol = mDepartureSpline.find_col( aT );
            }
        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::v( const uint & aIndex, const real & aT, const real & aP )
        {
            this->update_component_parameters( aT, aP );
            return mComponentV( aIndex );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::hdep( const uint & aIndex, const real & aT, const real & aP )
        {
            this->update_component_parameters( aT, aP );

            Matrix< real > & tData = mDepartureCoefficients( aIndex );

            return mComponentHDEP( aIndex ) -
                ( (   ( tData( 0, mComponentCol )   * aT
                      + tData( 1, mComponentCol ) ) * aT
                      + tData( 2, mComponentCol ) ) * aT
                      + tData( 3, mComponentCol ) );

        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::cpdep( const uint & aIndex, const real & aT, const real & aP )
        {
            this->update_component_parameters( aT, aP );

            Matrix< real > & tData = mDepartureCoefficients( aIndex );

            return mComponentCPDEP( aIndex ) -
                      ( ( ( 3.0 * tData( 0, mComponentCol )   * aT
                          + 2.0 * tData( 1, mComponentCol ) ) * aT
                          +       tData( 2, mComponentCol ) ) );

        }

//----------------------------------------------------------------------------

        real
        EoS_Cubic::chi( const real & aV )
        {
            if ( mDepartureV != aV )
            {
                mDepartureValue = std::log( ( aV - mB( 1 ) )/ ( aV -  mB( 2 ) )) / mB( 3 );
                mDepartureV = aV;
            }
            return mDepartureValue;
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dchidT( const real & aT, const real & aP, const real & aV  )
        {
            return -this->alpha( aT, aP ) * aV / ( ( aV - mB( 1 ) ) * ( aV - mB( 2 ) ) );
        }

//------------------------------------------------------------------------------

        real
        EoS_Cubic::dchidp( const real & aT, const real & aP, const real & aV  )
        {

            real tdVdP =  - mParent.kappa( aT, aP ) * aV ;

            return -tdVdP / ( ( aV - mB( 1 ) ) * ( aV - mB( 2 ) ) );
        }

//------------------------------------------------------------------------------
    }
}