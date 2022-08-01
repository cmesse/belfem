//
// Created by Christian Messe on 01.05.20.
//
#include "cl_GM_EoS_TableGas.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        EoS_TableGas::EoS_TableGas( TableGas & aParent ) :
                EoS( aParent ),
                mTable( aParent.lookup_table()),
                mIndexM( mTable.field_index( "M" )),
                mPimin( mTable.min( 1 ) + 1.0 ),
                mPimax( mTable.max( 1 ) - 1.0 ),
                mPmin( std::pow( 10,  mPimin*0.001 ) ),
                mPmax( std::pow( 10,  mPimax*0.001 ) )
        {
            // wirte a value into stateval p for iteration
            mStatevals.set( BELFEM_STATEVAL_T, gastables::gTref );
            mStatevals.set( BELFEM_STATEVAL_P, gastables::gPref );

        }

//----------------------------------------------------------------------------

        void
        EoS_TableGas::remix()
        {
            /* do nothing */
        }

//----------------------------------------------------------------------------
// Special
//----------------------------------------------------------------------------

        real
        EoS_TableGas::pi( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );

            if ( !mStatevals.test( BELFEM_STATEVAL_PI ))
            {
                mStatevals.set( BELFEM_STATEVAL_PI,
                        std::max( std::min( std::log10( aP ) * 1000.0, mPimax ), mPimin ) );
            }

            return mStatevals.get( BELFEM_STATEVAL_PI );
        }

//----------------------------------------------------------------------------

        real
        EoS_TableGas::dpidp( const real aT, const real aP )
        {
            return std::min( std::max( aP, mPmin ), mPmax ) * mdpscale ;
        }

//----------------------------------------------------------------------------

        const real &
        EoS_TableGas::M( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );

            if ( !mStatevals.test( BELFEM_STATEVAL_M ))
            {
                real tM = mTable.compute_value( mIndexM, aT, this->pi( aT, aP ));

                // set the M-value
                mStatevals.set( BELFEM_STATEVAL_M, tM );

                // set the R-value
                mStatevals.set( BELFEM_STATEVAL_R, constant::Rm / tM );

            }

            return mStatevals.get( BELFEM_STATEVAL_M );
        }

//----------------------------------------------------------------------------

        real
        EoS_TableGas::dMdT( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );

            if ( !mStatevals.test( BELFEM_STATEVAL_DMDT ))
            {
                // compute derivative
                Vector< real > tDM(
                        mTable.compute_derivative( mIndexM, aT,
                                                   this->pi( aT, aP )));

                // scale derivative to p
                tDM( 1 ) *= mdpscale / aP;

                // set the value
                mStatevals.set( BELFEM_STATEVAL_DMDT, tDM( 0 ));
                mStatevals.set( BELFEM_STATEVAL_DMDP, tDM( 1 ));

            }

            return mStatevals.get( BELFEM_STATEVAL_DMDT );
        }

//----------------------------------------------------------------------------

        real
        EoS_TableGas::dMdp( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );

            if ( !mStatevals.test( BELFEM_STATEVAL_DMDP ))
            {
                // compute derivative
                Vector< real > tDM(
                        mTable.compute_derivative( mIndexM, aT,
                                                   this->pi( aT, aP )));

                // scale derivative to p
                tDM( 1 ) *= mdpscale / aP;

                // set the value
                mStatevals.set( BELFEM_STATEVAL_DMDT, tDM( 0 ));
                mStatevals.set( BELFEM_STATEVAL_DMDP, tDM( 1 ));

            }

            return mStatevals.get( BELFEM_STATEVAL_DMDP );
        }

//----------------------------------------------------------------------------
// Thermodynamic State
//----------------------------------------------------------------------------

        real
        EoS_TableGas::p( const real aT, const real aV )
        {
            // relaxation factor
            const real tOmega0 = 0.9;
            real tOmega1;
            real tOmega ;

            // constant for function
            real tC = constant::Rm * aT / aV;

            real tF = BELFEM_REAL_MAX;
            real tdF;

            uint tCount = 0;

            // get a value to start the iteration
            real aP = mStatevals.get( BELFEM_STATEVAL_P );

            while ( std::abs( tF ) > 1e-7 )
            {
                // molar mass
                real tM = this->M( aT, aP );

                // function
                tF = tC / tM - aP;

                // derivative
                tdF = tC * this->dMdp( aT, aP )/( tM * tM ) - 1.0;

                tOmega1 = 0.9 * std::abs( ( aP - mPmin ) * tdF / tF  );

                if( tOmega1 < BELFEM_EPSILON )
                {
                    return mPmin ;
                }

                tOmega = tOmega1 < tOmega0 ? tOmega1 : tOmega0 ;

                // perform newton step
                aP -= tOmega * tF / tdF;

                BELFEM_ERROR( tCount++ < 1000,
                             "too many iterations in EoS_TableGas::p for T = %f, v= %f ",
                             ( float ) aT,
                             ( float ) aV );

            }

            // return pressure value
            return aP;
        }

//----------------------------------------------------------------------------

        real
        EoS_TableGas::v( const real aT, const real aP )
        {
            return mParent.R( aT, aP ) * aT / aP;
        }

//----------------------------------------------------------------------------

        real
        EoS_TableGas::T( const real aP, const real aV )
        {
            // relaxation factor
            const real tOmega = 0.9;

            // constant for function
            real tC = constant::Rm / aV;

            real tF = BELFEM_REAL_MAX;
            real tM;
            real tdF;

            uint tCount = 0;

            real aT = mStatevals.get( BELFEM_STATEVAL_T );

            while ( std::abs( tF ) < 1e-7 )
            {
                // compute M
                tM = this->M( aT, aP );


                // function
                tF = tC * aT / tM - aP;

                // derivative
                tdF = tC * ( tM - aT * this->dMdT( aT, aP )) / ( tM * tM );

                // perform newton step
                aT -= tOmega * tF / tdF;

                BELFEM_ERROR( tCount++ < 1000,
                             "too many iterations in EoS_TableGas::T for p = %f, v= %f ",
                             ( float ) aP,
                             ( float ) aV );

            }

            // return temperature value
            return aT;
        }

//------------------------------------------------------------------------------
// State Derivatives
//------------------------------------------------------------------------------

        real
        EoS_TableGas::dpdT( const real aT, const real aV )
        {
            real tP = this->p( aT, aV );

            return this->beta( aT, tP ) * tP;
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::d2pdT2( const real aT, const real aV )
        {
            real tP = this->p( aT, aV );
            real tM = this->M( aT, tP );
            real tdMdT = this->dMdT( aT, tP );

            Vector< real > tD2M( mTable.compute_second_derivative(
                    mIndexM, aT, this->pi( aT, tP )));

            real td2MdT2 = tD2M( 0 );

            return -constant::Rm / ( aV * tM * tM * tM ) *
                   ( 2.0 * tdMdT * ( tM - aT * tdMdT )
                     + tM * aT * td2MdT2 );
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::dpdv( const real aT, const real aV )
        {
            return -1.0 / ( this->kappa( aT, this->p( aT, aV )) * aV );
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::d2pdv2( const real aT, const real aV )
        {
            BELFEM_ERROR( false, "d2pdv2 is not implemented for TableGas" );
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::dvdT( const real aT, const real aV )
        {

            return this->alpha( aT, this->p( aT, aV )) * aV;
        }

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

        real
        EoS_TableGas::alpha( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );

            if ( !mStatevals.test( BELFEM_STATEVAL_ALPHA ))
            {


                real tM = this->M( aT, aP );

                real tV = this->v( aT, aP );

                real tdVdT = constant::Rm * ( tM - aT * this->dMdT( aT, aP )) /
                             ( tM * tM * aP );

                mStatevals.set( BELFEM_STATEVAL_ALPHA, tdVdT / tV );

            }

            return mStatevals.get( BELFEM_STATEVAL_ALPHA );
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::beta( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );
            if ( !mStatevals.test( BELFEM_STATEVAL_BETA ))
            {
                mStatevals.set( BELFEM_STATEVAL_BETA, this->alpha( aT, aP ) / ( this->kappa( aT, aP ) * aP ));
            }

            return mStatevals.get( BELFEM_STATEVAL_BETA );
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::kappa( const real aT, const real aP )
        {
            mStatevals.update_Tp( aT, aP );

            if ( !mStatevals.test( BELFEM_STATEVAL_KAPPA ))
            {
                real tM = this->M( aT, aP );

                real tV = this->v( aT, aP );

                real tdVdp = -constant::Rm * ( tM + aP * this->dMdp( aT, aP )) /
                             ( tM * tM * aP * aP );

                mStatevals.set( BELFEM_STATEVAL_KAPPA, -tdVdp / tV );
            }

            return mStatevals.get( BELFEM_STATEVAL_KAPPA );
        }

//------------------------------------------------------------------------------

        void
        EoS_TableGas::eval_critical_point( real & aT, real & aP, real & aV )
        {
            // an ideal gas has no critical point!
            aT = BELFEM_QUIET_NAN;
            aP = BELFEM_QUIET_NAN;
            aV = BELFEM_QUIET_NAN;
        }
//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

        real
        EoS_TableGas::v( const uint aIndex, const real aT, const real aP )
        {
            BELFEM_ERROR( false, "illegal function call in TableGas EoS");
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::hdep( const uint aIndex, const real aT, const real aP )
        {
            BELFEM_ERROR( false, "illegal function call in TableGas EoS");
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        real
        EoS_TableGas::cpdep( const uint aIndex, const real aT, const real aP )
        {
            BELFEM_ERROR( false, "illegal function call in TableGas EoS");
            return BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

    }

}