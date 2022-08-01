//
// Created by Christian Messe on 30.04.20.
//

#include "cl_TableGas.hpp"
#include "assert.hpp"
#include "cl_GM_EoS_TableGas.hpp"

namespace belfem
{
//---------------------------------------------------------------------------

    TableGas::TableGas( const string aPath)
    {
        this->create_table( aPath );
        this->init_parameters();
        this->create_eos();
    }

//---------------------------------------------------------------------------

    TableGas::~TableGas()
    {
        // destroyed by parent
        // delete mEoS ;

        // destroy the lookop table
        delete mTable ;
    }

//---------------------------------------------------------------------------

    void
    TableGas::create_table( const string aPath )
    {
        if( aPath == "" )
        {
            // create default table
            mTable = new bspline::LookupTable(
                    gastables::data_path() + "/hotair.hdf5" );
        }
        else
        {
            // create table from file
            mTable = new bspline::LookupTable( aPath );
        }
    }

//---------------------------------------------------------------------------

    void
    TableGas::init_parameters()
    {
        // a tablegas does not have any components
        mNumberOfComponents = 0;

        mIndexM = mTable->field_index( "M" );
        mIndexH = mTable->field_index( "h" );
        mIndexS = mTable->field_index( "s" );
        mIndexMu = mTable->field_index( "mu" );
        mIndexLambda = mTable->field_index( "lambda" );
        mIndexHd = mTable->field_index( "hd" );
    }

//---------------------------------------------------------------------------

    void
    TableGas::create_eos()
    {
        // TableGas is always an ideal gas
        mEoS = new gasmodels::EoS_TableGas( *this );
        mGasModel = GasModel::IDGAS ;
    }

//---------------------------------------------------------------------------

    void
    TableGas::remix( const Vector<real> & aMolarFractions,
                     bool aRemixHeat,
                     bool aRemixTransport )
    {
        BELFEM_ERROR( false, "Calling 'remix()' is forbidden for a TableGas" );
    }

//------------------------------------------------------------------------------

    void
    TableGas::remix_mass( const Vector<real> & aMassFractions,
                bool aRemixHeat,
                bool aRemixTransport )
    {
        BELFEM_ERROR( false, "Calling 'remix_mass()' is forbidden for a TableGas" );
    }

//------------------------------------------------------------------------------

    void
    TableGas::reset_mixture()
    {
        BELFEM_ERROR( false, "Calling 'reset_mixture()' is forbidden for a TableGas" );
    }

//------------------------------------------------------------------------------

    real
    TableGas::pi( const real aT, const real aP )
    {
        return mEoS->pi( aT, aP );
    }

//------------------------------------------------------------------------------
// MIXTURE
//------------------------------------------------------------------------------

    const real &
    TableGas::M( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_M ) )
        {
            real tM = mTable->compute_value( mIndexM, aT, this->pi( aT, aP ) ) ;

            // set the M-value
            mStatevals.set( BELFEM_STATEVAL_M, tM );

            // set the R-value
            mStatevals.set( BELFEM_STATEVAL_R, constant::Rm / tM );

        }

        return mStatevals.get( BELFEM_STATEVAL_M );
    }

//------------------------------------------------------------------------------

    const real &
    TableGas::R( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_R ) )
        {
            real tM = mTable->compute_value( mIndexM, aT, this->pi( aT, aP ) ) ;

            // set the M-value
            mStatevals.set( BELFEM_STATEVAL_M, tM );

            // set the R-value
            mStatevals.set( BELFEM_STATEVAL_R, constant::Rm / tM );

        }

        return mStatevals.get( BELFEM_STATEVAL_R );
    }

//------------------------------------------------------------------------------
// Caloric Properties
//------------------------------------------------------------------------------

    real
    TableGas::cp( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_CP ) )
        {

            Vector< real > tDH(
                    mTable->compute_derivative( mIndexH,
                            aT, this->pi( aT, aP ) ) );

            // set the value
            mStatevals.set( BELFEM_STATEVAL_CP, tDH( 0 ) ) ;
            mStatevals.set( BELFEM_STATEVAL_DHDP, tDH( 1 ) * mEoS->dpidp( aT, aP ) );
        }

        return mStatevals.get( BELFEM_STATEVAL_CP );
    }

//------------------------------------------------------------------------------

    real
    TableGas::cv( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_CV ) )
        {
            // set the value
            mStatevals.set( BELFEM_STATEVAL_CV,
                    this->cp( aT, aP ) - this->R( aT, aP ) );

        }

        return mStatevals.get( BELFEM_STATEVAL_CV );
    }

//------------------------------------------------------------------------------

    real
    TableGas::gamma( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_GAMMA ) )
        {
            real tCp = this->cp( aT, aP );

            // set the value
            mStatevals.set( BELFEM_STATEVAL_GAMMA,
                            tCp / ( tCp - this->R( aT, aP ) ) );

        }

        return mStatevals.get( BELFEM_STATEVAL_GAMMA );
    }

//------------------------------------------------------------------------------

    real
    TableGas::c( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_C ) )
        {
            real tCp = this->cp( aT, aP );
            real tR  = this->R( aT, aP );

            // set the value
            mStatevals.set( BELFEM_STATEVAL_C,
                            std::sqrt( tCp/( tCp - tR ) * tR * aT ) );

        }

        return mStatevals.get( BELFEM_STATEVAL_C );
    }

//------------------------------------------------------------------------------

    real
    TableGas::h( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_H ) )
        {
            // set the value
            mStatevals.set( BELFEM_STATEVAL_H,
                            mTable->compute_value( mIndexH, aT,
                                    this->pi( aT, aP ) ) ) ;

        }

        return mStatevals.get( BELFEM_STATEVAL_H );
    }

//------------------------------------------------------------------------------

    real
    TableGas::s( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_S ) )
        {

            // set the value
            mStatevals.set( BELFEM_STATEVAL_S,
                            mTable->compute_value( mIndexS, aT,
                                                   this->pi( aT, aP ) ) ) ;

        }

        return mStatevals.get( BELFEM_STATEVAL_S );
    }

//------------------------------------------------------------------------------

    real
    TableGas::dsdT( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_DSDT ) )
        {
            // compute derivative
            Vector< real > tDS(
                    mTable->compute_derivative( mIndexS, aT,
                                                this->pi( aT, aP ) ) );

            // scale derivative to p
            tDS( 1 ) *= mdpscale / aP ;

            // set the value
            mStatevals.set( BELFEM_STATEVAL_DSDT, tDS( 0 ) ) ;
            mStatevals.set( BELFEM_STATEVAL_DSDP, tDS( 1 ) ) ;

        }

        return mStatevals.get( BELFEM_STATEVAL_DSDT );
    }

//------------------------------------------------------------------------------

    real
    TableGas::dsdp( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_DSDP ) )
        {
            // compute derivative
            Vector< real > tDS(
                    mTable->compute_derivative( mIndexS, aT,
                                                this->pi( aT, aP ) ) );

            // scale derivative to p
            tDS( 1 ) *= mdpscale / aP ;

            // set the value
            mStatevals.set( BELFEM_STATEVAL_DSDT, tDS( 0 ) ) ;
            mStatevals.set( BELFEM_STATEVAL_DSDP, tDS( 1 ) ) ;

        }

        return mStatevals.get( BELFEM_STATEVAL_DSDP );
    }

//------------------------------------------------------------------------------

    real
    TableGas::dcpdT( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_DCPDT ) )
        {
            // compute derivative
            Vector< real > tDS(
                    mTable->compute_second_derivative( mIndexH, aT,
                                                this->pi( aT, aP ) ) );

            // set the value
            mStatevals.set( BELFEM_STATEVAL_DCPDT, tDS( 0 ) ) ;

        }

        return mStatevals.get( BELFEM_STATEVAL_DCPDT );
    }

//------------------------------------------------------------------------------
// Transport Properties
//------------------------------------------------------------------------------

    real
    TableGas::mu( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_MU ) )
        {
            // set the value
            mStatevals.set( BELFEM_STATEVAL_MU,
                            mTable->compute_value( mIndexMu, aT,
                                                   this->pi( aT, aP ) ) );

        }

        return mStatevals.get( BELFEM_STATEVAL_MU );
    }

//------------------------------------------------------------------------------

    real
    TableGas::lambda( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_LAMBDA ) )
        {
            // set the value
            mStatevals.set( BELFEM_STATEVAL_LAMBDA,
                            mTable->compute_value( mIndexLambda, aT,
                                                   this->pi( aT, aP ) ) );

        }

        return mStatevals.get( BELFEM_STATEVAL_LAMBDA );
    }

//------------------------------------------------------------------------------

    real
    TableGas::Pr( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_PR ) )
        {
            // set the value
            mStatevals.set( BELFEM_STATEVAL_PR,
                this->mu( aT, aP ) * this->cp( aT, aP ) / this->lambda( aT, aP ) );
        }

        return mStatevals.get( BELFEM_STATEVAL_PR );
    }

//------------------------------------------------------------------------------

    // enthalpy derivative to pressure ( needed for total temperature )
    real
    TableGas::dhdp( const real aT, const real aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( ! mStatevals.test( BELFEM_STATEVAL_CP ) )
        {

            Vector< real > tDH(
                    mTable->compute_derivative( mIndexH,
                                                aT, this->pi( aT, aP ) ) );

            // set the value
            mStatevals.set( BELFEM_STATEVAL_CP, tDH( 0 ) ) ;
            mStatevals.set( BELFEM_STATEVAL_DHDP, tDH( 1 ) * mEoS->dpidp( aT, aP ) );
        }

        return mStatevals.get( BELFEM_STATEVAL_DHDP );
    }

//------------------------------------------------------------------------------

    // dissociation enthalpy
    real
    TableGas::hd( const real aT, const real aP )
    {
        return mTable->compute_value( mIndexHd, aT,
                                      this->pi( aT, aP ) );
    }

//------------------------------------------------------------------------------
}
