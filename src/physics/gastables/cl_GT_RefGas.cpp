//
// Created by Christian Messe on 2019-08-19.
//
#include "constants.hpp"
#include "fn_linspace.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"

#include "GT_globals.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_HeatPolyCustom.hpp"
#include "cl_GT_HeatPolyGlue.hpp"
#include "cl_GT_HeatPolyEmpty.hpp"
#include "cl_GT_TransportPolyEmpty.hpp"
#include "cl_GT_TransportPolyCustom.hpp"
#include "cl_GT_TransportPolyGlue.hpp"
#include "cl_GT_ComparisonObjects.hpp"
#include "fn_GT_create_glue_poly.hpp"
#include "fn_create_fifth_order_beam_poly.hpp"

#include "fn_GT_idgas_mu.hpp"
#include "fn_GT_idgas_lambda.hpp"
#include "fn_GT_is_noble.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        RefGas::RefGas( const string & aLabel ) :
            mAmNoble( gastables::is_noble( aLabel )),
            mLabel( mData.label() ),
            mHeatSpline( gNumberOfSplinePoints, 0.0, gTmax ),
            mViscositySpline( gNumberOfSplinePoints, 0.0, gTmax ),
            mConductivitySpline( gNumberOfSplinePoints, 0.0, gTmax )
        {

            mData.set_label( aLabel );

            this->set_mode( RefGasMode::POLY );
        }

//------------------------------------------------------------------------------

        RefGas::~RefGas()
        {
            this->delete_heat_polys();
            this->delete_transport_polys();
        }

//------------------------------------------------------------------------------

        void
        RefGas::add_component( const string & aLabel, const real & aValue )
        {
            mData.mElements.push( aLabel );
            mData.mComposition[ aLabel ] = aValue;
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_molar_mass( const real aMolarMass )
        {
            mData.set_molar_mass( aMolarMass );
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_reference_formation_enthalpy( const real aDeltaHf )
        {
            mData.set_formation_enthalpy( aDeltaHf );
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_reference_enthalpy( const real aHref )
        {
            mData.set_reference_enthalpy( aHref );
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_liquid_flag()
        {
            mLiquidFlag = true;
        }

//------------------------------------------------------------------------------

        void
        RefGas::unset_liquid_flag()
        {
            mLiquidFlag = false;
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_cryo_thermo_flag()
        {
            mHaveCryoThermo = true;
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_cryo_transport_flag()
        {
            mHaveCryoTransport = true;
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_component_flag()
        {
            mHaveComponents = true;
        }

//------------------------------------------------------------------------------

        void
        RefGas::add_heat_poly( HeatPoly * aHeatPoly )
        {
            mHeatPolys.push( aHeatPoly );
            mHaveThermo = true;
        }

//------------------------------------------------------------------------------

        void
        RefGas::add_transport_poly( TransportPoly * aTransportPoly )
        {
            switch( aTransportPoly->type() )
            {
                case( TransportPolyType::VISCOSITY ):
                {
                    mHaveViscosity = true;
                    mViscosityPolys.push( aTransportPoly );
                    break;
                }
                case( TransportPolyType::CONDUCTIVITY ):
                {
                    mHaveConductivity = true;
                    mConductivityPolys.push( aTransportPoly );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false,
                            "something went terribly wrong while trying to add polynomial to %s",
                            mLabel.c_str() );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        RefGas::delete_heat_polys()
        {
            // delete heat polynomials
            for( auto tPoly : mHeatPolys )
            {
                delete tPoly;
            }

            // clear the container
            mHeatPolys.clear();

            // reset the flag
            mHaveThermo = false;
        }

//------------------------------------------------------------------------------

        void
        RefGas::delete_transport_polys()
        {
            // delete viscosity polynomoals
            for( auto tPoly : mViscosityPolys )
            {
                delete tPoly;
            }

            // delete transport polynomoals
            for( auto tPoly : mConductivityPolys )
            {
                delete tPoly;
            }

            // clear the containers
            mViscosityPolys.clear();
            mConductivityPolys.clear();

            // reset the flag
            mHaveConductivity = false;
            mHaveViscosity = false;
        }


//------------------------------------------------------------------------------

        void
        RefGas::finalize_thermo()
        {
            if( mHaveThermo  )
            {

                // change temperature range of first poly if cryodata exist
                if( this->has_cryo_thermo() )
                {
                    mHeatPolys( 1 )->set_T_min( mHeatPolys( 0 )->T_max() );
                }


                // epsilon environment for temperature
                const real tEps = 1e-3;

                // find the poly that contains T_ref
                HeatPoly * tRefPoly = nullptr;
                for ( HeatPoly * tPoly : mHeatPolys )
                {
                    if( tPoly->T_min()-tEps <= gTref && gTref <= tPoly->T_max()+tEps )
                    {
                        tRefPoly = tPoly;
                        break;
                    }
                }

                // enthalpy at Tref
                real tHref = mData.Href()  + mData.Hf();

                // entropy at Tref
                real tSref = tRefPoly->S( gTref );

                // number of polynomialy that are read from the table
                uint tNumberOfOriginalPolynomials = mHeatPolys.size();

                // fix reference points ( first run )
                this->fix_reference_points( 1, tNumberOfOriginalPolynomials );

                // create connecting polynomials
                this->create_glue_polys_heat( tNumberOfOriginalPolynomials );

                // extrapolate cold polynomial
                this->create_cryo_poly_heat();

                // extrapolate hot polynomial
                this->create_hot_poly_heat( tNumberOfOriginalPolynomials );

                // sort the polynomials
                sort( mHeatPolys, fHeatPoly );

                // fix first poly manually
                mHeatPolys(0 )->set_enthalpy_constant( 0.0 );
                mHeatPolys(0 )->set_entropy_constant( 0.0 );

                // make polynomials contunuous
                this->fix_reference_points( 1, mHeatPolys.size() );

                // get offsets from reference poly
                real tDeltaH = tHref-tRefPoly->H( gTref );
                real tDeltaS = tSref-tRefPoly->S( gTref );

                // shift polynomials
                for ( HeatPoly * tPoly : mHeatPolys )
                {
                    tPoly->set_enthalpy_constant( tPoly->enthalpy_constant() + tDeltaH );
                    tPoly->set_entropy_constant( tPoly->entropy_constant() + tDeltaS );
                }

                // remember entropy
                mData.set_reference_enthalpy( tHref );
                mData.set_reference_entropy( tSref );

            }
            else
            {
                mHeatPolys.push( new HeatPolyEmpty() );
            }
        }

//------------------------------------------------------------------------------

        void
        RefGas::finalize_transport()
        {
            if ( this->has_viscosity() )
            {
                if( this->has_cryo_transport() )
                {
                    mViscosityPolys( 1 )->set_T_min( mViscosityPolys( 0 )->T_max());
                }

                uint tN = mViscosityPolys.size();
                this->create_glue_polys_transport( mViscosityPolys, tN );
                this->create_hot_poly_transport( mViscosityPolys, tN );
                this->create_cryo_poly_transport( mViscosityPolys );

                // sort the polynomials
                sort( mViscosityPolys, fTransportPoly );
            }
            else
            {
                mViscosityPolys.push( new TransportPolyEmpty( TransportPolyType::VISCOSITY ) );
            }

            if( this->has_conductivity() )
            {
                uint tN = mConductivityPolys.size();

                if( this->has_cryo_transport() )
                {
                    mConductivityPolys( 1 )->set_T_min( mConductivityPolys( 0 )->T_max());
                }

                this->create_glue_polys_transport( mConductivityPolys, tN );
                this->create_hot_poly_transport( mConductivityPolys, tN );
                this->create_cryo_poly_transport( mConductivityPolys );

                // sort the polynomials
                sort( mConductivityPolys, fTransportPoly );
            }
            else
            {
                mConductivityPolys.push( new TransportPolyEmpty(TransportPolyType::CONDUCTIVITY ) );
            }
        }

//------------------------------------------------------------------------------

        void
        RefGas::finalize()
        {
            if( ! mFinalizedFlag )
            {
                this->finalize_thermo();
                this->finalize_transport();

                mFinalizedFlag = true;
            }
        }

//------------------------------------------------------------------------------

        void
        RefGas::create_glue_polys_heat( const uint & aNumberOfOriginalPolynomials )
        {

            const uint tN = 100;
            Vector< real > tValues( tN );
            Vector< real > tSteps( tN );

            // loop over all original polynomials
            for ( uint k=1; k<aNumberOfOriginalPolynomials; ++k )
            {
                real tTmid = mHeatPolys( k )->T_min();
                real tDeltaT =( tTmid < 2000 ) ? 5.0 : 10.0;

                real tDeltaTmax = ( tTmid < 500 ) ? 50 : 250;

                while( true )
                {
                    real tTmin = tTmid - tDeltaT;
                    real tTmax = tTmid + tDeltaT;

                    // get boundary conditions on the left side
                    Vector<real> tRHS( 5 );

                    // fill RHS
                    tRHS( 0 ) = mHeatPolys( k - 1 )->Cp( tTmin );
                    tRHS( 1 ) = mHeatPolys( k - 1 )->dCpdT( tTmin );
                    tRHS( 2 ) = 0.5 * ( mHeatPolys( k - 1 )->Cp( tTmid ) + mHeatPolys( k )->Cp( tTmid ));
                    tRHS( 3 ) = mHeatPolys( k )->Cp( tTmax );
                    tRHS( 4 ) = mHeatPolys( k )->dCpdT( tTmax );

                    // scale conditions
                    tRHS /= constant::Rm;

                    // Coefficients for the new cas
                    Vector<real> tCoefficients( 5 );

                    // create the polynomials
                    create_glue_poly( tTmid, tDeltaT, tRHS, tCoefficients );

                    // create a new polynomial
                    HeatPolyGlue * tPoly = new HeatPolyGlue( tTmin, tTmax, 0.0, 0.0, tCoefficients );

                    // create steps
                    linspace( tTmin, tTmax, tN, tSteps );

                    // evaluate steps
                    for( uint i=0; i<tN; ++i )
                    {
                        tValues( i ) = tPoly->d2CpdT2( tSteps( i ) );
                    }

                    if ( ( min( tValues ) < 0 && max( tValues ) < 0 ) || ( min( tValues ) > 0 && max( tValues ) > 0 ) || tTmid >= 500 )
                    {
                        // adjust boundaries
                        mHeatPolys( k - 1 )->set_T_max( tTmin );
                        mHeatPolys( k )->set_T_min( tTmax );
                        mHeatPolys.push( tPoly );
                        break;
                    }
                    else
                    {
                        delete tPoly;
                        tDeltaT += 1.0;
                    }

                    BELFEM_ERROR( tDeltaT <= tDeltaTmax,
                        "something went terribly wrong in RefGas::create_glue_polys_heat()" );

                }
            }
        }

//------------------------------------------------------------------------------

        void
        RefGas::create_cryo_poly_transport( Cell< TransportPoly * > & aPolys )
        {
            real tT = aPolys( 0 )->T_min();


            // pick the initial poly
            TransportPoly * tPoly = aPolys( 0 );

            real tF     = tPoly->rawpoly( tT );
            real tdFdT  = tPoly->drawpoly( tT );


            Vector< real > tExponents = { 1.0, 0.0 };
            Vector< real > tCoefficients( 2 );
            tCoefficients( 0 ) = tdFdT;
            tCoefficients( 1 ) =  tF - tdFdT * tT;

            aPolys.push( new TransportPolyCustom(
                    tPoly->type(),
                    0.0,
                    tT,
                    tCoefficients,
                    tExponents ) );
        }


//-------------------------------------------------------------------------------

        void
        RefGas::create_glue_polys_transport(
                Cell< TransportPoly * > & aPolys,
                const uint & aNumberOfOriginalPolynomials )
        {
            const uint tN = 100;
            Vector< real > tValues( tN );
            Vector< real > tSteps( tN );

            // vector with supporting values
            Vector<real> tF( 4 );

            Vector< real > tCoefficients( 6 );
            Vector< real > tExponents = { 5, 4, 3, 2, 1, 0 };

            // loop over all original polynomials
            for ( uint k=1; k<aNumberOfOriginalPolynomials; ++k )
            {
                TransportPoly * tPoly0 = aPolys( k-1 );
                TransportPoly * tPoly1 = aPolys( k );

                real tTmid =  tPoly1->T_min();

                real tTmin;
                real tTmax;

                real tDeltaTmax = ( tTmid < 500 ) ? 50 : 250;
                real tDeltaT =  tTmid < 4000 ? 5 : 1500;


                while( true )
                {
                    tTmin = tTmid - tDeltaT;
                    tTmax = tTmid + tDeltaT;

                    // calculate coefficients
                    create_fifth_order_beam_poly(
                            tTmin,
                            tPoly0->rawpoly( tTmin ),
                            tPoly0->drawpoly( tTmin ),
                            tPoly0->ddrawpoly( tTmin ),
                            tTmax,
                            tPoly1->rawpoly( tTmax ),
                            tPoly1->drawpoly( tTmax ),
                            tPoly1->ddrawpoly( tTmax ),
                            tCoefficients );

                    // create steps
                    linspace( tTmin, tTmax, tN, tSteps );

                    // create a new polynomial
                    TransportPoly * tPoly = new TransportPolyCustom(
                            tPoly0->type(),
                            tTmin,
                            tTmax,
                            tCoefficients,
                            tExponents );

                    // test polynomial
                    for ( uint k = 0; k < tN; ++k )
                    {
                        tValues( k ) = tPoly->ddrawpoly( tSteps( k ) );
                    }

                    // delete me
                    if (tDeltaT == tDeltaTmax )
                    {
                        for ( uint k = 0; k < tN; ++k )
                        {
                            tValues( k ) = std::exp( tPoly->rawpoly( tSteps( k ) ) );
                        }
                        tValues.print("Values");
                        tSteps.print("Steps");
                    }

                    // check curvature condition
                    if (   ( min( tValues ) < 0.0 && max( tValues ) < 0.0 )
                        || ( min( tValues ) > 0.0 && max( tValues ) > 0.0 ) || tTmid > 4000.0 )
                    {
                        aPolys.push( tPoly );
                        break;
                    }
                    else
                    {
                        delete tPoly;
                        tDeltaT += 1.0;
                    }





                    BELFEM_ERROR( tDeltaT <= tDeltaTmax,
                                 "something went terribly wrong in RefGas::create_glue_polys_transport()" );
                }

                // change boundaries of polynomuals
                tPoly0->set_T_max( tTmin );
                tPoly1->set_T_min( tTmax );
            }
        }

//-------------------------------------------------------------------------------

        void
        RefGas::create_hot_poly_transport( Cell< TransportPoly * > & aPolys,
                                   const uint & aNumberOfOriginalPolynomials )
        {
            // last polynomial
            TransportPoly * tLastPoly = aPolys( aNumberOfOriginalPolynomials- 1 );

            real tT =  tLastPoly->T_max();

            if( tT < gTmax )
            {
                Vector< real > tExponents = { 1.0, 0.0 };
                Vector< real > tCoefficients( 2 );

                // create linear extrapolation
                tCoefficients( 0 ) = tLastPoly->drawpoly( tT );
                tCoefficients( 1 ) = tLastPoly->rawpoly( tT ) -  tCoefficients( 0 ) * tT;

                // create the new polynomial
                aPolys.push( new TransportPolyCustom( tLastPoly->type(), tT, gTmax, tCoefficients, tExponents ) );
            }
        }

//-------------------------------------------------------------------------------
        void
        RefGas::create_cryo_poly_heat()
        {
            // beware that these polynomials are for numeric stability only. They are not realistic.

            Vector< real > tExponents = { 3.0, 2.0, 0.0 };
            Vector< real > tCoefficients( 3 );

            real tT        = mHeatPolys( 0 )->T_min();
            real tCp       = mHeatPolys( 0 )->Cp( tT );
            real tH        = mHeatPolys( 0 )->H( tT ) + mData.Href() - mHeatPolys( 0 )->H( gTref );
            real tdCpdT    = mHeatPolys( 0 )->dCpdT( tT );

            // this polynomial is tangential at T=0 and hits H, cp and dcpdT at T_min
            tCoefficients( 0 ) = ( 4.0 *( ( tdCpdT*tT - 3.0*tCp) *tT + 3.0*tH ) )/( 3.0*std::pow( tT, 4 ) * constant::Rm );
            tCoefficients( 1 ) = -( 3.0 *( ( tdCpdT*tT - 4.0*tCp) *tT + 4.0*tH ) )/( 2.0*std::pow( tT, 3 ) * constant::Rm );
            tCoefficients( 2 ) = ( ( tdCpdT*tT - 6.0*tCp ) *tT + 12.0*tH )/( 6.0*tT*constant::Rm );

            // now we make sure that cp never gets negative
            if( tCoefficients( 2 ) < 0 )
            {
                // use other parameters instead to enforce cp = 0 at 0 K
                tExponents( 0 ) = 4.0;
                tExponents( 1 ) = 3.0;
                tExponents( 2 ) = 2.0;

                tCoefficients( 0 ) =  (2.5*( ( tdCpdT*tT - 6.0 * tCp) *tT + 12.0 * tH))/( std::pow( tT, 5 ) * constant::Rm );
                tCoefficients( 1 ) = -(4.0*( ( tdCpdT*tT - 7.0 * tCp) *tT + 15.0 * tH))/( std::pow( tT, 4 ) * constant::Rm );
                tCoefficients( 2 ) =  (1.5*( ( tdCpdT*tT - 8.0 * tCp) *tT + 20.0 * tH))/( std::pow( tT, 3 ) * constant::Rm );


            }

            // create new polynomial
            HeatPolyCustom * tPoly = new HeatPolyCustom( 0.0, tT, 0.0, 0.0, tCoefficients, tExponents );

            // find minimal Cp
            real tT1 = 0.0;
            real tT2 = tT;
            real tTold = 1000;
            real tF1 = tPoly->dCpdT( tT1 );
            real tF;
            uint tCount = 0;

            // find minmal Cp in order to proof that Cp > 0 all the time
            while ( std::abs( tTold - tT ) > 1.0e-3 )
            {
                // shift T
                tTold = tT;

                // new beta
                tT = 0.5 * ( tT1 + tT2 );

                // call beta function
                tF = tPoly->dCpdT( tT );

                // test result
                if( tF1 * tF > 0 )
                {
                    tT1 = tT;
                    tF1 = tF;
                }
                else
                {
                    tT2 = tT;
                }

                // increment counter
                BELFEM_ERROR( tCount++ < 1000, "too many iterations" );
            }

            BELFEM_ERROR( tPoly->Cp( tT ) > 0, "Negative specific heat for gas %s detected. This is a serious bug.", mLabel.c_str() );
            // add polynomial to list
            mHeatPolys.push( tPoly );
        }

//-------------------------------------------------------------------------------

        void
        RefGas::create_hot_poly_heat( const uint & aNumberOfOriginalPolynomials )
        {
            // last polynomial
            HeatPoly * tLastPoly = mHeatPolys( aNumberOfOriginalPolynomials- 1 );

            real tT =  tLastPoly->T_max();

            if( tT < gTmax )
            {
                // Cps are linearly extrapolated
                Vector< real > tExponents = { 1.0, 0.0 };
                Vector< real > tCoefficients( 2 );
                tCoefficients( 0 ) = tLastPoly->dCpdT( tT ) / constant::Rm;
                tCoefficients( 1 ) = tLastPoly->Cp( tT ) / constant::Rm -  tCoefficients( 0 ) * tT;

                mHeatPolys.push( new HeatPolyCustom( tT, gTmax, 0.0, 0.0, tCoefficients, tExponents ) );
            }
        }

//-------------------------------------------------------------------------------

        void
        RefGas::fix_reference_points( const uint & aStart, const uint & aEnd )
        {
            // reset constant ( first run )
            for( uint k=aStart; k<aEnd; ++k )
            {
                mHeatPolys( k )->set_enthalpy_constant( 0.0 );
                mHeatPolys( k )->set_entropy_constant( 0.0 );

                // get initial temperature
                real tT = mHeatPolys( k )->T_min();

                // fix enthalpy
                real tValue = mHeatPolys( k-1 )->H( tT ) - mHeatPolys( k )->H( tT );
                mHeatPolys( k )->set_enthalpy_constant( tValue );

                // fix entropy
                tValue = mHeatPolys( k-1 )->S( tT ) - mHeatPolys( k )->S( tT );
                mHeatPolys( k )->set_entropy_constant( tValue );
            }
        }

//-------------------------------------------------------------------------------

        HeatPoly *
        RefGas::find_heat_poly( const real & aT )
        {
            for( HeatPoly * aPoly : mHeatPolys )
            {
                if ( aPoly->T_min() <= aT && aT <= aPoly->T_max() )
                {
                    return aPoly;
                }
            }

            BELFEM_ERROR( false,
            "Temperature T=%f out of bounds for gas %s.",
                           ( float ) aT,
                           mLabel.c_str() );

            return nullptr;
        }

//-------------------------------------------------------------------------------

        TransportPoly *
        RefGas::find_viscosity_poly( const real & aT )
        {

            for( TransportPoly * aPoly : mViscosityPolys )
            {
                if ( aPoly->T_min() <= aT && aT <= aPoly->T_max() )
                {
                    return aPoly;
                }
            }

            BELFEM_ERROR( false,
                         "Temperature T=%f out of bounds for gas %s.",
                         ( float ) aT,
                         mLabel.c_str() );

            return nullptr;
        }

//-------------------------------------------------------------------------------

        TransportPoly *
        RefGas::find_conductivity_poly( const real & aT )
        {
            for( TransportPoly * aPoly : mConductivityPolys )
            {
                if ( aPoly->T_min() <= aT && aT <= aPoly->T_max() )
                {
                    return aPoly;
                }
            }

            BELFEM_ERROR( false,
                         "Temperature T=%f out of bounds for gas %s.",
                         ( float ) aT,
                         mLabel.c_str() );

            return nullptr;
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_mode( const RefGasMode & aMode )
        {
            mMode = aMode;
            if( aMode == RefGasMode::POLY )
            {
                mFunctionCp           = & RefGas::poly_Cp;
                mFunctiondCpdT        = & RefGas::poly_dCpdT;
                mFunctiond2CpdT2      = & RefGas::poly_d2CpdT2;
                mFunctionH            = & RefGas::poly_H;
                mFunctionS            = & RefGas::poly_S;
                mFunctiondSdT         = & RefGas::poly_dSdT;
                mFunctionMu           = & RefGas::poly_Mu;
                mFunctiondMudT        = & RefGas::poly_dMudT;
                mFunctiond2MudT2      = & RefGas::poly_d2MudT2;
                mFunctionLambda       = & RefGas::poly_Lambda;
                mFunctiondLambdadT    = & RefGas::poly_dLambdadT;
                mFunctiond2LambdadT2  = & RefGas::poly_d2LambdadT2;
            }
            else if ( aMode == RefGasMode::SPLINE )
            {
                if ( this->has_thermo() )
                {
                    mFunctionCp           = & RefGas::spline_Cp;
                    mFunctiondCpdT        = & RefGas::spline_dCpdT;
                    mFunctiond2CpdT2      = & RefGas::poly_d2CpdT2; // <-- this is intended
                    mFunctionH            = & RefGas::spline_H;
                    mFunctionS            = & RefGas::spline_S;
                    mFunctiondSdT         = & RefGas::spline_dSdT;
                }
                else
                {
                    mFunctionCp           = & RefGas::zero;
                    mFunctiondCpdT        = & RefGas::zero;
                    mFunctiond2CpdT2      = & RefGas::zero;
                    mFunctionH            = & RefGas::zero;
                    mFunctionS            = & RefGas::zero;
                    mFunctiondSdT         = & RefGas::zero;
                }

                if ( this->has_viscosity() )
                {
                    mFunctionMu           = & RefGas::spline_Mu;
                    mFunctiondMudT        = & RefGas::spline_dMudT;
                    mFunctiond2MudT2      = & RefGas::spline_d2MudT2;
                }
                else
                {
                    mFunctionMu           = & RefGas::zero;
                    mFunctiondMudT        = & RefGas::zero;
                    mFunctiond2MudT2      = & RefGas::zero;
                }

                if ( this->has_conductivity() )
                {
                    mFunctionLambda       = & RefGas::spline_Lambda;
                    mFunctiondLambdadT    = & RefGas::spline_dLambdadT;
                    mFunctiond2LambdadT2  = & RefGas::spline_d2LambdadT2;
                }
                else
                {
                    mFunctionLambda       = & RefGas::zero;
                    mFunctiondLambdadT    = & RefGas::zero;
                    mFunctiond2LambdadT2  = & RefGas::zero;
                }
            }
        }

//------------------------------------------------------------------------------

        real
        RefGas::zero( const real & aT )
        {
            return 0.0;
        }
//------------------------------------------------------------------------------

        real
        RefGas::H_ref() const
        {
            return  mData.Href();
        }

//------------------------------------------------------------------------------

        real
        RefGas::h_ref() const
        {
            return  mData.Href() * mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::Cp( const real & aT )
        {
            return ( this->*mFunctionCp ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::H( const real & aT )
        {
            return ( this->*mFunctionH ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::S( const real & aT )
        {
            return ( this->*mFunctionS ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dSdT( const real & aT )
        {
            return ( this->*mFunctiondSdT ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dCpdT( const real & aT )
        {
            return ( this->*mFunctiondCpdT ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::d2CpdT2( const real & aT )
        {
            return ( this->*mFunctiond2CpdT2 ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::mu( const real & aT )
        {
            return ( this->*mFunctionMu ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dmudT( const real & aT )
        {
            return ( this->*mFunctiondMudT ) ( aT );
        }


//------------------------------------------------------------------------------

        real
        RefGas::d2mudT2( const real & aT )
        {
            return ( this->*mFunctiond2MudT2 ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::lambda( const real & aT )
        {
            return ( this->*mFunctionLambda ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dlambdadT( const real & aT )
        {
            return ( this->*mFunctiondLambdadT ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::d2lambdadT2( const real & aT )
        {
            return ( this->*mFunctiond2LambdadT2 ) ( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_Cp( const real & aT )
        {
           return this->find_heat_poly( aT )->Cp( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_H( const real & aT )
        {
            return this->find_heat_poly( aT )->H( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_S( const real & aT )
        {
            return this->find_heat_poly( aT )->S( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dSdT( const real & aT )
        {
            return this->find_heat_poly( aT )->dSdT( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dCpdT( const real & aT )
        {
            return this->find_heat_poly( aT )->dCpdT( aT );
        }
//------------------------------------------------------------------------------

        real
        RefGas::poly_d2CpdT2( const real & aT )
        {
            return this->find_heat_poly( aT )->d2CpdT2( aT );
        }

//-----------------------------------------------------------------------------

        real
        RefGas::poly_Mu( const real & aT )
        {
            return this->find_viscosity_poly( aT )->eval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dMudT( const real & aT )
        {
            return this->find_viscosity_poly( aT )->deval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_d2MudT2( const real & aT )
        {
            return this->find_viscosity_poly( aT )->ddeval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_Lambda( const real & aT )
        {
            return this->find_conductivity_poly( aT )->eval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dLambdadT( const real & aT )
        {
            return this->find_conductivity_poly( aT )->deval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_d2LambdadT2( const real & aT )
        {
            return this->find_conductivity_poly( aT )->ddeval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_Cp( const real & aT )
        {
            return mHeatSpline.deval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dCpdT( const real & aT )
        {
            return mHeatSpline.ddeval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_H( const real & aT )
        {
            return mHeatSpline.eval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_S( const real & aT )
        {
            return mHeatSpline.entropy( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dSdT( const real & aT )
        {
            return mHeatSpline.dentropy( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_Mu( const real & aT )
        {
            return mViscositySpline.eval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dMudT( const real & aT )
        {
            return mViscositySpline.deval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_d2MudT2( const real & aT )
        {
            return mViscositySpline.ddeval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_Lambda( const real & aT )
        {
            return mConductivitySpline.eval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dLambdadT( const real & aT )
        {
            return mConductivitySpline.deval( aT );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_d2LambdadT2( const real & aT )
        {
            return mConductivitySpline.ddeval( aT );
        }

//------------------------------------------------------------------------------

        void
        RefGas::create_splines(
                const Vector< real > & aT,
                      SpMatrix       & aHelpMatrix )
        {
            uint tNumberOfSamples = aT.length();

            Vector< real > tValues( tNumberOfSamples );

            this->set_mode( RefGasMode::POLY  );

            // Heats
            if ( this->has_thermo() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = this->H( aT( k ) );
                }

                mHeatSpline.initialize( aT, tValues, aHelpMatrix, gTref, this->S( gTref ) );
            }

            // Viscosity
            if ( this->has_viscosity() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = this->mu( aT( k ) );
                }

                mViscositySpline.initialize( aT, tValues, aHelpMatrix );
            }
            else if ( this->has_thermo() && mData.has_crit() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = idgas_mu( this, aT( k ) );
                }

                mViscositySpline.initialize( aT, tValues, aHelpMatrix );
                mHaveViscosity = true;
            }
            // Conductiviity
            if( this->has_conductivity() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = this->lambda( aT ( k ) );
                }
                mConductivitySpline.initialize(aT, tValues, aHelpMatrix );
            }
            else if ( this->has_thermo() && mData.has_crit() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = idgas_lambda( this, aT( k ) );
                }

                mConductivitySpline.initialize( aT, tValues, aHelpMatrix );
                mHaveConductivity = true;
            }

            this->set_mode(  RefGasMode::SPLINE );
        }

//------------------------------------------------------------------------------

        real
        RefGas::cp( const real & aT )
        {
            return this->Cp( aT ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::h( const real & aT )
        {
            return this->H( aT ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::s( const real & aT )
        {
            return this->S( aT ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::dcpdT( const real & aT )
        {
            return this->dCpdT( aT ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::d2cpdT2( const real & aT )
        {
            return this->d2CpdT2( aT ) / mData.M();
        }

//------------------------------------------------------------------------------
    }
}