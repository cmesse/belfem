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

#include <cmath>

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
        RefGas::add_component( const string & aLabel, const real aValue )
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

                BELFEM_ERROR( tRefPoly != nullptr,
                        "No heat polynomial of %s contains the reference temperature %f K",
                        mLabel.c_str(), ( double ) gTref );

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

            // A junction counts as already smooth below this, relative to cp.
            //
            // It isolates one case: where a fitted low temperature interval hands
            // over to the tabulated data, which it is constrained to meet in value
            // and slope, and does to 1e-10 or better. Everything else stays above
            // the threshold and still goes through the search below — a junction
            // between two tabulated intervals, whose slope kink runs from 4e-7 up
            // to 1e-2, and also a junction between two fitted intervals, which is
            // constrained the same way but only agrees to about 4e-6 because
            // evaluating either seven term polynomial there costs several digits
            // to cancellation.
            const real tSmooth = 1.0e-7;

            // loop over all original polynomials
            for ( uint k=1; k<aNumberOfOriginalPolynomials; ++k )
            {
                real tTmid = mHeatPolys( k )->T_min();

                // A glue polynomial repairs that kink, so it is only wanted where
                // there is one. Across a junction that is already smooth, the
                // search below looks for a sign change in what is round off, never
                // finds one, and runs into the error at the bottom of the loop.
                real tCp = mHeatPolys( k )->Cp( tTmid );

                if ( std::abs( mHeatPolys( k-1 )->Cp( tTmid ) - tCp )
                             < tSmooth * std::abs( tCp )
                  && std::abs( mHeatPolys( k-1 )->dCpdT( tTmid )
                             - mHeatPolys( k )->dCpdT( tTmid ) )
                             < tSmooth * std::abs( tCp ) / tTmid )
                {
                    continue;
                }

                real tDeltaT =( tTmid < 2000 ) ? 5.0 : 10.0;

                real tDeltaTmax = ( tTmid < 500 ) ? 50 : 250;

                Vector<real> tRHS( 5 );
                Vector<real> tCoefficients( 5 );

                while( true )
                {
                    real tTmin = tTmid - tDeltaT;
                    real tTmax = tTmid + tDeltaT;

                    // fill RHS
                    tRHS( 0 ) = mHeatPolys( k - 1 )->Cp( tTmin );
                    tRHS( 1 ) = mHeatPolys( k - 1 )->dCpdT( tTmin );
                    tRHS( 2 ) = 0.5 * ( mHeatPolys( k - 1 )->Cp( tTmid ) + mHeatPolys( k )->Cp( tTmid ));
                    tRHS( 3 ) = mHeatPolys( k )->Cp( tTmax );
                    tRHS( 4 ) = mHeatPolys( k )->dCpdT( tTmax );

                    // scale conditions
                    tRHS /= constant::Rm;

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

            /* A junction counts as already smooth below this, relative to the
             * property itself. Same role as the threshold in
             * create_glue_polys_heat, and the same case: a fitted low
             * temperature interval handing over to the tabulated data, which it
             * is constrained to meet in value and slope. Most species of the
             * shipped table carry such an interval, so without this test the
             * search below hunts a sign change in round off, never finds one -
             * the curvature condition excludes zero by construction - and runs
             * into the error at the bottom of the loop. */
            const real tSmooth = 1.0e-7;

            // loop over all original polynomials
            for ( uint k=1; k<aNumberOfOriginalPolynomials; ++k )
            {
                TransportPoly * tPoly0 = aPolys( k-1 );
                TransportPoly * tPoly1 = aPolys( k );

                real tTmid =  tPoly1->T_min();

                // a glue polynomial repairs a kink, so it is only wanted where
                // there is one
                real tValue = tPoly1->eval( tTmid );

                if ( std::abs( tPoly0->eval( tTmid ) - tValue )
                             < tSmooth * std::abs( tValue )
                  && std::abs( tPoly0->deval( tTmid ) - tPoly1->deval( tTmid ) )
                             < tSmooth * std::abs( tValue ) / tTmid )
                {
                    continue;
                }

                real tTmin;
                real tTmax;

                real tDeltaTmax = ( tTmid < 500 ) ? 50 : 250;

                /* The glue must not wiggle: phase 0 accepts only a window over
                 * which its curvature keeps one sign, which is the historical
                 * condition. Where the two branches genuinely curve in opposite
                 * directions - the low temperature conductivity fits carry such
                 * inflections - no such window exists at ANY width, since the
                 * glue matches both end curvatures. Phase 1 therefore accepts
                 * the geometric minimum of exactly one curvature sign change,
                 * which still rejects wiggly glues. Phase 1 only runs when
                 * phase 0 is exhausted, so every junction the old condition
                 * accepted is glued identically. */
                bool tSuccess = false ;

                for( uint tPhase = 0; tPhase < 2 && ! tSuccess; ++tPhase )
                {
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
                        for ( uint i = 0; i < tN; ++i )
                        {
                            tValues( i ) = tPoly->ddrawpoly( tSteps( i ) );
                        }

                        bool tAccept ;

                        if( tPhase == 0 )
                        {
                            // curvature keeps one sign over the whole window
                            tAccept = ( min( tValues ) < 0.0 && max( tValues ) < 0.0 )
                                   || ( min( tValues ) > 0.0 && max( tValues ) > 0.0 )
                                   || tTmid >= 4000.0 ;
                        }
                        else
                        {
                            // exactly one sign change: the inflection the
                            // branches demand, and nothing more
                            uint tNumSignChanges = 0 ;

                            for ( uint i = 1; i < tN; ++i )
                            {
                                if ( tValues( i-1 ) * tValues( i ) < 0.0 )
                                {
                                    ++tNumSignChanges ;
                                }
                            }

                            tAccept = tNumSignChanges == 1 ;
                        }

                        if ( tAccept )
                        {
                            aPolys.push( tPoly );
                            tSuccess = true ;
                            break;
                        }
                        else
                        {
                            delete tPoly;
                            tDeltaT += 1.0;
                        }

                        if ( tDeltaT > tDeltaTmax )
                        {
                            // this phase is exhausted
                            break;
                        }
                    }
                }

                BELFEM_ERROR( tSuccess,
                             "something went terribly wrong in RefGas::create_glue_polys_transport()" );

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

        const HeatPoly *
        RefGas::find_heat_poly( const real T ) const
        {
            for( HeatPoly * aPoly : mHeatPolys )
            {
                if ( aPoly->T_min() <= T && T <= aPoly->T_max() )
                {
                    return aPoly;
                }
            }

            BELFEM_ERROR( false,
            "Temperature T=%f out of bounds for gas %s.",
                           ( float ) T,
                           mLabel.c_str() );

            return nullptr;
        }

//-------------------------------------------------------------------------------

        const TransportPoly *
        RefGas::find_viscosity_poly( const real T ) const
        {

            for( TransportPoly * aPoly : mViscosityPolys )
            {
                if ( aPoly->T_min() <= T && T <= aPoly->T_max() )
                {
                    return aPoly;
                }
            }

            BELFEM_ERROR( false,
                         "Temperature T=%f out of bounds for gas %s.",
                         ( float ) T,
                         mLabel.c_str() );

            return nullptr;
        }

//-------------------------------------------------------------------------------

        const TransportPoly *
        RefGas::find_conductivity_poly( const real T ) const
        {
            for( TransportPoly * aPoly : mConductivityPolys )
            {
                if ( aPoly->T_min() <= T && T <= aPoly->T_max() )
                {
                    return aPoly;
                }
            }

            BELFEM_ERROR( false,
                         "Temperature T=%f out of bounds for gas %s.",
                         ( float ) T,
                         mLabel.c_str() );

            return nullptr;
        }

//------------------------------------------------------------------------------

        void
        RefGas::set_mode( const RefGasMode & aMode )
        {
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
            else
            {
                BELFEM_ERROR( false, "Invalid RefGasMode for gas %s", mLabel.c_str() );
            }
        }

//------------------------------------------------------------------------------

        real
        RefGas::zero( const real T ) const
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
            return  mData.Href() / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::Cp( const real T ) const
        {
            return ( this->*mFunctionCp ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::H( const real T ) const
        {
            return ( this->*mFunctionH ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::S( const real T ) const
        {
            return ( this->*mFunctionS ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dSdT( const real T ) const
        {
            return ( this->*mFunctiondSdT ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dCpdT( const real T ) const
        {
            return ( this->*mFunctiondCpdT ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::d2CpdT2( const real T ) const
        {
            return ( this->*mFunctiond2CpdT2 ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::mu( const real T ) const
        {
            return ( this->*mFunctionMu ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dmudT( const real T ) const
        {
            return ( this->*mFunctiondMudT ) ( T );
        }


//------------------------------------------------------------------------------

        real
        RefGas::d2mudT2( const real T ) const
        {
            return ( this->*mFunctiond2MudT2 ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::lambda( const real T ) const
        {
            return ( this->*mFunctionLambda ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::dlambdadT( const real T ) const
        {
            return ( this->*mFunctiondLambdadT ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::d2lambdadT2( const real T ) const
        {
            return ( this->*mFunctiond2LambdadT2 ) ( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_Cp( const real T ) const
        {
           return this->find_heat_poly( T )->Cp( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_H( const real T ) const
        {
            return this->find_heat_poly( T )->H( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_S( const real T ) const
        {
            return this->find_heat_poly( T )->S( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dSdT( const real T ) const
        {
            return this->find_heat_poly( T )->dSdT( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dCpdT( const real T ) const
        {
            return this->find_heat_poly( T )->dCpdT( T );
        }
//------------------------------------------------------------------------------

        real
        RefGas::poly_d2CpdT2( const real T ) const
        {
            return this->find_heat_poly( T )->d2CpdT2( T );
        }

//-----------------------------------------------------------------------------

        real
        RefGas::poly_Mu( const real T ) const
        {
            return this->find_viscosity_poly( T )->eval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dMudT( const real T ) const
        {
            return this->find_viscosity_poly( T )->deval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_d2MudT2( const real T ) const
        {
            return this->find_viscosity_poly( T )->ddeval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_Lambda( const real T ) const
        {
            return this->find_conductivity_poly( T )->eval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_dLambdadT( const real T ) const
        {
            return this->find_conductivity_poly( T )->deval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::poly_d2LambdadT2( const real T ) const
        {
            return this->find_conductivity_poly( T )->ddeval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_Cp( const real T ) const
        {
            return mHeatSpline.deval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dCpdT( const real T ) const
        {
            return mHeatSpline.ddeval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_H( const real T ) const
        {
            return mHeatSpline.eval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_S( const real T ) const
        {
            return mHeatSpline.entropy( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dSdT( const real T ) const
        {
            return mHeatSpline.dentropy( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_Mu( const real T ) const
        {
            return mViscositySpline.eval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dMudT( const real T ) const
        {
            return mViscositySpline.deval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_d2MudT2( const real T ) const
        {
            return mViscositySpline.ddeval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_Lambda( const real T ) const
        {
            return mConductivitySpline.eval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_dLambdadT( const real T ) const
        {
            return mConductivitySpline.deval( T );
        }

//------------------------------------------------------------------------------

        real
        RefGas::spline_d2LambdadT2( const real T ) const
        {
            return mConductivitySpline.ddeval( T );
        }

//------------------------------------------------------------------------------

        void
        RefGas::create_splines(
                const Vector< real > & T,
                      SpMatrix       & aHelpMatrix )
        {
            this->fix_switches();

            uint tNumberOfSamples = T.length();

            Vector< real > tValues( tNumberOfSamples );

            this->set_mode( RefGasMode::POLY );

            // Heats
            if ( this->has_thermo() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = this->H( T( k ) );
                }

                mHeatSpline.initialize( T, tValues, aHelpMatrix,
                spline::SplineBC::NoCurvature,spline::SplineBC::NoCurvature,0.0,0.0,
                gTref, this->S( gTref ) );
            }

            // Viscosity
            if ( this->has_viscosity() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = this->mu( T( k ) );
                }

                mViscositySpline.initialize( T, tValues, aHelpMatrix );
            }
            else if ( this->has_thermo() && mData.has_crit() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = idgas_mu( this, T( k ) );
                }

                mViscositySpline.initialize( T, tValues, aHelpMatrix );
                mHaveViscosity = true;
            }
            // Conductiviity
            if( this->has_conductivity() )
            {
                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = this->lambda( T ( k ) );
                }
                mConductivitySpline.initialize(T, tValues, aHelpMatrix );
            }
            else if ( this->has_thermo() && mData.has_crit() )
            {
                // idgas_lambda needs the viscosity, and this object is still
                // evaluating from polynomials. For a species whose viscosity was
                // synthesized just above there is no viscosity polynomial, only
                // the empty placeholder, so mu() has to read the spline that the
                // block above has already filled. Both branches leave it filled.
                this->set_mode( RefGasMode::SPLINE );

                for ( uint k = 0; k < tNumberOfSamples; ++k )
                {
                    tValues( k ) = idgas_lambda( this, T( k ) );
                }

                mConductivitySpline.initialize( T, tValues, aHelpMatrix );
                mHaveConductivity = true;
            }

            this->set_mode(  RefGasMode::SPLINE );
        }

//------------------------------------------------------------------------------

        real
        RefGas::cp( const real T ) const
        {
            return this->Cp( T ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::h( const real T ) const
        {
            return this->H( T ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::s( const real T ) const
        {
            return this->S( T ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::dcpdT( const real T ) const
        {
            return this->dCpdT( T ) / mData.M();
        }

//------------------------------------------------------------------------------

        real
        RefGas::d2cpdT2( const real T ) const
        {
            return this->d2CpdT2( T ) / mData.M();
        }

//------------------------------------------------------------------------------

        void
        RefGas::fix_switches()
        {
            if( this->has_viscosity() )
            {

                mHaveViscosity = mViscosityPolys( 0 )->kind() != TransportPolyKind::EMPTY ;
            }
            if( this->has_conductivity() )
            {
                mHaveConductivity = mConductivityPolys( 0 )->kind() != TransportPolyKind::EMPTY ;
            }
        }

//------------------------------------------------------------------------------
    }
}