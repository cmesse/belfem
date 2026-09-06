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

#include <cctype>

#include "cl_NgspiceCircuitFactory.hpp"
#include "cl_ElectricalCircuit.hpp"
#include "cl_SourceFunction.hpp"
#include "fn_spice_number.hpp"
#include "assert.hpp"
#include "stringtools.hpp"
#include "constants.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------
//  construction / destruction
//-----------------------------------------------------------------------------

        NgspiceCircuitFactory::NgspiceCircuitFactory( const NetlistParser & aNetlist ) :
            mSource( aNetlist.source() )
        {
            this->build( aNetlist );
        }

//-----------------------------------------------------------------------------

        NgspiceCircuitFactory::NgspiceCircuitFactory( const string & aPath ) :
            mSource( aPath )
        {
            NetlistParser tNetlist( aPath );
            this->build( tNetlist );
        }

//-----------------------------------------------------------------------------

        NgspiceCircuitFactory::~NgspiceCircuitFactory() = default;
        // produce-and-hand-over: an unclaimed circuit dies with the
        // unique_ptr -- including when a refusal throws mid-build

//-----------------------------------------------------------------------------

        ElectricalCircuit *
        NgspiceCircuitFactory::circuit()
        {
            BELFEM_ERROR( mCircuit != nullptr,
                          "circuit() called twice on the netlist factory for %s "
                          "-- the first caller owns the circuit",
                          mSource.c_str() );
            return mCircuit.release();
        }

//-----------------------------------------------------------------------------
//  the build sequence
//-----------------------------------------------------------------------------

        void
        NgspiceCircuitFactory::build( const NetlistParser & aNetlist )
        {
            this->build_node_map( aNetlist );

            Map< string, uint > tOrders;
            this->collect_orders( aNetlist, tOrders );

            mCircuit.reset( new ElectricalCircuit( mNumberOfNodes ) );

            // duplicate instance labels would make the current-output
            // lookup ambiguous
            Map< string, index_t > tSeenLabels;
            auto tCheckLabel = [ this, &tSeenLabels ]( const string & aLabel,
                                                       const index_t aLine )
            {
                BELFEM_ERROR( !tSeenLabels.key_exists( aLabel ),
                              "%s:%lu: duplicate instance name '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aLine,
                              aLabel.c_str() );
                tSeenLabels[ aLabel ] = aLine;
            };

            // create components in netlist line order -- elements and
            // directives merged -- so the unknown-current layout is frozen
            // by the file ( O9 )
            const Cell< NetlistElement > & tElements = aNetlist.elements();
            const Cell< NetlistDirective > & tDirectives = aNetlist.directives();

            index_t tE = 0;
            index_t tD = 0;
            while ( tE < tElements.size() || tD < tDirectives.size() )
            {
                const bool tTakeElement =
                        tD >= tDirectives.size() ||
                        ( tE < tElements.size() &&
                          tElements( tE ).mLine < tDirectives( tD ).mLine );

                if ( tTakeElement )
                {
                    tCheckLabel( tElements( tE ).mName, tElements( tE ).mLine );
                    this->create_element( tElements( tE ), aNetlist, tOrders );
                    ++tE;
                }
                else
                {
                    const NetlistDirective & tDirective = tDirectives( tD );
                    if ( tDirective.mKind != "order" )
                    {
                        BELFEM_ERROR( tDirective.mArgs.size() > 0,
                                      "%s:%lu: directive '%s' needs a component name",
                                      mSource.c_str(),
                                      ( long unsigned int ) tDirective.mLine,
                                      tDirective.mKind.c_str() );
                        tCheckLabel( tDirective.mArgs( 0 ), tDirective.mLine );
                    }
                    this->create_from_directive( tDirective );
                    ++tD;
                }
            }
        }

//-----------------------------------------------------------------------------
//  node map ( O9 packing, O5 lookup )
//-----------------------------------------------------------------------------

        bool
        NgspiceCircuitFactory::is_ground( const string & aName )
        {
            return aName == "0" || aName == "gnd";
        }

//-----------------------------------------------------------------------------

        void
        NgspiceCircuitFactory::build_node_map( const NetlistParser & aNetlist )
        {
            const Cell< string > & tNames = aNetlist.node_names();

            bool tHaveGround = false;
            index_t tCount = 0;

            for ( index_t k = 0; k < tNames.size(); ++k )
            {
                if ( is_ground( tNames( k ) ) )
                {
                    tHaveGround = true;
                }
                else if ( !mNodeMap.key_exists( tNames( k ) ) )
                {
                    mNodeMap[ tNames( k ) ] = tCount++;
                }
            }

            // manual, "Ground node": every circuit has to have one
            BELFEM_ERROR( tHaveGround,
                          "netlist %s never references ground ( node 0 or gnd ) "
                          "-- the circuit would be singular",
                          mSource.c_str() );

            // BELFEM convention: ground is the LAST index
            mNodeMap[ "0" ] = tCount;
            mNumberOfNodes = tCount + 1;

            // O9: the frozen packing, visible in the log ( rank 0 only --
            // the factory is constructed on every rank )
            if ( comm_rank() != 0 )
            {
                return;
            }
            message( InfoLevel::Verbose,
                     "    netlist %s: %u nodes, ground = index %lu",
                     mSource.c_str(),
                     ( unsigned int ) mNumberOfNodes,
                     ( long unsigned int ) tCount );
            for ( index_t k = 0; k < tNames.size(); ++k )
            {
                if ( !is_ground( tNames( k ) ) )
                {
                    message( InfoLevel::Verbose,
                             "        node %-12s -> %lu",
                             tNames( k ).c_str(),
                             ( long unsigned int ) mNodeMap( tNames( k ) ) );
                }
            }
        }

//-----------------------------------------------------------------------------

        index_t
        NgspiceCircuitFactory::node_index( const string & aName ) const
        {
            const string tFolded = string_to_lower( aName );
            const string tCanonical = is_ground( tFolded ) ? string( "0" ) : tFolded;
            BELFEM_ERROR( mNodeMap.key_exists( tCanonical ),
                          "netlist %s has no node named '%s'",
                          mSource.c_str(),
                          aName.c_str() );
            return mNodeMap( tCanonical );
        }

//-----------------------------------------------------------------------------

        index_t
        NgspiceCircuitFactory::node_index_checked( const string & aName,
                                                   const string & aCard,
                                                   const index_t aLine ) const
        {
            const string tFolded = string_to_lower( aName );
            const string tCanonical = is_ground( tFolded ) ? string( "0" ) : tFolded;
            BELFEM_ERROR( mNodeMap.key_exists( tCanonical ),
                          "%s:%lu: unknown node '%s' in '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aName.c_str(),
                          aCard.c_str() );
            return mNodeMap( tCanonical );
        }

//-----------------------------------------------------------------------------
//  order directives ( O8 )
//-----------------------------------------------------------------------------

        void
        NgspiceCircuitFactory::collect_orders( const NetlistParser & aNetlist,
                                               Map< string, uint > & aOrders ) const
        {
            const Cell< NetlistDirective > & tDirectives = aNetlist.directives();

            for ( index_t k = 0; k < tDirectives.size(); ++k )
            {
                const NetlistDirective & tDirective = tDirectives( k );
                if ( tDirective.mKind != "order" )
                {
                    continue;
                }

                BELFEM_ERROR( tDirective.mArgs.size() == 2 &&
                              tDirective.mKwargs.size() == 0,
                              "%s:%lu: the order directive is "
                              "'* belfem: order <instance> <n>': '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) tDirective.mLine,
                              tDirective.mCard.c_str() );

                const string & tInstance = tDirective.mArgs( 0 );
                const string & tDigits   = tDirective.mArgs( 1 );

                // a single digit 1..6 ( the BDF kernel refuses orders
                // above 6 ) -- no multi-digit parse, no wraparound
                BELFEM_ERROR( tDigits.length() == 1 &&
                              tDigits[ 0 ] >= '1' && tDigits[ 0 ] <= '6',
                              "%s:%lu: BDF order must be a single digit in 1..6, got '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) tDirective.mLine,
                              tDigits.c_str() );
                const uint tOrder = static_cast< uint >( tDigits[ 0 ] - '0' );

                // the instance must be an L or C card of this netlist
                bool tFound = false;
                for ( index_t e = 0; e < aNetlist.elements().size(); ++e )
                {
                    const NetlistElement & tElement = aNetlist.elements()( e );
                    if ( tElement.mName == tInstance )
                    {
                        BELFEM_ERROR( tElement.mType == 'l' || tElement.mType == 'c',
                                      "%s:%lu: order directive targets '%s', which is "
                                      "not an inductor or capacitor",
                                      mSource.c_str(),
                                      ( long unsigned int ) tDirective.mLine,
                                      tInstance.c_str() );
                        tFound = true;
                        break;
                    }
                }
                BELFEM_ERROR( tFound,
                              "%s:%lu: order directive targets unknown instance '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) tDirective.mLine,
                              tInstance.c_str() );

                BELFEM_ERROR( !aOrders.key_exists( tInstance ),
                              "%s:%lu: duplicate order directive for '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) tDirective.mLine,
                              tInstance.c_str() );

                aOrders[ tInstance ] = tOrder;
            }
        }

//-----------------------------------------------------------------------------
//  element cards
//-----------------------------------------------------------------------------

        real
        NgspiceCircuitFactory::rcl_value( const NetlistElement & aElement,
                                          const char * aKeyword ) const
        {
            const bool tHavePositional = aElement.mValues.size() > 0;
            const bool tHaveKeyword    = aElement.mKwargs.key_exists( aKeyword );

            BELFEM_ERROR( aElement.mValues.size() <= 1,
                          "%s:%lu: too many value fields on '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aElement.mLine,
                          aElement.mCard.c_str() );

            BELFEM_ERROR( tHavePositional != tHaveKeyword,
                          "%s:%lu: expected exactly one value ( positional or %s= ) "
                          "on '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aElement.mLine,
                          aKeyword,
                          aElement.mCard.c_str() );

            // every other keyword is refused -- ic= and friends are
            // value-changing ( plan §12.2 )
            for ( auto tIterator = aElement.mKwargs.begin();
                  tIterator != aElement.mKwargs.end();
                  ++tIterator )
            {
                BELFEM_ERROR( tIterator->first == aKeyword,
                              "%s:%lu: parameter '%s=' is not supported in v1: '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aElement.mLine,
                              tIterator->first.c_str(),
                              aElement.mCard.c_str() );
            }

            const string & tToken = tHavePositional ?
                    aElement.mValues( 0 ) : aElement.mKwargs( aKeyword );

            const real tValue = spice_number_to_si( tToken );

            BELFEM_ERROR( tValue != 0.0,
                          "%s:%lu: zero value on '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aElement.mLine,
                          aElement.mCard.c_str() );

            return tValue;
        }

//-----------------------------------------------------------------------------

        SourceFunction *
        NgspiceCircuitFactory::source_function( const NetlistElement & aElement ) const
        {
            BELFEM_ERROR( aElement.mKwargs.size() == 0,
                          "%s:%lu: keyword parameters are not supported on "
                          "source cards in v1: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aElement.mLine,
                          aElement.mCard.c_str() );

            const Cell< string > & tValues = aElement.mValues;

            BELFEM_ERROR( tValues.size() > 0,
                          "%s:%lu: source card without a value: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aElement.mLine,
                          aElement.mCard.c_str() );

            // bare value or "DC <value>" -> constant
            if ( tValues.size() == 1 ||
                 ( tValues.size() == 2 && tValues( 0 ) == "dc" ) )
            {
                const string & tToken =
                        tValues.size() == 1 ? tValues( 0 ) : tValues( 1 );
                const real tAmplitude = spice_number_to_si( tToken );

                SourceFunction * tFunction = new SourceFunction();
                tFunction->set_constant( tAmplitude );
                return tFunction;
            }

            if ( tValues( 0 ) == "sin" )
            {
                // SIN(VO VA FREQ [TD [THETA [PHASE]]]), manual §"transient"
                BELFEM_ERROR( tValues.size() >= 4 && tValues.size() <= 7,
                              "%s:%lu: SIN takes 3 to 6 arguments: '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aElement.mLine,
                              aElement.mCard.c_str() );

                const real tOffset    = spice_number_to_si( tValues( 1 ) );
                const real tAmplitude = spice_number_to_si( tValues( 2 ) );
                const real tFrequency = spice_number_to_si( tValues( 3 ) );
                const real tDelay     = tValues.size() > 4 ?
                        spice_number_to_si( tValues( 4 ) ) : 0.0;
                const real tTheta     = tValues.size() > 5 ?
                        spice_number_to_si( tValues( 5 ) ) : 0.0;
                const real tPhaseDeg  = tValues.size() > 6 ?
                        spice_number_to_si( tValues( 6 ) ) : 0.0;

                // BELFEM sine has no offset, delay or damping -- refusing
                // beats silently changing every sample ( plan §12.2 )
                BELFEM_ERROR( tOffset == 0.0 && tDelay == 0.0 && tTheta == 0.0,
                              "%s:%lu: SIN offset, delay and damping must be zero "
                              "( BELFEM sine cannot represent them ): '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aElement.mLine,
                              aElement.mCard.c_str() );

                BELFEM_ERROR( tFrequency > 0.0,
                              "%s:%lu: SIN frequency must be positive: '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aElement.mLine,
                              aElement.mCard.c_str() );

                // SPICE phase is in degrees; set_periodic takes radians
                SourceFunction * tFunction = new SourceFunction();
                tFunction->set_periodic( SourceFunctionType::Sine,
                                         tAmplitude,
                                         1.0 / tFrequency,
                                         tPhaseDeg * constant::pi / 180.0 );
                return tFunction;
            }

            BELFEM_ERROR( false,
                          "%s:%lu: source function '%s' is not supported in v1 "
                          "( PULSE and PWL come with plan Phase 6 ): '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aElement.mLine,
                          tValues( 0 ).c_str(),
                          aElement.mCard.c_str() );
            return nullptr;  // unreachable
        }

//-----------------------------------------------------------------------------

        void
        NgspiceCircuitFactory::create_element( const NetlistElement & aElement,
                                               const NetlistParser & aNetlist,
                                               const Map< string, uint > & aOrders )
        {
            const index_t tPlus = this->node_index_checked(
                    aElement.mNodes( 0 ), aElement.mCard, aElement.mLine );
            const index_t tMinus = this->node_index_checked(
                    aElement.mNodes( 1 ), aElement.mCard, aElement.mLine );

            const uint tOrder = aOrders.key_exists( aElement.mName ) ?
                    aOrders( aElement.mName ) : 1;

            switch ( aElement.mType )
            {
                case 'r' :
                {
                    mCircuit->create_resistor(
                            this->rcl_value( aElement, "r" ),
                            tPlus, tMinus, aElement.mName );
                    break;
                }
                case 'c' :
                {
                    mCircuit->create_capacitor(
                            this->rcl_value( aElement, "c" ),
                            tOrder, tPlus, tMinus, aElement.mName );
                    break;
                }
                case 'l' :
                {
                    mCircuit->create_inductor(
                            this->rcl_value( aElement, "l" ),
                            tOrder, tPlus, tMinus, aElement.mName );
                    break;
                }
                case 'v' :
                {
                    mCircuit->create_voltage_source(
                            this->source_function( aElement ),
                            tPlus, tMinus, aElement.mName );
                    break;
                }
                case 'i' :
                {
                    mCircuit->create_current_source(
                            this->source_function( aElement ),
                            tPlus, tMinus, aElement.mName );
                    break;
                }
                case 'd' :
                {
                    BELFEM_ERROR( aElement.mValues.size() == 1 &&
                                  aElement.mKwargs.size() == 0,
                                  "%s:%lu: a diode card is 'Dxxx n+ n- <model>' "
                                  "( no area factor, no parameters ): '%s'",
                                  mSource.c_str(),
                                  ( long unsigned int ) aElement.mLine,
                                  aElement.mCard.c_str() );

                    const string & tModelName = aElement.mValues( 0 );
                    const NetlistModel * tModel = nullptr;
                    for ( index_t k = 0; k < aNetlist.models().size(); ++k )
                    {
                        if ( aNetlist.models()( k ).mName == tModelName )
                        {
                            tModel = &aNetlist.models()( k );
                            break;
                        }
                    }
                    BELFEM_ERROR( tModel != nullptr,
                                  "%s:%lu: no .model card named '%s' for '%s'",
                                  mSource.c_str(),
                                  ( long unsigned int ) aElement.mLine,
                                  tModelName.c_str(),
                                  aElement.mCard.c_str() );
                    BELFEM_ERROR( tModel->mType == "d",
                                  "%s:%lu: model '%s' has type '%s', a diode "
                                  "needs type d",
                                  mSource.c_str(),
                                  ( long unsigned int ) aElement.mLine,
                                  tModelName.c_str(),
                                  tModel->mType.c_str() );

                    real tIs = 1.0e-14;  // SPICE default saturation current
                    real tN  = 1.0;      // emission coefficient
                    for ( auto tIterator = tModel->mKwargs.begin();
                          tIterator != tModel->mKwargs.end();
                          ++tIterator )
                    {
                        if ( tIterator->first == "is" )
                        {
                            tIs = spice_number_to_si( tIterator->second );
                        }
                        else if ( tIterator->first == "n" )
                        {
                            tN = spice_number_to_si( tIterator->second );
                        }
                        else
                        {
                            BELFEM_ERROR( false,
                                          "%s:%lu: model parameter '%s=' is not "
                                          "supported in v1 ( only is= and n= ): '%s'",
                                          mSource.c_str(),
                                          ( long unsigned int ) tModel->mLine,
                                          tIterator->first.c_str(),
                                          tModel->mCard.c_str() );
                        }
                    }

                    // the Shockley law divides by Vt = n * kT/q and
                    // scales by Is -- nonpositive values are user error,
                    // not a device
                    BELFEM_ERROR( tIs > 0.0 && tN > 0.0,
                                  "%s:%lu: diode model '%s' needs is > 0 and n > 0",
                                  mSource.c_str(),
                                  ( long unsigned int ) aElement.mLine,
                                  tModelName.c_str() );

                    // Vt = n * kT/q at room temperature ( plan §3 )
                    mCircuit->create_diode( tIs, tN * 0.026,
                                            tPlus, tMinus, aElement.mName );
                    break;
                }
                default :
                {
                    // the parser only emits r/c/l/v/i/d
                    BELFEM_ERROR( false,
                                  "%s:%lu: unexpected element type '%c'",
                                  mSource.c_str(),
                                  ( long unsigned int ) aElement.mLine,
                                  aElement.mType );
                }
            }
        }

//-----------------------------------------------------------------------------
//  extension directives
//-----------------------------------------------------------------------------

        void
        NgspiceCircuitFactory::create_from_directive( const NetlistDirective & aDirective )
        {
            if ( aDirective.mKind == "order" )
            {
                // consumed by collect_orders()
                return;
            }

            // shared helpers for the component directives
            auto tRequire = [ this, &aDirective ]( const char * aKey ) -> const string &
            {
                BELFEM_ERROR( aDirective.mKwargs.key_exists( aKey ),
                              "%s:%lu: directive '%s' needs '%s=': '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aDirective.mLine,
                              aDirective.mKind.c_str(),
                              aKey,
                              aDirective.mCard.c_str() );
                return aDirective.mKwargs( aKey );
            };

            auto tCheckKeys = [ this, &aDirective ]( const Cell< string > & aAllowed )
            {
                for ( auto tIterator = aDirective.mKwargs.begin();
                      tIterator != aDirective.mKwargs.end();
                      ++tIterator )
                {
                    bool tKnown = false;
                    for ( index_t k = 0; k < aAllowed.size(); ++k )
                    {
                        if ( tIterator->first == aAllowed( k ) )
                        {
                            tKnown = true;
                            break;
                        }
                    }
                    BELFEM_ERROR( tKnown,
                                  "%s:%lu: unknown key '%s=' on directive '%s': '%s'",
                                  mSource.c_str(),
                                  ( long unsigned int ) aDirective.mLine,
                                  tIterator->first.c_str(),
                                  aDirective.mKind.c_str(),
                                  aDirective.mCard.c_str() );
                }
            };

            BELFEM_ERROR( aDirective.mArgs.size() == 1,
                          "%s:%lu: directive '%s' takes exactly one component "
                          "name: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aDirective.mLine,
                          aDirective.mKind.c_str(),
                          aDirective.mCard.c_str() );

            const string & tLabel = aDirective.mArgs( 0 );

            if ( aDirective.mKind == "superconductor" )
            {
                tCheckKeys( { "n+", "n-", "ic", "n", "ec", "length" } );
                const index_t tPlus = this->node_index_checked(
                        tRequire( "n+" ), aDirective.mCard, aDirective.mLine );
                const index_t tMinus = this->node_index_checked(
                        tRequire( "n-" ), aDirective.mCard, aDirective.mLine );

                const real tIc     = spice_number_to_si( tRequire( "ic" ) );
                const real tExpN   = spice_number_to_si( tRequire( "n" ) );
                const real tEc     = spice_number_to_si( tRequire( "ec" ) );
                const real tLength = spice_number_to_si( tRequire( "length" ) );

                // the E-J power law divides by all four
                BELFEM_ERROR( tIc > 0.0 && tExpN > 0.0 && tEc > 0.0 && tLength > 0.0,
                              "%s:%lu: superconductor '%s' needs Ic, n, Ec and "
                              "length > 0",
                              mSource.c_str(),
                              ( long unsigned int ) aDirective.mLine,
                              tLabel.c_str() );

                mCircuit->create_superconductor(
                        tIc, tExpN, tEc, tLength, tPlus, tMinus, tLabel );
            }
            else if ( aDirective.mKind == "switch" )
            {
                tCheckKeys( { "n+", "n-", "state", "t_switch" } );
                const index_t tPlus = this->node_index_checked(
                        tRequire( "n+" ), aDirective.mCard, aDirective.mLine );
                const index_t tMinus = this->node_index_checked(
                        tRequire( "n-" ), aDirective.mCard, aDirective.mLine );

                const string & tState = tRequire( "state" );
                BELFEM_ERROR( tState == "open" || tState == "closed",
                              "%s:%lu: switch state must be open or closed, "
                              "got '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aDirective.mLine,
                              tState.c_str() );

                const real tSwitchTime = spice_number_to_si( tRequire( "t_switch" ) );
                BELFEM_ERROR( tSwitchTime >= 0.0,
                              "%s:%lu: switch '%s' needs t_switch >= 0",
                              mSource.c_str(),
                              ( long unsigned int ) aDirective.mLine,
                              tLabel.c_str() );

                mCircuit->create_switch(
                        tState == "closed", tSwitchTime, tPlus, tMinus, tLabel );
            }
            else
            {
                BELFEM_ERROR( false,
                              "%s:%lu: unknown '* belfem:' directive '%s' "
                              "( v1 knows superconductor, switch, order ): '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aDirective.mLine,
                              aDirective.mKind.c_str(),
                              aDirective.mCard.c_str() );
            }
        }

//-----------------------------------------------------------------------------
    }
}
