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

#include <variant>
#include "cl_MaxwellBoundaryConditionFactory.hpp"
#include "cl_Vector.hpp"
#include "constants.hpp"
#include "fn_check_unit.hpp"

namespace belfem
{
    namespace fem
    {
//-------------------------------------------------------------------------------

        MaxwellBoundaryConditionFactory::MaxwellBoundaryConditionFactory( const input::Section * aSection, const uint aNumberOfDimensions ) :
            mNumDimensions( aNumberOfDimensions )
        {
            BELFEM_ERROR( aSection->type()=="boundary conditions" || aSection->parent()->type()=="boundary conditions" ,"Maxwell Boundary Condition Factory should be defined with a Boundary Condition section" ) ;
            uint n = aSection->num_sections();

            // pre-count the global-name bases so repeated sections get a
            // section ordinal on EVERY occurrence ( "background_1",
            // "background_2" ), never only on the second: a bare base sitting
            // next to an ordinal-suffixed one reads as an inconsistent naming
            // scheme, and the ordinal is what tells two same-typed sections
            // apart. Terminal sections are not counted: they publish no block
            // global ( see has_block_global and the labelling loop below )
            for ( uint d=0; d<n; ++d )
            {
                const input::Section * tSection = aSection->section( d ) ;
                if ( has_block_global( boundary_condition_type( tSection->type() ) ) )
                {
                    string tBase = tSection->label().size() > 0 ?
                            tSection->label() : tSection->type() ;
                    mGlobalNameCount[ tBase ] =
                        mGlobalNameCount.key_exists( tBase ) ?
                            mGlobalNameCount( tBase ) + 1 : 1 ;
                }
            }

            uint tCount = 0;
            for ( uint d=0; d<n; ++d )
            {

                uint tNumCondition = 1; //number of boundary conditions associated with this section

                const input::Section * tSection = aSection->section( d ) ;

                Vector < id_t > tDomains ;
                string tUnits ;

                //Set the boundary condition type
                BoundaryConditionType tType = boundary_condition_type( tSection->type() ) ;

                // only needed if current BC is applied to thin shell

                //Read the domains on which the condition is applied and set the domains
                switch ( tType )
                {
                    case BoundaryConditionType::Bearing :
                    case BoundaryConditionType::Gauge :
                    {
                        BELFEM_ERROR(tSection->key_exists("nodes"), "Nodes undefined for gauge/bearing") ;
                        tSection->get_ids("nodes",tDomains) ;
                        mPhysicalBoundaryConditions.push(new PhysicalBoundaryCondition( tCount ));
                        mPhysicalBoundaryConditions( tCount )->set_type(tType) ;
                        mPhysicalBoundaryConditions(tCount++)->set_domains(tDomains) ;
                        break ;
                    }
                    case BoundaryConditionType::Current :
                    case BoundaryConditionType::Voltage :
                    {

                        string tKey ;
                        bool tIsThinShell = false ;
                        if (tSection->key_exists( "input terminal" ))
                        {
                            tKey = "input terminal" ;
                        }
                        else if (tSection->key_exists( "input terminals" ))
                        {
                            tKey = "input terminals" ;
                        }
                        else if (tSection->key_exists( "input curve" ) )
                        {
                            tKey = "input curve" ;
                            tIsThinShell = true ;
                        }
                        else if (tSection->key_exists( "input curves" ))
                        {
                            tKey = "input curves" ;
                            tIsThinShell = true ;
                        }
                        else
                        {
                            BELFEM_ERROR( false , "Input terminals or curves undefined for current boundary condition" ) ;
                        }


                        //Look for output terminals
                        bool tOutputExists = true ;
                        string tKey2 ;
                        if (tSection->key_exists( "output terminal" ))
                        {
                            tKey2 = "output terminal" ;
                        }
                        else if (tSection->key_exists( "output terminals" ))
                        {
                            tKey2 = "output terminals" ;
                        }
                        else if (tSection->key_exists( "output curve" ) )
                        {
                            tKey2 = "output curve" ;
                        }
                        else if (tSection->key_exists( "output curves" ))
                        {
                            tKey2 = "output curves" ;
                        }
                        else
                        {
                            tOutputExists = false ;
                        }

                        //Separate the groups defined with []. Each group becomes one boundary condition
                        Cell<Cell < id_t > >tDomainGroupsIn ;
                        tSection->get_id_groups(tKey,tDomainGroupsIn) ;
                        real tScale = 1.0 ;
                        Cell<Cell < id_t > >tDomainGroupsOut ;
                        if (tOutputExists)
                        {
                            tSection->get_id_groups(tKey2,tDomainGroupsOut) ;
                            BELFEM_ERROR(tDomainGroupsIn.size() == tDomainGroupsOut.size() ,
                                        "The current/voltage boundary condition should contain the same number of input and output terminals") ;
                        }
                        else //If the output terminals are not defined, the problem is 2-D, we read only the input terminals instead
                        {
                            tSection->get_id_groups(tKey,tDomainGroupsOut) ;
                            if (tType == BoundaryConditionType::Voltage)
                            {
                                BELFEM_ERROR( tSection->key_exists("length"), "Undefined length for voltage boundary condition in 2-D" ) ;
                                string tUnitDef = tSection->get_units( "length" );
                                value tValue = unit_to_si( tUnitDef );
                                BELFEM_ERROR(
                                            check_unit( tValue, "m" ),
                                            "required unit for the amplitude of the boundary condition is: m");

                                real tLength = tSection->get_value( "length","m" ).first ;

                                // same guard as the circuit terminal-pair path:
                                // zero gives an infinite scale, negative flips
                                // the sign of the imposed voltage
                                BELFEM_ERROR( tLength > 0.0,
                                              "length for the 2-D voltage boundary condition must be positive, but is %f m",
                                              tLength );

                                tScale *= 1.0/tLength ;

                            }
                        }


                        for(uint i = 0; i < tDomainGroupsIn.size(); ++i)
                        {
                            //Grouping the input and output terminals all into tDomains
                            tDomains.set_size(tDomainGroupsIn(i).size()+tDomainGroupsOut(i).size(),0) ;
                            for(uint j = 0 ; j < tDomainGroupsIn(i).size();++j)
                            {
                                tDomains(j) = tDomainGroupsIn(i)(j) ;
                            }
                            for(uint j = 0 ; j < tDomainGroupsOut(i).size();++j)
                            {
                                tDomains(tDomainGroupsIn(i).size()+j) = tDomainGroupsOut(i)(j) ;
                            }

                            mPhysicalBoundaryConditions.push(new PhysicalBoundaryCondition( tCount ));
                            mPhysicalBoundaryConditions( tCount )->set_type(tType) ;
                            mPhysicalBoundaryConditions( tCount )->set_scale( tScale ) ;
                            mPhysicalBoundaryConditions(tCount++)->set_domains(tDomains, tIsThinShell ) ;

                        }
                        tNumCondition = tDomainGroupsIn.size();
                        break ;
                    }
                    case BoundaryConditionType::Background :
                    case BoundaryConditionType::Neumann:
                    case BoundaryConditionType::Dirichlet:
                    {
                        BELFEM_ERROR( tSection->key_exists( "sideset" ) xor tSection->key_exists( "sidesets" ), "Sidesets undefined for boundary condition" ) ;
                        tSection->get_ids( tSection->key_exists( "sideset" ) ? "sideset" : "sidesets",tDomains ) ;
                        mPhysicalBoundaryConditions.push(new PhysicalBoundaryCondition( tCount ));
                        mPhysicalBoundaryConditions( tCount )->set_type(tType) ;
                        mPhysicalBoundaryConditions(tCount++)->set_domains(tDomains) ;
                        break ;
                    }
                    default:
                    {
                        // no boundary condition was pushed for this section,
                        // but the loop below still runs tNumCondition times and
                        // indexes mPhysicalBoundaryConditions( tCount - i ) --
                        // which would reach into a previously created condition,
                        // or underflow the unsigned index when this is the
                        // first one. Refuse the type instead of corrupting a
                        // neighbor ( this is how "background dirichlet" used
                        // to behave )
                        BELFEM_ERROR( false,
                                      "Unsupported boundary condition type: %s",
                                      tSection->type().c_str() );
                        break ;
                    }

                }

                // Name the mesh global this section publishes ( written by
                // impose_bc ): the section header label wins
                // ( "background : outer { }" — existing grammar ), else the
                // section type. If the base name repeats in the deck, EVERY
                // occurrence carries its section ordinal ( pre-counted
                // above ). A section publishes ONE global however many groups
                // it has — see the labelling loop below. Bearing imposes a
                // bare node constraint with no evaluated scalar — no global;
                // Gauge publishes its imposed potential. Current and voltage
                // sections publish no block global either: Controller::save_IV
                // writes one I/U pair per condition, and takes the section
                // label from the condition to name it
                string tGlobalBase ;
                if ( has_block_global( tType ) )
                {
                    tGlobalBase = tSection->label().size() > 0 ?
                            tSection->label() : tSection->type() ;

                    if ( mGlobalNameCount( tGlobalBase ) > 1 )
                    {
                        uint tOcc = mGlobalNameSeen.key_exists( tGlobalBase ) ?
                                mGlobalNameSeen( tGlobalBase ) + 1 : 1 ;
                        mGlobalNameSeen[ tGlobalBase ] = tOcc ;

                        tGlobalBase += "_" + std::to_string( tOcc );
                    }
                }

                for (uint i = 1; i <= tNumCondition; ++i)
                {
                    if ( tType == BoundaryConditionType::Current
                      || tType == BoundaryConditionType::Voltage )
                    {
                        // terminal conditions: EVERY member carries the section
                        // label ( empty if the section has none ), because
                        // Controller::save_IV names each condition's I/U pair
                        // after it and adds a running suffix where a label
                        // repeats — the groups of one section, or two sections
                        // sharing a label. No block global exists for these, so
                        // update_global() no-ops on the label
                        mPhysicalBoundaryConditions( tCount - i )->set_label( tSection->label() );
                    }
                    else if ( tGlobalBase.size() > 0 && i == tNumCondition )
                    {
                        // Only ONE member of a group carries the label, and so only
                        // one mesh global is published per section: every condition
                        // of a section shares that section's value function, so
                        // per-member ordinals would write copies of one number. i
                        // runs the conditions in reverse creation order, so the
                        // first-created member is the one reached at
                        // i == tNumCondition; the others keep an empty label, which
                        // update_global() already no-ops on. Caveat accepted
                        // deliberately: a "userdefined" section builds one
                        // SourceFunction per member from the same ( file, label ),
                        // so a source library that is impure in those arguments
                        // would show only the first member's value
                        mPhysicalBoundaryConditions( tCount - i )->set_label( tGlobalBase );
                    }

                    //Set the direction for the background field
                    if ( tType == BoundaryConditionType::Background )
                    {
                        BELFEM_ERROR( tSection->key_exists( "direction" ), "Undefined direction for the background field" ) ;
                        Vector < real > tDirection ;
                        tSection->get_reals( "Direction", tDirection ) ;
                        mPhysicalBoundaryConditions( tCount-i )->set_direction( tDirection( 0 ),tDirection( 1 ),tDirection( 2 ) ) ;
                    }

                    //Set the required units in order to define the function
                    switch ( tType )
                    {
                        case BoundaryConditionType::Current :
                        {
                            tUnits = "A" ;
                            break ;
                        }
                        case (BoundaryConditionType::Voltage) :
                        {
                            //tUnits = tDim == 3 ? "V" : "V/m" ;
                            tUnits = "V" ;
                            break ;
                        }
                        case BoundaryConditionType::Background :
                        case BoundaryConditionType::BackgroundDirichlet :
                        {
                            tUnits = "A/m" ;

                            // A background field is naturally written as a flux
                            // density, so a Tesla-family amplitude is accepted and
                            // read as B. The condition is imposed on a far-field
                            // boundary, where B = mu0*H, so H = B/mu0 = nu0*B. The
                            // factor rides on the condition's scale, which impose_bc
                            // applies once ( mValue = compute( t ) * mScale ).
                            //
                            // The token is matched by NAME rather than by dimension.
                            // This was originally FORCED: unit_to_si used to code the
                            // whole Tesla family with VOLT's exponents, so a dimension
                            // test accepted "1 V" as one tesla. That defect was fixed
                            // on 2026-08-30 ( stringtools.cpp, the "magnetic
                            // density" block ), and check_unit( x, "T" ) is now sound.
                            //
                            // The whitelist is KEPT because it is what was audited,
                            // not because it is still required. Switching to
                            // check_unit would be a small behaviour change, not a
                            // refactor: it would additionally accept compound
                            // spellings such as "V*s/m^2" and "kg/(s^2*A)". Make that
                            // change deliberately, with its own round -- do not slip
                            // it in as a cleanup.
                            //
                            // userdefined carries its dimension on a separate `units`
                            // key and never reads `amplitude`, so the key to inspect
                            // follows the waveform, not a precedence order.
                            // Gated on `type` as a whole: with no `type` key the
                            // waveform parser below reads NOTHING, so converting
                            // here would flag a condition that never acquires a
                            // source function -- and would raise the ferro warning
                            // for a section that is going to fail anyway
                            const string tAmplitudeKey = ! tSection->key_exists( "type" ) ?
                                    "" :
                                    tSection->get_string( "type" ) == "userdefined" ?
                                            "units" : "amplitude" ;

                            if ( tAmplitudeKey.size() > 0
                                 && tSection->key_exists( tAmplitudeKey ) )
                            {
                                // `units` holds the token itself, `amplitude` a
                                // number followed by one. Fold the micro sign the
                                // same way unit_to_si does ( stringtools.cpp:367 ),
                                // or "1 µT" would be converted by the unit layer
                                // and missed by this whitelist
                                const string tToken = search_and_replace(
                                        tAmplitudeKey == "units" ?
                                                tSection->get_string( tAmplitudeKey ) :
                                                tSection->get_units( tAmplitudeKey ),
                                        "µ", "mu" ) ;

                                if (    tToken == "T"   || tToken == "mT"
                                     || tToken == "muT" || tToken == "kT"
                                     || tToken == "MT"  || tToken == "G" )
                                {
                                    tUnits = "T" ;

                                    mPhysicalBoundaryConditions( tCount-i )
                                            ->scale() *= constant::nu0 ;

                                    mPhysicalBoundaryConditions( tCount-i )
                                            ->set_amplitude_is_flux_density( true ) ;
                                }
                            }
                            break ;
                        }
                        default :
                        {
                            break ;
                        }
                    }

                    //Define the function on the boundary
                    SourceFunction * tFunction = new SourceFunction() ;
                    switch ( tType )
                    {
                        case BoundaryConditionType::Bearing :
                        {
                            tSection->get_ids( "nodes",tDomains ) ;
                            break ;
                        }
                        case BoundaryConditionType::Current :
                        case (BoundaryConditionType::Voltage) :
                        case BoundaryConditionType::Gauge :
                        case BoundaryConditionType::Background :
                        case BoundaryConditionType::BackgroundDirichlet :
                        case BoundaryConditionType::Neumann :
                        case BoundaryConditionType::Dirichlet :
                        {
                            // 'expk' was parsed and stored for years without ever reaching the
                            // waveform. Refusing it beats dropping the read in silence, which is
                            // what this check exists to close. This sits BEFORE the type switch
                            // on purpose: the key was only ever read on the sigmoid branch, so a
                            // sigmoid-only check would leave it silently ignored on every other
                            // source shape, and on a section with no 'type' at all
                            BELFEM_ERROR( ! tSection->key_exists( "expk" ),
                                "'expk' is no longer accepted, and it never affected the waveform: "
                                "function_sigmoid takes its rate from 'fuzzyness' and 'period' alone. "
                                "DELETE the key to keep the results you get today. Do NOT translate it "
                                "-- if you meant it to set the slope, the equivalent is "
                                "fuzzyness = 1/(1+expk) for expk > 0, which CHANGES the waveform." ) ;

                            if ( tSection->key_exists( "type" ) )
                            {
                                const string tFType =  tSection->get_string( "type" ) ;
                                if ( tFType == "ramp" )
                                {
                                    //First check if the required data is there
                                    BELFEM_ERROR( tSection->key_exists( "amplitude" ), "Amplitude not defined in the boundary condition" ) ;
                                    BELFEM_ERROR( tSection->key_exists( "period" ), "Period not defined in the boundary condition" ) ;
                                    BELFEM_ERROR( tSection->key_exists( "offset" ), "Offset not defined in the boundary condition" ) ;

                                    //Then check the units
                                    string tUnitDef = tSection->get_units( "amplitude" );
                                    value tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, tUnits ),
                                            "required unit for the amplitude of the boundary condition is: %s",
                                            tUnits.c_str() );
                                    mPhysicalBoundaryConditions( tCount-i )->set_units( tValue.second ) ;

                                    tUnitDef = tSection->get_units( "period" );
                                    tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, "s" ),
                                            "required unit for the period of the boundary condition is: s" );

                                    tUnitDef = tSection->get_units( "offset" );
                                    tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, "s" ),
                                            "required unit for the offset of the boundary condition is: s" );


                                    //Then set the function
                                    tFunction->set_ramp( tSection->get_value( "amplitude",tUnits ).first,
                                                         tSection->get_value( "period","s" ).first,
                                                         tSection->get_value( "offset","s" ).first ) ;
                                }
                                else if ( tFType == "sigmoid" )
                                {
                                    //First check if the required data is there
                                    BELFEM_ERROR( tSection->key_exists( "amplitude" ), "Amplitude not defined in the boundary condition" ) ;
                                    BELFEM_ERROR( tSection->key_exists( "period" ), "Period not defined in the boundary condition" ) ;
                                    BELFEM_ERROR( tSection->key_exists( "offset" ), "Offset not defined in the boundary condition" ) ;

                                    //Then check the units
                                    string tUnitDef = tSection->get_units( "amplitude" );
                                    value tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, tUnits ),
                                            "required unit for the amplitude of the boundary condition is: %s",
                                            tUnits.c_str() );
                                    mPhysicalBoundaryConditions( tCount-i )->set_units( tValue.second ) ;

                                    tUnitDef = tSection->get_units( "period" );
                                    tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, "s" ),
                                            "required unit for the period of the boundary condition is: s" );

                                    tUnitDef = tSection->get_units( "offset" );
                                    tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, "s" ),
                                            "required unit for the offset of the boundary condition is: s" );

                                    //Then set the function
                                    tFunction->set_sigmoid( tSection->get_value( "amplitude",tUnits ).first,
                                                            tSection->get_value( "period","s" ).first,
                                                            tSection->get_value( "offset","s" ).first,
                                                            tSection->key_exists( "fuzzyness" )?tSection->get_real( "fuzzyness" ):0.01 ) ;
                                }
                                else if ( tFType == "sine" || tFType == "square" ||
                                          tFType == "triangle" || tFType == "sawtooth" )
                                {

                                    //First check if the required data is there
                                    BELFEM_ERROR( tSection->key_exists( "amplitude" ), "Amplitude not defined in the boundary condition" ) ;

                                    //Then check the units
                                    string tUnitDef = tSection->get_units( "amplitude" );
                                    value tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, tUnits ),
                                            "required unit for the amplitude of the boundary condition is: %s",
                                            tUnits.c_str() );
                                    mPhysicalBoundaryConditions( tCount-i )->set_units( tValue.second ) ;

                                    real tPeriod = 0.0;
                                    if ( tSection->key_exists( "period" ) )
                                    {
                                        tUnitDef = tSection->get_units( "period" );
                                        tValue = unit_to_si( tUnitDef );
                                        BELFEM_ERROR(
                                                check_unit( tValue, "s" ),
                                                "required unit for the  of the boundary condition is: s" ) ;

                                        tPeriod = tSection->get_value( "period","s" ).first ;
                                    }
                                    else if ( tSection->key_exists( "frequency" ) )
                                    {
                                        tUnitDef = tSection->get_units( "frequency" );
                                        tValue = unit_to_si( tUnitDef );
                                        BELFEM_ERROR(
                                                check_unit( tValue, "Hz" ),
                                                "required unit for the  of the boundary condition is: Hz" ) ;

                                        tPeriod = 1.0/( tSection->get_value( "frequency","Hz" ).first ) ;
                                    }
                                    else
                                    {
                                        BELFEM_ERROR( false,"Period or frequency undefined" ) ;
                                    }

                                    //Then set the function
                                    tFunction->set_periodic( boundary_condition_function_type( tFType ),
                                                             tSection->get_value( "amplitude",tUnits ).first,
                                                             tPeriod,
                                                             tSection->key_exists( "phase" )?tSection->get_value( "phase","rad" ).first:0.0 ) ;
                                }
                                else if ( tFType == "constant" )
                                {

                                    //First check if the required data is there
                                    BELFEM_ERROR( tSection->key_exists( "amplitude" ), "Amplitude not defined in the boundary condition" ) ;

                                    //Then check the units
                                    string tUnitDef = tSection->get_units( "amplitude" );
                                    value tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, tUnits ),
                                            "required unit for the amplitude of the boundary condition is: %s",
                                            tUnits.c_str() );
                                    mPhysicalBoundaryConditions( tCount-i )->set_units( tValue.second ) ;

                                    //Then set the function
                                    tFunction->set_constant( tSection->get_value( "amplitude",tUnits ).first ) ;
                                }
                                else if ( tFType == "userdefined" )
                                {

                                    //First check if the required data is there
                                    BELFEM_ERROR( tSection->key_exists( "file" ), "File not defined for user-defined function in the boundary condition" ) ;
                                    BELFEM_ERROR( tSection->key_exists( "label" ), "Label not defined for user-defined function in the boundary condition" ) ;
                                    BELFEM_ERROR( tSection->key_exists( "units" ), "Units not defined for user-defined function in the boundary condition" ) ;

                                    string tUnitDef = tSection->get_string( "units" );
                                    value tValue = unit_to_si( tUnitDef );
                                    BELFEM_ERROR(
                                            check_unit( tValue, tUnits ),
                                            "required unit for the amplitude of the boundary condition is: %s",
                                            tUnits.c_str() );
                                    mPhysicalBoundaryConditions( tCount-i )->set_units( tValue.second ) ;
                                    mPhysicalBoundaryConditions( tCount-i )->scale()*=tValue.first ;

                                    //Then set the function
                                    tFunction->read_user_defined(tSection->get_string( "file" ), tSection->get_string( "label" )) ;
                                }
                                else
                                {
                                    BELFEM_ERROR( false, "Undefined boundary condition type" ) ;
                                }
                            }
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Undefined Boundary Condition" );
                        }

                    }
                    mPhysicalBoundaryConditions( tCount-i )->set_function( tFunction ) ;
                }


            }
        }

//-----------------------------------------------------------------------------

        Cell < PhysicalBoundaryCondition * > &
        MaxwellBoundaryConditionFactory::boundary_conditions()
        {
            return mPhysicalBoundaryConditions ;
        }

//-----------------------------------------------------------------------------

        Cell < PhysicalBoundaryCondition * >
        MaxwellBoundaryConditionFactory::current_boundary_conditions()
        {
            Cell< PhysicalBoundaryCondition * > aCurrentBCs = Cell< PhysicalBoundaryCondition * >();
            for (PhysicalBoundaryCondition * tBC: mPhysicalBoundaryConditions)
            {
                if (tBC->type() == BoundaryConditionType::Current || tBC->type() == BoundaryConditionType::CircuitCurrent)
                {
                    aCurrentBCs.push( tBC );
                }
            }
            return aCurrentBCs;
        }


//-------------------------------------------------------------------------------

        Cell < PhysicalBoundaryCondition * >
        MaxwellBoundaryConditionFactory::voltage_boundary_conditions()
        {
            Cell < PhysicalBoundaryCondition * > aVoltageBCs = Cell < PhysicalBoundaryCondition * >();
            for(PhysicalBoundaryCondition * tBC : mPhysicalBoundaryConditions)
            {
                if (tBC->type() == BoundaryConditionType::Voltage)
                {
                    aVoltageBCs.push(tBC);
                }
            }
            return aVoltageBCs ;
        }

//-------------------------------------------------------------------------------


        void
        MaxwellBoundaryConditionFactory::set_fields( DofManager * aField )
        {
            for ( PhysicalBoundaryCondition * tBC : mPhysicalBoundaryConditions )
            {
                tBC->set_field( aField ) ;
            }
        }

//-------------------------------------------------------------------------------
    }
}
