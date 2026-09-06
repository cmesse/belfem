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

#include "commtools.hpp"
#include "cl_ElectricalCircuitFactory.hpp"

#include "fn_check_unit.hpp"
#include "cl_NetlistParser.hpp"
#include "cl_NgspiceCircuitFactory.hpp"
#include "stringtools.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        ElectricalCircuitFactory::ElectricalCircuitFactory( const string & aInputFile, Cell< fem::PhysicalBoundaryCondition * > & aPhysicalBoundaryConditions ) :
                mCommRank( comm_rank() ),
                mInputFile( new InputFile( aInputFile ) ),
                mPhysicalBoundaryConditions(aPhysicalBoundaryConditions)
        {
            if (mInputFile->section_exists("circuit") )
            {
                // hybrid path ( plan decision 6b ): the netlist holds the
                // lumped topology, the deck holds terminal pairs and output
                if ( mInputFile->section("circuit")->key_exists("file") )
                {
                    this->read_circuit_from_netlist( mInputFile->section("circuit") ) ;
                }
                else
                {
                    BELFEM_ERROR(mInputFile->section("circuit")->key_exists("number of nodes"),
                                 "The electrical circuit should be defined with a number of nodes");

                    mCircuit.reset( new ElectricalCircuit(mInputFile->section("circuit")->get_int("number of nodes")) );
                    this->read_circuit( mInputFile->section("circuit") ) ;
                }
            }
        }


//----------------------------------------------------------------------------

        ElectricalCircuitFactory::~ElectricalCircuitFactory() = default;
        // the unique_ptr members also clean up when a refusal throws
        // mid-construction ( the dtor of a half-built object never runs )

//----------------------------------------------------------------------------

        void
        ElectricalCircuitFactory::read_circuit( const input::Section * aSection )
        {
            //First, read the topology of the circuit
            BELFEM_ERROR(aSection->section("topology"), "No topology defined for the circuit") ;
            const input::Section * tTopologySection = aSection->section("topology") ;

            uint tCount = mPhysicalBoundaryConditions.size() ;

            for( id_t i = 0 ; i < tTopologySection->num_sections(); ++i)
            {

                const input::Section * tComponentSection = tTopologySection->section(i) ;
                ComponentType tComponentType = component_type(tComponentSection->key()) ;

                //Read the label if it exists
                string tLabel = "" ;
                if(tComponentSection->key_exists("label"))
                {
                    tLabel = tComponentSection->get_string("label") ;
                }

                //Read the nodes of the component
                BELFEM_ERROR(tComponentSection->key_exists("node +"), "Every electrical component in the circuit must have a node +") ;
                BELFEM_ERROR(tComponentSection->key_exists("node -"), "Every electrical component in the circuit must have a node -") ;
                Vector< id_t > tNodes =  Vector< id_t >(2,0);
                tNodes(0) = tComponentSection->get_int("node +") ;
                tNodes(1) = tComponentSection->get_int("node -") ;

                BELFEM_ERROR(tNodes(0) < mCircuit->number_of_nodes(), "Invalid node + index for %s in circuit",
                             to_string(tComponentType).c_str()) ;

                BELFEM_ERROR(tNodes(1) < mCircuit->number_of_nodes(), "Invalid node - index for %s in circuit",
                             to_string(tComponentType).c_str()) ;

                //Insert the correct component in the circuit
                switch (tComponentType)
                {
                    case ComponentType::RESISTOR :
                    {
                        BELFEM_ERROR(tComponentSection->key_exists("value"), "Undefined value of the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;

                        string tUnits = "Ohm";
                        string tUnitDef = tComponentSection->get_units( "value" );
                        value tValue = unit_to_si( tUnitDef );
                        BELFEM_ERROR(
                                check_unit( tValue, tUnits ),
                                "required unit for the %s is: %s",
                                to_string(tComponentType).c_str(), tUnits.c_str() );

                        mCircuit->create_resistor(tComponentSection->get_value( "value",tUnits ).first,tNodes(0),tNodes(1),tLabel) ;
                        break ;
                    }
                    case ComponentType::CAPACITOR :
                    {
                        BELFEM_ERROR(tComponentSection->key_exists("value"), "Undefined value of the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;

                        string tUnits = "F";
                        string tUnitDef = tComponentSection->get_units( "value" );
                        value tValue = unit_to_si( tUnitDef );
                        BELFEM_ERROR(
                                check_unit( tValue, tUnits ),
                                "required unit for the %s is: %s",
                                to_string(tComponentType).c_str(), tUnits.c_str() );

                        uint tOrder = 1;
                        if (tComponentSection->key_exists("order")) tOrder = tComponentSection->get_int("order") ;
                        mCircuit->create_capacitor(tComponentSection->get_value( "value",tUnits ).first,tOrder,tNodes(0),tNodes(1),tLabel) ;
                        break ;
                    }
                    case ComponentType::INDUCTOR :
                    {
                        BELFEM_ERROR(tComponentSection->key_exists("value"), "Undefined value of the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;

                        string tUnits = "H";
                        string tUnitDef = tComponentSection->get_units( "value" );
                        value tValue = unit_to_si( tUnitDef );
                        BELFEM_ERROR(
                                check_unit( tValue, tUnits ),
                                "required unit for the %s is: %s",
                                to_string(tComponentType).c_str(), tUnits.c_str() );

                        uint tOrder = 1;
                        if (tComponentSection->key_exists("order")) tOrder = tComponentSection->get_int("order") ;
                        mCircuit->create_inductor(tComponentSection->get_value( "value",tUnits ).first,tOrder,tNodes(0),tNodes(1),tLabel) ;
                        break ;
                    }
                    case ComponentType::VOLTAGESOURCE :
                    case ComponentType::CURRENTSOURCE :
                    {
                        string tUnits ;
                        if (tComponentType == ComponentType::VOLTAGESOURCE) tUnits = "V" ;
                        else if (tComponentType == ComponentType::CURRENTSOURCE) tUnits = "A" ;
                        // # CHANGE THE DEFINITION OF THE VOLTAGE SOURCE TO ADD A BOUNDARY CONDITION FUNCTION AS AN ATTRIBUTE

                        // 'expk' was parsed and stored for years without ever reaching the
                        // waveform. Refusing it beats dropping the read in silence, which is
                        // what this check exists to close. This sits BEFORE the type switch
                        // on purpose: the key was only ever read on the sigmoid branch, so a
                        // sigmoid-only check would leave it silently ignored on every other
                        // source shape, and on a section with no 'type' at all
                        BELFEM_ERROR( ! tComponentSection->key_exists( "expk" ),
                            "'expk' is no longer accepted, and it never affected the waveform: "
                            "function_sigmoid takes its rate from 'fuzzyness' and 'period' alone. "
                            "DELETE the key to keep the results you get today. Do NOT translate it "
                            "-- if you meant it to set the slope, the equivalent is "
                            "fuzzyness = 1/(1+expk) for expk > 0, which CHANGES the waveform." ) ;

                        BELFEM_ERROR(tComponentSection->key_exists("type"),"Source type not defined in the circuit") ;

                        const string tFType =  tComponentSection->get_string( "type" ) ;
                        SourceFunction * tFunction = new SourceFunction() ;

                        if ( tFType == "ramp" )
                        {
                            //First check if the required data is there
                            BELFEM_ERROR( tComponentSection->key_exists( "amplitude" ), "Amplitude not defined for the source in the circuit" ) ;
                            BELFEM_ERROR( tComponentSection->key_exists( "period" ), "Period not defined for the source in the circuit" ) ;
                            BELFEM_ERROR( tComponentSection->key_exists( "offset" ), "Offset not defined for the source in the circuit" ) ;

                            //Then check the units
                            string tUnitDef = tComponentSection->get_units( "amplitude" );
                            value tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, tUnits ),
                                    "required unit for the amplitude of the source in the circuit is: %s",
                                    tUnits.c_str() );

                            tUnitDef = tComponentSection->get_units( "period" );
                            tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, "s" ),
                                    "required unit for the period of the source in the circuit is: s" );

                            tUnitDef = tComponentSection->get_units( "offset" );
                            tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, "s" ),
                                    "required unit for the offset of the source in the circuit is: s" );


                            //Then set the function
                            tFunction->set_ramp( tComponentSection->get_value( "amplitude",tUnits ).first,
                                                 tComponentSection->get_value( "period","s" ).first,
                                                 tComponentSection->get_value( "offset","s" ).first ) ;
                        }
                        else if ( tFType == "sigmoid" )
                        {
                            //First check if the required data is there
                            BELFEM_ERROR( tComponentSection->key_exists( "amplitude" ), "Amplitude not defined for the source in the circuit" ) ;
                            BELFEM_ERROR( tComponentSection->key_exists( "period" ), "Period not defined for the source in the circuit" ) ;
                            BELFEM_ERROR( tComponentSection->key_exists( "offset" ), "Offset not defined for the source in the circuit" ) ;

                            //Then check the units
                            string tUnitDef = tComponentSection->get_units( "amplitude" );
                            value tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, tUnits ),
                                    "required unit for the amplitude of the source in the circuit is: %s",
                                    tUnits.c_str() );

                            tUnitDef = tComponentSection->get_units( "period" );
                            tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, "s" ),
                                    "required unit for the period of the source in the circuit is: s" );

                            tUnitDef = tComponentSection->get_units( "offset" );
                            tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, "s" ),
                                    "required unit for the offset of the source in the circuit is: s" );

                            //Then set the function
                            tFunction->set_sigmoid( tComponentSection->get_value( "amplitude",tUnits ).first,
                                                    tComponentSection->get_value( "period","s" ).first,
                                                    tComponentSection->get_value( "offset","s" ).first,
                                                    tComponentSection->key_exists( "fuzzyness" )?tComponentSection->get_real( "fuzzyness" ):0.01 ) ;
                        }
                        else if ( tFType == "sine" || tFType == "square" ||
                                  tFType == "triangle" || tFType == "sawtooth" )
                        {

                            //First check if the required data is there
                            BELFEM_ERROR( tComponentSection->key_exists( "amplitude" ), "Amplitude not defined for the source in the circuit" ) ;

                            //Then check the units
                            string tUnitDef = tComponentSection->get_units( "amplitude" );
                            value tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, tUnits ),
                                    "required unit for the amplitude of the source in the circuit is: %s",
                                    tUnits.c_str() );

                            real tPeriod = 0.0;
                            if ( tComponentSection->key_exists( "period" ) )
                            {
                                tUnitDef = tComponentSection->get_units( "period" );
                                tValue = unit_to_si( tUnitDef );
                                BELFEM_ERROR(
                                        check_unit( tValue, "s" ),
                                        "required unit for the period of the source in the circuit is: s" ) ;

                                tPeriod = tComponentSection->get_value( "period","s" ).first ;
                            }
                            else if ( tComponentSection->key_exists( "frequency" ) )
                            {
                                tUnitDef = tComponentSection->get_units( "frequency" );
                                tValue = unit_to_si( tUnitDef );
                                BELFEM_ERROR(
                                        check_unit( tValue, "Hz" ),
                                        "required unit for the frequency of the source in the circuit is: Hz" ) ;

                                tPeriod = 1.0/( tComponentSection->get_value( "frequency","Hz" ).first ) ;
                            }
                            else
                            {
                                BELFEM_ERROR( false,"Period or frequency undefined for the source in the circuit" ) ;
                            }

                            //Then set the function
                            tFunction->set_periodic( boundary_condition_function_type( tFType ),
                                                     tComponentSection->get_value( "amplitude",tUnits ).first,
                                                     tPeriod,
                                                     tComponentSection->key_exists( "phase" )?tComponentSection->get_value( "phase","rad" ).first:0.0 ) ;
                        }
                        else if ( tFType == "constant" )
                        {

                            //First check if the required data is there
                            BELFEM_ERROR( tComponentSection->key_exists( "amplitude" ), "Amplitude not defined for the source in the circuit" ) ;

                            //Then check the units
                            string tUnitDef = tComponentSection->get_units( "amplitude" );
                            value tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, tUnits ),
                                    "required unit for the amplitude of the source in the circuit is: %s",
                                    tUnits.c_str() );

                            //Then set the function
                            tFunction->set_constant( tComponentSection->get_value( "amplitude",tUnits ).first ) ;
                        }
                        else if ( tFType == "userdefined" )
                        {

                            //First check if the required data is there
                            BELFEM_ERROR( tComponentSection->key_exists( "file" ), "File not defined for user-defined function in the boundary condition" ) ;
                            BELFEM_ERROR( tComponentSection->key_exists( "label" ), "Label not defined for user-defined function in the boundary condition" ) ;
                            BELFEM_ERROR( tComponentSection->key_exists( "units" ), "Units not defined for user-defined function in the boundary condition" ) ;

                            string tUnitDef = tComponentSection->get_string( "units" );
                            value tValue = unit_to_si( tUnitDef );
                            BELFEM_ERROR(
                                    check_unit( tValue, tUnits ),
                                    "required unit for the amplitude of the boundary condition is: %s",
                                    tUnits.c_str() );
                            mPhysicalBoundaryConditions( tCount-i )->set_units( tValue.second ) ;
                            mPhysicalBoundaryConditions( tCount-i )->scale()*=tValue.first ;

                            //Then set the function
                            tFunction->read_user_defined(tComponentSection->get_string( "file" ), tComponentSection->get_string( "label" )) ;
                        }
                        else
                        {
                            BELFEM_ERROR( false, "Undefined source type in the circuit" ) ;
                        }
                        if (tComponentType == ComponentType::VOLTAGESOURCE) mCircuit->create_voltage_source(tFunction,tNodes(0),tNodes(1),tLabel) ;
                        else if (tComponentType == ComponentType::CURRENTSOURCE) mCircuit->create_current_source(tFunction,tNodes(0),tNodes(1),tLabel) ;

                        break ;
                    }
                    case ComponentType::SWITCH:
                    {
                        //First get the initial state of the switch
                        BELFEM_ERROR(tComponentSection->key_exists("initial state"), "Undefined initial state of the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;

                        bool tIsClosed = false;
                        string tState = tComponentSection->get_string("initial state") ;
                        if (tState == "closed")
                        {
                            tIsClosed = true ;
                        }
                        else if (tState == "open")
                        {
                            tIsClosed = false ;
                        }
                        else
                        {
                            BELFEM_ERROR(false,"The switch state is either open or closed") ;
                        }

                        //Then get the switch time
                        BELFEM_ERROR(tComponentSection->key_exists("switch time"), "Undefined switch time for the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;
                        string tUnits = "s";
                        string tUnitDef = tComponentSection->get_units( "switch time" );
                        value tValue = unit_to_si( tUnitDef );
                        BELFEM_ERROR(
                                check_unit( tValue, tUnits ),
                                "required unit for the switch time is: %s",
                                tUnits.c_str() );

                        mCircuit->create_switch(tIsClosed,tComponentSection->get_value( "switch time",tUnits ).first,tNodes(0),tNodes(1),tLabel) ;
                        break ;
                    }
                    case ComponentType::DIODE:
                    {
                        real tIs = 0.1e-3 ; //(A) standard value for a diode
                        real tVt = 0.026 ; //(V) standard value for a diode

                        //Overwrite the values of Is and Vt if the were defined by the user
                        if (tComponentSection->key_exists("Is"))
                        {
                            string tUnitsIs = "A";
                            string tUnitDefIs = tComponentSection->get_units( "Is" );
                            value tValueIs = unit_to_si( tUnitDefIs );
                            BELFEM_ERROR(
                                    check_unit( tValueIs, tUnitsIs ),
                                    "required unit for the %s is: %s",
                                    to_string(tComponentType).c_str(), tUnitsIs.c_str() );

                            tIs = tComponentSection->get_value( "Is",tUnitsIs ).first ;
                        }

                        if (tComponentSection->key_exists("Vt"))
                        {
                            string tUnitsVt = "V";
                            string tUnitDefVt = tComponentSection->get_units( "Vt" );
                            value tValueVt = unit_to_si( tUnitDefVt );
                            BELFEM_ERROR(
                                    check_unit( tValueVt, tUnitsVt ),
                                    "required unit for the %s is: %s",
                                    to_string(tComponentType).c_str(), tUnitsVt.c_str() );

                            tVt = tComponentSection->get_value( "Vt",tUnitsVt ).first ;
                        }

                        mCircuit->create_diode(tIs, tVt,tNodes(0),tNodes(1),tLabel) ;
                        break ;
                    }
                    case ComponentType::SUPERCONDUCTOR:
                    {

                        real tIc ;
                        real tn ;
                        real tEc ;
                        real tLength ;

                        //Get the Ic value
                        BELFEM_ERROR(tComponentSection->key_exists("Ic"), "Undefined Ic in the circuit");
                        string tUnitsIc = "A";
                        string tUnitDefIc = tComponentSection->get_units( "Ic" );
                        value tValueIc = unit_to_si( tUnitDefIc );
                        BELFEM_ERROR(
                                check_unit( tValueIc, tUnitsIc ),
                                "required unit for Ic is: %s",
                                tUnitsIc.c_str() );

                        tIc = tComponentSection->get_value( "Ic",tUnitsIc ).first ;

                        //Get the n value
                        BELFEM_ERROR(tComponentSection->key_exists("n"), "Undefined n for the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;
                        string tUnitsn = "";
                        string tUnitDefn = tComponentSection->get_units( "n" );
                        value tValuen = unit_to_si( tUnitDefn );
                        BELFEM_ERROR(
                                check_unit( tValuen, tUnitsn ),
                                "required unit for n is: %s",
                                tUnitsn.c_str() );

                        tn = tComponentSection->get_value( "n",tUnitsn ).first ;

                        //Get the Ec value
                        BELFEM_ERROR(tComponentSection->key_exists("Ec"), "Undefined Ec for the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;
                        string tUnitsEc = "V/m";
                        string tUnitDefEc = tComponentSection->get_units( "Ec" );
                        value tValueEc = unit_to_si( tUnitDefEc );
                        BELFEM_ERROR(
                                check_unit( tValueEc, tUnitsEc ),
                                "required unit for Ec is: %s",
                                tUnitsEc.c_str() );

                        tEc = tComponentSection->get_value( "Ec",tUnitsEc ).first ;

                        //Get the length value
                        BELFEM_ERROR(tComponentSection->key_exists("length"), "Undefined length for the %s in the circuit",
                                     to_string(tComponentType).c_str()) ;
                        string tUnitsLength = "m";
                        string tUnitDefLength = tComponentSection->get_units( "length" );
                        value tValueLength = unit_to_si( tUnitDefLength );
                        BELFEM_ERROR(
                                check_unit( tValueLength, tUnitsLength ),
                                "required unit for length is: %s",
                                tUnitsLength.c_str() );

                        tLength = tComponentSection->get_value( "length",tUnitsLength ).first ;

                        mCircuit->create_superconductor(tIc, tn, tEc, tLength, tNodes(0),tNodes(1),tLabel) ;
                        break ;
                    }
                    case ComponentType::TERMINALPAIR:
                    {
                        this->read_terminal_pair( tComponentSection, i,
                                                  tNodes(0), tNodes(1),
                                                  tLabel, tCount );
                        break ;
                    }
                    default:
                    {
                        break ;
                    }
                }
            }

            //Then, read the output
            this->read_output( aSection, nullptr );
        }

//----------------------------------------------------------------------------

        void
        ElectricalCircuitFactory::read_terminal_pair(
                const input::Section * aComponentSection,
                const id_t aSectionIndex,
                const index_t aNodePlus,
                const index_t aNodeMinus,
                const string & aLabel,
                uint & aCount )
        {
            // local aliases keep the body identical to the pre-refactor
            // switch case it was lifted from
            const input::Section * tComponentSection = aComponentSection;
            const id_t i = aSectionIndex;
            const string & tLabel = aLabel;
            uint & tCount = aCount;
                    {

                        //Create the terminal pair in the circuit
                        mCircuit->create_terminal_pair(aNodePlus, aNodeMinus,tLabel);

                        // The terminal pair represents a pair of terminals in the FEM problem.
                        // Therefore, they are associated with boundary conditions that we create here
                        string tKey ;
                        bool tIsThinShell = false ;
                        if (tComponentSection->key_exists( "input terminal" ))
                        {
                            tKey = "input terminal" ;
                        }
                        else if (tComponentSection->key_exists( "input terminals" ))
                        {
                            tKey = "input terminals" ;
                        }
                        else if (tComponentSection->key_exists( "input curve" ) )
                        {
                            tKey = "input curve" ;
                            tIsThinShell = true ;
                        }
                        else if (tComponentSection->key_exists( "input curves" ))
                        {
                            tKey = "input curves" ;
                            tIsThinShell = true ;
                        }
                        else
                        {
                            BELFEM_ERROR( false , "Input terminals or curves undefined for current boundary condition" ) ;
                        }

                        Vector < id_t > tDomains ;
                        fem::BoundaryConditionType tType = fem::BoundaryConditionType::CircuitVoltage;


                        //Look for output terminals
                        bool tOutputExists = true ;
                        string tKey2 ;
                        if (tComponentSection->key_exists( "output terminal" ))
                        {
                            tKey2 = "output terminal" ;
                        }
                        else if (tComponentSection->key_exists( "output terminals" ))
                        {
                            tKey2 = "output terminals" ;
                        }
                        else if (tComponentSection->key_exists( "output curve" ) )
                        {
                            tKey2 = "output curve" ;
                        }
                        else if (tComponentSection->key_exists( "output curves" ))
                        {
                            tKey2 = "output curves" ;
                        }
                        else
                        {
                            tOutputExists = false ;
                        }

                        //Separate the groups defined with []. Each group becomes one boundary condition
                        Cell<Cell < id_t > >tDomainGroupsIn ;
                        tComponentSection->get_id_groups(tKey,tDomainGroupsIn) ;

                        // FROZEN GRAMMAR ( 2026-08-15 ): a terminal
                        // pair creates ONE lumped component
                        // ( create_terminal_pair above runs once ), but this
                        // walk creates one boundary condition PER bracket
                        // group — and the controller pairs conditions with
                        // components by POSITION, so every group past the
                        // first would read beyond the component list: a
                        // silent index overrun. Until an N:1 pairing is
                        // designed, more than one group is refused here, on
                        // every rank ( the factory ctor is not rank-guarded )
                        {
                            string tBase = tLabel.size() > 0 ?
                                    tLabel : "terminalpair_" + std::to_string( i + 1 ) ;
                            BELFEM_ERROR( tDomainGroupsIn.size() == 1,
                                "terminal pair '%s' defines %u bracket groups, but the circuit coupling\n"
                                "supports exactly one boundary condition per terminal pair ( the controller\n"
                                "pairs conditions with components by position ). Write the terminals as ONE\n"
                                "bracketed group, e.g. [1,2] — an unbracketed comma list creates one group\n"
                                "per id. For several pairs, define several components.",
                                tBase.c_str(), ( unsigned int ) tDomainGroupsIn.size() );
                        }

                        real tScale = 1.0 ;
                        Cell<Cell < id_t > >tDomainGroupsOut ;
                        if (tOutputExists)
                        {
                            tComponentSection->get_id_groups(tKey2,tDomainGroupsOut) ;
                            BELFEM_ERROR(tDomainGroupsIn.size() == tDomainGroupsOut.size() ,
                                        "The current/voltage boundary condition should contain the same number of input and output terminals") ;
                        }
                        else //If the output terminals are not defined, the problem is 2-D, we read only the input terminals instead
                        {
                            tComponentSection->get_id_groups(tKey,tDomainGroupsOut) ;
                            // a terminal pair is always CircuitVoltage ( set
                            // above and never reassigned ), so the old test
                            // against plain Voltage never fired and the 2-D
                            // scale silently stayed 1.0
                            if (tType == fem::BoundaryConditionType::CircuitVoltage)
                            {
                                BELFEM_ERROR( tComponentSection->key_exists("length"), "Undefined length for voltage boundary condition in 2-D" ) ;
                                string tUnitDef = tComponentSection->get_units( "length" );
                                value tValue = unit_to_si( tUnitDef );
                                BELFEM_ERROR(
                                            check_unit( tValue, "m" ),
                                            "required unit for the amplitude of the boundary condition is: m");

                                real tLength = tComponentSection->get_value( "length","m" ).first ;

                                // we are about to divide by this: zero gives an
                                // infinite scale and a negative value silently
                                // flips the sign of the imposed voltage
                                BELFEM_ERROR( tLength > 0.0,
                                              "length for the 2-D voltage boundary condition must be positive, but is %f m",
                                              tLength );

                                tScale *= 1.0/tLength ;

                            }
                        }


                        for(uint k = 0; k < tDomainGroupsIn.size(); ++k)
                        {
                            //Grouping the input and output terminals all into tDomains
                            tDomains.set_size(tDomainGroupsIn(k).size()+tDomainGroupsOut(k).size(),0) ;
                            for(uint j = 0 ; j < tDomainGroupsIn(k).size();++j)
                            {
                                tDomains(j) = tDomainGroupsIn(k)(j) ;
                            }
                            for(uint j = 0 ; j < tDomainGroupsOut(k).size();++j)
                            {
                                tDomains(tDomainGroupsIn(k).size()+j) = tDomainGroupsOut(k)(j) ;
                            }

                            mPhysicalBoundaryConditions.push(new fem::PhysicalBoundaryCondition( tCount ));
                            mPhysicalBoundaryConditions( tCount )->set_type(tType) ;
                            mPhysicalBoundaryConditions( tCount )->set_scale( tScale ) ;

                            // global-variable name: the component label
                            // ( existing circuit grammar ), else a running
                            // terminal-pair name. One group per pair
                            // ( enforced above ), so no group suffix
                            {
                                string tBase = tLabel.size() > 0 ?
                                        tLabel : "terminalpair_" + std::to_string( i + 1 ) ;
                                mPhysicalBoundaryConditions( tCount )->set_label( tBase );
                            }

                            mPhysicalBoundaryConditions(tCount++)->set_domains(tDomains, tIsThinShell ) ;

                        }
                    }
        }

//----------------------------------------------------------------------------

        void
        ElectricalCircuitFactory::read_output( const input::Section * aSection,
                                               const NgspiceCircuitFactory * aNetlist )
        {
            if (aSection->section_exists("output"))
            {
                const input::Section * tOutputSection = aSection->section("output") ;

                //first set the output file name
                BELFEM_ERROR(tOutputSection->key_exists("file"),"Undefined file name for circuit outputs") ;
                mCircuit->set_file( tOutputSection->get_string("file")) ;

                // Then list the components and nodes for current/voltage output
                if (tOutputSection->key_exists("currents"))
                {
                    Cell< string > tCurrents = string_to_words(
                            search_and_replace( search_and_replace(
                                    tOutputSection->get_string("currents"), " ","" ),","," "));

                    // netlist labels are case-folded instance names ( O5 );
                    // fold the deck side so the lookup cannot miss on case
                    if ( aNetlist != nullptr )
                    {
                        for ( index_t k = 0; k < tCurrents.size(); ++k )
                        {
                            tCurrents( k ) = string_to_lower( tCurrents( k ) );
                        }
                    }
                    mCircuit->set_output_currents(tCurrents) ;

                }

                if (tOutputSection->key_exists("voltages"))
                {
                    Vector< id_t > tVoltages ;
                    if ( aNetlist == nullptr )
                    {
                        tOutputSection->get_ids("voltages", tVoltages) ;
                    }
                    else
                    {
                        // hybrid grammar ( O5 ): voltages are netlist node
                        // NAMES, resolved through the same map as the
                        // terminal pairs -- never raw indices, which would
                        // be ambiguous against the ground remap
                        Cell< string > tNames = string_to_words(
                                search_and_replace( search_and_replace(
                                        tOutputSection->get_string("voltages"), " ","" ),","," "));
                        tVoltages.set_size( tNames.size(), 0 );
                        for ( index_t k = 0; k < tNames.size(); ++k )
                        {
                            tVoltages( k ) = aNetlist->node_index( tNames( k ) );

                            // the CSV header prints the internal index --
                            // give the user the decoder ring
                            if ( comm_rank() == 0 )
                            {
                                message( InfoLevel::Default,
                                         "    circuit output: node %-12s -> column V%lu",
                                         tNames( k ).c_str(),
                                         ( long unsigned int ) tVoltages( k ) );
                            }
                        }
                    }
                    mCircuit->set_output_voltages(tVoltages) ;
                }

                mCircuit->init_output_file() ;

            }
        }

//----------------------------------------------------------------------------

        void
        ElectricalCircuitFactory::read_circuit_from_netlist( const input::Section * aSection )
        {
            // the node count is derived from the netlist -- restating it
            // invites a silent mismatch
            BELFEM_ERROR( !aSection->key_exists("number of nodes"),
                          "circuit: 'number of nodes' must not be set when a netlist "
                          "'file' is given -- the count is derived from the netlist" );

            const string tPath = aSection->get_string("file");

            NetlistParser tNetlist( tPath );
            NgspiceCircuitFactory tNetlistFactory( tNetlist );

            // O7: the deck's solver{timestep{}} owns time integration; a
            // .tran card is stored by the parser but deliberately unused
            if ( tNetlist.controls().size() > 0 && comm_rank() == 0 )
            {
                message( InfoLevel::Default,
                         "    circuit: ignoring .tran in %s -- solver{timestep{}} "
                         "owns the time integration",
                         tPath.c_str() );
            }

            // the lumped circuit comes from the netlist; this factory owns it
            mCircuit.reset( tNetlistFactory.circuit() );

            uint tCount = mPhysicalBoundaryConditions.size() ;

            // the deck may only contribute terminal pairs -- lumped elements
            // live in the netlist ( O6 )
            if ( aSection->section_exists("topology") )
            {
                const input::Section * tTopologySection = aSection->section("topology") ;

                for( id_t i = 0 ; i < tTopologySection->num_sections(); ++i)
                {
                    const input::Section * tComponentSection = tTopologySection->section(i) ;

                    BELFEM_ERROR( component_type(tComponentSection->key())
                                  == ComponentType::TERMINALPAIR,
                                  "circuit topology section '%s': with a netlist file, "
                                  "the deck may only define terminal pairs -- lumped "
                                  "elements belong in %s",
                                  tComponentSection->key().c_str(),
                                  tPath.c_str() );

                    // hybrid mode is case-insensitive throughout ( O5,
                    // one fold policy ): the label folds like the netlist
                    // instance names, so the currents lookup cannot miss
                    string tLabel = "" ;
                    if(tComponentSection->key_exists("label"))
                    {
                        tLabel = string_to_lower(
                                tComponentSection->get_string("label") );
                    }

                    // a label shared with a netlist instance ( or another
                    // pair ) would make the current output ambiguous --
                    // every match writes a column under one header
                    if ( tLabel.size() > 0 )
                    {
                        for ( uint k = 0; k < mCircuit->number_of_components(); ++k )
                        {
                            BELFEM_ERROR( mCircuit->component( k )->get_label() != tLabel,
                                          "terminal pair label '%s' is already used by "
                                          "another component ( netlist instance or pair )",
                                          tLabel.c_str() );
                        }
                    }

                    // hybrid grammar ( O5 ): node references are netlist
                    // node NAMES ( "0"/"gnd" included ), resolved through
                    // the netlist's map -- never raw indices
                    BELFEM_ERROR(tComponentSection->key_exists("node +"), "Every electrical component in the circuit must have a node +") ;
                    BELFEM_ERROR(tComponentSection->key_exists("node -"), "Every electrical component in the circuit must have a node -") ;

                    const index_t tNodePlus = tNetlistFactory.node_index(
                            tComponentSection->get_string("node +") );
                    const index_t tNodeMinus = tNetlistFactory.node_index(
                            tComponentSection->get_string("node -") );

                    this->read_terminal_pair( tComponentSection, i,
                                              tNodePlus, tNodeMinus,
                                              tLabel, tCount );
                }
            }

            this->read_output( aSection, &tNetlistFactory );
        }

//----------------------------------------------------------------------------

        Circuit *
        ElectricalCircuitFactory::circuit()
        {
            // old contract: the factory keeps ownership for its lifetime
            return mCircuit.get() ;
        }

//-----------------------------------------------------------------------------
    }

}
