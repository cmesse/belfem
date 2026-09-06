//
// Created by gregorygiard on 10/24/25.
//

#include <variant>
#include "cl_ThermalBoundaryConditionFactory.hpp"
#include "cl_Vector.hpp"
#include "fn_check_unit.hpp"

namespace belfem
{
    namespace fem
    {

//-------------------------------------------------------------------------------

        ThermalBoundaryConditionFactory::ThermalBoundaryConditionFactory()
        {

        }
//-------------------------------------------------------------------------------

        ThermalBoundaryConditionFactory::ThermalBoundaryConditionFactory( const input::Section * aSection )
        {
            BELFEM_ERROR( aSection->type()=="boundary conditions" || aSection->parent()->type()=="boundary conditions" ,"Maxwell Boundary Condition Factory should be defined with a Boundary Condition section" ) ;
            uint n = aSection->num_sections();

            // pre-count global-name bases, cf. the Maxwell twin: repeated
            // bases carry their section ordinal on every occurrence
            for ( uint d=0; d<n; ++d )
            {
                const input::Section * tSection = aSection->section( d ) ;
                if ( boundary_condition_type( tSection->type() )
                    != BoundaryConditionType::Bearing )
                {
                    string tBase = tSection->label().size() > 0 ?
                            tSection->label() : "thermal_" + tSection->type() ;
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
                    case BoundaryConditionType::Dirichlet:
                    {
                        BELFEM_ERROR( tSection->key_exists( "sideset" ) xor tSection->key_exists( "sidesets" ), "Sidesets undefined for boundary condition" ) ;
                        tSection->get_ids( tSection->key_exists( "sideset" ) ? "sideset" : "sidesets",tDomains ) ;
                        mPhysicalBoundaryConditions.push(new PhysicalBoundaryCondition( tCount ));
                        mPhysicalBoundaryConditions( tCount )->set_type(tType) ;
                        mPhysicalBoundaryConditions(tCount++)->set_domains(tDomains) ;
                        break ;
                    }
                    case BoundaryConditionType::Neumann:
                    {
                        BELFEM_ERROR( false, "Neumann boundary conditions not implemented yet for thermal" ) ;
                    }
                    default:
                    {
                        // no condition was pushed for this section, but the
                        // loops below index mPhysicalBoundaryConditions(
                        // tCount - i ) — which would relabel a previous
                        // section's condition, or underflow on the first.
                        // Refuse the type instead of corrupting a neighbor
                        // ( same hardening as the Maxwell twin, 2026-08-11 )
                        BELFEM_ERROR( false,
                            "Unsupported thermal boundary condition type: %s",
                            tSection->type().c_str() );
                        break ;
                    }

                }

                // global-variable name, cf. the Maxwell twin: header label
                // wins, else "thermal_" + type ( the prefix keeps thermal
                // and magnetic globals apart on the SHARED mesh ); a base
                // that repeats carries its section ordinal on every
                // occurrence. Bearing publishes none
                string tGlobalBase ;
                if ( tType != BoundaryConditionType::Bearing )
                {
                    tGlobalBase = tSection->label().size() > 0 ?
                            tSection->label() : "thermal_" + tSection->type() ;

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
                    if ( tGlobalBase.size() > 0 )
                    {
                        mPhysicalBoundaryConditions( tCount - i )->set_label(
                            tNumCondition > 1 ?
                                tGlobalBase + "_" + std::to_string( tNumCondition - i + 1 )
                              : tGlobalBase );
                    }

                    //Set the required units in order to define the function
                    switch ( tType )
                    {
                        case BoundaryConditionType::Gauge :
                        case BoundaryConditionType::Dirichlet :
                        {
                            // a gauge imposes a Dirichlet temperature on its
                            // nodes ( cl_FEM_PhysicalBoundaryCondition.cpp,
                            // Gauge -> bearing()->impose_dirichlet ), so its
                            // amplitude is a temperature like Dirichlet's
                            tUnits = "K" ;
                            break ;
                        }
                        case BoundaryConditionType::Neumann :
                        {
                            // flux density: W/m^2 in 2-D and 3-D alike. The
                            // dimensionality lives in the integration measure,
                            // not in the prescribed value
                            tUnits = "W/m^2" ;
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
                        case BoundaryConditionType::Gauge :
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
                            BELFEM_ERROR( false, "Undefined Boundary Condition for thermal problem" );
                        }

                    }
                    mPhysicalBoundaryConditions( tCount-i )->set_function( tFunction ) ;
                }


            }
        }

//-----------------------------------------------------------------------------

        Cell < PhysicalBoundaryCondition * > &
        ThermalBoundaryConditionFactory::boundary_conditions()
        {
            return mPhysicalBoundaryConditions ;
        }

//-----------------------------------------------------------------------------

        void
        ThermalBoundaryConditionFactory::set_fields( DofManager * aField )
        {
            for ( PhysicalBoundaryCondition * tBC : mPhysicalBoundaryConditions )
            {
                tBC->set_field( aField ) ;
            }
        }

//-------------------------------------------------------------------------------

    }
}