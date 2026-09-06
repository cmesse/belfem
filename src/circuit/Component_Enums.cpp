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

#include "stringtools.hpp"
#include "assert.hpp"
#include "Component_Enums.hpp"

namespace belfem
{
    namespace electronics
    {
//------------------------------------------------------------------------------

        std::string
        to_string( const ComponentType aComponentType )
        {
            switch ( aComponentType )
            {
                case ( ComponentType::RESISTOR ) :
                {
                    return "resistor";
                }
                case ( ComponentType::CAPACITOR ) :
                {
                    return "capacitor";
                }
                case ( ComponentType::INDUCTOR ) :
                {
                    return "inductor";
                }
                case ( ComponentType::SWITCH ) :
                {
                    return "switch";
                }
                case ( ComponentType::VOLTAGESOURCE ) :
                {
                    return "voltage source";
                }
                case ( ComponentType::CURRENTSOURCE ) :
                {
                    return "current source";
                }
                case ( ComponentType::DIODE ) :
                {
                    return "diode";
                }
                case ( ComponentType::SUPERCONDUCTOR ) :
                {
                    return "superconductor";
                }
                case ( ComponentType::TERMINALPAIR ) :
                {
                    return "terminal pair";
                }
                default:
                {
                    return "undefined";
                }
            }
        }

//------------------------------------------------------------------------------

        ComponentType
        component_type( const string & aString )
        {
            string tString = string_to_lower( aString );

            if ( tString == "resistor" )
            {
                return ComponentType::RESISTOR;
            }
            else if ( tString == "capacitor" )
            {
                return ComponentType::CAPACITOR;
            }
            else if ( tString == "inductor" )
            {
                return ComponentType::INDUCTOR;
            }
            else if ( tString == "switch" )
            {
                return ComponentType::SWITCH;
            }
            else if ( tString == "voltage source" )
            {
                return ComponentType::VOLTAGESOURCE;
            }
            else if ( tString == "current source" )
            {
                return ComponentType::CURRENTSOURCE;
            }
            else if ( tString == "diode" )
            {
                return ComponentType::DIODE;
            }
            else if ( tString == "superconductor" )
            {
                return ComponentType::SUPERCONDUCTOR;
            }
            else if ( tString == "terminal pair" )
            {
                return ComponentType::TERMINALPAIR;
            }
            else
            {
                BELFEM_ERROR( false, "Unknown Electrical Component Type: %s", aString.c_str());
                return ComponentType::UNDEFINED;
            }
        }

//------------------------------------------------------------------------------
    }
}
