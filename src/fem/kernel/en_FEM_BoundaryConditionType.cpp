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

#include "en_FEM_BoundaryConditionType.hpp"
#include "stringtools.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        string
        to_string( const BoundaryConditionType aBoundaryConditiontype )
        {
            switch (aBoundaryConditiontype)
            {
                case (BoundaryConditionType::Neumann) :
                {
                    return "Neumann";
                }
                case (BoundaryConditionType::Dirichlet) :
                {
                    return "Dirichlet";
                }
                case (BoundaryConditionType::Bearing) :
                {
                    return "Bearing";
                }
                case (BoundaryConditionType::Gauge) :
                {
                    return "Gauge";
                }
                case (BoundaryConditionType::Current) :
                case (BoundaryConditionType::CircuitCurrent) :
                {
                    return "Current";
                }
                case (BoundaryConditionType::Voltage) :
                case (BoundaryConditionType::CircuitVoltage) :
                {
                    return "Voltage";
                }
                case (BoundaryConditionType::Background) :
                {
                    return "Background";
                }
                case (BoundaryConditionType::BackgroundDirichlet) :
                {
                    return "Background";
                }
                default:
                {
                    return "Undefined";
                }
            }
        }

//------------------------------------------------------------------------------

        BoundaryConditionType
        boundary_condition_type( const string & aString )
        {
            string tString = string_to_lower( aString );

            if (tString == "neumann")
            {
                return BoundaryConditionType::Neumann;
            }
            else if (tString == "dirichlet")
            {
                return BoundaryConditionType::Dirichlet;
            }
            else if (tString == "bearing" )
            {
                return BoundaryConditionType::Bearing;
            }
            else if (tString == "gauge")
            {
                return BoundaryConditionType::Gauge;
            }
            else if (tString == "current")
            {
                return BoundaryConditionType::Current;
            }
            else if (tString == "voltage")
            {
                return BoundaryConditionType::Voltage;
            }
            else if (tString == "background")
            {
                return BoundaryConditionType::Background;
            }
            else if (tString == "background dirichlet")
            {
                return BoundaryConditionType::BackgroundDirichlet;
            }
            else
            {
                BELFEM_ERROR( false, "Unknown Boundary Condition Type: %s", aString.c_str());
                return BoundaryConditionType::UNDEFINED;
            }
        }

//-----------------------------------------------------------------------------

        bool
        has_block_global( const BoundaryConditionType aBoundaryConditionType )
        {
            switch ( aBoundaryConditionType )
            {
                case ( BoundaryConditionType::Bearing ) :
                case ( BoundaryConditionType::Current ) :
                case ( BoundaryConditionType::Voltage ) :
                case ( BoundaryConditionType::CircuitCurrent ) :
                case ( BoundaryConditionType::CircuitVoltage ) :
                {
                    return false ;
                }
                default :
                {
                    return true ;
                }
            }
        }

//-----------------------------------------------------------------------------
    }
}
