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

#include "en_DomainType.hpp"
#include "stringtools.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    string
    to_string( const DomainType aDomainType )
    {
        switch ( aDomainType )
        {
            case DomainType::Conductor :
            {
                return "Conductor";
            }
            case DomainType::Coil :
            {
                return "Coil";
            }
            case DomainType::Ferro :
            {
                return "Ferro";
            }
            case DomainType::Air :
            {
                return "Air";
            }
            case DomainType::Buffer :
            {
                return "Buffer";
            }
            case DomainType::InterfaceCondAir :
            case DomainType::InterfaceFerroAir :
            case DomainType::InterfaceCondFerro :
            {
                return "Interface";
            }
            case DomainType::AirSymmetry :
            case DomainType::BufferSymmetry:
            case DomainType::ConductorSymmetry :
            case DomainType::FerroSymmetry :
            {
                return "Symmetry";
            }
            case DomainType::AirAntiSymmetry :
            case DomainType::BufferAntiSymmetry :
            case DomainType::ConductorAntiSymmetry :
            case DomainType::FerroAntiSymmetry :
            {
                return "AntiSymmetry";
            }
            case DomainType::AirPeriodic :
            case DomainType::BufferPeriodic :
            case DomainType::ConductorPeriodic :
            case DomainType::FerroPeriodic :
            {
                return "Periodic";
            }
            case DomainType::Cut :
            {
                return "Cut";
            }
            case DomainType::BackgroundField :
            {
                return "BackgroundField";
            }
            case DomainType::ThinShell :
            {
                return "Tape";
            }
            case DomainType::Curve :
            {
                return "Curve";
            }
            case DomainType::LeftCoating :
            case DomainType::RightCoating :
            {
                return "Coating";
            }
            default:
            {
                return "Undefined";
            }
        }
    }

//------------------------------------------------------------------------------

    DomainType
    domain_type( const string & aString )
    {
        string tString = string_to_lower( aString );

        if ( tString == "conductor" ||  tString == "superconductor"  )
        {
            return DomainType::Conductor;
        }
        else if ( tString == "buffer")
        {
            return DomainType::Buffer;
        }
        else if ( tString == "coil" )
        {
            return DomainType::Coil;
        }
        else if ( tString == "ferro" || tString == "iron" )
        {
            return DomainType::Ferro;
        }
        else if ( tString == "air" || tString == "vacuum" || tString == "void" )
        {
            return DomainType::Air;
        }
        // NOTE: "cut" is deliberately NOT accepted as an input section type.
        // DomainType::Cut is still live for cuts generated internally from the
        // cohomology, but no parser path ever reads a topology { cut { } }
        // section, so accepting the string here only produced a section that
        // parsed and then did nothing. Let it fail loudly below instead.
        else if ( tString == "background field" )
        {
            return DomainType::BackgroundField;
        }
        else if ( tString == "air symmetry" )
        {
            return DomainType::AirSymmetry;
        }
        else if ( tString == "buffer symmetry" )
        {
            return DomainType::BufferSymmetry;
        }
        else if ( tString == "ferro symmetry" )
        {
            return DomainType::FerroSymmetry;
        }
        else if ( tString == "conductor symmetry" )
        {
            return DomainType::ConductorSymmetry;
        }
        else if ( tString == "air antisymmetry" )
        {
            return DomainType::AirAntiSymmetry;
        }
        else if ( tString == "buffer antisymmetry" )
        {
            return DomainType::BufferAntiSymmetry;
        }
        else if ( tString == "ferro antisymmetry" )
        {
            return DomainType::FerroAntiSymmetry;
        }
        else if ( tString == "conductor antisymmetry" )
        {
            return DomainType::ConductorAntiSymmetry;
        }
        else if ( tString == "thinshell" || tString == "tape" || tString == "shell"  )
        {
            return DomainType::ThinShell;
        }
        else if ( tString == "curve" || tString == "curves"  )
        {
            return DomainType::Curve;
        }
        else if ( tString == "periodic" || tString == "periodic"  )
        {
            return DomainType::Periodic;
        }
        else
        {
            BELFEM_ERROR( false, "Unknown Domain Type: %s", aString.c_str());
            return DomainType::UNDEFINED;
        }
    }
}