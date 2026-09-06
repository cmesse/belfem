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

#ifndef EN_DOMAINTYPE_HPP
#define EN_DOMAINTYPE_HPP

#include "typedefs.hpp"

namespace belfem
{
    enum class DomainType
    {
        Default = 0,           // default setting for blocks and sidesets

        // the domain types for Coil, Air, Ferro and Conductor must not be changed
        // they control the order in which master and slave priority is given for facets
        // must not be bigger than 255

        Air                    =  1,  // Maxwell Specific
        Buffer                 =  2,  // Maxwell Specific
	    Ferro                  =  3,  // Maxwell Specific
        Coil                   =  4,  // Maxwell Specific
        Conductor              =  5,  // Maxwell Specific
        LeftCoating          =  6,  // Maxwell Specific
        RightCoating         =  7,  // Maxwell Specific
        Cut                    =  8,  // Maxwell Specific

        ThinShell              =  9,  // Maxwell Specific

        Symmetry               = 10,  // Maxwell Specific
        AntiSymmetry           = 11,  // Maxwell Specific
        Periodic               = 12,  // Maxwell Specific

        AirSymmetry            = 13,  // Maxwell Specific
        AirAntiSymmetry        = 14,  // Maxwell Specific
        AirPeriodic            = 15,  // Maxwell Specific

        BufferSymmetry         = 16,  // Maxwell Specific
        BufferAntiSymmetry     = 17,  // Maxwell Specific
        BufferPeriodic         = 18,  // Maxwell Specific

        ConductorSymmetry      = 19,  // Maxwell Specific
        ConductorAntiSymmetry  = 20,  // Maxwell Specific
        ConductorPeriodic      = 21,  // Maxwell Specific

        FerroSymmetry          = 22,  // Maxwell Specific
        FerroAntiSymmetry      = 23,  // Maxwell Specific
        FerroPeriodic          = 24,  // Maxwell Specific

        InterfaceCondAir       = 25,   // Maxwell Specific
        InterfaceCondFerro     = 26,   // Maxwell Specific
        InterfaceFerroAir      = 27,   // Maxwell Specific

        InterfaceAirCoil       = 28,   // Maxwell Specific
        InterfaceFerroCoil     = 29,   // Maxwell Specific
        InterfaceTsCond        = 30,
        InterfaceTsConnector   = 31,
        BackgroundField        = 32,   // Maxwell Specific

        Terminal               = 33,   // Maxwell Specific

        Curve                  = 34,  // for alpha BC on wetted surface

        Inactive               = 35,   // Maxwell Specific
        GeometryOnly           = 36,
        Ghost                  = 37,   // Maxwell Specific
        EnrichedInterface      = 38,   // Maxwell Specific

        Dirichlet              = 39,  // for dot T BC on wetted surface
        Neumann                = 40,  // for dot Q BC on wetted surface
        ThermalAlpha           = 41,  // for alpha BC on wetted surface
        UNDEFINED              = 42
    };

    string
    to_string( const DomainType aDomainType ) ;

    DomainType
    domain_type( const string & aString );
}
#endif //EN_DOMAINTYPE_HPP
