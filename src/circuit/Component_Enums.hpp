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

#ifndef BELFEM_COMPONENT_ENUMS_HPP
#define BELFEM_COMPONENT_ENUMS_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace electronics
    {
        enum class ComponentType
        {
            RESISTOR = 0,
            CAPACITOR = 1,
            INDUCTOR = 2,
            VOLTAGESOURCE = 3,
            CURRENTSOURCE = 4,
            SWITCH = 5,
            DIODE = 6,
            SUPERCONDUCTOR = 7,
            TERMINALPAIR = 8,
            UNDEFINED = 9
        };

        string
        to_string( const ComponentType aComponentType );

//------------------------------------------------------------------------------

        ComponentType
        component_type( const string & aString );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_COMPONENT_ENUMS_HPP
