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

#ifndef BELFEM_EN_FEM_BOUNDARYCONDITIONTYPE_HPP
#define BELFEM_EN_FEM_BOUNDARYCONDITIONTYPE_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace fem
    {
//-----------------------------------------------------------------------------

        enum class BoundaryConditionType
        {
            Neumann,
            Dirichlet,
            Bearing,
            Gauge,


            Current, //Maxwell Specific
            Voltage, //Maxwell Specific
            CircuitCurrent, //Maxwell Specific
            CircuitVoltage, //Maxwell Specific
            Background, //Maxwell Specific, impose a background field weakly
            BackgroundDirichlet, //Maxwell Specific, impose a background field through phi


            UNDEFINED
        };

//------------------------------------------------------------------------------

        string
        to_string( const BoundaryConditionType aBoundaryConditiontype );

//------------------------------------------------------------------------------

        BoundaryConditionType
        boundary_condition_type( const string & aString );

//------------------------------------------------------------------------------

        /**
         * true if a condition of this type publishes its imposed value as a
         * mesh global named after its deck block. False for Bearing, which
         * imposes a bare node constraint with no evaluated scalar, and for
         * the four terminal types: their value AND its response are written
         * by Controller::save_IV as one I/U pair per condition, so a block
         * global would duplicate the I half under a second name
         */
        bool
        has_block_global( const BoundaryConditionType aBoundaryConditionType );

//------------------------------------------------------------------------------

    }
}

#endif //BELFEM_EN_FEM_BOUNDARYCONDITIONTYPE_HPP
