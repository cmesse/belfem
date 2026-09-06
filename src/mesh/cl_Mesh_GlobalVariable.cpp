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

#include "cl_Mesh_GlobalVariable.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        GlobalVariable::GlobalVariable(
                const string & aLabel,
                const id_t   & aID,
                const real     aValue ) :
                mLabel( aLabel ),
                mID( aID ),
                mValue( aValue )
        {
        }

//------------------------------------------------------------------------------
    }
}