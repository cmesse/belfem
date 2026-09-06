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

#include "cl_Cell.hpp"
#include "assert.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace cell
    {
       void
       error_out_of_bounds( const index_t aI, const index_t aN )
        {
            const string tMessage = sprint( "Index %lu out of bounds for cell ( expect < %lu )",
                (long unsigned int) aI, (long unsigned int) aN );
            BELFEM_ERROR( false, tMessage.c_str() );
        }
    }
}

