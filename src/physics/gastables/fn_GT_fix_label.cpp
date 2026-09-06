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

#include "fn_GT_fix_label.hpp"
#include "stringtools.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        string
        fix_label( const string & aLabel )
        {
            if( aLabel == "C2H3,vinyl" )
            {
                return "C2H3";
            }
            else
            {
                return search_and_replace(
                        search_and_replace( search_and_replace( aLabel,
                        "AR", "Ar" ),
                        "AL", "Al" ),
                        "CL", "Cl" );
            }
        }

//------------------------------------------------------------------------------
    }
}