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

#include "fn_GT_fix_capitals.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        string
        fix_capitals( const string & aLabel )
        {
            uint tN = aLabel.length();
            if ( tN == 0 )
            {
                return aLabel;
            }
            else if ( tN == 1 )
            {
                return string_to_upper( aLabel );
            }
            else
            {
                string tResult = string_to_upper( aLabel.substr( 0, 1 ) );
                tResult += string_to_lower( aLabel.substr( 1, tN - 1 ) );
                return tResult;
            }
        }

//------------------------------------------------------------------------------
    }
}