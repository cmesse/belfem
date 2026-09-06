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


#ifndef BELFEM_FN_CHECK_UNIT_HPP
#define BELFEM_FN_CHECK_UNIT_HPP

#include "typedefs.hpp"
#include "stringtools.hpp"

namespace belfem
{
    /**
     * check if a value has the correct unit
     */
    inline bool
    check_unit( const value & aValue, const string & aUnit )
    {
        return aValue.second == unit_to_si( aUnit ).second ;
    }
}
#endif //BELFEM_FN_CHECK_UNIT_HPP
