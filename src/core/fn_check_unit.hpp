//
// Created by christian on 9/29/21.
//

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
