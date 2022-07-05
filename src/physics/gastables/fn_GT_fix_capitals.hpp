//
// Created by Christian Messe on 02.09.19.
//

#ifndef BELFEM_FN_GT_FIX_CAPITALS_HPP
#define BELFEM_FN_GT_FIX_CAPITALS_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        /**
         * fix capitalization of components, eg XE --> Xe
         */
        string
        fix_capitals(  const string & aLabel  );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_GT_FIX_CAPITALS_HPP
