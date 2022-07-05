//
// Created by Christian Messe on 02.09.19.
//

#ifndef BELFEM_FN_GT_FIX_LABEL_HPP
#define BELFEM_FN_GT_FIX_LABEL_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
         *
         * fix illogical naming in CEA, eg AL --> Al, CL --> Cl
         */
        string
        fix_label( const string & aLabel );

    }
}

#endif //BELFEM_FN_GT_FIX_LABEL_HPP
