//
// Created by Christian Messe on 02.09.19.
//

#ifndef BELFEM_FN_GT_IS_NOBLE
#define BELFEM_FN_GT_IS_NOBLE

#include "typedefs.hpp"

namespace belfem
{
    namespace gastables
    {
        bool
        is_noble( const string & aLabel )
        {
            return aLabel == "He" || aLabel == "Ne" ||
                   aLabel == "Ar" || aLabel == "Kr" ||
                   aLabel == "Xe" || aLabel == "Rn";
        }
    }
}
#endif // BELFEM_FN_GT_IS_NOBLE
