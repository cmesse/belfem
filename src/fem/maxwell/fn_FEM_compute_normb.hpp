//
// Created by christian on 11/1/21.
//

#ifndef BELFEM_FN_FEM_COMPUTE_NORMB_HPP
#define BELFEM_FN_FEM_COMPUTE_NORMB_HPP

#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
        void
        compute_normb( DofManager * aField, const bool aBiotSavartFlag=false );

    }
}
#endif //BELFEM_FN_FEM_COMPUTE_NORMB_HPP
