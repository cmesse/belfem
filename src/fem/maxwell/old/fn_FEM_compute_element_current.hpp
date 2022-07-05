//
// Created by christian on 10/14/21.
//

#ifndef BELFEM_FN_FEM_COMPUTE_ELEMENT_CURRENT_HPP
#define BELFEM_FN_FEM_COMPUTE_ELEMENT_CURRENT_HPP
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
//----------------------------------------------------------------------------

        void
        compute_element_current( DofManager * aField );

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_FN_FEM_COMPUTE_ELEMENT_CURRENT_HPP