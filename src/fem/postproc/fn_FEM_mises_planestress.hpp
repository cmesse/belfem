//
// Created by Christian Messe on 18.11.19.
//

#ifndef BELFEM_FN_FEM_MISES_PLANESTRESS_HPP
#define BELFEM_FN_FEM_MISES_PLANESTRESS_HPP

namespace belfem
{
    namespace fem
    {
        class Field;
        class DofManagerBase ;

//------------------------------------------------------------------------------

        void
        mises_planestress( DofManagerBase * aField );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_FEM_MISES_PLANESTRESS_HPP
