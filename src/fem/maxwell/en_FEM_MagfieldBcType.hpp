//
// Created by christian on 10/29/21.
//

#ifndef BELFEM_EN_FEM_MAGFIELDBCTYPE_HPP
#define BELFEM_EN_FEM_MAGFIELDBCTYPE_HPP

namespace belfem
{
    namespace fem
    {
        enum class MagfieldBcType
        {
            Wave          = 0,
            Symmetry      = 1,  // = farfield
            AntiSymmetry  = 2,
            Bearing       = 3,
            UNDEFINED     = 4
        };
    }
}
#endif //BELFEM_EN_FEM_MAGFIELDBCTYPE_HPP
