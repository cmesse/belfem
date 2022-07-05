//
// Created by Christian Messe on 23.01.21.
//

#ifndef BELFEM_CL_MATERIAL_INCONEL750X_HPP
#define BELFEM_CL_MATERIAL_INCONEL750X_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class Inconel750X : public IsotropicMaterial
        {
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Inconel750X() ;

//----------------------------------------------------------------------------

            ~Inconel750X() = default ;

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_MATERIAL_INCONEL750X_HPP
