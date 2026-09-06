//
// Created by grgia@ge.polymtl.ca on 22/03/25.
//

#ifndef BELFEM_MT_MAXWELL_BACKGROUND_HPP
#define BELFEM_MT_MAXWELL_BACKGROUND_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_TimestepMatrices.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            background_phi( Calculator * aCalc, TimestepMatrices * aMatrices );

        }
    }
}

#endif //BELFEM_MT_MAXWELL_BACKGROUND_HPP
