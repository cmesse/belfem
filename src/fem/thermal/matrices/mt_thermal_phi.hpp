//
// Created by gregorygiard on 10/29/25.
//

#ifndef BELFEM_MT_THERMAL_PHI_HPP
#define BELFEM_MT_THERMAL_PHI_HPP

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
        void
        T_phi( Calculator * aCalc, TimestepMatrices * aMatrices );
    }
}

#endif //BELFEM_MT_THERMAL_PHI_HPP