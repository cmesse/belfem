//
// Created by christian on 9/20/24.
//

#ifndef MT_MAXWELL_SYMMETRY_HPP
#define MT_MAXWELL_SYMMETRY_HPP

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
            symmetry_phi_2d( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            symmetry_phi_3d( Calculator * aCalc, TimestepMatrices * aMatrices);

            void
            h_symmetry_2d( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            h_symmetry_3d( Calculator * aCalc, TimestepMatrices * aMatrices );

        }
    }
}

#endif //MT_MAXWELL_SYMMETRY_HPP
