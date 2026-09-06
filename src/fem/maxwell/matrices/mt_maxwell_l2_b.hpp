//
// Created by christian on 12/4/24.
//

#ifndef MT_MAXWELL_L2_B_HPP
#define MT_MAXWELL_L2_B_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            l2_b2d( Calculator * aCalc, Matrix< real > & aK,  Vector< real > & aF );

        }
    }
}

#endif //MT_MAXWELL_L2_B_HPP
