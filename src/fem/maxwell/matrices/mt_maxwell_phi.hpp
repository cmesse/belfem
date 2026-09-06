//
// Created by christian on 9/19/24.
//

#ifndef MT_MAXWELL_PHI_HPP
#define MT_MAXWELL_PHI_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_TimestepMatrices.hpp"


namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            phi( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            phi_ferro_picard( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            phi_ferro_newton( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            phi_tri3( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            phi_tet4( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            phi_tri6_tet10( Calculator * aCalc, TimestepMatrices * aMatrices );
        }
    }
}

#endif //MT_MAXWELL_PHI_HPP
