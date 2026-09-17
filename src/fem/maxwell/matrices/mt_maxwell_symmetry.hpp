/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

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
