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

#ifndef BELFEM_MT_THERMAL_H_HPP
#define BELFEM_MT_THERMAL_H_HPP

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
        /**
         * Collapsed thermal Picard kernel: all per-point material math is
         * delegated to aCalc->maxwell() ( thermal-side instance; the Maxwell
         * peer element is linked automatically by Calculator::link ), whose
         * constructor-time dispatch covers the bulk/thin-shell, metal/alloy/
         * HTS, defect, and piecewise axes. Dispatch in
         * IWG_MaxwellThermal::link_to_group() is by domain type only.
         *
         * The Newton companion T_h_newton adds the consistent tangent
         * blocks ( dcp/dT with qhist contraction, the dlambda/dT
         * mixed-operator outer product, and the drho/dT quench feedback —
         * B and beta are frozen per thermal solve, not thermal dofs );
         * derivation in thermal_matrices_cleanup_and_newton_plan.md §3.
         */
        void
        T_h_picard( Calculator * aCalc, TimestepMatrices * aMatrices );

        void
        T_h_newton( Calculator * aCalc, TimestepMatrices * aMatrices );

    }
}

#endif //BELFEM_MT_THERMAL_H_HPP
