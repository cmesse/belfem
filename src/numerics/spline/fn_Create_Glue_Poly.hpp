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

#ifndef BELFEM_FN_CREATE_GLUE_POLY_HPP
#define BELFEM_FN_CREATE_GLUE_POLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    /**
     * consider an element such as
     *
     *
     * f(x)
     *  ^
     *  |  (X-dX)  (X)  (X+dX)
     *  . o ---- o ---- o    --> X
     *    0      1      2
     *
     * This function creates a polynomial of the shape
     *
     * f = a*x^4 + b*x^3 + c*x^2 + d*x + e
     *
     * @param aX
     * @param aDeltaX
     * @param aF = { f0, dfdx0, f1, f2, dfdx2 }
     * @param aC = { a, b, c, d, e }
     */
    void
    create_glue_poly(
            const         real   & aX,
            const         real   & aDeltaX,
            const Vector< real > & aF,
                  Vector< real > & aC );

}
#endif //BELFEM_FN_CREATE_GLUE_POLY_HPP
