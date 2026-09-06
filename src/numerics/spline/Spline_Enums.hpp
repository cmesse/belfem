/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_SPLINE_ENUMS_HPP
#define BELFEM_SPLINE_ENUMS_HPP

namespace belfem
{
    namespace spline
    {
        enum class SplineBC {
            NoCurvature,  // Natural spline: second derivative = 0 at endpoint
            Parabolic,    // Parabolic runout: first/last polynomial is quadratic (cubic coeff = 0)
            Tangent       // Clamped: prescribed first derivative at endpoint
        };

        enum class ExtraMode
        {
            None,
            Entropy,
            Integral
        };

    }
}
#endif // BELFEM_SPLINE_ENUMS_HPP
