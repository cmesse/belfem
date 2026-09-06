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

#ifndef BELFEM_FN_INTPOINTS_GAUSS_PYRA1_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_PYRA1_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"


namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------
        inline void
        gauss_pyra1(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {

            aPoints.set_size( 3, 1, 0.0 );

            aPoints( 2, 0 ) = 0.25 ;
            aWeights.set_size( 1, 4.0/3.0 );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_PYRA1_HPP
