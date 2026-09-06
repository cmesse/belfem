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

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI3_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI3_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        inline void
        gauss_tri3(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            // source : 10.1002/nme.1620070316
            aPoints.set_size( 3, 3 );

            aPoints( 0, 0 ) = 4.0;
            aPoints( 1, 0 ) = 1.0;
            aPoints( 2, 0 ) = 1.0;

            aPoints( 0, 1 ) = 1.0;
            aPoints( 1, 1 ) = 4.0;
            aPoints( 2, 1 ) = 1.0;

            aPoints( 0, 2 ) = 1.0;
            aPoints( 1, 2 ) = 1.0;
            aPoints( 2, 2 ) = 4.0;

            aPoints /= 6.0;

            aWeights.set_size( 3, 1.0 );

            aWeights/= 6.0;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_TRI3_HPP
