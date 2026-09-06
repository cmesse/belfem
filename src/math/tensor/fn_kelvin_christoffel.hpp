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

#ifndef BELFEM_FN_KELVIN_CHRISTOFFEL_HPP
#define BELFEM_FN_KELVIN_CHRISTOFFEL_HPP


#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Tensor.hpp"
#include "fn_TR_kelvin_christoffel.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    template < typename T >
    inline void
    kelvin_christoffel(  const Tensor< T > & aA, const Vector< T > & aB, Matrix< T > & aC )
    {
        BELFEM_ASSERT( aA.is_3333(),
                      "operating tensor must be 3x3x3x3" );

        BELFEM_ASSERT( aB.length() == 3,
                              "argument vector must be allocated as 3x1" );

        BELFEM_ASSERT(    aC.n_rows() == 3
                      && aC.n_cols() == 3,
                      "target matrix must be allocated as 3x3" );

        // canonical mixed signature ( tensor data pointer, vector/matrix
        // refs ) — the only form both backends provide; under Blaze the
        // matrices are padded, so raw matrix pointers must not be used
        tensor::kelvin_christoffel( aA.data(), aB, aC );
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_FN_KELVIN_CHRISTOFFEL_HPP
