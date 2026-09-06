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

#ifndef BELFEM_FN_ROTATE_HPP
#define BELFEM_FN_ROTATE_HPP

#include "cl_Matrix.hpp"
#include "fn_TR_rotate42.hpp"
#include "cl_Tensor.hpp"


namespace belfem
{
//----------------------------------------------------------------------------

    template < typename T >
    void
    rotate(  const Tensor< T > & aB, const Matrix< T > & aR, Tensor< T > & aA )
    {
        BELFEM_ASSERT( aB.is_3333(), "Tensor A must be of type 3x3x3x3" );
        BELFEM_ASSERT( aR.n_rows() == 3 && aR.n_cols() == 3,
                      "Matrix B must be of size 3x3" );
        BELFEM_ASSERT( aA.is_3333(), "Tensor C must be of type 3x3x3x3" );

        tensor::rotate42( aB.data(), aR, aA.data() );
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_FN_ROTATE_HPP
