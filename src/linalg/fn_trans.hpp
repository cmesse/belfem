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

/**
 * @file
 * @brief Matrix transpose.
 * @ingroup grp_linalg
 *
 * This header dispatches directly, because Armadillo and Blaze expose the operation
 * under different names (`arma::strans` and `blaze::trans`). The overload accepts a
 * generic expression type, so a caller can transpose an unevaluated backend expression
 * as well as a materialised Matrix.
 *
 * The result is the plain transpose, never the conjugate transpose. For the real matrices
 * BELFEM uses almost everywhere the two coincide; for the complex `Matrix<T>`
 * instantiations the LAPACK wrappers accept, `trans` does not conjugate.
 */

#ifndef BELFEM_FN_TRANS_HPP
#define BELFEM_FN_TRANS_HPP

#ifdef BELFEM_ARMADILLO
#include "armadillo.hpp"
#elif  BELFEM_BLAZE
#include <blaze/math/functors/Trans.h>
#endif

#include "cl_Matrix.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

#ifdef BELFEM_ARMADILLO
    template < typename ET >
    auto
    trans( ET & aMatrix )
        -> decltype( arma::strans( aMatrix ) )
    {
        return arma::strans( aMatrix );
    }

#elif  BELFEM_BLAZE
    template < typename ET >
    auto
    trans( ET & aMatrix )
        -> decltype( blaze::trans( aMatrix ) )
    {
        return blaze::trans( aMatrix );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename ET >
    auto
    trans( const ET & aMatrix )
        -> decltype( blaze::trans( aMatrix ) )
    {
        return blaze::trans( aMatrix );
    }

#endif
//------------------------------------------------------------------------------

    template < typename T >
    auto
    trans( Matrix < T > & aMatrix )
        -> decltype( trans( aMatrix.matrix_data() ) )
    {
        return  trans( aMatrix.matrix_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    auto
    trans( const Matrix < T > & aMatrix )
        -> decltype( trans( aMatrix.matrix_data() ) )
    {
        return  trans( aMatrix.matrix_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_TRANS_HPP
