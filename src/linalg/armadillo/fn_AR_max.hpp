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

#ifndef BELFEM_FN_AR_MAX_HPP
#define BELFEM_FN_AR_MAX_HPP

#include "cl_AR_Vector.hpp"
#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template <typename T >
    T
    max( const Vector<T> & aVector )
    {
        return aVector.vector_data().max();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    T
    max( const Matrix<T> & aMatrix )
    {
        return aMatrix.matrix_data().max();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <typename T >
    T
    max( const arma::Mat <T> & aMatrix )
    {
        return aMatrix.max();
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_MAX_HPP
