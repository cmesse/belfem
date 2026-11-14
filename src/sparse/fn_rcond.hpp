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

#ifndef BELFEM_FN_RCOND_HPP
#define BELFEM_FN_RCOND_HPP

#include "cl_SpMatrix.hpp"

namespace belfem
{
    /**
     * rough reciprocal condition estimate of a symmetric positive definite
     * matrix, from CHOLMOD's Cholesky factor ( cholmod_rcond ).
     *
     * CSC needs the upper triangle, CSR the lower one; storing both always
     * works. A matrix that is not SPD, not square, not CSC/CSR, or a build
     * with 64-bit int_t is refused with an error. The indexing
     * base is switched to C++ for the call and restored afterwards.
     * Needs SuiteSparse and OpenMP.
     */
    real
    rcond( SpMatrix & aMatrix );
}
#endif //BELFEM_FN_RCOND_HPP
