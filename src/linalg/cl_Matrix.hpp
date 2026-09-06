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

#ifndef BELFEM_CL_MATRIX_HPP
#define BELFEM_CL_MATRIX_HPP

#include "cl_Vector.hpp"

// include implementation
#ifdef BELFEM_ARMADILLO
#include "cl_AR_Matrix.hpp"
#include "op_AR_MatrixEqualEqual.hpp"
#elif BELFEM_BLAZE
#include "cl_BZ_Matrix.hpp"
#include "op_BZ_MatrixEqualEqual.hpp"
#endif

#include "op_MatrixPlus.hpp"
#include "op_MatrixMinus.hpp"
#include "op_MatrixTimes.hpp"
#endif //BELFEM_CL_MATRIX_HPP
