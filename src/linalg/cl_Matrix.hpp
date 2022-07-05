//
// Created by Christian Messe on 03.09.19.
//

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
