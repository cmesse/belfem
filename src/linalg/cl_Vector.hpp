//
// Created by Christian Messe on 03.09.19.
//

#ifndef BELFEM_CL_VECTOR_HPP
#define BELFEM_CL_VECTOR_HPP

// include implementation
#ifdef BELFEM_ARMADILLO
#include "cl_AR_Vector.hpp"
#include "op_AR_VectorEqualEqual.hpp"
#include "op_AR_VectorElementwiseMultiplication.hpp"
#elif BELFEM_BLAZE
#include "cl_BZ_Vector.hpp"
#include "op_BZ_VectorEqualEqual.hpp"
#include "op_BZ_VectorElementwiseMultiplication.hpp"
#endif

// include operators
#include "op_VectorPlus.hpp"
#include "op_VectorMinus.hpp"
#include "op_VectorTimes.hpp"
#include "op_VectorDivide.hpp"
#endif
