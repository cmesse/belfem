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
 * @brief Largest entry of a vector, matrix or column view.
 * @ingroup grp_linalg
 *
 * The input must not be empty.
 *
 */

#ifndef BELFEM_FN_MAX_HPP
#define BELFEM_FN_MAX_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_max.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_max.hpp"
#endif


#endif //BELFEM_FN_MAX_HPP
