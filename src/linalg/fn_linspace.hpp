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
 * @brief Evenly spaced values over an interval.
 * @ingroup grp_linalg
 *
 * Produces exactly @p aN values including both endpoints, so @p aN must be at least 2 --
 * the spacing is computed as the interval divided by `aN - 1`.
 *
 */

#ifndef BELFEM_FN_LINSPACE_HPP
#define BELFEM_FN_LINSPACE_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_linspace.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_linspace.hpp"
#endif

#endif //BELFEM_FN_LINSPACE_HPP
