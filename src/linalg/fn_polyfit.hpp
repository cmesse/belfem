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
 * @brief Least-squares polynomial fit.
 * @ingroup grp_linalg
 *
 * Preconditions the caller must meet: @p aX and @p aY have the same length, that length
 * exceeds the requested degree, and **the first and last entry of @p aX differ**. Both
 * backends scale the abscissa by `( n - 1 ) / ( aX(n-1) - aX(0) )` to improve the
 * conditioning of the fit and then scale the coefficients back, so equal endpoints are a
 * division by zero rather than a diagnostic.
 *
 * The coefficient vector is resized to degree + 1.
 *
 * **Coefficient order is descending**: `aCoeffs(0)` is the coefficient of the highest
 * power and the last entry is the constant term, so a vector of length n+1 describes a
 * polynomial of degree n. This is the convention Armadillo and MATLAB use, and it is the
 * opposite of the ascending order some other libraries take. Getting it backwards
 * produces a plausible-looking wrong answer rather than an error.
 */

#ifndef BELFEM_FN_POLYFIT_HPP
#define BELFEM_FN_POLYFIT_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_polyfit.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_polyfit.hpp"
#endif

#endif //BELFEM_FN_POLYFIT_HPP
