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

#ifndef BELFEM_INTPOINTS_HPP
#define BELFEM_INTPOINTS_HPP

#ifdef __cplusplus
extern"C"
{
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    intpoints_lobatto( int * n, double * w, double * xi );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    intpoints_gauss( int * n, double * w, double * xi );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef __cplusplus
}
#endif

#endif //BELFEM_INTPOINTS_HPP
