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

#ifndef BELFEM_PARDISOTOOLS_HPP
#define BELFEM_PARDISOTOOLS_HPP

#include "typedefs.hpp"

#ifdef __cplusplus
extern"C" {
#endif

    belfem::int_t
    pardisotools_initialize_parameters( const belfem::int_t * aParameters );

    belfem::int_t
    pardisotools_symbolic_factorization(
            const belfem::int_t  & aN,
            const belfem::int_t  & aNNZ,
            const belfem::int_t  & aNRHS,
            const belfem::int_t  * aPointers,
            const belfem::int_t  * aIndices,
            const double         * aValues );

    belfem::int_t
    pardisotools_solve(  const belfem::int_t    &  aN,
                         const belfem::int_t    & aNNZ,
                         const belfem::int_t    & aNRHS,
                         const belfem::int_t    * aPointers,
                         const belfem::int_t    * aIndices,
                         const double * aValues,
                         double       * aLHS,
                         const double * aRHS,
                         belfem::int_t          * aInfo   );

    belfem::int_t
    pardisotools_free() ;

    double
    pardisotools_get_determinant() ;

#ifdef __cplusplus
}
#endif

#endif //BELFEM_PARDISOTOOLS_HPP
