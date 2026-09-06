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

#ifndef BELFEM_FN_SIGN_HPP
#define BELFEM_FN_SIGN_HPP

namespace belfem
{
    template < typename T >
    inline T sign( const T aX )
    {
        return ( aX > 0) ? 1 : ( ( aX < 0) ? -1 : 0);
    }

}

#endif //BELFEM_FN_SIGN_HPP
