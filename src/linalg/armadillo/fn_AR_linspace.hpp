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

#ifndef BELFEM_FN_AR_LINSPACE_HPP
#define BELFEM_FN_AR_LINSPACE_HPP

#include "cl_AR_Vector.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    void
    linspace( const T              & aStart,
              const T              & aEnd,
              const belfem::size_t  & aN,
              Vector< T >          & aValues)
    {
        aValues = arma::linspace< arma::Mat< T > >( aStart, aEnd, aN );
    }

//------------------------------------------------------------------------------

    template< typename T >
    auto
    linspace( const T              & aStart,
              const T              & aEnd,
              const belfem::size_t  & aN)
        ->decltype( arma::linspace< arma::Mat< T > >( aStart, aEnd, aN ) )
    {
        return arma::linspace< arma::Mat< T > >( aStart, aEnd, aN );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AR_LINSPACE_HPP
