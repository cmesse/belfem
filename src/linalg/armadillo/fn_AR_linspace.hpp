//
// Created by Christian Messe on 04.09.19.
//

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
