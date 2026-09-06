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

#ifndef BELFEM_FN_AR_DOT_HPP
#define BELFEM_FN_AR_DOT_HPP

#include "cl_AR_Vector.hpp"
#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    auto
    dot( const Vector< T > & aA, const Vector< T > & aB )
        -> decltype( arma::dot( aA.vector_data() , aB.vector_data() ) )
    {
        return arma::dot( aA.vector_data() , aB.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    dot( const Matrix< T > & aA, const Vector< T > & aB )
    -> decltype( arma::dot( aA.matrix_data() , aB.vector_data() ) )
    {
        return arma::dot( aA.matrix_data() , aB.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    dot( const Vector< T > & aA, const Matrix< T > & aB )
    -> decltype( arma::dot( aA.vector_data() , aB.matrix_data() ) )
    {
        return arma::dot( aA.vector_data() , aB.matrix_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET >
    auto
    dot( const ET & aA, const ET & aB )
        -> decltype( arma::dot( aA , aB ) )
    {
        return arma::dot( aA , aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T, typename ET >
    auto
     dot( Vector< T > & aA,  const ET & aB)
        ->decltype( arma::dot( aA.vector_data(), aB) )
    {
        return arma::dot( aA.vector_data(), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename T >
    auto
    dot( const ET & aA,  const Vector< T > & aB )
        ->decltype( arma::dot( aA, aB.vector_data() ) )
    {
        return arma::dot( aA, aB.vector_data()  );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_DOT_HPP
