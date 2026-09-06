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


#ifndef BELFEM_OP_AR_MATRIXEQUALEQUAL_HPP
#define BELFEM_OP_AR_MATRIXEQUALEQUAL_HPP

#include "cl_AR_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    template< typename T >
    bool
    operator==( const Matrix <T> & aA,
                const Matrix <T> & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                    arma::accu( aA.matrix_data() == aB.matrix_data() ) )
                   ==  ( std::size_t ) aA.n_rows()*aA.n_cols() ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    bool
    operator==( const T & aA,
                const Matrix <T> & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                    arma::accu( aA == aB.matrix_data() ) )
                   ==  ( std::size_t ) aA.n_rows()*aA.n_cols() ;
        return aA == aB.matrix_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    bool
    operator==( const Matrix <T> & aA,
                const T & aB )
    {
        return ( std::size_t ) arma::as_scalar(
                    arma::accu( aA.matrix_data() == aB ) )
                   ==  ( std::size_t ) aA.n_rows()*aA.n_cols() ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_OP_AR_MATRIXEQUALEQUAL_HPP
