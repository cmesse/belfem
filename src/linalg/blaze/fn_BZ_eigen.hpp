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

#ifndef BELFEM_FN_BZ_EIGEN_HPP
#define BELFEM_FN_BZ_EIGEN_HPP
#include <blaze/util/typetraits/IsComplex.h>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_BZ_Vector.hpp"
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    inline int_t
    eigen( const Matrix< real > &  aMatrix,
           Vector< real > & aValues,
           const bool aAbortOnComplex = true )
    {
        BELFEM_ASSERT( aMatrix.n_cols() == aMatrix.n_rows(),
                      "Matrix must be quadratic");

        blaze::DynamicVector<blaze::complex<real>, blaze::columnVector > tValues ;
        blaze::eigen( aMatrix.matrix_data(), tValues );

        size_t tN = aMatrix.n_cols() ;
        aValues.set_size( aMatrix.n_cols() );

        int_t tNumComplex = 0 ;

        for( size_t k=0; k<tN; ++k )
        {
            if( std::abs( std::imag( tValues[ k ] ) ) > BELFEM_EPSILON )
            {
                BELFEM_ERROR( ! aAbortOnComplex,
                    "eigen(): eigenvalue %u of this matrix is complex ( %g %+g i ), and a real "
                    "Vector cannot hold it. Use eigen_sym() if the matrix is symmetric, or pass "
                    "aAbortOnComplex = false to receive a count and NaN entries instead.",
                    ( unsigned int ) k,
                    ( double ) std::real( tValues[ k ] ),
                    ( double ) std::imag( tValues[ k ] ) );

                aValues( k ) = BELFEM_QUIET_NAN ;
                ++tNumComplex ;
            }
            else
            {
                aValues( k ) = std::real( tValues[ k ] );
            }
        }

        return tNumComplex ;
   }

//------------------------------------------------------------------------------

    inline void
    eigen_sym( const Matrix< real > & aMatrix,
               Vector< real > & aValues )
    {
        BELFEM_ASSERT( aMatrix.n_cols() == aMatrix.n_rows(),
                      "Matrix must be quadratic");

        aValues.set_size( aMatrix.n_cols() );

        // blaze::syev overwrites the matrix it is given, so it gets a copy
        Matrix< real > tWork( aMatrix.matrix_data() );

        // 'N' = eigenvalues only. 'U' = read the upper triangle, matching what
        // arma::eig_sym does, so both backends agree even when the caller passes a
        // matrix that is not actually symmetric.
        blaze::syev( tWork.matrix_data(), aValues.vector_data(), 'N', 'U' );
   }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BZ_EIGEN_HPP
