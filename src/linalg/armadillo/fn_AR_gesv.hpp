//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_AR_GESV_HPP
#define BELFEM_FN_AR_GESV_HPP

#include "assert.hpp"
#include "cl_AR_Vector.hpp"
#include "cl_AR_Matrix.hpp"

#ifdef __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------

void
sgesv_( int    * N,
        int    * NRHS,
        float  * A,
        int    * LDA,
        int    * IPIV,
        float  * B,
        int    * LDB,
        int    * INFO );

//------------------------------------------------------------------------------

void
dgesv_( int    * N,
        int    * NRHS,
        double * A,
        int    * LDA,
        int    * IPIV,
        double * B,
        int    * LDB,
        int    * INFO );

//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    inline void
    gesv(
            int    * N,
            int    * NRHS,
            T      * A,
            int    * LDA,
            int    * IPIV,
            T      * B,
            int    * LDB,
            int    * INFO )
    {
        BELFEM_ERROR( false, "gesv not implemented for selected data type" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline void
    gesv(
            int    *N,
            int    * NRHS,
            float  * A,
            int    * LDA,
            int    * IPIV,
            float  * B,
            int    * LDB,
            int    * INFO )
    {
        sgesv_( N, NRHS, A, LDA, IPIV, B, LDB, INFO );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline void
    gesv(
            int    * N,
            int    * NRHS,
            double * A,
            int    * LDA,
            int    * IPIV,
            double * B,
            int    * LDB,
            int    * INFO )
    {
        dgesv_( N, NRHS, A, LDA, IPIV, B, LDB, INFO );
    }

//------------------------------------------------------------------------------

    template< typename T >
    inline void
    gesv( Matrix< T > & aA, Vector< T > & aX, Vector< int > & aP )
    {
        BELFEM_ASSERT( aA.n_rows() == aX.length(),
                      "Number of rows of matrix does not match." );
        BELFEM_ASSERT( aA.n_cols() == aX.length(),
                      "Number of cols of matrix does not match." );
        BELFEM_ASSERT( aP.length() >= aX.length(),
                      "Length of pivot vector does not match." );

        // size of matrix
        int tN = aX.length();

        int tNRHS = 1;

        // error code
        int tStatus = 0;

        // call lapack
        gesv(   &tN,
                &tNRHS,
                aA.data(),
                &tN,
                aP.data(),
                aX.data(),
                &tN,
                &tStatus );

        BELFEM_ERROR( tStatus == 0,
                   "LAPACK gesv has thrown an error: %i", tStatus );
    }

//------------------------------------------------------------------------------

    template< typename T >
    inline void
    gesv( Matrix< T > & aA, Matrix< T > & aX, Vector< int > & aP )
    {
        BELFEM_ASSERT( aA.n_rows() == aX.n_rows(),
                      "Number of rows of matrix does not match." );
        BELFEM_ASSERT( aA.n_cols() == aX.n_rows(),
                      "Number of cols of matrix does not match." );
        BELFEM_ASSERT( aP.length() == aX.n_rows(),
                      "Length of Pivot Vector does not match." );

        // size of matrix
        int tN    = aX.n_rows();

        int tNRHS = aX.n_cols();

        // error code
        int tStatus = 0;

        // call lapack
        gesv(   &tN,
                &tNRHS,
                aA.data(),
                &tN,
                aP.data(),
                aX.data(),
                &tN,
                &tStatus );

        BELFEM_ERROR( tStatus == 0,
                   "LAPACK gesv has thrown an error: %i", tStatus );
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_AR_GESV_HPP
