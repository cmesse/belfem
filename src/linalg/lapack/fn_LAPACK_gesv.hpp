//
// Created by Christian Messe on 18.12.20.
//

#ifndef BELFEM_FN_LAPACK_GESV_HPP
#define BELFEM_FN_LAPACK_GESV_HPP

#include "assert.hpp"

//------------------------------------------------------------------------------
namespace belfem
{
    namespace lapack
    {
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------

        template< typename T >
        inline void
        gesv(   const int    & N,
                const int    & NRHS,
                      T      * A,
                const int    & LDA,
                      int    * IPIV,
                      T      * B,
                const int    & LDB,
                      int    & INFO )
        {
            BELFEM_ERROR( false, "gesv not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gesv(   const int    & N,
                const int    & NRHS,
                      float  * A,
                const int    & LDA,
                      int    * IPIV,
                      float  * B,
                const int    & LDB,
                      int    & INFO )
        {
            sgesv_( const_cast< int * >( &N ),
                    const_cast< int * >( &NRHS ),
                                          A,
                    const_cast< int * >( &LDA ),
                                          IPIV,
                                          B,
                    const_cast< int * >( &LDB ),
                                         &INFO );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gesv(   const int    & N,
                const int    & NRHS,
                      double * A,
                const int    & LDA,
                      int    * IPIV,
                      double * B,
                const int    & LDB,
                      int    & INFO )
        {
            dgesv_( const_cast< int * >( &N ),
                    const_cast< int * >( &NRHS ),
                                          A,
                    const_cast< int * >( &LDA ),
                                          IPIV,
                                          B,
                    const_cast< int * >( &LDB ),
                                         &INFO );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
} /* end namespace belfem */

#endif //BELFEM_FN_LAPACK_GESV_HPP
