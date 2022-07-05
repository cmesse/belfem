//
// Created by Christian Messe on 18.12.20.
//

#ifndef BELFEM_FN_LAPACK_GEMM_HPP
#define BELFEM_FN_LAPACK_GEMM_HPP

//------------------------------------------------------------------------------
namespace belfem
{
    namespace lapack
    {

#ifdef __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------

            void
            sgemm_(  char   * TRANSA,
                     char   * TRANSB,
                     int    * M,
                     int    * N,
                     int    * K,
                     float  * ALPHA,
                     float  * A,
                     int    * LDA,
                     float  * B,
                     int    * LDB,
                     float  * BETA,
                     float  * C,
                     int    * LDC
                      );

//------------------------------------------------------------------------------

            void
            dgemm_(  char   * TRANSA,
                     char   * TRANSB,
                     int    * M,
                     int    * N,
                     int    * K,
                     double * ALPHA,
                     double * A,
                     int    * LDA,
                     double * B,
                     int    * LDB,
                     double * BETA,
                     double * C,
                     int    * LDC );

//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------------------

        template< typename T >
        inline void
        gemm(    const char   & TRANSA,
                 const char   & TRANSB,
                 const int    & M,
                 const int    & N,
                 const int    & K,
                 const T      & ALPHA,
                 const T      * A,
                 const int    & LDA,
                 const T      * B,
                 const int    & LDB,
                 const T      & BETA,
                       T      * C,
                 const int    & LDC )
        {
            BELFEM_ERROR( false, "gemm not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gemm( const char  & TRANSA,
              const char  & TRANSB,
              const int   & M,
              const int   & N,
              const int   & K,
              const float & ALPHA,
              const float * A,
              const int   & LDA,
              const float * B,
              const int   & LDB,
              const float & BETA,
                    float * C,
              const int   & LDC )
        {
            sgemm_( const_cast< char * >  ( &TRANSA ),
                    const_cast< char * >  ( &TRANSB ),
                    const_cast< int * >   ( &M ),
                    const_cast< int * >   ( &N ),
                    const_cast< int * >   ( &K ),
                    const_cast< float * > ( &ALPHA ),
                    const_cast< float * > ( A ),
                    const_cast< int * >   ( &LDA ),
                    const_cast< float * > ( B ),
                    const_cast< int * >   ( &LDB ),
                    const_cast< float * > ( &BETA ),
                                            C,
                    const_cast< int * >   ( &LDC ) );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gemm( const char   & TRANSA,
              const char   & TRANSB,
              const int    & M,
              const int    & N,
              const int    & K,
              const double & ALPHA,
              const double * A,
              const int    & LDA,
              const double * B,
              const int    & LDB,
              const double & BETA,
                    double * C,
              const int    & LDC )
        {
            dgemm_( const_cast< char * >  ( &TRANSA ),
                    const_cast< char * >  ( &TRANSB ),
                    const_cast< int * >   ( &M ),
                    const_cast< int * >   ( &N ),
                    const_cast< int * >   ( &K ),
                    const_cast< double * >( &ALPHA ),
                    const_cast< double * >( A ),
                    const_cast< int * >   ( &LDA ),
                    const_cast< double * >( B ),
                    const_cast< int * >   ( &LDB ),
                    const_cast< double * >( &BETA ),
                                            C,
                    const_cast< int * >   ( &LDC ) );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
} /* end namespace belfem */
//------------------------------------------------------------------------------
#endif //BELFEM_FN_LAPACK_GEMM_HPP
