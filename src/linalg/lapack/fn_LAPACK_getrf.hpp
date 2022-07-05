//
// Created by Christian Messe on 18.12.20.
//

#ifndef BELFEM_FN_LAPACK_GETRF_HPP
#define BELFEM_FN_LAPACK_GETRF_HPP

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
    sgetrf_( int    * M,
             int    * N,
             float  * A,
             int    * LDA,
             int    * IPIV,
             int    * INFO );

//------------------------------------------------------------------------------

    void
    dgetrf_( int    * M,
             int    * N,
             double * A,
             int    * LDA,
             int    * IPIV,
             int    * INFO );

//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------------------

        template< typename T >
        inline void
        getrf( const int    & M,
               const int    & N,
                     T      * A,
               const int    & LDA,
                     int    * IPIV,
                     int    & INFO )
        {
            BELFEM_ERROR( false, "getri not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getrf( const int    & M,
               const int    & N,
                     float  * A,
               const int    & LDA,
                     int    * IPIV,
                     int    & INFO )
        {
            sgetrf_( const_cast< int * >( &M ),
                     const_cast< int * >( &N ),
                                           A,
                     const_cast< int * >( &LDA ),
                                           IPIV,
                                           &INFO );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getrf( const int    & M,
               const int    & N,
                     double * A,
               const int    & LDA,
                     int    * IPIV,
                     int    & INFO )
        {
            dgetrf_( const_cast< int * >( &M ),
                     const_cast< int * >( &N ),
                                           A,
                     const_cast< int * >( &LDA ),
                                           IPIV,
                                          &INFO );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
} /* end namespace belfem */

#endif //BELFEM_FN_LAPACK_GETRF_HPP
