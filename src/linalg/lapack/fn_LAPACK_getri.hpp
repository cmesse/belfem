//
// Created by Christian Messe on 18.12.20.
//

#ifndef BELFEM_FN_LAPACK_GETRI_HPP
#define BELFEM_FN_LAPACK_GETRI_HPP

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
            sgetri_( int    * N,
                     float  * A,
                     int    * LDA,
                     int    * IPIV,
                     float  * WORK,
                     int    * LWORK,
                     int    * INFO );

//------------------------------------------------------------------------------

            void
            dgetri_( int    * N,
                     double * A,
                     int    * LDA,
                     int    * IPIV,
                     double * WORK,
                     int    * LWORK,
                     int    * INFO );

//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------------------

        template< typename T >
        inline void
        getri(  const int     & N,
                      T       * A,
                const int     & LDA,
                      int     * IPIV,
                      T       * WORK,
                const int     & LWORK,
                      int     & INFO )
        {
            BELFEM_ERROR( false, "getri not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getri(  const int     & N,
                      float   * A,
                const int     & LDA,
                      int     * IPIV,
                      float   * WORK,
                const int     & LWORK,
                      int     & INFO )
        {
            sgetri_( const_cast< int * > ( &N ),
                                            A,
                     const_cast< int * > ( &LDA ),
                                            IPIV,
                                            WORK,
                     const_cast< int * > ( &LWORK ),
                                            &INFO );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getri(  const int     & N,
                      double  * A,
                const int     & LDA,
                      int     * IPIV,
                      double  * WORK,
                const int     & LWORK,
                      int     & INFO )
        {
            dgetri_( const_cast< int * > ( &N ),
                                            A,
                     const_cast< int * > ( &LDA ),
                                            IPIV,
                                            WORK,
                     const_cast< int * > ( &LWORK ),
                                           &INFO );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
} /* end namespace belfem */
#endif //BELFEM_FN_LAPACK_GETRI_HPP
