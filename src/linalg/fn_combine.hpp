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

/**
 * @file
 * @brief Concatenates two, three or four vectors into one.
 * @ingroup grp_linalg
 *
 * The output is the **last** argument. It is resized to the combined length and receives
 * the inputs in argument order. It must not alias any of them.
 *
 */

#ifndef FN_COMBINE_HPP
#define FN_COMBINE_HPP

#include <cstring> // for std::memcpy
#include "cl_Vector.hpp"

namespace belfem
{
    template< typename T >
    void
    combine( const Vector< T > & aA, const Vector< T > & aB, Vector< T > & aC   )
    {
#ifdef BELFEM_ARMADILLO
        aC.vector_data() = arma::join_cols( aA.vector_data() , aB.vector_data() );
#else
        aC.set_size( aA.length() + aB.length() );

        index_t tCount = 0 ;

        for ( T tVal : aA )
        {
            aC( tCount++ ) = tVal ;
        }
        for ( T tVal : aB )
        {
            aC( tCount++ ) = tVal;
        }
#endif
    }

    template< typename T >
void
combine( const Vector< T > & aA, const Vector< T > & aB, const Vector< T > & aC, Vector< T > & aD   )
    {
#ifdef BELFEM_ARMADILLO
        aD.vector_data() = arma::join_cols( arma::join_cols( aA.vector_data(), aB.vector_data() ), aC.vector_data() );
#else
        aD.set_size( aA.length() + aB.length() + aC.length() );
        index_t tCount = 0 ;
        for ( T tVal : aA )
        {
            aD( tCount++ ) = tVal ;
        }
        for ( T tVal : aB )
        {
            aD( tCount++ ) = tVal;
        }
        for ( T tVal : aC )
        {
            aD( tCount++ ) = tVal;
        }
#endif
    }

   template< typename T >
   void
   combine( const Vector< T > & aA, const Vector< T > & aB, const Vector< T > & aC, const Vector< T > & aD, Vector< T > & aE )
    {
#ifdef BELFEM_ARMADILLO
        aE.vector_data() = arma::join_cols( arma::join_cols( arma::join_cols( aA.vector_data(), aB.vector_data() ), aC.vector_data() ), aD.vector_data() );
#else
        aE.set_size( aA.length() + aB.length() + aC.length() + aD.length() );
        index_t tCount = 0 ;
        for ( T tVal : aA )
        {
            aE( tCount++ ) = tVal ;
        }
        for ( T tVal : aB )
        {
            aE( tCount++ ) = tVal;
        }
        for ( T tVal : aC )
        {
            aE( tCount++ ) = tVal;
        }
        for ( T tVal : aD )
        {
            aE( tCount++ ) = tVal;
        }
#endif
    }
}

#endif //FN_COMBINE_HPP
