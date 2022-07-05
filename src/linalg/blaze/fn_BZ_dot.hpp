//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_DOT_HPP
#define BELFEM_FN_BZ_DOT_HPP

#include "cl_BZ_Vector.hpp"
#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    auto
    dot( const Vector< T > & aA , const Vector< T > & aB )
        -> decltype( blaze::dot( aA.vector_data() , aB.vector_data() ) )
    {
        return blaze::dot( aA.vector_data() , aB.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    dot( const Matrix< T > & aA , const Vector< T > & aB )
    -> decltype( blaze::dot( aA.matrix_data() , aB.vector_data() ) )
    {
        return blaze::dot( aA.matrix_data() , aB.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    auto
    dot( const Vector< T > & aA , const Matrix< T > & aB )
    -> decltype( blaze::dot( aA.vector_data() , aB.matrix_data() ) )
    {
        return blaze::dot( aA.vector_data() , aB.matrix_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET >
    auto
    dot( const ET & aA, const ET & aB )
        -> decltype( blaze::dot( aA , aB ) )
    {
        return blaze::dot( aA , aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    T
    dot( const blaze::Rows<blaze::DynamicMatrix<T, true>, false, true, false> & aA,
         const blaze::Column<blaze::DynamicMatrix<T, true>, true, true, false> & aB )
    {
        auto aValue = aA * aB ;
        return aValue[ 0 ];
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T >
    T
    dot( const blaze::Rows<blaze::DynamicMatrix<T, true>, false, true, false> & aA,
            const Vector< T > & aB )
    {
        auto aValue = aA * aB.vector_data() ;
        return aValue[ 0 ];
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T, typename ET >
    auto
    dot( Vector< T > & aA,  const ET & aB)
        ->decltype( blaze::dot( aA.vector_data(), aB) )
    {
        return blaze::dot( aA.vector_data(), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename T >
    auto
    dot( const ET & aA,  const Vector< T > & aB)
        ->decltype( blaze::dot( aA, aB.vector_data() ) )
    {
        return blaze::dot(aA, aB.vector_data());
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_DOT_HPP
