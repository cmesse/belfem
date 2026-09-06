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
//   BELFEM wrapper types: Vector<T> and Matrix<T>
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

//------------------------------------------------------------------------------
//   two Blaze expression types (Row, Column, DynamicVector, etc.)
//------------------------------------------------------------------------------

    template< typename ET1, typename ET2 >
    auto
    dot( const ET1 & aA, const ET2 & aB )
        -> decltype( blaze::dot( aA , aB ) )
    {
        return blaze::dot( aA , aB );
    }

//------------------------------------------------------------------------------
//   single-row Rows matrix view: extract row vector via blaze::row()
//------------------------------------------------------------------------------

    template< typename MT, bool SO, bool SF, bool DF, typename ET >
    auto
    dot( const blaze::Rows< MT, SO, SF, DF > & aA, const ET & aB )
        -> decltype( blaze::dot( blaze::row( aA, 0UL ), aB ) )
    {
        return blaze::dot( blaze::row( aA, 0UL ), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename MT, bool SO, bool SF, bool DF >
    auto
    dot( const ET & aA, const blaze::Rows< MT, SO, SF, DF > & aB )
        -> decltype( blaze::dot( aA, blaze::row( aB, 0UL ) ) )
    {
        return blaze::dot( aA, blaze::row( aB, 0UL ) );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename MT1, bool SO1, bool SF1, bool DF1,
              typename MT2, bool SO2, bool SF2, bool DF2 >
    auto
    dot( const blaze::Rows< MT1, SO1, SF1, DF1 > & aA,
         const blaze::Rows< MT2, SO2, SF2, DF2 > & aB )
        -> decltype( blaze::dot( blaze::row( aA, 0UL ), blaze::row( aB, 0UL ) ) )
    {
        return blaze::dot( blaze::row( aA, 0UL ), blaze::row( aB, 0UL ) );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename MT, bool SO, bool SF, bool DF, typename T >
    T
    dot( const blaze::Rows< MT, SO, SF, DF > & aA, const Vector< T > & aB )
    {
        return blaze::dot( blaze::row( aA, 0UL ), aB.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T, typename MT, bool SO, bool SF, bool DF >
    T
    dot( const Vector< T > & aA, const blaze::Rows< MT, SO, SF, DF > & aB )
    {
        return blaze::dot( aA.vector_data(), blaze::row( aB, 0UL ) );
    }

//------------------------------------------------------------------------------
//   single-column Columns matrix view: extract column vector via blaze::column()
//------------------------------------------------------------------------------

    template< typename MT, bool SO, bool SF, bool DF, typename ET >
    auto
    dot( const blaze::Columns< MT, SO, SF, DF > & aA, const ET & aB )
        -> decltype( blaze::dot( blaze::column( aA, 0UL ), aB ) )
    {
        return blaze::dot( blaze::column( aA, 0UL ), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename MT, bool SO, bool SF, bool DF >
    auto
    dot( const ET & aA, const blaze::Columns< MT, SO, SF, DF > & aB )
        -> decltype( blaze::dot( aA, blaze::column( aB, 0UL ) ) )
    {
        return blaze::dot( aA, blaze::column( aB, 0UL ) );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename MT1, bool SO1, bool SF1, bool DF1,
              typename MT2, bool SO2, bool SF2, bool DF2 >
    auto
    dot( const blaze::Columns< MT1, SO1, SF1, DF1 > & aA,
         const blaze::Columns< MT2, SO2, SF2, DF2 > & aB )
        -> decltype( blaze::dot( blaze::column( aA, 0UL ), blaze::column( aB, 0UL ) ) )
    {
        return blaze::dot( blaze::column( aA, 0UL ), blaze::column( aB, 0UL ) );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename MT, bool SO, bool SF, bool DF, typename T >
    T
    dot( const blaze::Columns< MT, SO, SF, DF > & aA, const Vector< T > & aB )
    {
        return blaze::dot( blaze::column( aA, 0UL ), aB.vector_data() );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename T, typename MT, bool SO, bool SF, bool DF >
    T
    dot( const Vector< T > & aA, const blaze::Columns< MT, SO, SF, DF > & aB )
    {
        return blaze::dot( aA.vector_data(), blaze::column( aB, 0UL ) );
    }

//------------------------------------------------------------------------------
//   cross-combinations: Rows × Columns (resolve ambiguity)
//------------------------------------------------------------------------------

    template< typename MT1, bool SO1, bool SF1, bool DF1,
              typename MT2, bool SO2, bool SF2, bool DF2 >
    auto
    dot( const blaze::Rows< MT1, SO1, SF1, DF1 > & aA,
         const blaze::Columns< MT2, SO2, SF2, DF2 > & aB )
        -> decltype( blaze::dot( blaze::row( aA, 0UL ), blaze::column( aB, 0UL ) ) )
    {
        return blaze::dot( blaze::row( aA, 0UL ), blaze::column( aB, 0UL ) );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename MT1, bool SO1, bool SF1, bool DF1,
              typename MT2, bool SO2, bool SF2, bool DF2 >
    auto
    dot( const blaze::Columns< MT1, SO1, SF1, DF1 > & aA,
         const blaze::Rows< MT2, SO2, SF2, DF2 > & aB )
        -> decltype( blaze::dot( blaze::column( aA, 0UL ), blaze::row( aB, 0UL ) ) )
    {
        return blaze::dot( blaze::column( aA, 0UL ), blaze::row( aB, 0UL ) );
    }

//------------------------------------------------------------------------------
//   BELFEM Vector<T> with Blaze expression type
//------------------------------------------------------------------------------

    template< typename T, typename ET >
    auto
    dot( const Vector< T > & aA, const ET & aB )
        -> decltype( blaze::dot( aA.vector_data(), aB ) )
    {
        return blaze::dot( aA.vector_data(), aB );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename ET, typename T >
    auto
    dot( const ET & aA, const Vector< T > & aB )
        -> decltype( blaze::dot( aA, aB.vector_data() ) )
    {
        return blaze::dot( aA, aB.vector_data() );
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BZ_DOT_HPP
