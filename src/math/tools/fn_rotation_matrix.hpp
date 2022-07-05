//
// Created by Christian Messe on 26.10.19.
//

#ifndef BELFEM_FN_ROTATION_MATRIX_HPP
#define BELFEM_FN_ROTATION_MATRIX_HPP

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * rotate around an axis with an angle
     *
     * @tparam T
     * @param aAxis
     * @param aAngle
     * @param aMatrix
     */
    template < typename T >
    void
    rotation_matrix(
            const Vector< T >    & aAxis,
            const T              & aAngle,
                  Matrix< real > & aMatrix )
    {
        BELFEM_ASSERT( aAxis.length() == 3,
                      "Input vector must be of length 3." );

        BELFEM_ASSERT( aMatrix.n_rows() == 3 && aMatrix.n_cols() == 3,
            "Input matrix must be a 3x3 matrix" );

        T c = std::cos( aAngle );
        T s = std::sin( aAngle );
        T d = 1.0 - c;

        T * tData = aMatrix.data() ;

#ifdef BELFEM_ARMADILLO
        tData[ 0 ] = aAxis( 0 ) * aAxis ( 0 ) * d + c;
        tData[ 1 ] = aAxis( 1 ) * aAxis ( 0 ) * d + aAxis( 2 ) * s;
        tData[ 2 ] = aAxis( 2 ) * aAxis ( 0 ) * d - aAxis( 1 ) * s;

        tData[ 3 ] = aAxis( 0 ) * aAxis ( 1 ) * d - aAxis( 2 ) * s;
        tData[ 4 ] = aAxis( 1 ) * aAxis ( 1 ) * d + c;
        tData[ 5 ] = aAxis( 2 ) * aAxis ( 1 ) * d + aAxis( 0 ) * s;

        tData[ 6 ] = aAxis( 0 ) * aAxis ( 2 ) * d + aAxis( 1 ) * s;
        tData[ 7 ] = aAxis( 1 ) * aAxis ( 2 ) * d - aAxis( 0 ) * s;
        tData[ 8 ] = aAxis( 2 ) * aAxis ( 2 ) * d + c;
#elif  BELFEM_BLAZE
        tData[  0 ] = aAxis( 0 ) * aAxis ( 0 ) * d + c;
        tData[  1 ] = aAxis( 1 ) * aAxis ( 0 ) * d + aAxis( 2 ) * s;
        tData[  2 ] = aAxis( 2 ) * aAxis ( 0 ) * d - aAxis( 1 ) * s;

        tData[  4 ] = aAxis( 0 ) * aAxis ( 1 ) * d - aAxis( 2 ) * s;
        tData[  5 ] = aAxis( 1 ) * aAxis ( 1 ) * d + c;
        tData[  6 ] = aAxis( 2 ) * aAxis ( 1 ) * d + aAxis( 0 ) * s;

        tData[  8 ] = aAxis( 0 ) * aAxis ( 2 ) * d + aAxis( 1 ) * s;
        tData[  9 ] = aAxis( 1 ) * aAxis ( 2 ) * d - aAxis( 0 ) * s;
        tData[ 10 ] = aAxis( 2 ) * aAxis ( 2 ) * d + c;
#endif

    }

//------------------------------------------------------------------------------

    template < typename T >
    void
    rotation_matrix_strip(
            const Vector< T >    & aAxis,
            const T              & aAngle,
            Matrix< real > & aMatrix )
    {
        BELFEM_ASSERT( aAxis.length() == 3,
                      "Input vector must be of length 3." );

        BELFEM_ASSERT( aMatrix.n_rows() == 2 && aMatrix.n_cols() == 3,
                      "Input matrix must be a 2x3 matrix" );

        T c = std::cos( aAngle );
        T s = std::sin( aAngle );
        T d = 1.0 - c;

        aMatrix( 0, 0 ) = aAxis( 0 ) * aAxis ( 0 ) * d + c;
        aMatrix( 1, 0 ) = aAxis( 1 ) * aAxis ( 0 ) * d + aAxis( 2 ) * s;
        //aMatrix( 2, 0 ) = aAxis( 2 ) * aAxis ( 0 ) * d - aAxis( 1 ) * s;

        aMatrix( 0, 1 ) = aAxis( 0 ) * aAxis ( 1 ) * d - aAxis( 2 ) * s;
        aMatrix( 1, 1 ) = aAxis( 1 ) * aAxis ( 1 ) * d + c;
        //aMatrix( 2, 1 ) = aAxis( 2 ) * aAxis ( 1 ) * d + aAxis( 0 ) * s;

        aMatrix( 0, 2 ) = aAxis( 0 ) * aAxis ( 2 ) * d + aAxis( 1 ) * s;
        aMatrix( 1, 2 ) = aAxis( 1 ) * aAxis ( 2 ) * d - aAxis( 0 ) * s;
        //aMatrix( 2, 2 ) = aAxis( 2 ) * aAxis ( 2 ) * d + c;
    }

//------------------------------------------------------------------------------

    /**
     * rotate using euler angles
     * project from body to global
     * according to DIN 9300
     *
     * @param aYaw
     * @param aPitch
     * @param aRoll
     * @param aMatrix
     */
    template < typename T >
    void
    rotation_matrix(
            const T        & aYaw,
            const T        & aPitch,
            const T        & aRoll,
            Matrix< real > & aMatrix )
    {
        BELFEM_ASSERT( aMatrix.n_rows() == 3 && aMatrix.n_cols() == 3,
                      "Input matrix must be a 3x3 matrix" );

        const real sa = std::sin( aRoll );
        const real ca = std::cos( aRoll );

        const real sb = std::sin( aPitch );
        const real cb = std::cos( aPitch );

        const real sc = std::sin( aYaw );
        const real cc = std::cos( aYaw );

        T * tData = aMatrix.data() ;
#ifdef BELFEM_ARMADILLO
        tData[ 0 ] = ca*cb ;
        tData[ 1 ] = cb*sa ;
        tData[ 2 ] = -sb ;

        tData[ 3 ] = ca*sb*sc - cc*sa ;
        tData[ 4 ] = ca*cc + sa*sb*sc ;
        tData[ 5 ] = cb*sc ;

        tData[ 6 ] = sa*sc + ca*cc*sb ;
        tData[ 7 ] = cc*sa*sb - ca*sc ;
        tData[ 8 ] = cb*cc;
#elif  BELFEM_BLAZE
        tData[  0 ] = ca*cb ;
        tData[  1 ] = cb*sa ;
        tData[  2 ] = -sb ;

        tData[  4 ] = ca*sb*sc - cc*sa ;
        tData[  5 ] = ca*cc + sa*sb*sc ;
        tData[  6 ] = cb*sc ;

        tData[  8 ] = sa*sc + ca*cc*sb ;
        tData[  9 ] = cc*sa*sb - ca*sc ;
        tData[ 10 ] = cb*cc;
#endif

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_ROTATION_MATRIX_HPP
