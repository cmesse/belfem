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

#ifndef CL_QUATERNION_HPP
#define CL_QUATERNION_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>

#include "assert.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    /**
     * @brief Scalar-first quaternion value type for 3D rotations.
     *
     * @ingroup grp_math_quaternion
     * @see @ref math_quaternion_quaternion_usage_guide
     */
    template < typename T >
    class Quaternion
    {
        //! the four components (w, x, y, z), stored inline. The quaternion is
        //! a plain value type: trivially copyable, no heap allocation, safe to
        //! memcpy / MPI-transfer. (Earlier versions optionally aliased external
        //! memory via a pointer; that introduced copy/move hazards for owners
        //! and bought nothing for a four-element value — so it was removed.)
        T mData[ 4 ] ;

    public:
//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------

        Quaternion() :
            mData{ T ( 0 ), T ( 0 ), T ( 0 ), T ( 0 ) }
        {
        }

        Quaternion ( T a, T b, T c, T d ) :
            mData{ a, b, c, d }
        {
        }

        Quaternion ( const Vector < T > & aVector )
        {
            BELFEM_ASSERT( aVector.length() == 3, "Must pass a vector of length 3" );
            mData [ 0 ] = T ( 0 );
            mData [ 1 ] = aVector( 0 );
            mData [ 2 ] = aVector( 1 );
            mData [ 3 ] = aVector( 2 );
        }

        /**
         * @brief Construct a rotation quaternion from axis and angle
         * @param aAxis  The rotation axis (will be normalized internally)
         * @param aAngle The rotation angle in radians
         */
        Quaternion ( const Vector < T > & aAxis, const T aAngle )
        {
            BELFEM_ASSERT( aAxis.length() == 3, "Axis must be a vector of length 3" );

            // normalize the axis
            T tNorm = std::sqrt( aAxis( 0 ) * aAxis( 0 ) + aAxis( 1 ) * aAxis( 1 ) + aAxis( 2 ) * aAxis( 2 ) );
            BELFEM_ASSERT( tNorm > BELFEM_EPSILON, "Axis cannot be zero vector" );

            T tHalfAngle = aAngle * T ( 0.5 );
            T tSinHalf = std::sin( tHalfAngle ) / tNorm;

            mData [ 0 ] = std::cos( tHalfAngle );
            mData [ 1 ] = aAxis( 0 ) * tSinHalf;
            mData [ 2 ] = aAxis( 1 ) * tSinHalf;
            mData [ 3 ] = aAxis( 2 ) * tSinHalf;
        }

        // Copy/move construction, copy/move assignment, and destruction are all
        // trivial (rule of zero): the only data member is a four-element value
        // array. This is what makes Quaternion safe to copy and MPI-transfer.

        static Quaternion identity()
        {
            return { T ( 1 ), T ( 0 ), T ( 0 ), T ( 0 ) };
        }

//------------------------------------------------------------------------------
// Accessors
//------------------------------------------------------------------------------

        T * data()
        {
            return mData;
        }

        const T * data() const
        {
            return mData;
        }

        T & a()
        {
            return mData [ 0 ];
        }

        const T & a() const
        {
            return mData [ 0 ];
        }

        T & b()
        {
            return mData [ 1 ];
        }

        const T & b() const
        {
            return mData [ 1 ];
        }

        T & c()
        {
            return mData [ 2 ];
        }

        const T & c() const
        {
            return mData [ 2 ];
        }

        T & d()
        {
            return mData [ 3 ];
        }

        const T & d() const
        {
            return mData [ 3 ];
        }

//------------------------------------------------------------------------------
// Iterators
//------------------------------------------------------------------------------

        T * begin()
        {
            return mData;
        }

        const T * begin() const
        {
            return mData;
        }

        T * end()
        {
            return mData + 4;
        }

        const T * end() const
        {
            return mData + 4;
        }

//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

        Quaternion & operator= ( std::initializer_list < T > initList )
        {
            BELFEM_ASSERT ( initList.size() == 4, "Initializer list must be of length 4" );
            size_t tCount = 0;
            for ( const T & value : initList )
            {
                mData [ tCount++ ] = value;
            }
            return *this;
        }

        Quaternion & operator+= ( const Quaternion & aOther )
        {
            mData [ 0 ] += aOther.mData [ 0 ];
            mData [ 1 ] += aOther.mData [ 1 ];
            mData [ 2 ] += aOther.mData [ 2 ];
            mData [ 3 ] += aOther.mData [ 3 ];
            return *this;
        }

        Quaternion & operator-= ( const Quaternion & aOther )
        {
            mData [ 0 ] -= aOther.mData [ 0 ];
            mData [ 1 ] -= aOther.mData [ 1 ];
            mData [ 2 ] -= aOther.mData [ 2 ];
            mData [ 3 ] -= aOther.mData [ 3 ];
            return *this;
        }

        Quaternion & operator*= ( T aScalar )
        {
            mData [ 0 ] *= aScalar;
            mData [ 1 ] *= aScalar;
            mData [ 2 ] *= aScalar;
            mData [ 3 ] *= aScalar;
            return *this;
        }

        Quaternion & operator/= ( T aScalar )
        {
            BELFEM_ERROR ( std::abs ( aScalar ) > BELFEM_EPSILON, "Division by zero" );
            mData [ 0 ] /= aScalar;
            mData [ 1 ] /= aScalar;
            mData [ 2 ] /= aScalar;
            mData [ 3 ] /= aScalar;
            return *this;
        }

        Quaternion & operator*= ( const Quaternion & aOther )
        {
            T a1 = mData [ 0 ];
            T b1 = mData [ 1 ];
            T c1 = mData [ 2 ];
            T d1 = mData [ 3 ];

            T a2 = aOther.mData [ 0 ];
            T b2 = aOther.mData [ 1 ];
            T c2 = aOther.mData [ 2 ];
            T d2 = aOther.mData [ 3 ];

            mData [ 0 ] = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2;
            mData [ 1 ] = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2;
            mData [ 2 ] = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2;
            mData [ 3 ] = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2;
            return *this;
        }

//------------------------------------------------------------------------------
// Special functions
//------------------------------------------------------------------------------

        T norm() const
        {
            return std::sqrt ( mData [ 0 ] * mData [ 0 ] + mData [ 1 ] * mData [ 1 ] + mData [ 2 ] * mData [ 2 ] + mData [ 3 ] * mData [ 3 ] );
        }

        Quaternion< T >
        conj() const
        {
            return { mData [ 0 ], -mData [ 1 ], -mData [ 2 ], -mData [ 3 ] };
        }

        Quaternion< T >
        inv() const
        {
            T tSqNorm = mData [ 0 ] * mData [ 0 ] + mData [ 1 ] * mData [ 1 ] + mData [ 2 ] * mData [ 2 ] + mData [ 3 ] * mData [ 3 ];
            BELFEM_ERROR ( tSqNorm > BELFEM_EPSILON, "Zero quaternion has no inverse" );
            return conj() / tSqNorm;
        }

        Quaternion< T > &
        normalize()
        {
            T tNorm = norm();
            BELFEM_ERROR ( tNorm > BELFEM_EPSILON, "Cannot normalize zero quaternion" );
            *this /= tNorm;
            return *this;
        }

        /**
         * @brief Rotate a 3D vector by this unit quaternion
         * @param aVector The vector to rotate (must have length 3)
         * @return The rotated vector
         *
         * Uses the formula: v' = q * v * q* (conjugate)
         * Requires a unit quaternion.
         */
        Vector < T >
        rotate( const Vector < T > & aVector ) const
        {
            BELFEM_ASSERT( std::abs( this->norm() - T( 1 ) ) < T( 100 ) * BELFEM_EPSILON,
                           "rotate() requires a unit quaternion" );
            BELFEM_ASSERT( aVector.length() == 3, "Vector must have length 3" );

            // embed vector as pure quaternion (0, vx, vy, vz)
            Quaternion < T > tV ( T ( 0 ), aVector( 0 ), aVector( 1 ), aVector( 2 ) );

            // compute q * v * q* (unit quaternion, so conjugate = inverse)
            Quaternion < T > tResult = ( *this ) * tV * this->conj();

            // extract the vector part
            Vector < T > tOut( 3 );
            tOut( 0 ) = tResult.b();
            tOut( 1 ) = tResult.c();
            tOut( 2 ) = tResult.d();
            return tOut;
        }
    };

    // Non-member operators
    template < typename T >
    Quaternion < T > operator+ ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        Quaternion < T > aResult = aLHS;
        aResult += aRHS;
        return aResult;
    }

    template < typename T >
    Quaternion < T > operator- ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        Quaternion < T > aResult = aLHS;
        aResult -= aRHS;
        return aResult;
    }

    template < typename T >
    Quaternion < T > operator* ( const Quaternion < T > & aLHS, T aScalar )
    {
        Quaternion < T > aResult = aLHS;
        aResult *= aScalar;
        return aResult;
    }

    template < typename T >
    Quaternion < T > operator* ( T aScalar, const Quaternion < T > & aRHS )
    {
        return aRHS * aScalar;
    }

    template < typename T >
    Quaternion < T > operator/ ( const Quaternion < T > & aLHS, T aScalar )
    {
        Quaternion < T > aResult = aLHS;
        aResult /= aScalar;
        return aResult;
    }

    template < typename T >
    Quaternion < T > operator* ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        Quaternion < T > aResult = aLHS;
        aResult *= aRHS;
        return aResult;
    }

    // Component-wise comparison. Note: for unit quaternions representing
    // rotations, q and -q encode the same rotation. This operator does NOT
    // account for that — it compares quaternions as 4-vectors.
    template < typename T >
    bool operator== ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        if ( std::abs ( aLHS.a() - aRHS.a() ) > BELFEM_EPSILON )
            return false;
        if ( std::abs ( aLHS.b() - aRHS.b() ) > BELFEM_EPSILON )
            return false;
        if ( std::abs ( aLHS.c() - aRHS.c() ) > BELFEM_EPSILON )
            return false;
        if ( std::abs ( aLHS.d() - aRHS.d() ) > BELFEM_EPSILON )
            return false;
        return true;
    }

    template < typename T >
    bool operator!= ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        return ! ( aLHS == aRHS );
    }

    template < typename T >
    T dot ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        return aLHS.a() * aRHS.a() + aLHS.b() * aRHS.b() + aLHS.c() * aRHS.c() + aLHS.d() * aRHS.d();
    }

    template < typename T >
    Quaternion < T > cross ( const Quaternion < T > & aLHS, const Quaternion < T > & aRHS )
    {
        return { T ( 0 ), aLHS.c() * aRHS.d() - aLHS.d() * aRHS.c(), aLHS.d() * aRHS.b() - aLHS.b() * aRHS.d(), aLHS.b() * aRHS.c() - aLHS.c() * aRHS.b() };
    }

    /**
     * @brief Spherical linear interpolation between two quaternions
     * @param aQ1 The starting quaternion (t=0)
     * @param aQ2 The ending quaternion (t=1)
     * @param aT  Interpolation parameter in [0, 1]
     * @return Interpolated quaternion
     *
     * Interpolates along the shortest arc on the 4D unit sphere.
     * Both input quaternions should be unit quaternions for meaningful results.
     */
    template < typename T >
    Quaternion < T > slerp ( const Quaternion < T > & aQ1, const Quaternion < T > & aQ2, T aT )
    {
        BELFEM_ASSERT( std::abs( aQ1.norm() - T( 1 ) ) < T( 100 ) * BELFEM_EPSILON,
                       "slerp: aQ1 must be a unit quaternion" );
        BELFEM_ASSERT( std::abs( aQ2.norm() - T( 1 ) ) < T( 100 ) * BELFEM_EPSILON,
                       "slerp: aQ2 must be a unit quaternion" );

        // compute cosine of angle between quaternions
        T tCosTheta = dot( aQ1, aQ2 );

        // if dot product is negative, negate one quaternion to take the shorter path
        Quaternion < T > tQ2 = aQ2;
        if ( tCosTheta < T ( 0 ) )
        {
            tQ2 = aQ2 * T ( -1 );
            tCosTheta = -tCosTheta;
        }

        // if quaternions are very close, use linear interpolation to avoid division by zero
        if ( tCosTheta > T ( 1 ) - BELFEM_EPSILON )
        {
            Quaternion < T > tResult = aQ1 * ( T ( 1 ) - aT ) + tQ2 * aT;
            tResult.normalize();
            return tResult;
        }

        // standard slerp formula
        T tTheta = std::acos( tCosTheta );
        T tSinTheta = std::sin( tTheta );

        T tW1 = std::sin( ( T ( 1 ) - aT ) * tTheta ) / tSinTheta;
        T tW2 = std::sin( aT * tTheta ) / tSinTheta;

        return aQ1 * tW1 + tQ2 * tW2;
    }
}
#endif //CL_QUATERNION_HPP
