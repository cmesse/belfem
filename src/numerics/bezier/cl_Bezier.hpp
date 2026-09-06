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

#ifndef BELFEM_CL_BEZIER_HPP
#define BELFEM_CL_BEZIER_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    enum class BezierType
    {
        Horizontal,
        Vertical
    };

//------------------------------------------------------------------------------

    class Bezier
    {
        // X-coordinates
        Vector< real > mX ;

        // Y-Coordinates
        Vector< real > mY ;

        // Work vector
        mutable Vector< real > mWork ;


//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Bezier( const real aX0=0.0, const real aY0=0.0, const real adYdx0=0.0,
                const real aX1=1.0, const real aY1=1.0, const real adYdx1=0.0,
                const BezierType aType = BezierType::Horizontal );

        Bezier( const Vector< real > & aX,
                const Vector< real > & aY,
                const BezierType aType = BezierType::Horizontal );

//------------------------------------------------------------------------------

        ~Bezier() = default;

//------------------------------------------------------------------------------

        /**
         * expose the x-coordinates of the basis
         */
         Vector< real > &
         basis_x();

         const Vector< real > &
         basis_x() const ;

//------------------------------------------------------------------------------

        /**
         * expose the y-coordinates of the basis
         */
        Vector< real > &
        basis_y();

        const Vector< real > &
        basis_y() const ;

//------------------------------------------------------------------------------

        /**
         * The parameter coordinate goes from -1 <= xi <= 1.
         * This funciton computes the spatial x-coordinate
         * @param  aXi
         * @return aX
         */
        real
        x_by_xi( const real aXi ) const ;

//------------------------------------------------------------------------------

        /**
         * The parameter coordinate goes from -1 <= xi <= 1.
         * This funciton computes the spatial y-coordinate
         * @param  aXi
         * @return aY
         */
        real
        y_by_xi( const real aXi ) const ;

//------------------------------------------------------------------------------

       /**
         * Inverts the function x_by_xi. Queries outside [ mX(0), mX(3) ] saturate to xi = -1 / xi = +1.
         * @param  aX
         * @return aXi
         */
        real
        xi_by_x( const real aX ) const ;

//------------------------------------------------------------------------------

        /**
          * Inverts the function y_by_xi. Queries outside [ mY(0), mY(3) ] saturate to xi = -1 / xi = +1.
          * @param  aY
          * @return aXi
          */
        real
        xi_by_y( const real aY ) const ;

//------------------------------------------------------------------------------

        /**
         * x-coordinate as funciton of the y-coordinate
         * @param aY
         * @return
         */
        real
        x( const real aY ) const ;

//------------------------------------------------------------------------------

        /**
         * y-coordinate as funciton of the x-coordinate
         * @param aX
         * @return
         */
        real
        y( const real aX ) const ;

//------------------------------------------------------------------------------

        /**
         * point as function of parameter coordinate
         * @param aXi
         * @param aX
         * @param aY
         */
        void
        point( const real aXi, real & aX, real & aY ) const ;

//------------------------------------------------------------------------------

        /**
         * tangent vector as function of parameter coordinate
         * @param aXi
         * @param adXdXi
         * @param adYdXi
         */
        void
        dpoint( const real aXi, real & adXdXi, real & adYdXi ) const ;

//------------------------------------------------------------------------------

        /**
         * curvature vector as function of parameter coordinate
         * @param aXi
         * @param ad2XdXi2
         * @param ad2YdXi2
         */
        void
        ddpoint( const real aXi, real & ad2XdXi2, real & ad2YdXi2 ) const ;

//------------------------------------------------------------------------------

        /**
         * aberrancy vector as function of parameter coordinate
         * @param ad3XdXi3
         * @param ad3YdXi3
         */
        void
        dddpoint( real & ad3XdXi3, real & ad3YdXi3 ) const ;


//------------------------------------------------------------------------------

        /**
         * derivative of Y with respect to X
         * @param aX
         * @return aY
         */
        real
        dydx( const real aX ) const ;

//------------------------------------------------------------------------------

        /**
         * derivative of X with respect to Y
         * @param aY
         * @return aX
         */
        real
        dxdy( const real aY ) const ;

//------------------------------------------------------------------------------

        /**
         * second derivative of Y with respect to X
         * @param aX
         * @return aY
         */
        real
        d2ydx2( const real aX ) const ;

//------------------------------------------------------------------------------

        /**
         * second derivative of X with respect to Y
         * @param aY
         * @return second derivative of X
         */
        real
        d2xdy2( const real aY ) const ;

//------------------------------------------------------------------------------

        /**
        * third derivative of Y with respect to X
        * @param aX
        * @return aY
        */
        real
        d3ydx3( const real aX ) const ;

//------------------------------------------------------------------------------

        /**
        * third derivative of X with respect to Y
        * @param aY
        * @return aX
        */
        real
        d3xdy3( const real aY ) const ;

//------------------------------------------------------------------------------

        /**
         * computes the length of the curve
         */
         real
         compute_length( const uint aNumIntegrationPoints=21 ) const ;

//------------------------------------------------------------------------------

        /**
         * computes the length of the curve, but provide points and weights
         * note that we need the doubles here, not the reals
         */
        real
        compute_length(
                const Vector< double > & aW,
                const Vector< double > & aXi ) const ;

//------------------------------------------------------------------------------


        void
        compute_basis_xwise(
                const real aX0, const real aY0,const real adYdx0,
                const real aX1, const real aY1, const real adYdx1 ) ;

//------------------------------------------------------------------------------

        void
        compute_basis_ywise(
                const real aX0, const real aY0,const real adYdx0,
                const real aX1, const real aY1, const real adYdx1 ) ;
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

        // the Bezier type has a simple shape function
        void
        compute_N( const real aXi ) const ;

//------------------------------------------------------------------------------

        // first derivative of the shape function
        void
        compute_dNdXi( const real aXi ) const ;

//------------------------------------------------------------------------------

        // second derivative of the shape functionon
        void
        compute_d2NdXi2( const real aXi ) const ;

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------

    inline Vector< real > &
    Bezier::basis_x()
    {
        return mX ;
    }

    inline const Vector< real > &
    Bezier::basis_x() const
    {
        return mX ;
    }


//------------------------------------------------------------------------------

    inline Vector< real > &
    Bezier::basis_y()
    {
        return mY ;
    }

    inline const Vector< real > &
    Bezier::basis_y() const
    {
        return mY ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_BEZIER_HPP
