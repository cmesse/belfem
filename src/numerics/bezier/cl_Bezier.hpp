//
// Created by Christian Messe on 30.11.20.
//

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
        Vector< real > mWork ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Bezier( const real & aX0, const real & aY0, const real & adYdx0,
                const real & aX1, const real & aY1, const real & adYdx1,
                const BezierType aType = BezierType::Horizontal );

//------------------------------------------------------------------------------

        ~Bezier() = default;

//------------------------------------------------------------------------------

        /**
         * expose the x-coordinates of the basis
         */
         Vector< real > &
         basis_x();

//------------------------------------------------------------------------------

        /**
         * expose the y-coordinates of the basis
         */
        Vector< real > &
        basis_y();

//------------------------------------------------------------------------------

        /**
         * The parameter coordinate goes from -1 <= xi <= 1.
         * This funciton computes the spatial x-coordinate
         * @param  aXi
         * @return aX
         */
        real
        x_by_xi( const real & aXi );

//------------------------------------------------------------------------------

        /**
         * The parameter coordinate goes from -1 <= xi <= 1.
         * This funciton computes the spatial y-coordinate
         * @param  aXi
         * @return aY
         */
        real
        y_by_xi( const real & aXi );

//------------------------------------------------------------------------------

       /**
         * Inverts the function x_by_xi. Xi must be within the defined range
         * @param  aX
         * @return aXi
         */
        real
        xi_by_x( const real & aX );

//------------------------------------------------------------------------------

        /**
          * Inverts the function y_by_xi. Xi must be within the defined range
          * @param  aY
          * @return aXi
          */
        real
        xi_by_y( const real & aY );

//------------------------------------------------------------------------------

        /**
         * x-coordinate as funciton of the y-coordinate
         * @param aY
         * @return
         */
        real
        x( const real & aY );

//------------------------------------------------------------------------------

        /**
         * y-coordinate as funciton of the x-coordinate
         * @param aX
         * @return
         */
        real
        y( const real & aX );

//------------------------------------------------------------------------------

        /**
         * point as function of parameter coordinate
         * @param aXi
         * @param aX
         * @param aY
         */
        void
        point( const real & aXi, real & aX, real & aY );

//------------------------------------------------------------------------------

        /**
         * tangent vector as function of parameter coordinate
         * @param aXi
         * @param aX
         * @param aY
         */
        void
        dpoint( const real & aXi, real & adXdXi, real & adYdXi );

//------------------------------------------------------------------------------

        /**
         * curvature vector as function of parameter coordinate
         * @param aXi
         * @param aX
         * @param aY
         */
        void
        ddpoint( const real & aXi, real & ad2XdXi2, real & ad2YdXi2 );

//------------------------------------------------------------------------------

        /**
         * derivative of Y with respect to X
         * @param aX
         * @return aY
         */
        real
        dydx( const real & aX );

//------------------------------------------------------------------------------

        /**
         * derivative of X with respect to Y
         * @param aY
         * @return aX
         */
        real
        dxdy( const real & aY );

//------------------------------------------------------------------------------

        /**
         * second derivative of Y with respect to X
         * @param aX
         * @return aY
         */
        real
        d2ydx2( const real & aX );

//------------------------------------------------------------------------------

        /**
         * second derivative of X with respect to Y
         * @param aX
         * @return aY
         */
        real
        d2xdy2( const real & aY );

//------------------------------------------------------------------------------

        /**
         * computes the length of the curve
         */
         real
         compute_length( const uint aNumIntegrationPoints=21 );

//------------------------------------------------------------------------------

        /**
         * computes the length of the curve, but provide points and weights
         * note that we need the doubles here, not the reals
         */
        real
        compute_length(
                const Vector< double > & aW,
                const Vector< double > & aXi );

//------------------------------------------------------------------------------


        void
        compute_basis_xwise(
                const real & aX0, const real & aY0,const real & adYdx0,
                const real & aX1, const real & aY1, const real & adYdx1 );

//------------------------------------------------------------------------------

        void
        compute_basis_ywise(
                const real & aX0, const real & aY0,const real & adYdx0,
                const real & aX1, const real & aY1, const real & adYdx1 );
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

        // the Bezier type has a simple shape function
        void
        compute_N( const real & aXi );

//------------------------------------------------------------------------------

        // first derivative of the shape function
        void
        compute_dNdXi( const real & aXi );

//------------------------------------------------------------------------------

        // second derivative of the shape functionon
        void
        compute_d2NdXi2( const real & aXi );

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------

    inline Vector< real > &
    Bezier::basis_x()
    {
        return mX ;
    }

//------------------------------------------------------------------------------

    inline Vector< real > &
    Bezier::basis_y()
    {
        return mY ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_BEZIER_HPP
