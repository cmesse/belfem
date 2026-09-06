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

#include "cl_Bezier.hpp"
#include "assert.h"
#include "fn_dot.hpp"
#include "intpoints.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    Bezier::Bezier(
            const real aX0, const real aY0, const real adYdx0,
            const real aX1, const real aY1, const real adYdx1,
            const BezierType aType )
    {
        // allocate memory
        mX.set_size( 4 );
        mY.set_size( 4 );
        mWork.set_size( 4 );

        // compute the coordiantes for the basis
        switch( aType )
        {
            case( BezierType::Horizontal ) :
            {
                this->compute_basis_xwise( aX0, aY0, adYdx0,
                                           aX1, aY1, adYdx1 );
                break ;
            }
            case( BezierType::Vertical ) :
            {
                this->compute_basis_ywise( aX0, aY0, adYdx0,
                                           aX1, aY1, adYdx1 );
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid Bezier Type" );
            }
        }
    }

    Bezier::Bezier( const Vector< real > & aX,
               const Vector< real > & aY,
               const BezierType aType )
    {
        mX = aX ;
        mY = aY ;
        mWork.set_size( 4 );
    }

//------------------------------------------------------------------------------

    real
    Bezier::x_by_xi( const real aXi ) const
    {
        // compute the shape function
        this->compute_N( aXi );

        // return the interpolated value
        return dot( mWork, mX );
    }

//------------------------------------------------------------------------------

    real
    Bezier::y_by_xi( const real aXi ) const
    {
        // compute the shape function
        this->compute_N( aXi );

        // return the interpolated value
        return dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    real
    Bezier::xi_by_x( const real aX ) const
    {
        // The curve is only defined for -1 <= xi <= 1, and a Bezier interpolates
        // its first and last control point, so mX( 0 ) and mX( 3 ) are the ends
        // of the x-range. Searching outside that range is not merely inaccurate:
        // if the outer control-point spans are lopsided, the cubic folds back on
        // itself just past the end, the bisection below sees a second root, and
        // it can return a xi from the wrong branch without any sign of trouble.
        // Queries outside the range therefore saturate at the end of the curve.
        if ( aX <= mX( 0 ) ) return -1.0 ;
        if ( aX >= mX( 3 ) ) return  1.0 ;

        real tXi0 = -1.0;
        real tF0 = this->x_by_xi( tXi0 ) - aX ;

        real tXi1 =  1.0;
        //real tF1 = this->x_by_xi( tXi1 ) - aX ;

        real aXi = BELFEM_QUIET_NAN ;
        real tF = BELFEM_REAL_MAX ;

        uint tCount = 0 ;

        // start bisection
        while( std::abs( tF ) > 1e-12 )
        {
            // aXi = tXi0 - tF0 * ( tXi1 - tXi0 ) / ( tF1 - tF0 );
            aXi = 0.5 * ( tXi0 + tXi1 );

            tF = this->x_by_xi( aXi ) - aX ;

            if ( tF * tF0 > 0.0 )
            {
                tXi0 = aXi ;
                tF0 = tF ;
            }
            else
            {
                tXi1 = aXi ;
                //tF1 = tF ;
            }

            BELFEM_ERROR( tCount++ < 100, "Too many iterations.");
        }

        return aXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::xi_by_y( const real aY ) const
    {
        // same saturation as in xi_by_x, see the comment there
        if ( aY <= mY( 0 ) ) return -1.0 ;
        if ( aY >= mY( 3 ) ) return  1.0 ;

        real tXi0 = -1.0;
        real tF0 = this->y_by_xi( tXi0 ) - aY ;

        real tXi1 =  1.0;
        //real tF1 = this->y_by_xi( tXi1 ) - aY ;

        real aXi = BELFEM_QUIET_NAN ;
        real tF = BELFEM_REAL_MAX ;

        uint tCount = 0 ;

        // start bisection
        while( std::abs( tF ) > 1e-12 )
        {
            //aXi = tXi0 - tF0 * ( tXi1 - tXi0 ) / ( tF1 - tF0 );
            aXi = 0.5 * ( tXi0 + tXi1 );
            tF = this->y_by_xi( aXi ) - aY ;

            if ( tF * tF0 > 0.0 )
            {
                tXi0 = aXi ;
                tF0 = tF ;
            }
            else
            {
                tXi1 = aXi ;
                //tF1 = tF ;
            }

            BELFEM_ERROR( tCount++ < 100, "Too many iterations.");
        }

        return aXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::x( const real aY ) const
    {
        return this->x_by_xi( this->xi_by_y( aY ) );
    }

//------------------------------------------------------------------------------

    real
    Bezier::y( const real aX ) const
    {
        return this->y_by_xi( this->xi_by_x( aX ) );
    }

//------------------------------------------------------------------------------

    void
    Bezier::point( const real aXi, real & aX, real & aY ) const
    {
        this->compute_N( aXi );

        aX = dot( mWork, mX );
        aY = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    void
    Bezier::dpoint( const real aXi, real & adXdXi, real & adYdXi ) const
    {
        this->compute_dNdXi( aXi );

        adXdXi = dot( mWork, mX );
        adYdXi = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    void
    Bezier::ddpoint( const real aXi, real & ad2XdXi2, real & ad2YdXi2 ) const
    {
        this->compute_d2NdXi2( aXi );

        ad2XdXi2 = dot( mWork, mX );
        ad2YdXi2 = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    void
    Bezier::dddpoint( real & ad3XdXi3, real & ad3YdXi3 ) const
    {
        // d3/dxi3 of a cubic Bezier is constant. The basis is defined on
        // xi in [-1,1], hence the 1/8 relative to the [0,1] Bernstein form
        mWork = { -0.75, 2.25, -2.25, 0.75 };

        ad3XdXi3 = dot( mWork, mX );
        ad3YdXi3 = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    real
    Bezier::dydx( const real aX ) const
    {
        real tdXdXi ;
        real tdYdXi ;
        this->dpoint( this->xi_by_x( aX ), tdXdXi, tdYdXi ) ;
        return tdYdXi / tdXdXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::dxdy( const real aY ) const
    {
        real tdXdXi ;
        real tdYdXi ;
        this->dpoint( this->xi_by_y( aY ), tdXdXi, tdYdXi ) ;

        return tdXdXi / tdYdXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::d2ydx2( const real aX ) const
    {
        real xi = this->xi_by_x( aX );

        real dxdXi, dydXi;
        real d2xdXi2, d2ydXi2;

        this->dpoint( xi, dxdXi, dydXi );
        this->ddpoint( xi, d2xdXi2, d2ydXi2 );

        // Avoid division by zero (though xi_by_x should prevent this)
        if ( std::abs( dxdXi ) < BELFEM_EPSILON )
        {
            return 0.0;
        }

        real denom = dxdXi * dxdXi * dxdXi;  // (dx/dξ)^3

        return ( d2ydXi2 * dxdXi - dydXi * d2xdXi2 ) / denom;
    }

//------------------------------------------------------------------------------

    real
    Bezier::d2xdy2( const real aY ) const
    {
        real xi = this->xi_by_y( aY );
        real dxdXi, dydXi;
        real d2xdXi2, d2ydXi2;

        this->dpoint( xi, dxdXi, dydXi );
        this->ddpoint( xi, d2xdXi2, d2ydXi2 );

        real denom = dydXi * dydXi * dydXi;

        return ( d2xdXi2 * dydXi - dxdXi * d2ydXi2 ) / denom;
    }

//------------------------------------------------------------------------------

    /*
    * third derivative of Y with respect to X (documented in cl_Bezier.hpp)
    */
    real
    Bezier::d3ydx3( const real aX ) const
    {
        real xi = this->xi_by_x( aX );

        real dxdXi, dydXi;
        real d2xdXi2, d2ydXi2;
        real d3xdXi3, d3ydXi3 ;

        this->dpoint( xi, dxdXi, dydXi );
        this->ddpoint( xi, d2xdXi2, d2ydXi2 );
        this->dddpoint( d3xdXi3, d3ydXi3 );

        if ( std::abs( dxdXi ) < BELFEM_EPSILON )
        {
            return 0.0;
        }

        // d3y/dx3 = [ y''' x'^2 - 3 y'' x' x'' + y' ( 3 x''^2 - x' x''' ) ] / x'^5
        // with primes denoting derivatives with respect to xi
        real denom = dxdXi * dxdXi ;
        denom *= denom ;
        denom *= dxdXi ;

        return (   d3ydXi3 * dxdXi * dxdXi
             - 3. * d2ydXi2 * dxdXi * d2xdXi2
             + dydXi * ( 3. * d2xdXi2 * d2xdXi2 - dxdXi * d3xdXi3 ) ) / denom ;
    }

//------------------------------------------------------------------------------

    /*
    * third derivative of X with respect to Y (documented in cl_Bezier.hpp)
    */
    real
    Bezier::d3xdy3( const real aY ) const
    {
        real xi = this->xi_by_y( aY );

        real dxdXi, dydXi;
        real d2xdXi2, d2ydXi2;
        real d3xdXi3, d3ydXi3 ;

        this->dpoint( xi, dxdXi, dydXi );
        this->ddpoint( xi, d2xdXi2, d2ydXi2 );
        this->dddpoint( d3xdXi3, d3ydXi3 );

        if ( std::abs( dydXi ) < BELFEM_EPSILON )
        {
            return 0.0;
        }

        // d3x/dy3 = [ x''' y'^2 - 3 x'' y' y'' + x' ( 3 y''^2 - y' y''' ) ] / y'^5
        // with primes denoting derivatives with respect to xi
        real denom = dydXi * dydXi ;
        denom *= denom ;
        denom *= dydXi ;

        return (   d3xdXi3 * dydXi * dydXi
             - 3. * d2xdXi2 * dydXi * d2ydXi2
             + dxdXi * ( 3. * d2ydXi2 * d2ydXi2 - dydXi * d3ydXi3 ) ) / denom ;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_basis_xwise(
            const real aX0, const real aY0,const real adYdx0,
            const real aX1, const real aY1, const real adYdx1 )
    {
        real tDX = ( aX1 - aX0 ) / 3. ;

        mX( 0 ) = aX0 ;
        mX( 1 ) = aX0 + tDX ;
        mX( 2 ) = aX1 - tDX ;
        mX( 3 ) = aX1 ;

        mY( 0 ) = aY0 ;
        mY( 1 ) = aY0 + tDX * adYdx0 ;
        mY( 2 ) = aY1 - tDX * adYdx1 ;
        mY( 3 ) = aY1 ;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_basis_ywise(
            const real aX0, const real aY0,const real adYdx0,
            const real aX1, const real aY1, const real adYdx1 )
    {
        real tDY = ( aY1 - aY0 ) / 3. ;

        mX( 0 ) = aX0 ;
        mX( 1 ) = aX0 - tDY * adYdx0 ;
        mX( 2 ) = aX1 + tDY * adYdx1 ;
        mX( 3 ) = aX1 ;

        mY( 0 ) = aY0 ;
        mY( 1 ) = aY0 + tDY ;
        mY( 2 ) = aY1 - tDY ;
        mY( 3 ) = aY1 ;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_N( const real aXi ) const
    {
        mWork( 0 ) = std::pow(1 - aXi, 3);
        mWork( 1 ) = 3. * ((aXi * (aXi - 1.) - 1.) * aXi + 1.);
        mWork( 2 ) = 3. * ((1. - (aXi + 1.) * aXi) * aXi + 1.);
        mWork( 3 ) = ((aXi * (3. + aXi) + 3.) * aXi + 1.);

        mWork *= 0.125;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_dNdXi( const real aXi ) const
    {
        mWork( 0 ) = - (1. - aXi) * (1. - aXi);
        mWork( 1 ) =  ((3. * aXi - 2.) * aXi - 1.);
        mWork( 2 ) =    1. - aXi * (2. + 3. * aXi);
        mWork( 3 ) =   (1. + aXi) * (1. + aXi);

        mWork *= 0.375;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_d2NdXi2( const real aXi ) const
    {
        mWork( 0 ) =  1. - aXi;
        mWork( 1 ) =  3. * aXi - 1.;
        mWork( 2 ) = -3. * aXi - 1.;
        mWork( 3 ) =  aXi + 1.;

        mWork *= 0.75;
    }

//------------------------------------------------------------------------------

    /**
     * computes the length of the curve
     */
    real
    Bezier::compute_length( const uint aNumIntegrationPoints ) const
    {
        // number of points
        int tN = aNumIntegrationPoints ;

        // integration weights
        Vector< double > tW( aNumIntegrationPoints );

        // integration points
        Vector< double > tXi( aNumIntegrationPoints );

        // compute the integration points
        intpoints_gauss(
                &tN,
                tW.data(),
                tXi.data() );

        return this->compute_length( tW, tXi );
    }

//------------------------------------------------------------------------------

    /**
     * computes the length of the curve, but provide points and weights
     */
    real
    Bezier::compute_length(
            const Vector< double > & aW,
            const Vector< double > & aXi ) const
    {
        // check input
        BELFEM_ASSERT( aW.length() == aXi.length(),
                      "Length of weights and points does not match ( %u vs %u )",
                      ( unsigned int ) aW.length(),
                      ( unsigned int ) aXi.length() );

        uint tN = aW.length() ;

        // integrate
        real aLength = 0.0 ;
        real tdXdXi ;
        real tdYdXi ;
        for( uint k=0; k<tN; ++k )
        {
            this->dpoint( aXi( k ), tdXdXi, tdYdXi ) ;

            aLength += std::sqrt( tdXdXi * tdXdXi + tdYdXi * tdYdXi );
        }

        return aLength ;
    }
//------------------------------------------------------------------------------
}